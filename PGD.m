function [x,timer,error_t]=PGD(xs,K,n1,n2,nd,r,p,tol,max_iter,opt,stepsize,x_star)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs 
%          xs: observtions.
%           K: observation indices
%          n1: 1st dim of Hankel matrix
%          n2: 2st dim of Hankel matrix
%          nd:  length of the 1-dim signal
%           r: model order of the spectrally sparse signal
%           p: sample rate
%         tol: convergence tolerance
%    max_iter: maximum number of iterations.
%         opt: options for stepsize strategy. 1 (line search), 0 (fixed).
%    stepsize: initial stepsize.
%     x_star : true signal
%Outputs
%           x: recoverd signal 
%       timer: run time versus iterations
%     error_t: recovery error versur iterations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
proj=true;
T = max_iter;

thresh_up = 1e2;
%thresh_low = 1e-12;
% errors_GD = zeros(1, T);
omega = zeros(nd, 1);
omega(K,:) = 1 ;
w = zeros(nd, 1);
for k = 1:nd
    w(k) = min([k, n1+n2-k, n1, n2]); %length of skew-diagonals
end
%     kappa = kappa_list(i_kappa);
%     sigma_star = linspace(1, 1/kappa, r);
%     x_star = exp(2*pi*1i * (0:(nd-1))' * freq_seed) * sigma_star' / sqrt(n1*n2);
    omega(K,:) = 1 ;
     timer = zeros(T,1);
    error_t = zeros(T,1);
    xerr_storage = zeros(T,nd);
    y = xs;
    %% Spectral initialization via Lanczos algorithm
    [U0, Sigma0, V0] = svds(@(v,tflag)HankelVecMul(y/p, v, tflag, nd, n1, n2), [n1, n2], r);
    mu=sqrt(2*max(sum(abs(U0).^2,2)));
    sqrtsigma=sqrt(max(Sigma0,[],'all'));
    step = zeros(T,1);
    %% PGD
    L = U0*sqrt(Sigma0);
    R = V0*sqrt(Sigma0);
    lambda=1/4;    
    if proj
    L = projop(L,mu,sqrtsigma);
    R = projop(R,mu,sqrtsigma);
    end
    tic
    x = HankelProj(L, R, nd, w); %this operator is H*, not G*, but doesn't matter, as long as the following operator is H(because HH*=GG*). 
   
    for t = 1:T
        xerr_storage(t,:) =  x(1:nd)-x_star;
        z = (omega.*x-y)/p - x;
        [Lz, Rz] = HankelMul(z, L, R, nd, n1, n2);
        Lzv = lambda*L*(L'*L)+(1-lambda)*L*(R'*R) + Lz;
        Rzv = lambda*R*(R'*R)+(1-lambda)*R*(L'*L) + Rz;        
        x_last=x;
        if opt
        gd_nm2 = norm(Lzv,'fro')^2+norm(Rzv,'fro')^2;
        f_val_current = norm((omega.*x-y).*sqrt(w))^2/p;    
        k = 1;
        l = 0; % number of line search
        eta=1.5/Sigma0(1,1);
        while k
        L_update = L - eta*Lzv;
        R_update = R - eta*Rzv;        
        x = HankelProj(L_update, R_update, nd, w);
        f_val_update = norm((omega.*x-y).*sqrt(w))^2/p;               
            
        if f_val_update < f_val_current-1/10*eta*gd_nm2
            k = 0;
        else
            eta = 3*eta/4;
            l = l+1;
        end 
        if l > 5
            break;
        end
        end
        step(t) = 1.5*(3/4)^l;
        else
        eta=stepsize/Sigma0(1,1);   
        L_update = L - eta*Lzv;
        R_update = R - eta*Rzv;   
        x = HankelProj(L_update, R_update, nd, w);
       step(t) = stepsize;
        end
        L_plus = L_update;
        R_plus = R_update;
        
        if proj
            L_plus = projop(L_plus,mu,sqrtsigma);
            R_plus = projop(R_plus,mu,sqrtsigma);
        end
        L = L_plus;
        R = R_plus;     
        timer(t) = toc;   
        error_t(t) = norm(xerr_storage(t,:))/norm(x_star);
        upd_error = norm((x - x_last))/norm(x_last);

%         error=norm([L_plus]-[L],'fro');
        if ~isfinite(upd_error) || upd_error < tol || upd_error > thresh_up
            break;
        end        
        
      end
function x = HankelProj(L, R, nd, w)
% Hankel projection H[x]=H[LR']
    n = 2^nextpow2(nd);
    x = sum(ifft(fft(L, n) .* fft(conj(R), n)), 2);
    x = x(1:nd)./w;
end

function [Lx, Rx] = HankelMul(x, L, R, nd, n1, n2)
% Hankel multiplication Lx=H[x]R, Rx=H[x]'L
    n = 2^nextpow2(nd);
    Lx = ifft(bsxfun(@times, fft(flip(R), n), fft(x, n)));
    Lx = Lx(n2:nd, :);  
    Rx = ifft(bsxfun(@times, fft(flip(L), n), fft(conj(x), n)));
    Rx = Rx(n1:nd, :);
end


function z = HankelVecMul(x, v, tflag, nd, n1, n2)
% Hankel multiplication z=H[x]v if tflag='notransp'; z=H[x]'v if tflag='transp'
    n = 2^nextpow2(nd);
    if strcmp(tflag, 'notransp')
        z = ifft(fft(flip(v), n) .* fft(x, n));
        z = z(n2:nd);
    else
        z = ifft(fft(flip(v), n) .* fft(conj(x), n));
        z = z(n1:nd);
    end
end


function [Z] = projop(Z,mu,sqrtsig)
% 2,infty incoherence projection Z=P_C(Z);
incoh=mu*sqrtsig;
Z_row_norm=sqrt(sum(abs(Z).^2,2));
bigrow_Z=Z_row_norm>incoh;

Z(bigrow_Z,:) = bsxfun(@times,Z(bigrow_Z,:),(incoh ./ Z_row_norm(bigrow_Z))); 
end

  

end