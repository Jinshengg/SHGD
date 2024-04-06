function [x,timer,error_t]=SHGD(xs,K,nd,r,p,tol,max_iter,opt,stepsize,x_star)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs 
%          xs: observtions.
%           K: observation indices
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
% eta = 0.5;
thresh_up = 1e2;
omega = zeros(nd, 1);
 omega(K,:) = 1 ;
 timer = zeros(T,1);
 error_t = zeros(T,1);
 xerr_storage = zeros(T,nd);
 y=xs;
%% SHGD
    %% preprocess (zero-padding) to form a square Hankel matrix
    if mod(nd,2)==0
    n_s=(nd+2)/2;
    y_s=[y;0];
    omegas=[omega;0];
    else
    n_s=(nd+1)/2;
    y_s=y;
    omegas=[omega];
    end

    nds=2*n_s-1;
    ws = zeros(nds, 1);
for k = 1:nds
    ws(k) = min([k, n_s+n_s-k, n_s, n_s]); %length of skew-diagonals ¡ª¡ªsquare case
end
  %% Spectral initialization
  % SVD of the initial estimated Hankel matrix
    [U0, Sigma0, V0] = svds(@(v,tflag)HankelVecMul(y_s/p, v, tflag, nds, n_s, n_s), [n_s, n_s], r);
    % From SVD to Takagi factorization.
    % Refer to equation (5) and (7) in  "Singular value decompossition for the Takagi factorization of symmetric matrices", Alexander M. Chebotarev,  Alexander E. 
    % Apply fast computation version supposing singular values are different.
     U = zeros(n_s,r);
     for i=1:r
        U(:,i)=sqrt(U0(:,i)'*conj(V0(:,i)))*U0(:,i);
     end
     Z=U*sqrt(Sigma0);

    % Projection
    mu=sqrt(2*max(sum(abs(U0).^2,2)));
    sqrtsigma=sqrt(max(Sigma0,[],'all'));   
    if proj
    Z=projop( Z,mu,sqrtsigma);
    end
    step = zeros(T,1);
    
    tic
    % H^*(ZZ^T)
    x = HankelProj(Z, nds, ws); %this operator is H*, not G*, but doesn't matter, as long as the following operator is H ( HH*=GG*).
    for t = 1:T   
%         error = norm((x(1:nd,:) - x_star).*sqrt(w))/norm(x_star.*sqrt(w));
%         errors_SHGD(i_kappa, t) = error;   
        xerr_storage(t,:) =  x(1:nd)-x_star;
        z = (omegas.*x-y_s)/p - x;
        Zz = HankelMul_S(z, Z, nds, n_s);
        Zzv = Z*(Z.'*conj(Z)) + Zz;    %gradient  
        x_last=x;
        
       
        if opt 
        % backtracking ;
        gd_nm2 = norm(Zzv,'fro')^2;
        f_val_current = norm((omegas.*x-y_s).*sqrt(ws))^2/p; 
        k = 1;
        l = 0; % number of line search
        eta=1.5/Sigma0(1,1);
        while k
        Z_update = Z - eta*Zzv;  
        x = HankelProj(Z_update, nds, ws);
        f_val_update = norm((omegas.*x-y_s).*sqrt(ws))^2/p;                          
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
        % fixed stepsized
        eta=stepsize/Sigma0(1,1);   
        Z_update = Z - eta*Zzv;
        x = HankelProj(Z_update, nds, ws);
       step(t) = stepsize;
        end
        Z_plus = Z_update;
         
        % projection
        if proj
             Z_plus=projop(Z_plus,mu,sqrtsigma);
        end       
        Z = Z_plus;
        timer(t) = toc;   

        upd_error = norm((x - x_last))/norm(x_last);
        error_t(t) = norm(xerr_storage(t,:))/norm(x_star);
        if ~isfinite(upd_error)|| upd_error < tol ||upd_error>thresh_up
            break;
        end
    end
    


function x = HankelProj(Z, nd, w)
% Hankel projection H[x]=H[ZZ^T]
    n = 2^nextpow2(nd);
    x = sum(ifft(fft(Z, n) .* fft(Z, n)), 2);
    x = x(1:nd)./w;
end

function [Zx] = HankelMul_S(x, Z, nd, n)
% Hankel multiplication Zx=H[z]conj(Z)
    nn = 2^nextpow2(nd);
    Zx = ifft(bsxfun(@times, fft(flip(conj(Z)), nn), fft(x, nn)));
    Zx = Zx(n:nd, :);  

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