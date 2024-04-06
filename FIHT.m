function [x,timer,error_t]=FIHT(xs,K,nd,r,p,tol,max_iter,x_star)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs 
%          xs: observtions.
%           K: observation indices
%          nd:  length of the 1-dim signal
%           r: model order of the spectrally sparse signal
%           p: sample rate
%         tol: convergence tolerance
%    max_iter: maximum number of iterations.
%     x_star : true signal
%Outputs
%           x: recoverd signal 
%       timer: run time versus iterations
%     error_t: recovery error versur iterations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = max_iter;
thresh_up = 2; 
y = xs;
 timer = zeros(T,1);
 error_t = zeros(T,1);
 xerr_storage = zeros(T,nd);
omega = zeros(nd, 1);
omega(K,:) = 1 ;
% w = zeros(nd, 1);
%%  preprocess to form square Hankel matrix

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
U = zeros(n_s,r);
Ucs = zeros(2*r,r);


    %% Spectral initialization via Lanczos algorithm
    [U0, Sigma0, V0] = svds(@(v,tflag)HankelVecMul(y_s/p, v, tflag, nds, n_s, n_s), [n_s, n_s], r);

     for i=1:r
        U(:,i)=sqrt(U0(:,i)'*conj(V0(:,i)))*U0(:,i);
     end
          

    mu=sqrt(2*max(sum(abs(U).^2,2)));      
    U = projop(U,mu,1); 
    Z=U*sqrt(Sigma0);  
    tic 
    x = HankelProj_2(Z,nds, ws);
    for t = 1:T
        xerr_storage(t,:) =  x(1:nd)-x_star;
        z =  x + (y_s-omegas.*x)/p;
        [Uz] = HankelMul_S(z,U, nds, n_s);
        C =  U'*Uz;
        X = Uz-U*C;
        [Q, R] = qr(X,0);
        Intermid = [C R.';R zeros(r,r) ];       
        [Uc, Sigmac, Vc] = svds(Intermid,r);    
        
        for i=1:r
        Ucs(:,i)=sqrt(Uc(:,i)'*conj(Vc(:,i)))*Uc(:,i);
        end      
        U_plus = [U Q]*Ucs;
%         U_plus = [U Q]*Uc;
        U = U_plus;
        U = projop(U,mu,1);
        Z = U*sqrt(Sigmac); 
        x_last = x; 
        x = HankelProj_2(Z,nds, ws);

        timer(t) = toc;   
        upd_error = norm((x - x_last))/norm(x_last);
        error_t(t) = norm(xerr_storage(t,:))/norm(x_star);
        if ~isfinite(upd_error) || upd_error > thresh_up || upd_error < tol
            break;
        end
        
    end
   

function x = HankelProj_2(Z, nd, w)
% Hankel projection H[x]=H[LR']
    n = 2^nextpow2(nd);
    x = sum(ifft(fft(Z, n) .* fft(Z, n)), 2);
    x = x(1:nd)./w;
end

function [Zx] = HankelMul_S(x, Z, nd, n)
% Hankel multiplication Zx=H[x]conj(Z)
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
% end
end

end