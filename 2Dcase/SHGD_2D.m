function [x ,timer,iter_count,error_t]=SHGD_2D(Xs,K,N1i,N2i,r,p,tol,max_iter,opt,stepsize,x_star)
% clear; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs 
%          Xs: observtions, 2-D matrix.
%           K: observation indices
%          N1i:  1st dim of 2D signal
%          N2i:  2st dim of 2D signal
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
omega = zeros(N1i*N2i,1);
omega(K,:) = 1 ;
Omega=reshape(omega,[N1i N2i]); 
iter_count = zeros(T,1); 
error_t = zeros(T,1);
xerr_storage = zeros(T,N1i,N2i);
timer = zeros(T,1);
thresh_up = 1e2;

%% preprocess : each dim to odd dim , then form a square 2-level block Hankel matrix

if mod(N1i,2)
    p1 = (N1i+1)/2;
    w1 = [1:p1 p1-1:-1:1]';
    N1 = N1i;
else
    %  row zero padding
    N1=N1i+1;
    p1 = (N1+1)/2;
    w1 = [1:p1 p1-1:-1:1]';
    Xs = [Xs;zeros(1,N2i)];
    Omega = [Omega;zeros(1,N2i)];
end

if mod(N2i,2)
    p2 = (N2i+1)/2;
    w2 = [1:p2 p2-1:-1:1]';
    N2 = N2i;
else
    % column zero padding
    N2=N2i+1;
    p2 = (N2+1)/2;
    w2 = [1:p2 p2-1:-1:1]';
    Xs = [Xs zeros(N1,1)];
    Omega = [Omega zeros(N1,1)];
end
y = Xs; %2D observed signal y


WW = kron(w2,w1);
WW = reshape(WW,N1,N2);
q1 = p1;
q2 = p2;

l1 = p1*p2; % 1st dim of constructed 2-level block square Hankel matrix
l2 = q1*q2; % 2st dim of constructed 2-level block square Hankel matrix

% indicies pre-computed for fhmvmultiply_2D to use
ind1 = zeros(l2,1);
for i = 1:q2
    ind1((i-1)*q1+1:i*q1) = (i-1)*N1+1:(i-1)*N1+q1;
end
ind2 = zeros(l1,1);
for i = 1:p2
    ind2((i-1)*p1+1:i*p1) = (q2+i-2)*N1+q1:(q2+i-1)*N1;
end
ind3 = zeros(l1,1);
for i = 1:p2
    ind3((i-1)*p1+1:i*p1) = (i-1)*N1+1:(i-1)*N1+p1;
end
ind4 = zeros(l2,1);
for i = 1:q2
    ind4((i-1)*q1+1:i*q1) = (p2+i-2)*N1+p1:(p2+i-1)*N1;
end

% thresh_up = 1e3;
   
    %% Spectral initialization via Lanczos algorithm
    [U0, Sigma0, V0] = svds(@(v,tflag)HankelVecMul_2D(y/p, v, tflag, p1,p2,q1,q2,ind1,ind2,ind3,ind4), [l1, l2], r);   
    U = zeros(l1,r); 
     for i=1:r
        U(:,i)=sqrt(U0(:,i)'*conj(V0(:,i)))*U0(:,i);
     end
     Z=U*sqrt(Sigma0);  
     mu=sqrt(2*max(sum(abs(U0).^2,2)));
     sqrtsigma=sqrt(max(Sigma0,[],'all'));
    
     if proj
        Z=projop(Z,mu,sqrtsigma);
     end 
     
    tic     
     x=HankelProj2D(Z, p1,p2,q1,q2,WW,r); %this operator is H*, not G*, but doesn't matter, as long as the following operator is H ( HH*=GG*).
    
    for t = 1:T

%% iteration time vs error 
        iter_count(t) = t;
        xerr_storage(t,:,:) =  x(1:N1i,1:N2i)-x_star;
        
        z = (Omega.*x-y)/p - x;  % 2D signal x y z
        
        [Zz] = HankelMul2D_S(z,Z, r,l1,q1,q2,ind1,ind2);       
        
        x_last=x;
       
        % backtracking for t=1;
        if opt 
          Zzv = Z*(Z.'*conj(Z)) + Zz;     
        gd_nm2 = norm(Zzv,'fro')^2;
        f_val_current = norm((Omega.*x-y).*sqrt(WW))^2/p;    
        k = 1;
        l = 0; % number of line search
        eta=1.5/Sigma0(1,1);
        while k
        Z_update = Z - eta*Zzv;  
        x = HankelProj2D(Z_update,  p1,p2,q1,q2,WW,r);
        f_val_update = norm((Omega.*x-y).*sqrt(WW))^2/p;               
            
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
        end
%         Z_plus = Z_update;    
        
         Z_plus = Z - eta*(Z*(Z.'*conj(Z)) + Zz);
         
        if proj
             Z_plus=projop( Z_plus,mu,sqrtsigma);
        end
        
        timer(t) = toc;    
        Z = Z_plus;    
        
        x = HankelProj2D(Z, p1,p2,q1,q2,WW,r);
       
        upd_error = norm((x - x_last))/norm(x_last);

        if ~isfinite(upd_error) || upd_error > thresh_up || upd_error < tol 
            break;
        end        
        
    end
    
    for t = 1:T
    error_t(t) = norm(reshape(xerr_storage(t,:,:),[N1i,N2i]),'fro')/norm(x_star,'fro');
    end
    

    
function x = HankelProj2D(Z, p1,p2,q1,q2,WW,r)
% Hankel projection H[x]=H[ZZ^T], p1 = q1, p2 = q2.
Z = Z(:,1:r);
n1 = p1+q1-1;
n2 = p2+q2-1;
x = zeros(n1,n2);
for tt = 1:r
    li = reshape(Z(:,tt),p1,p2);
%     ri = reshape(R(:,tt),q1,q2);
    x = x+conv_fft(li,li);
end
x = x./WW;
end

function [Zx] = HankelMul2D_S(x, Z,r, l1,q1,q2,ind1,ind2)
% Hankel multiplication Zx=H[z]conj(Z)
    Zx = zeros(l1,r);
    for tt = 1:r
        zi = Z(:,tt);
        Zx(:,tt) = fhmvmultiply_2D(x,conj(zi),q1,q2,ind1,ind2);
    end
end

function z = HankelVecMul_2D(X, v, tflag, p1,p2,q1,q2, ind1,ind2,ind3,ind4)  
% Hankel multiplication z=H[x]v if tflag='notransp'; z=H[x]'v if tflag='transp'
%     n = 2^nextpow2(nd);
    if strcmp(tflag, 'notransp')
        z = fhmvmultiply_2D(X,v,q1,q2,ind1,ind2);
    else
        z =  fhmvmultiply_2D(conj(X),v,p1,p2,ind3,ind4);
    end
end


function [Z] = projop(Z,mu,sqrtsig)
% 2,infty incoherence projection Z=P_C(Z);
incoh=mu*sqrtsig;
Z_row_norm=sqrt(sum(abs(Z).^2,2));
bigrow_Z=Z_row_norm>incoh;
% sprintf('proj');
Z(bigrow_Z,:) = bsxfun(@times,Z(bigrow_Z,:),(incoh ./ Z_row_norm(bigrow_Z))); 
end

  

end