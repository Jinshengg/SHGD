  function [xs,K,x_star,f,amplitude]=generate_signal_2D(m,nd,r,seperation,damp)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [xs,K,x_star]=generate_signal_2D(m,nd,r,seperation,damp)
% 
% Inputs:
% m     : num of observsyions  
% nd    : dim of 2-D signal, nd(1) and nd(2)
% r     : model order of 2-D spectral sparse signal
%seperation : frequency seperation£ºfalse (0) or true (1)
%damp: :dampling para: false or true
% Outputs:
% xs : observed 2D spectrally r-sparse signal.
% K : observed indices 
% x_star : groudtruth data matrix as column vector form
% f: ground truth 2D-frequency
% amplitude: amplitudes of ground truth 2D-frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate signals
    dynamic_range=10;
    c=exp(1i*2*pi*rand(r,1)).*(1+10.^(rand(r,1)*dynamic_range/20));
    amplitude = c;
%     K=randsample(nd,m);
      
    N1 = nd(1); s1 = 0:N1-1;
    N2 = nd(2); s2 = 0:N2-1;   
    xs=zeros(N1*N2,1);
    K = randsample(N1*N2,m); % sample without replacement
        switch seperation
            case {false,'false',0} % without separation 
                f1 = rand(r,1);
                f2 = rand(r,1);
            case {true,'true',1} % with separation
                d1 = 1.5/N1;
                E1 = 1-r*d1; % excess space
                if E1 < 0
                    error('Model order is too big, separation condition fails!')
                else
                    f1 = rand(r+1,1); % wrap around distance
                    fs = E1*f1(1:r)/sum(f1);
                    fs = d1*ones(r,1)+fs;
                    f1 = cumsum(fs);
                end
                d2 = 1.5/N2;
                E2 = 1-r*d2; % excess space
                if E2 < 0
                    error('Model order is too big, separation condition fails!')
                else
                    f2 = rand(r+1,1); % wrap around distance
                    fs = E2*f2(1:r)/sum(f2);
                    fs = d2*ones(r,1)+fs;
                    f2 = cumsum(fs);
                end
            otherwise
                error('Separation should either be true or false!')
        end
       switch damp        
            case {false,'false',0} % without damping
                x_star = exp(kron(ones(N2,1),s1')*(1i*2*pi*f1')+kron(s2',ones(N1,1))...
                    *(1i*2*pi*f2'))*c;    
            case {true,'true',1} % with damping
%                 d1 = round(N1/4)+round(N1/4)*rand(r,1);
%                 d1 = -1./d1;
%                 d2 = round(N2/4)+round(N2/4)*rand(r,1);
%                 d2 = -1./d2;
                d1 = 8+8*rand(r,1);
                d1 = -1./d1;
                d2 = 8+8*rand(r,1);
                d2 = -1./d2;
                x_star = exp(kron(ones(N2,1),s1')*(d1'+1i*2*pi*f1')+kron(s2',ones(N1,1))...
                    *(d2'+1i*2*pi*f2'))*c;
            otherwise
                error('Damping should either be true or false!')
        end
        
        f = [f1 f2];
        
%         if m < 1
%             m = round(N1*N2*m);
%         end
  
%     omega(K,:) = 1 ;
    xs(K) = x_star(K);
    
  end

    
        
     