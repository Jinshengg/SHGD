function [xs,K,x_star,f,amp]=generate_signal_1D(m,nd,r,seperation,damp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs 
%          m: number of observtions.
%          nd:  length of the 1-dim signal
%           r: model order of the spectrally sparse signal
%  seperation: frequencies seperation
%        damp: dampled signal

%Outputs
%         xs : observed signal
%           K: observation indices
%      x_star: true signal 
%           f: normlized true frequencies
%         amp: amplititude of the frequencies spikes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate signals
    dynamic_range=10;
    c=exp(1i*2*pi*rand(r,1)).*(1+10.^(rand(r,1)*dynamic_range/20));
    xs=zeros(nd,1);
    K=randsample(nd,m);
    
    switch seperation 
          case {false,'false',0} % without separation 
            freq_seed = randperm(nd, r)/nd;
          case {true,'true',1} % with separation
                d1 = 1.5/nd;
                E1 = 1-r*d1; % excess space
             if E1 < 0
                    error('Model order is too big, separation condition fails!')
             else
                    freq_seed = rand(1,r+1); % wrap around distance
                    fs = E1*freq_seed(1:r)/sum(freq_seed);
                    fs = d1*ones(1,r)+fs;
                    freq_seed = cumsum(fs);
             end
            otherwise
                error('Separation should either be true or false!')
    end
    
        switch damp %  damping
            case {false,'false',0} % without damping
               x_star = exp(2*pi*1i * (0:(nd-1))' * freq_seed) * c;  
            case {true,'true',1} % with damping
%                 d1 = round(nd/4)+round(nd/4)*rand(r,1);
                d1 = 16+16*rand(r,1);
                d1 = -1./d1;
                x_star = exp((0:(nd-1))'*(d1'+1i*2*pi*freq_seed))*c;
            otherwise
                error('Damping should either be true or false!')
        end
%     omega(K,:) = 1 ;
    xs(K) = x_star(K);
    f = freq_seed;
    amp = c;
end