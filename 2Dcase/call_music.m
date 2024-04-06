function [freq,amp,signal] = call_music(Y,M,N,isPlot)
% Editor: Yinchuan Li Date:2019.05.31

tol = 0.002;

% T_cp = simParam.T_cp;
% T = simParam.T;
% len = T_cp/T;
% c = 3e8;
% simParam.PRT = 1500e-6;
% simParam.delta_f = 1/10e-6;
% simParam.fc = 2e9;
% simParam.wavelength = c/simParam.fc;

c = 1;
simParam.PRT = 1;
simParam.delta_f = 1;
simParam.fc = 1;
simParam.wavelength = 1;

M2 = floor(M/2);
N2 = floor(N/2);

A=[];
for m2 = 1:M-M2+1
    for n2 = 1:N-N2+1
        Y_sub=Y(m2:1:m2+M2-1,n2:1:n2+N2-1);
        A=[A, vec(Y_sub)];
    end
end

[U,D,V]=svd(A);
d=diag(D);
E = sum(d>max(d)*tol);
if E<M2*N2
    Un = U(:,E+1:end);
else
    Un = U(:,E);
end
m=[0:1:M2-1]';
n=[0:1:N2-1]';

K=5e2;
f1=[0:1/K:1-1/K];
f2=[0:1/K:1-1/K];

for k1=1:K
    for k2 = 1:K
        F=exp(1i*2*pi*f1(k1)*m) * exp(1i*2*pi*f2(k2)*n).';
        f=vec(F);
        P(k1,k2)=abs((f'*Un)*(f'*Un)');
    end
end

[idx1,idx2]=find(abs(P)<10);
p1=[];
p2=[];
for k=1:length(idx1)
    p = [idx1(k)-1,idx1(k)-1,idx1(k)-1,...
        idx1(k),idx1(k),...
        idx1(k)+1,idx1(k)+1,idx1(k)+1;
        idx2(k)-1,idx2(k),idx2(k)+1,...
        idx2(k)-1,idx2(k)+1,...
        idx2(k)-1,idx2(k),idx2(k)+1];
    
    for i=1:length(p)
        if p(1,i)<1
            p(1,i) = K-p(1,i);
        end
        if p(1,i)>K
            p(1,i) = 1;
        end
        if p(2,i)<1
            p(2,i) = K-p(2,i);
        end
        if p(2,i)>K
            p(2,i) = 1;
        end
    end
    
    for i=1:length(p)
        P_around(i) = P(p(1,i),p(2,i));
    end
    
    if abs(P(idx1(k),idx2(k)))<min(abs(P_around))
        p1 = [p1,idx1(k)];
        p2 = [p2,idx2(k)];
    end
end

for k=1:length(p1)
    freq(1,k)=f1(p1(k));
    freq(2,k)=f2(p2(k));
end

m=[0:1:M-1]';
n=[0:1:N-1]';
C=[];
for k=1:length(p1)
    F=exp(1i*2*pi*freq(1,k)*m) * exp(1i*2*pi*freq(2,k)*n).';
    C=[C,vec(F)];
end

y=vec(Y);
if ~isempty(p1)
    amp = pinv(C)*y;
    signal = reshape(C*amp,M,N);
else
    amp = [];
    freq = [];
    signal = zeros(M,N);
end

switch isPlot
    case 1
        PRT = simParam.PRT;
        wavelength = simParam.wavelength;
        delta_f = simParam.delta_f;

        f1=[-0.5+1/K:1/K:0.5];
        velocity = f1 / PRT * wavelength;
        range = (1 - f2)/delta_f*c;
        spec = fftshift(1./P,1);
        spec = spec/max(max(spec));
        figure,imagesc(range,velocity,db(spec)/2,[-40,0]),axis xy,hold on;
        xlabel('Range [m]','fontsize',12),ylabel('Velocity [m/s]','fontsize',12);
    case 2
        PRT = simParam.PRT;
        wavelength = simParam.wavelength;
        delta_f = simParam.delta_f;

        f1=[-0.5+1/K:1/K:0.5];
        velocity = f1 / PRT * wavelength;
        range = (1 - f2)/delta_f*c;
        spec = fftshift(1./P,1);
        spec = spec/max(max(spec));
        figure,imagesc(range,velocity,db(spec)/2,[-40,0]),axis xy,hold on;
        xlabel('Range [m]','fontsize',12),ylabel('Velocity [m/s]','fontsize',12);
    otherwise
        %         error('No such plotting strategy.');
end