clear; close all; clc;
k = 5;
N1=2.^k-2;N2=2.^k-2;

r=4;
sigma = 0.1; % noise level
success=1e-3;
max_iter = 100;
tol = 1e-7;
opt = 0;
stepsize = 0.5;
% Monte=10;
seperation=false;
m = 100;    
p= m/(N1*N2);

%% undamped case
damp=false;
% generate signal
[xs,K,x_star,f,amp_star] = generate_signal_2D(m,[N1 N2],r,seperation,damp);
noise =   randn(m,1)+1i*randn(m,1);
xs_n = xs; 
xs_n (K) = xs(K) + norm(xs) * sigma * noise/norm(noise);% noise contamination
X_star = reshape(x_star,[N1 N2]);
Xs_n = reshape(xs_n,[N1 N2]);
% SHGD_2D
[X ,~,~,~] = SHGD_2D(Xs_n(1:N1,1:N2),K,N1,N2,r,p,tol,max_iter,opt,stepsize,X_star);
error_SHGD = norm((X(1:N1,1:N2) - X_star),'fro')/norm(X_star,'fro');
%% obtain location and amplitude of frequencies 

[freq,amp,signal]=call_music(X(1:N1,1:N2),N1,N2,3);
freq = freq';
[~,ind] = maxk(abs(amp),r);
stem3(f(:,1),f(:,2),abs(amp_star),'--+','Markersize',11.5,'color','k','LineWidth',3.5);
hold on;
stem3(freq(ind,1),freq(ind,2),abs(amp(ind)),':o','Markersize',11.5,'color','r','LineWidth',3.5);
xlabel('$f_1$','interpreter','latex');
ylabel('$f_2$','interpreter','latex');
zlabel('magnitude');
axis([0,1,0,1,0,4]);
lgd = legend('Original','SHGD');
set(lgd,'position',[0.45 0.78 0.25 0.14]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gca,'FontName','times new roman','FontSize',28,'Layer','top','LineWidth',4);
fig_name = 'figure/stable_rec_2D_withoutdamp';
myfig = gcf;
% print( myfig, fig_name, '-depsc' );
 %% damped case 
%  m = 100;    
% p= m/(N1*N2);
% damp = true;
% [xs,K,x_star,f,amp_star] = generate_signal_2D(m,[N1 N2],r,seperation,damp);
% noise =   randn(m,1)+1i*randn(m,1);
% xs_n = xs; 
% xs_n (K) = xs(K) + norm(xs) * sigma * noise/norm(noise);% noise contamination
% X_star = reshape(x_star,[N1 N2]);
% 
% Xs_n = reshape(xs_n,[N1 N2]);
% [X ,~,~,~] = SHGD_2D(Xs_n(1:N1,1:N2),K,N1,N2,r,p,tol,max_iter,opt,stepsize,X_star);
% error_SHGD = norm((X(1:N1,1:N2) - X_star),'fro')/norm(X_star,'fro');
% 
% 
% [f1,f2,amp]=ESPRIT2D_damped(X,N1+1,N2+1,r);
% figure;
% stem3(f(:,1),f(:,2),abs(amp_star),'--+','Markersize',11.5,'color','k','LineWidth',3.5);
% hold on;
% stem3(f1,f2,abs(amp),':o','Markersize',11.5,'color','r','LineWidth',3.5);
% xlabel('$f_1$','interpreter','latex');
% ylabel('$f_2$','interpreter','latex');
% zlabel('magnitude');
% lgd = legend('Original','SHGD');
% axis([0,1,0,1,0,4]);
% set(lgd,'position',[0.45 0.78 0.25 0.14]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperPosition', [0 0 8 6]);
% set(gca,'FontName','times new roman','FontSize',28,'Layer','top','LineWidth',4);
% fig_name = 'figure/stable_rec_2D_withdamp';

% myfig = gcf;
% print( myfig, fig_name, '-depsc' );