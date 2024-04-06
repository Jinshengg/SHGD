clear; close all; clc;
%% parameters 
n1 =128; n2=128; nd=n1+n2-1; % length of signal
r=8; % model order
success=1e-3;  % success rate
max_iter =600; % maximum iterations
tol = 1e-9; % convergence tolerance
seperation = false; % frequencies separation
damp = false;% damped exponential signal
opt = 0; % backtracking or not
stepsize = 0.5; % stepsize 
m = 80;    % number observations 
p = m/(nd); % sample ratio
%% generate 1D signal
[xs,K,x_star]=generate_signal_1D(m,nd,r,seperation,damp);
%% SHGD for recovery 
[x ,timer_SHGD,error_t] = SHGD(xs,K,nd,r,p,tol,max_iter,opt,stepsize,x_star); 

%% plot error versus time
itend = length(find(error_t~=0));
clrs = {[.5,0,.5], [1,.5,0], [1,0,0], [0,.5,0], [0,0,1]};
mks = {'o', 'x', 'p', 's', 'd'};
figure('Position', [0,0,800,600], 'DefaultAxesFontSize', 20);
lgd = {'SHGD'};

semilogy(timer_SHGD(1:4:itend),error_t(1:4:itend),'Color', 'b', 'Marker', mks{1}, 'MarkerSize', 8,'LineWidth',1);

legend(lgd, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 28);
grid on;
xlabel('Time (secs)');
ylabel('Relative error');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gca,'FontName','times new roman','FontSize',22,'Layer','top');
myfig = gcf;