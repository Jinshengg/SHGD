clear; close all; clc;
%% parameters
n1 =1024; n2=1024; nd=n1+n2-1;
r=150;
success=1e-3;
max_iter =600;
max_iter_FIHT = 600;
tol = 0;
seperation = false;
damp = false;
opt = 0;
stepsize = 0.75;
m=round(115*log(nd));    
p= m/(nd);
Monte= 20;
errgrid_num=50;

% [xs,K,x_star]=generate_signal_1D(m,nd,r,seperation,damp);


time_PGD = zeros(max_iter,Monte);
time_SHGD = zeros(max_iter,Monte);
time_FIHT= zeros(max_iter_FIHT,Monte);

iter_PGD = zeros(max_iter,Monte);
iter_SHGD = zeros(max_iter,Monte);
iter_FIHT= zeros(max_iter_FIHT,Monte);

error_PGDt = zeros(max_iter,Monte);
error_SHGDt = zeros(max_iter,Monte);
error_FIHTt= zeros(max_iter_FIHT,Monte);
% Different recovery accuracy levels: minimum: 5e-7; max: 0.5
rec_err_grid = linspace(log(5e-7),log(0.5),errgrid_num);
rec_err_grid = exp(rec_err_grid);

errSHGD_versus_time = zeros(errgrid_num,Monte);
errPGD_versus_time = zeros(errgrid_num,Monte);
errFIHT_versus_time = zeros(errgrid_num,Monte);



for i1=1:1:Monte
%% generate 1D signal
[xs,K,x_star,~,~] = generate_signal_1D(m,nd,r,seperation,damp);

%% SHGD
tic
[x ,timer_SHGD,error_t] = SHGD(xs,K,nd,r,p,tol,max_iter,opt,stepsize,x_star);
toc
% iter_SHGD(:,i1) = iter;
time_SHGD(:,i1) = timer_SHGD;
error_SHGDt(:,i1) = error_t;
%% PGD
tic
[x ,timer_PGD,error_t] = PGD(xs,K,n1,n2,nd,r,p,tol,max_iter,opt,stepsize,x_star);
toc
% iter_PGD(:,i1) = max(iter);
time_PGD (:,i1) = timer_PGD;
error_PGDt(:,i1) = error_t;
%% FIHT
tic
[x ,timer_FIHT,error_t] =FIHT(xs,K,nd,r,p,tol,max_iter_FIHT,x_star);
toc
% iter_FIHT(:,i1) = iter;
time_FIHT(:,i1) = timer_FIHT;
error_FIHTt(:,i1) = error_t;

% find the indices which meet minimum required accuracy 
ind_SHGD = min(find(error_SHGDt(:,i1)<rec_err_grid(1)));
ind_PGD = min(find(error_PGDt(:,i1)<rec_err_grid(1)));
ind_FIHT = min(find(error_FIHTt(:,i1)<rec_err_grid(1)));
% if any algorithm doesn't attach the minimum accuracy, discard this trial
if (length(ind_SHGD)==0||length(ind_PGD)==0||length(ind_FIHT)==0)
0
continue;
end
% construct different recovery accuracies versus time  
for i2= 1:length(rec_err_grid) 
% the indice to attach the required accuracy
ind_SHGD = min(find(error_SHGDt(:,i1)<rec_err_grid(i2)));
ind_PGD = min(find(error_PGDt(:,i1)<rec_err_grid(i2)));
ind_FIHT = min(find(error_FIHTt(:,i1)<rec_err_grid(i2)));
% corresponding recovery time
errSHGD_versus_time(i2,i1) =  time_SHGD(ind_SHGD,i1);
errPGD_versus_time(i2,i1) =  time_PGD(ind_PGD,i1);
errFIHT_versus_time(i2,i1) =  time_FIHT(ind_FIHT,i1);
end


end

%%  Recovery accuracies versus average time

errSHGD_versus_time_avg = zeros(length(rec_err_grid),1);
errPGD_versus_time_avg = zeros(length(rec_err_grid),1);
errFIHT_versus_time_avg =zeros(length(rec_err_grid),1);
for i2= 1:length(rec_err_grid) 
errSHGD_versus_time_avg(i2) = mean(errSHGD_versus_time(i2,find(errSHGD_versus_time(i2,:)~=0)));
errPGD_versus_time_avg(i2) = mean(errPGD_versus_time(i2,find(errPGD_versus_time(i2,:)~=0)));
errFIHT_versus_time_avg(i2) =  mean(errFIHT_versus_time(i2,find(errFIHT_versus_time(i2,:)~=0)));
end

%% plot
clrs = {[.5,0,.5], [1,.5,0], [1,0,0], [0,.5,0], [0,0,1]};
mks = {'o', 'x', 'p', 's', 'd'};
figure('Position', [0,0,800,600], 'DefaultAxesFontSize', 20);

lgd = {'SHGD','PGD','FIHT'};
semilogy(errSHGD_versus_time_avg(1:2:end),rec_err_grid(1:2:end),'Color', 'b', 'Marker', mks{1}, 'MarkerSize', 8,'LineWidth',1);
hold on;grid on;
semilogy(errPGD_versus_time_avg(1:2:end),rec_err_grid(1:2:end),'Color', 'r', 'Marker', 'x', 'MarkerSize', 8,'LineWidth',1);
hold on;grid on;
semilogy(errFIHT_versus_time_avg(1:2:end),rec_err_grid(1:2:end),'Color', clrs{2}, 'Marker', 'd', 'MarkerSize', 8,'LineWidth',1);
xlabel('Time (secs)');
ylabel('Relative error');
legend(lgd, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 28);
fig_name = 'Relative_Error_vs_Runtime_1D_withoutsep';

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 8 6]);
set(gca,'FontName','times new roman','FontSize',22,'Layer','top');
myfig = gcf;
axis([0 10  10^(-5) 1]);
data_name=strcat('timecompar_1Dwithoutsep',datestr(now,30),'.mat');
%  save(data_name,'errSHGD_versus_time_avg','errPGD_versus_time_avg','errFIHT_versus_time_avg');
% print( myfig, fig_name, '-depsc' );
