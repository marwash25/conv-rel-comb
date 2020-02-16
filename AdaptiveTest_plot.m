clc
clear all
close all

corr_level = 0; %0.5
d = 250;
k = 100;
n = 500;
tol = 1e-8;
corr_str = strrep(num2str(corr_level), '.', '');
tol_str =  strrep(num2str(tol), '.', ''); 
fload = sprintf('Results/AdaptiveTest_corr%s_d%d_k%d_n%d_tol%s',corr_str,d,k,n,tol_str);
load(fload)

% Choose which methods to run
l1 = 1;
l1_ad = 1;
glinf = 1;
glinf_ad = 1;
rml2 = 0;
rml2_ad = 0;
rl2 = 0;
rl2_ad = 0;
rlinf = 1;
rlinf_ad =  1;

x_axis_vec = sig_noise_vec;
lbd_NotAdaptive = logspace(-6,3,20); 
lbd_Adaptive = logspace(-6,3,20); 
alphas = [0.3];
n_vec = n;
lbd_NotAdaptive = sort(lbd_NotAdaptive,'descend');
lbd_Adaptive = sort(lbd_Adaptive,'descend');
%%
if 0
close all
%for ind_n = 1:numel(n_vec)
%    for ind_alpha = 1:numel(alphas);
iter = 1;
ind_x = 1;
ind_alpha = 1;

figure 
loglog(lbd_NotAdaptive ,EstErr_l1(ind_x,:,iter),'k-*','linewidth',2); hold on
loglog(lbd_Adaptive ,EstErr_l1_ad(ind_x,:,iter,ind_alpha),'r-*','linewidth',2); hold on
loglog(lbd_NotAdaptive ,EstErr_glinf(ind_x,:,iter),'b-*','linewidth',2); hold on
loglog(lbd_Adaptive ,EstErr_glinf_ad(ind_x,:,iter,ind_alpha),'g-*','linewidth',2); 
loglog(lbd_NotAdaptive ,EstErr_rlinf(ind_x,:,iter),'m-*','linewidth',2); hold on
loglog(lbd_Adaptive ,EstErr_rlinf_ad(ind_x,:,iter,ind_alpha),'c-*','linewidth',2); 
legend('l1','l1-ad','glinf','glinf-ad','rlinf','rlinf-ad')
title(['Estimation Error, \alpha = ',num2str(alphas(ind_alpha)),', n = ',num2str(n_vec(ind_x))])

figure 
semilogx(lbd_NotAdaptive ,SuppErr_l1(ind_x,:,iter),'k-*','linewidth',2); hold on
semilogx(lbd_Adaptive ,SuppErr_l1_ad(ind_x,:,iter,ind_alpha),'r-*','linewidth',2); hold on
semilogx(lbd_NotAdaptive ,SuppErr_glinf(ind_x,:,iter),'b-*','linewidth',2); hold on
semilogx(lbd_Adaptive ,SuppErr_glinf_ad(ind_x,:,iter,ind_alpha),'g-*','linewidth',2); 
semilogx(lbd_NotAdaptive ,SuppErr_rlinf(ind_x,:,iter),'m-*','linewidth',2); hold on
semilogx(lbd_Adaptive ,SuppErr_rlinf_ad(ind_x,:,iter,ind_alpha),'c-*','linewidth',2); 
legend('l1','l1-ad','glinf','glinf-ad','rlinf','rlinf-ad')
title(['Support Error, \alpha = ',num2str(alphas(ind_alpha)),', n = ',num2str(n_vec(ind_x))])

%    end
%end
end

%% Plot errors w.r.t best lbd choice

%average over monter carlo runs, min errors over all possible values of lbd

%EstErr_l1_ad = zeros(numel(x_axis_vec),numel(lbd_Adaptive),nits,numel(alphas));

for ind_alpha = 1:numel(alphas)

%% plot estimation error
figure 
count = 1;
if l1
EstErr_l1_min_mean = min(mean(EstErr_l1,3),[],2);
loglog(x_axis_vec ,EstErr_l1_min_mean,'b-*','linewidth',2); hold on
legendStr{count} = '$\ell_1$';
%legendStr{count} = 'Homogeneous Relaxation';
count = count + 1;
end

if l1_ad
EstErr_l1_ad_min_mean = min(mean(EstErr_l1_ad(:,:,:,ind_alpha),3),[],2);
loglog(x_axis_vec ,EstErr_l1_ad_min_mean,'b--*','linewidth',2); 
legendStr{count} = '$\ell_1$-adaptive';
count = count + 1;
end

if glinf
EstErr_glinf_min_mean = min(mean(EstErr_glinf,3),[],2);
loglog(x_axis_vec ,EstErr_glinf_min_mean,'g-*','linewidth',2); hold on
legendStr{count} = '$\ell_1 / \ell_\infty$';
count = count + 1;
end

if glinf_ad
EstErr_glinf_ad_min_mean = min(mean(EstErr_glinf_ad(:,:,:,ind_alpha),3),[],2);
loglog(x_axis_vec ,EstErr_glinf_ad_min_mean,'g--*','linewidth',2); 
legendStr{count} = '$\ell_1 / \ell_\infty$-adaptive';
count = count + 1;
end

if rml2
EstErr_rml2_min_mean = min(mean(EstErr_rml2,3),[],2);
loglog(x_axis_vec ,EstErr_rml2_min_mean,'r-*','linewidth',2); hold on
legendStr{count} = 'rml2';
count = count + 1;
end

if rml2_ad
EstErr_rml2_ad_min_mean = min(mean(EstErr_rml2_ad(:,:,:,ind_alpha),3),[],2);
loglog(x_axis_vec ,EstErr_rml2_ad_min_mean,'r--*','linewidth',2); 
legendStr{count} = 'rml2-ad';
count = count + 1;
end

if rl2
EstErr_rl2_min_mean = min(mean(EstErr_rl2,3),[],2);
loglog(x_axis_vec ,EstErr_rl2_min_mean,'r-*','linewidth',2); hold on
legendStr{count} = 'rl2';
count = count + 1;
end

if rl2_ad
EstErr_rl2_ad_min_mean = min(mean(EstErr_rl2_ad(:,:,:,ind_alpha),3),[],2);
loglog(x_axis_vec ,EstErr_rl2_ad_min_mean,'r--*','linewidth',2); 
legendStr{count} = 'rl2-ad';
count = count + 1;
end

if rlinf
EstErr_rlinf_min_mean = min(mean(EstErr_rlinf,3),[],2);
loglog(x_axis_vec ,EstErr_rlinf_min_mean,'r-*','linewidth',2); hold on
legendStr{count} = '$\Theta^{r}_\infty$';
%legendStr{count} = 'Non-Homogeneous Relaxation';
count = count + 1;
end

if rlinf_ad
EstErr_rlinf_ad_min_mean = min(mean(EstErr_rlinf_ad(:,:,:,ind_alpha),3),[],2);
loglog(x_axis_vec ,EstErr_rlinf_ad_min_mean,'r--*','linewidth',2); 
legendStr{count} = '$\Theta^{r}_\infty$-adaptive';
count = count + 1;
end

%h = legend(legendStr,'Location', 'Best');
%title(['Estimation Error, \alpha = ',num2str(alphas(ind_alpha))])
title(['Estimation Error, \rho = ',num2str(corr_level)])
%set(h,'Interpreter','latex')
set(gca,'fontsize',25)
shg
fname = sprintf('ObjErr_corr%s_d%d_k%d_n%d',corr_str,d,k,n);
%print(gcf,'-dpdf','-r600',fname)

%% plot support error
start_ind = 1;

figure 
if l1
SuppErr_l1_min_mean = min(mean(SuppErr_l1,3),[],2);
semilogx(x_axis_vec(start_ind:end),SuppErr_l1_min_mean(start_ind:end),'b-*','linewidth',2); hold on
end

if l1_ad
SuppErr_l1_ad_min_mean = min(mean(SuppErr_l1_ad(:,:,:,ind_alpha),3),[],2);
semilogx(x_axis_vec(start_ind:end) ,SuppErr_l1_ad_min_mean(start_ind:end),'b--*','linewidth',2); 
end

if glinf
SuppErr_glinf_min_mean = min(mean(SuppErr_glinf,3),[],2);
semilogx(x_axis_vec ,SuppErr_glinf_min_mean,'g-*','linewidth',2); hold on
end

if glinf_ad
SuppErr_glinf_ad_min_mean = min(mean(SuppErr_glinf_ad(:,:,:,ind_alpha),3),[],2);
semilogx(x_axis_vec ,SuppErr_glinf_ad_min_mean,'g--*','linewidth',2); 
end

if rml2
SuppErr_rml2_min_mean = min(mean(SuppErr_rml2,3),[],2);
semilogx(x_axis_vec ,SuppErr_rml2_min_mean,'b-*','linewidth',2); hold on
end

if rml2_ad
SuppErr_rml2_ad_min_mean = min(mean(SuppErr_rml2_ad(:,:,:,ind_alpha),3),[],2);
semilogx(x_axis_vec ,SuppErr_rml2_ad_min_mean,'b--*','linewidth',2); 
end

if rl2
SuppErr_rl2_min_mean = min(mean(SuppErr_rl2,3),[],2);
semilogx(x_axis_vec ,SuppErr_rl2_min_mean,'r-*','linewidth',2); hold on
end

if rl2_ad
SuppErr_rl2_ad_min_mean = min(mean(SuppErr_rl2_ad(:,:,:,ind_alpha),3),[],2);
semilogx(x_axis_vec ,SuppErr_rl2_ad_min_mean,'r--*','linewidth',2); 
end

if rlinf
SuppErr_rlinf_min_mean = min(mean(SuppErr_rlinf,3),[],2);
semilogx(x_axis_vec(start_ind:end) ,SuppErr_rlinf_min_mean(start_ind:end),'r-*','linewidth',2); hold on
end

if rlinf_ad
SuppErr_rlinf_ad_min_mean = min(mean(SuppErr_rlinf_ad(:,:,:,ind_alpha),3),[],2);
semilogx(x_axis_vec(start_ind:end) ,SuppErr_rlinf_ad_min_mean(start_ind:end),'r--*','linewidth',2); 
end

h = legend(legendStr);
%title(['Support Error, \alpha = ',num2str(alphas(ind_alpha))])
title(['Support Recovery Error, \rho = ',num2str(corr_level)])
set(h,'Interpreter','latex','Location', 'Best')
set(gca,'fontsize',25)
axis tight
shg
fname = sprintf('SuppErr_corr%s_d%d_k%d_n%d',corr_str,d,k,n);
%print(gcf,'-dpdf','-r600',fname)

end
