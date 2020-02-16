function AdaptiveTest(corr_level,n,tol,nruns)

% set parameters
d = 5;
k = 2;
%n = 500;

nlbds = 7;
nalphas = 1;

sig_noise_vec = logspace(-3,0,7);%[1e-5,1e-4,1e-3,1e-2,1e-1,1];
%sig_noise_vec  = sig_noise_vec(4:end);
correlated = (corr_level ~= 0);

EstErr_l1_ad = zeros(numel(sig_noise_vec),nlbds,nruns,nalphas);
SuppErr_l1_ad = zeros(numel(sig_noise_vec),nlbds,nruns,nalphas);
EstErr_l1 = zeros(numel(sig_noise_vec),nlbds,nruns,nalphas);
SuppErr_l1 = zeros(numel(sig_noise_vec),nlbds,nruns,nalphas);

EstErr_glinf_ad = zeros(numel(sig_noise_vec),nlbds,nruns,nalphas);
SuppErr_glinf_ad = zeros(numel(sig_noise_vec),nlbds,nruns,nalphas);
EstErr_glinf = zeros(numel(sig_noise_vec),nlbds,nruns);
SuppErr_glinf = zeros(numel(sig_noise_vec),nlbds,nruns);

EstErr_rlinf_ad = zeros(numel(sig_noise_vec),nlbds,nruns,nalphas);
SuppErr_rlinf_ad = zeros(numel(sig_noise_vec),nlbds,nruns,nalphas);
EstErr_rlinf = zeros(numel(sig_noise_vec),nlbds,nruns,nalphas);
SuppErr_rlinf = zeros(numel(sig_noise_vec),nlbds,nruns,nalphas);

parfor ind = 1:numel(sig_noise_vec) 
    !hostname 

    %%
    % set seed
    rng(1)
    sigma = sig_noise_vec(ind);
    fprintf('Running Adaptive Test with noise = %2.5f and correlation %1.1f \n\n',sigma,corr_level);
    
    [EstErr_l1_ad(ind,:,:,:),SuppErr_l1_ad(ind,:,:,:),EstErr_l1(ind,:,:,:),SuppErr_l1(ind,:,:,:),EstErr_glinf_ad(ind,:,:,:),SuppErr_glinf_ad(ind,:,:,:),EstErr_glinf(ind,:,:,:),SuppErr_glinf(ind,:,:,:),EstErr_rlinf_ad(ind,:,:,:),SuppErr_rlinf_ad(ind,:,:,:),EstErr_rlinf(ind,:,:,:),SuppErr_rlinf(ind,:,:,:)]= AdaptiveTestfct(d,n,k,sigma,nruns,correlated,corr_level,tol)
end

% remove dot from corr_level
corr_str = strrep(num2str(corr_level), '.', '');
tol_str = strrep(num2str(tol), '.', '');
fsave = sprintf('Results/AdaptiveTest_corr%s_d%d_k%d_n%d_tol%s',corr_str,d,k,n,tol_str);
save(fsave)
end