function [EstErr_l1_ad,SuppErr_l1_ad,EstErr_l1,SuppErr_l1,EstErr_glinf_ad,SuppErr_glinf_ad,EstErr_glinf,SuppErr_glinf,EstErr_rlinf_ad,SuppErr_rlinf_ad,EstErr_rlinf,SuppErr_rlinf]= AdaptiveTestfct(d,n,k,sigma,nits,correlated,corr_level,tol)

%% Parameters

constant = 1;
%tol = 1e-8;

alphas = [0.3];%0.1:0.1:0.9; 

n_vec = n; 
sig_noise_vec = sigma; 
lbd_NotAdaptive = logspace(-6,3,20); 
lbd_Adaptive = logspace(-6,3,20);

lbd_NotAdaptive = sort(lbd_NotAdaptive,'descend');
lbd_NotAdaptive = lbd_NotAdaptive(5:11);

lbd_Adaptive = sort(lbd_Adaptive,'descend');
lbd_Adaptive = lbd_Adaptive(5:11);
% FISTA parameters
param.AR = 1; %flag for Adaptive restart
param.maxit = 2e3;
param.tolx = 1e-9;
param.verbose = 0;
param_ad = param;

% Choose which methods to run
l1 = 1;
l1_ad = 1;
glinf = 1;
glinf_ad = 1;
rlinf = 1;
rlinf_ad = 1;

% Choose if we want to plot against measurements n or noise level sig
x_axis_vec = sig_noise_vec; 
xaxis_noise = 1;

EstErr_l1_ad = zeros(numel(x_axis_vec),numel(lbd_Adaptive),nits,numel(alphas));
SuppErr_l1_ad = zeros(numel(x_axis_vec),numel(lbd_Adaptive),nits,numel(alphas));
EstErr_l1 = zeros(numel(x_axis_vec),numel(lbd_NotAdaptive),nits);
SuppErr_l1 = zeros(numel(x_axis_vec),numel(lbd_NotAdaptive),nits);

EstErr_glinf_ad = zeros(numel(x_axis_vec),numel(lbd_Adaptive),nits,numel(alphas));
SuppErr_glinf_ad = zeros(numel(x_axis_vec),numel(lbd_Adaptive),nits,numel(alphas));
EstErr_glinf = zeros(numel(x_axis_vec),numel(lbd_NotAdaptive),nits);
SuppErr_glinf = zeros(numel(x_axis_vec),numel(lbd_NotAdaptive),nits);

EstErr_rlinf_ad = zeros(numel(x_axis_vec),numel(lbd_Adaptive),nits,numel(alphas));
SuppErr_rlinf_ad = zeros(numel(x_axis_vec),numel(lbd_Adaptive),nits,numel(alphas));
EstErr_rlinf = zeros(numel(x_axis_vec),numel(lbd_NotAdaptive),nits);
SuppErr_rlinf = zeros(numel(x_axis_vec),numel(lbd_NotAdaptive),nits);
%% Define interval groups

M = 2*d-1; % number of groups

triu_ones = triu(ones(d));
tril_ones = tril(ones(d));
B = sparse([triu_ones,tril_ones(:,2:end)]);

% define corresponding GL norms
GLinf_norm = @(x) sum(max(B.*repmat(abs(x),[1,M])));

% define corresponding SPAMS graph
graph.groups = sparse([triu_ones(1:end-1,:),zeros(d-1);zeros(1,M);zeros(d-1),tril_ones(1:end-1,:)]);
graph.groups_var = sparse([eye(d),[zeros(1,d-1);eye(d-1)]]);
graph.eta_g = ones(1,M);

% parameters for range 
I = zeros(d*(d+1)/2,d); % rows of M are all possible intervals I.
b_mod = zeros(d*(d+1)/2,1); % b = F(I) for all I.
b = zeros(d*(d+1)/2,1); % b = F(I) for all I.
r = 0;
for i = 1:d
    for j = i:d
       r = r+1;
       I(r,i:j) = 1;
       b_mod(r) = j-i + d;
       b(r) = j-i + 1;
    end
end

%% Do nits monte Carlo runs
for iter = 1:nits
    %% Generate Gaussian/Constant signal with interval support (simply take first k coeff)

    x = zeros(d,1);
    J = [ones(k,1);zeros(d-k,1)];
    if constant
        xJ = rand*ones(k,1);
        %xJ = ones(k,1);
    else
        xJ = randn(k,1);
        xJ = xJ/norm(xJ,inf); % normalize to have ||x||inf =1;
    end
    x(1:k) = xJ; 

    EstErr = @(xhat) norm(xhat - x,2);
    SuppErr = @(Jhat) norm(J - Jhat,1);


    for ind_n = 1:numel(n_vec)
        n  = n_vec(ind_n);

        if correlated
            c = corr_level; % ~correlation b/w cols
            L = triu(c + .01*randn(d),1);
            Q = L + L' + eye(d);
            min_eig  = min(eig(Q));
            if min_eig <0 %make sure Q is psd
                Q = (Q - (min_eig-1e-7)*eye(d))/(1 - (min_eig-1e-7));
            end
            R = chol(Q);
            A = 1/sqrt(n) * randn(n, d)*R;
        else
            A = 1/sqrt(n) * randn(n, d);
        end

        Ax = A*x;
        
        for ind_sig = 1:numel(sig_noise_vec)
            sig_noise = sig_noise_vec(ind_sig);
            
            if xaxis_noise
                ind_x = ind_sig;
            else
                ind_x = ind_n;
            end
            
            noise = sig_noise*randn(n,1);
            y = Ax+ noise;


            % FISTA parameters
            gradf = @(x)A'*(A*x-y); 
            fx = @(x) 0.5*norm(A*x-y,2)^2;
            param.Lips = normest(A)^2;
            %% Run adaptive methods 
            for ind_alpha = 1:numel(alphas)
                alpha = alphas(ind_alpha);

                % Use LS for Adaptive weights
                time_ls = tic;
                w = abs(A\y + 1e-7).^(alpha - 1); %add 1e-7 to avoid dividing by zero
                time_ls = toc(time_ls);

                % Define corresponding fx and A (if we can't add the weights explicitly to the prox) 
                % FISTA parameters
                A_ad = A.*repmat((w.^(-1))',[n,1]);
                gradf_ad = @(x)A_ad'*(A_ad*x-y); 
                fx_ad = @(x) 0.5*norm(A_ad*x-y,2)^2;
                param_ad.Lips = normest(A_ad)^2;

                x_l1_ad_prev = zeros(d,1);
                x_glinf_ad_prev = zeros(d,1);

                for ind_lbd = 1:numel(lbd_Adaptive)
                    lbd = lbd_Adaptive(ind_lbd); %start by largest alpha

                    %% Adaptive Lasso
                    if l1_ad
              
                    g_l1 = @(x) lbd*norm(w.*x,1);
                    proxg_l1 = @(x) SofTh(x, lbd/param.Lips,w);

                    param.x0 = x_l1_ad_prev;

                    [x_l1_ad, info.l1_ad] = FISTA(fx, gradf, g_l1, proxg_l1, param);

                    x_l1_ad_prev = x_l1_ad;

                    %=======================================
                    EstErr_l1_ad(ind_x,ind_lbd,iter,ind_alpha)= EstErr(x_l1_ad);
                    SuppErr_l1_ad(ind_x,ind_lbd,iter,ind_alpha)= SuppErr(abs(x_l1_ad)>tol);
                    end
                     %% Adaptive L1/Linf Overlapping Group Lasso      
                     if glinf_ad 
                            % GL prox parameters
                            param_GL_prox.lambda = lbd/param_ad.Lips; % regularization parameter
                            param_GL_prox.regul='graph';

                            g_glinf = @(x) lbd*GLinf_norm(x);

                            proxg_glinf = @(x) mexProximalGraph(x,graph, param_GL_prox);

                            param_ad.x0 = x_glinf_ad_prev;

                            [x_glinf_ad, info.glinf_ad] = FISTA(fx_ad, gradf_ad, g_glinf, proxg_glinf, param_ad);

                            x_glinf_ad_prev = x_glinf_ad;
                            x_glinf_ad = x_glinf_ad./w;

                            %=======================================
                            EstErr_glinf_ad(ind_x,ind_lbd,iter,ind_alpha)= EstErr(x_glinf_ad);
                            SuppErr_glinf_ad(ind_x,ind_lbd,iter,ind_alpha)= SuppErr(abs(x_glinf_ad)>tol);
                     end   
                    
                    %% Adaptive Linf-Convex Relaxation of (unmodified) range function
                    if rlinf_ad
                    %============ cvx  ================
                    cvx_begin quiet
                        cvx_precision best
                        variables x_rlinf_cvx_ad(d) alpha(r);

                        minimize (0.5*square_pos(norm(y - A*x_rlinf_cvx_ad, 2))+lbd*b'*alpha)
                        subject to
                            I'*alpha >= abs(w.*x_rlinf_cvx_ad)
                            alpha>= 0
                            sum(alpha) <= 1
                    cvx_end

                    EstErr_rlinf_ad(ind_x,ind_lbd,iter,ind_alpha)= EstErr(x_rlinf_cvx_ad);
                    SuppErr_rlinf_ad(ind_x,ind_lbd,iter,ind_alpha)= SuppErr(abs(x_rlinf_cvx_ad)>tol);
                    %=======================================
                    end
                end
            end
            %% Run non-adaptive methods 
            x_l1_prev = zeros(d,1);
            x_glinf_prev = zeros(d,1);

            for ind_lbd = 1:numel(lbd_NotAdaptive)
                lbd = lbd_NotAdaptive(ind_lbd);

                %% Lasso
                if l1
                    
                g_l1 = @(x) lbd*norm(x,1);
                proxg_l1 = @(x) SofTh(x, lbd/param.Lips,ones(d,1));

                param.x0 = x_l1_prev;

                time_l1 = tic;
                [x_l1, info.l1] = FISTA(fx, gradf, g_l1, proxg_l1, param);
                time_l1 = toc(time_l1);

                x_l1_prev = x_l1;

                %=======================================
                EstErr_l1(ind_x,ind_lbd,iter)= EstErr(x_l1);
                SuppErr_l1(ind_x,ind_lbd,iter)= SuppErr(abs(x_l1)>tol);
                end
                %% L1/Linf Overlapping Group Lasso      
                if glinf

                % GL prox parameters
                param_GL_prox.lambda = lbd/param.Lips; % regularization parameter
                param_GL_prox.regul='graph';

                g_glinf = @(x) lbd*GLinf_norm(x);
                proxg_glinf = @(x) mexProximalGraph(x,graph, param_GL_prox);

                param.x0 = x_glinf_prev;

                [x_glinf, info.glinf] = FISTA(fx, gradf, g_glinf, proxg_glinf, param);

                x_glinf_prev = x_glinf;

                %=======================================
                EstErr_glinf(ind_x,ind_lbd,iter,ind_alpha)= EstErr(x_glinf);
                SuppErr_glinf(ind_x,ind_lbd,iter,ind_alpha)= SuppErr(abs(x_glinf)>tol);
                end
          
              
                %%  Linf-Convex Relaxation of (unmodified) range function
                if rlinf
                %=============== cvx ================
                cvx_begin quiet
                    cvx_precision best
                    variables x_rlinf_cvx(d) z(r);

                    minimize (0.5*square_pos(norm(y - A*x_rlinf_cvx, 2))+lbd*b'*z)
                    subject to
                        I'*z >= abs(x_rlinf_cvx)
                        z>= 0
                        sum(z) <= 1
                cvx_end

                EstErr_rlinf(ind_x,ind_lbd,iter,ind_alpha)= EstErr(x_rlinf_cvx);
                SuppErr_rlinf(ind_x,ind_lbd,iter,ind_alpha)= SuppErr(abs(x_rlinf_cvx)>tol);
                %=======================================
                end
            end
        end
    end
end


end