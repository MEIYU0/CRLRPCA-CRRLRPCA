clc
clear all
alpha = 0.1;
r = 5;
d1 = 4000;	
d2 = 4000;	

U0_t = randn(d1, r) / sqrt(d1);
V0_t = randn(d2, r) / sqrt(d2);

idx = randperm(d1 * d2);
idx = idx(1:floor(alpha * d1 * d2));
X_star = U0_t * V0_t';
Y0_t = X_star(:);
s_range = mean(abs(Y0_t));
S0_t = rand(length(idx), 1);
S0_t = s_range * (2.0 * S0_t - 1.0);
Y0_t(idx) = Y0_t(idx) + S0_t;
Y0_t = reshape(Y0_t, [d1, d2]);
X_star = double(X_star);
Y = double(Y0_t);
X_star = sparse(X_star);
Y = sparse(Y);

% Load Mode
model_path = "./trained_models/alpha0.1.mat";
load(model_path); % load the trained parameters
zeta = ths * 1;              % thresholds (adaptive to n,r)
eta  = step;                                % step sizes

[~, ~, dist_C, time_C] = CLRPCA(Y, r, X_star, zeta, eta);

[~, ~, dist_U, time_U] = ULRPCA(Y, r, X_star, zeta, eta);
