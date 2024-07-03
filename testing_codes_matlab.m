clc
clear all
alpha = 0.15;
%% Load Data
r = 5;
% dis_L = zeros(1, 5);
% dis_C = zeros(1, 5);
% dis_U = zeros(1, 5);
% dis_R = zeros(1, 5);
% ti_L = zeros(1, 5);
ti_C = zeros(1, 20);
% ti_U = zeros(1, 5);
% ti_R = zeros(1, 5);
% ti_Ac = zeros(1, 10);
% ti_Sc = zeros(1, 5);
d1 = 4000;	
d2 = 4000;	

% data_path = "./data/n2000_r5_alpha0.1.mat";
% data_path = "./alpha=0.2.mat";
% load(data_path); 
% X_star = sparse(X_star);
% Y = sparse(Y);
for i = 1:20
%     alpha = 0.05 * i;
    % % d1 = 2000 * i;
    % % d2 = 2000 * i;
    % % 初始化 U0_t 和 V0_t
    U0_t = randn(d1, r) / sqrt(d1);
    V0_t = randn(d2, r) / sqrt(d2);

    % 生成随机排列并选择前 alpha * d1 * d2 个元素的索引
    idx = randperm(d1 * d2);
    idx = idx(1:floor(alpha * d1 * d2));

    % 计算 Y0_t
    X_star = U0_t * V0_t';

    % 将 Y0_t 转换为列向量
    Y0_t = X_star(:);

    % 计算 s_range
    s_range = mean(abs(Y0_t));

    % 生成 S0_t
    S0_t = rand(length(idx), 1);
    S0_t = s_range * (2.0 * S0_t - 1.0);

    % 更新 Y0_t 的指定索引
    Y0_t(idx) = Y0_t(idx) + S0_t;
    % 
    % % 将 Y0_t 重塑为 d1 x d2 矩阵
    Y0_t = reshape(Y0_t, [d1, d2]);
    X_star = double(X_star);
    Y = double(Y0_t);
    X_star = sparse(X_star);
    Y = sparse(Y);
    % save('test.mat', 'U0_t', 'V0_t', 'Y0_t');
    % Load Mode
    % model_path = "./归一化1/LRPCA_alpha0.1.mat";
    model_path = "./trained_models/lrpcanet_alpha0.1.mat";
    % model_path = "./LRPCA_alpha0.3.mat";
    load(model_path); % load the trained parameters
    zeta = ths * 1;              % thresholds (adaptive to n,r)
    eta  = step;                                % step sizes
    
    % Running LRPCA
    % [~, T] = size(zeta);
    % time_limit = 20; % 20秒
        % [L1, S1, dist_Acc, time_Acc] = AccAltProj( Y, r, X_star, '' );
        % ti_Ac(i) = time_Acc;
    % end

    % [dist_SCGD, time_SCGD] = ScaledGD(Y, 0.01, X_star, r);
    %
    [X,L, dist_LRPCA, time_LRPCA] = LearnedRPCA(Y,r,X_star,zeta,eta); % x_star is only used for printing
    % load('zetanew.mat');
   % load('zeta36.mat')
   % zeta = zeta * 1;
     % zeta = ths * 1; 
    % [~, ~, dist_C, time_C] = CLRPCA(Y, r, X_star, zeta, eta);
 
    % [~, ~, dist_U, time_U] = ULRPCA(Y, r, X_star, zeta, eta);
    % ti_C(i) = time_U;
end
    % [~, ~, dist_R, time_R] = RLRPCA(Y, r, X_star, zeta, eta);
    % dis_L(i) = dist_LRPCA;
    % dis_C(i) = dist_C;
    % dis_U(i) = dist_U;
    % dis_R(i) = dist_R;
    % ti_Ac(i) = time_Acc;
    % ti_Sc(i) = time_SCGD;
    ti_L(i) = time_LRPCA;
    ti_C(i) = time_C;
    ti_U(i) = time_U;
    % ti_R(i) = time_R;
% end
% x = 2000:2000:10000;
% plot(x, ti_L, '-o', x, ti_C, '-o', x, ti_U, '-o', x, ti_R, '-o');
plot(time_LRPCA, dist_LRPCA, '-o', time_C, dist_C, '-o', time_U, dist_U, '-o', time_Acc, dist_Acc, '-o', time_SCGD, dist_SCGD, '-o');
% plot(1:T-1, dist_LRPCA, '-o', 1:T-1, dist_C, '-o', 1:T-1, dist_U, '-o', 1:T-1, dist_R, '-o');
xlabel('Time(sec)');
% ylabel('Time(sec)')
ylabel('$\Vert L-L_*\Vert / \Vert L_* \Vert$', 'Interpreter', 'latex');
% title('Distortion vs. Time');
% xlim([0, time_limit]); % 设置横坐标范围为0到20秒
ylim([1e-4, 1]);
set(gca, 'YScale', 'log');  % 设置纵坐标为对数刻度
legend('LRPCA', 'Col-random', 'Col-Row-random', 'AccAltProj', 'ScaledGD');
% axes('Position',[0.7, 0.3, 0.2, 0.3]);
% plot(x, ti_C, '-*', x, ti_U, '-s', x, ti_R, '-^');
% axis([7000, 8000, 6, 9])