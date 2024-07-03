function [X,S,dist, time_counter] = CLRPCA(Y,r,X_star, zeta, eta)
    [~, T] = size(zeta);
    % T=25;
    % dist_X_history = zeros(1, T-1);
    % time_history = zeros(1, T-1);
    time_counter = 0;

    % zeta_init = max(abs(Y(:)));
    % Initialization
    tStart = tic;
    % Y1 = full(Y);
    S = Thre(Y, zeta(1));
    X0 = full(Y - S);
    [C, U, R] = CUR(X0, r);
    L = C;
    R = U * R;
    % R = R';
    X = L*R;
    % S_old = 1;
    [m, n] = size(L);
    % Li = Y - S;
    % time_counter = time_counter + toc(tStart);
    % X_old = X_star;
    % main loop
    % fprintf("===============LRPCA logs=============\n");
    % loop_start = tic;
   for t = 1:(T-1)
       % 开始单次循环计时
        % loop_start = tic;
        % zeta = 0.7 * zeta_init;
        
        S = Thre(Y - X, zeta(t+1));
        Li = Y - S;
        [C, U, R] = CUR(Li, r);
        L = C;
        R = U * R;
        R = R - 0.8*(L/(L'*L + eps('double')*eye(n)))'*(X-Li);
        X = L * R;
        % L = Li(:, cols);
        % U = L(rows, :);
        % [Uu,Su,Vu] = svd(U);
        % d = diag(Su);
        % % Su = 1./Su(1:r);
        % % U = Vu(:,1:r)*Su*(Uu(:,1:r))';
        % U = Vu(:, 1:r) * diag(1./d(1:r)) * Uu(:, 1:r)';
        % R = Li(rows, :);
        % R = U * R;
        % R = R';
        % R = R' - eta(t+1)*(X+S-Y)'*L/(L'*L + eps('double')*eye(n));

        % time_counter = time_counter + toc(tStart);
        % X = L * R;
        dist = norm(X - X_star, 'fro')/norm(X_star, 'fro');
        % dist = norm(Y - X - S, "fro")/norm(Y, "fro");
        % dist_X = norm(X - X_old, 'fro')/norm(X_old, 'fro');
        % dist_S = norm(S - S_old, 'fro')/norm(S_old, 'fro');
        % if dist_S < 1e-3 && dist_X  < 1e-3
        %     break;
        % end
        % X_old = X;
        % S_old = S;
        dist_X_history(t) = dist;
        % time_counter = time_counter + toc(loop_start);
        % time_history(t) = time_counter; % 更新时间历史
        % time_counter = time_counter + toc(tStart);
        % time_history(t) = time_counter; % 记录累计时间
        if dist < 1e-4
            break
        end
    end
    % X = L*R';
    % time_counter = time_counter + toc(loop_start);
    time_counter = toc(tStart);
    % time_counter = time_counter + toc(tStart);
    % dist = norm(X - X_star, 'fro')/norm(X_star, 'fro');
    fprintf("k: %d Time: %f Err: %e\n", t , time_counter, dist);
end