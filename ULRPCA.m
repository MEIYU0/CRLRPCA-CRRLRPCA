function [X,S, dist_X_history, time_counter] = ULRPCA(Y,r,X_star, zeta, eta)
    [~, T] = size(zeta);
    % T = 25;
    % dist_X_history = zeros(1, T-1);
    % time_history = zeros(1, T-1);
    time_counter = 0;

    % Initialization
    tStart = tic;
    X0 = full(Y - Thre(Y, zeta(1)));
   
    [C, U, R] = CUR(X0, r); 
    X = C * U * R;
    % S_old = 1;
    [m, n] = size(U);
    % time_counter = time_counter + toc(tStart);

    % X_old = X_star;
    % main loop
    % fprintf("===============LRPCA logs=============\n");
   for t = 1:(T-1)
       % loop_start = tic;
        S = Thre(Y - X, zeta(t+1));
        Li = Y - S;
        [C, U, R] = CUR(Li, r);
        % U = U - eta(t+1) * (C' * C + eps('double') * eye(n))^(-1) * C' * (X-Li) * R' / (R * R' + eps('double') * eye(n));
        U = U - 0.8 * (C/(C' * C + eps('double') * eye(m)))' * (X-Li) * R' / (R * R' + eps('double') * eye(n));
        X = C * U *R;
        % time_counter = time_counter + toc(tStart);
        % X = L * R';
        % X = C * U * R;
        dist = norm(X - X_star, 'fro')/norm(X_star, 'fro');
        % dist = norm(Y - X - S, "fro")/norm(Y, "fro");
        % dist_X = norm(X - X_old, 'fro')/norm(X_old, 'fro');
        % dist_S = norm(S - S_old, 'fro')/norm(S_old, 'fro');
        % if dist_S < 1e-3 && dist_X <1e-3
        %     break;
        % end
        % X_old = X;
        % S_old = S;
        dist_X_history(t) = dist;
        % time_counter = time_counter + toc(loop_start);
        % time_history(t) = time_counter; % 更新时间历史
        if dist < 1e-4
            break
        end
   end
    % X = L*R';
    % X = C * U * R;
    % dist = norm(X - X_star, 'fro')/norm(X_star, 'fro');
    time_counter = toc(tStart);
    fprintf("k: %d Time: %f Err: %e\n", t , time_counter, dist);
end