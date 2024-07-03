function [X,S,dist, time_counter] = CLRPCA(Y,r,X_star, zeta, eta)
    [~, T] = size(zeta);
    time_counter = 0;

    % Initialization
    tStart = tic;
    S = Thre(Y, zeta(1));
    X0 = full(Y - S);
    [C, U, R] = CUR(X0, r);
    L = C;
    R = U * R;
    X = L*R;
    
    % main loop
   for t = 1:(T-1)
       % 开始单次循环计时
        loop_start = tic;      
        S = Thre(Y - X, zeta(t+1));
        Li = Y - S;
        [C, U, R] = CUR(Li, r);
        L = C;
        R = U * R;
        R = R - eta(t+1)*(L/(L'*L + eps('double')*eye(n)))'*(X-Li);
        X = L * R;
        dist = norm(X - X_star, 'fro')/norm(X_star, 'fro');
        dist_X_history(t) = dist;
        if dist < 1e-4
            break
        end
    end
    time_counter = toc(tStart);
    fprintf("k: %d Time: %f Err: %e\n", t , time_counter, dist);
end
