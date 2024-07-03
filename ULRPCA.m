function [X,S, dist_X_history, time_counter] = ULRPCA(Y,r,X_star, zeta, eta)
    [~, T] = size(zeta);
    time_counter = 0;

    % Initialization
    tStart = tic;
    X0 = full(Y - Thre(Y, zeta(1)));
   
    [C, U, R] = CUR(X0, r); 
    X = C * U * R;
    [m, n] = size(U);
    
    % main loop
   for t = 1:(T-1)
       % loop_start = tic;
        S = Thre(Y - X, zeta(t+1));
        Li = Y - S;
        [C, U, R] = CUR(Li, r);
        U = U - eta(t+1) * (C/(C' * C + eps('double') * eye(m)))' * (X-Li) * R' / (R * R' + eps('double') * eye(n));
        X = C * U *R;
        dist = norm(X - X_star, 'fro')/norm(X_star, 'fro');
        if dist < 1e-4
            break
        end
   end
    time_counter = toc(tStart);
    fprintf("k: %d Time: %f Err: %e\n", t , time_counter, dist);
end
