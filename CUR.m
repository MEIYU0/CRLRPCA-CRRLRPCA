function [C, U, R, rows, cols] = CUR(A, r)
con = 1;
[m, n] = size(A);
siz_col = ceil(con * r*log(n));
siz_row = ceil(con * r*log(m));

% A_sq = A .* A;
% sum_A_sq = sum(A_sq(:));
% sum_A_sq_0 = sum(A_sq, 1) ;    % 按列相加
% sum_A_sq_1 = sum(A_sq, 2);  % 按行相加
% % % %
% P_x_c = sum_A_sq_0 / sum_A_sq;  % 抽取每一列的概率
% P_x_r = sum_A_sq_1 / sum_A_sq;  % 抽取每一行的概率
% % % % 
% [~, rows] = sort(sum_A_sq_1, 'descend');
% rows = rows(1:siz_col);
% [~, cols] = sort(sum_A_sq_0, 'descend');
% cols = cols(1:siz_row);
rows = randsample(m, siz_row);
cols = randsample(n, siz_col);
% rows = randsample(m, siz_row, false, P_x_r);
% cols = randsample(n, siz_col, false, P_x_c);
% 
% rows = ssample(P_x_r, siz_row);
% cols = ssample(P_x_c, siz_col);
% rows = randperm(m, siz_row);
% cols = randperm(n, siz_col);
% % rows = unique(rows);
% cols = unique(cols);

C = A(:, cols);
R = A(rows, :);
W = C(rows, :);
[Uu,Su,Vu] = svds(W);
d = diag(Su);
% Su = 1./Su(1:r);
% U = Vu(:,1:r)*Su*(Uu(:,1:r))';
U = Vu(:, 1:r) * diag(1./d(1:r)) * Uu(:, 1:r)';
% C = C * Vu(:, 1:r);
% U = diag(1./d(1:r));
% R = Uu(:, 1:r)' * R;

% U = pinv(W);
% U = 1 ./ W;
end

