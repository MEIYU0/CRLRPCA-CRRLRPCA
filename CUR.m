function [C, U, R, rows, cols] = CUR(A, r)
con = 1;
[m, n] = size(A);
siz_col = ceil(con * r*log(n));
siz_row = ceil(con * r*log(m));

rows = randsample(m, siz_row);
cols = randsample(n, siz_col);

C = A(:, cols);
R = A(rows, :);
W = C(rows, :);
[Uu,Su,Vu] = svds(W);
d = diag(Su);
U = Vu(:, 1:r) * diag(1./d(1:r)) * Uu(:, 1:r)';

end
