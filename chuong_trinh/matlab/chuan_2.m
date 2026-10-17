function [n, k] = chuan_2(A, tol, max_lap)
%CHUAN_2  Chuẩn 2 (chuẩn phổ) = sqrt(lambda_max(A'*A)).
%   lambda_max tìm bằng PHƯƠNG PHÁP LŨY THỪA:
%       y = M*x,  x = y/||y||,  lambda = x'*M*x
%   Dừng khi |lambda_k - lambda_(k-1)| <= tol*lambda_k.
%   [n, k] = chuan_2(A) trả thêm số vòng lặp k.
if nargin < 2, tol = 1e-12; end
if nargin < 3, max_lap = 10000; end
M = A' * A;                          % đối xứng, xác định không âm
m = size(M, 1);
x = ones(m, 1) / sqrt(m);
if norm(M * x) == 0                  % vector đầu nằm trong nhân của M
    rng(0); x = randn(m, 1); x = x / norm(x);
end
lam_cu = 0;
for k = 1:max_lap
    y = M * x;
    ny = sqrt(y' * y);
    if ny == 0, n = 0; return; end   % A là ma trận không
    x = y / ny;
    lam = x' * M * x;
    if abs(lam - lam_cu) <= tol * abs(lam), break; end
    lam_cu = lam;
end
n = sqrt(lam);
end
