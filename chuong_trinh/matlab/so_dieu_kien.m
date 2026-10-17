function k = so_dieu_kien(A, p)
%SO_DIEU_KIEN  kappa_p(A) = ||A||_p * ||A^-1||_p,  p = 1, 2, inf hoặc 'fro'
if nargin < 2, p = 2; end
Ai = inv(A);
if ischar(p) || isstring(p)
    k = chuan_F(A) * chuan_F(Ai);
elseif p == 1
    k = chuan_1(A) * chuan_1(Ai);
elseif p == 2
    k = chuan_2(A) * chuan_2(Ai);
else
    k = chuan_vc(A) * chuan_vc(Ai);
end
end
