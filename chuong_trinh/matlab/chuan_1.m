function n = chuan_1(A)
%CHUAN_1  Chuẩn 1 (chuẩn cột): ||A||_1 = max_j sum_i |a_ij|
n = max(sum(abs(A), 1));
end
