function n = chuan_vc(A)
%CHUAN_VC  Chuẩn vô cùng (chuẩn hàng): ||A||_inf = max_i sum_j |a_ij|
n = max(sum(abs(A), 2));
end
