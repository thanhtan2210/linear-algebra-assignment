function n = chuan_F(A)
%CHUAN_F  Chuẩn Frobenius: ||A||_F = sqrt(sum_ij |a_ij|^2)
n = sqrt(sum(abs(A(:)).^2));
end
