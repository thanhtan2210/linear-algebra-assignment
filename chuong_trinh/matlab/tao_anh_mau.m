function A = tao_anh_mau(N)
%TAO_ANH_MAU  Ảnh thử nghiệm tổng hợp (cùng công thức với bai4_nen_anh_svd.py).
if nargin < 1, N = 256; end
[x, y] = meshgrid((0:N-1) / (N-1));
A = 90 + 80 * x;                                         % nền chuyển sắc
A((x - 0.30).^2 + (y - 0.35).^2 < 0.18^2) = 235;         % hình tròn sáng
A((x - 0.72).^2 + (y - 0.70).^2 < 0.14^2) = 30;          % hình tròn tối
vung = (x > 0.55) & (x < 0.92) & (y > 0.08) & (y < 0.40);
A(vung) = 128 + 100 * sin(2 * pi * 10 * (x(vung) + y(vung)));   % vùng sọc chéo
A(y > x + 0.45) = 190;                                   % tam giác góc dưới trái
A = round(A);
end
