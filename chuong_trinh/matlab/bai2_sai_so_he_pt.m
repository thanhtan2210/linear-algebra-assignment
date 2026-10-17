%% BÀI TOÁN 2 - Đánh giá sai số và độ ổn định của hệ Ax = b (mục 3.4, 3.5, 4.2)
% Nhiễu vế phải:  (1/k)||db||/||b|| <= ||dx||/||x|| <= k ||db||/||b||,  k = ||A||.||A^-1||
% Nhiễu ma trận:  ||dx||/||x|| <= k/(1 - k||dA||/||A||) * ||dA||/||A||
% Phần A: mạch điện 3 vòng (Kirchhoff) - điều kiện TỐT
% Phần B: xấp xỉ đa thức bình phương tối thiểu -> ma trận Hilbert - điều kiện RẤT XẤU
clc; clear; rng(2026);

%% PHẦN A: MẠCH ĐIỆN 3 VÒNG  R*I = E
fprintf('%s\nPHẦN A: MẠCH ĐIỆN 3 VÒNG - phương pháp dòng mắt lưới  R.I = E\n', repmat('=', 1, 76));
% R1=10, R2=5, R3=15 (riêng từng vòng); R12=5 (chung vòng 1-2), R23=10 (chung vòng 2-3)
R = [15 -5 0; -5 20 -10; 0 -10 25];
E = [12; 0; 6];                          % nguồn 12V ở vòng 1, 6V ở vòng 3
I = R \ E;
disp('  R ='); disp(R);
fprintf('  E = [%g %g %g] (V)\n  Dòng điện mắt lưới I = [%.5f %.5f %.5f] (A)\n', E, I);

fprintf('\n  A1) Nguồn đo có sai số +-1%%  (nhiễu vế phải b)\n');
fprintf('  %-16s%9s%15s%15s%12s%12s  Đúng?\n', 'Chuẩn', 'k(R)', '||dE||/||E||', '||dI||/||I||', 'Cận dưới', 'Cận trên');
ten_p = {'Chuẩn 1', 'Chuẩn 2', 'Chuẩn vô cùng'};
ds_p = {1, 2, inf};
for t = 1:3
    p = ds_p{t};
    dE = E .* (2 * rand(3, 1) - 1) * 0.01;
    dI = R \ (E + dE) - I;
    k = so_dieu_kien(R, p);
    eb = norm(dE, p) / norm(E, p);
    ex = norm(dI, p) / norm(I, p);
    kq = 'SAI'; if eb / k <= ex && ex <= k * eb, kq = 'ĐÚNG'; end
    fprintf('  %-16s%9.4f%15.3e%15.3e%12.3e%12.3e  %s\n', ten_p{t}, k, eb, ex, eb / k, k * eb, kq);
end

fprintf('\n  A2) Điện trở có dung sai +-5%%  (nhiễu ma trận A), dùng chuẩn vô cùng\n');
dR = R .* (2 * rand(3) - 1) * 0.05;
dI = (R + dR) \ E - I;
k = so_dieu_kien(R, inf);
eA = chuan_vc(dR) / chuan_vc(R);
ex = norm(dI, inf) / norm(I, inf);
can = k / (1 - k * eA) * eA;
kq = 'SAI'; if ex <= can, kq = 'ĐÚNG'; end
fprintf('    ||dR||/||R|| = %.4f, sai số thực tế ||dI||/||I|| = %.4f, cận trên = %.4f -> %s\n', eA, ex, can, kq);

fprintf('\n  A3) Trường hợp XẤU NHẤT (chuẩn 2): E theo hướng u_1, nhiễu theo hướng u_n (SVD)\n');
[U, S, ~] = svd(R);
s = diag(S);
E_xau = U(:, 1);  dE_xau = 1e-3 * U(:, end);
he_so = (norm(R \ dE_xau) / norm(R \ E_xau)) / (norm(dE_xau) / norm(E_xau));
fprintf('    Hệ số khuếch đại sai số = %.4f = k_2(R) = %.4f -> cận trên đạt được\n', he_so, s(1) / s(end));
k2 = so_dieu_kien(R, 2);
fprintf('\n  Đánh giá: k_2(R) = %.3f -> %s\n', k2, danh_gia_on_dinh(k2));

%% PHẦN B: MA TRẬN HILBERT (xấp xỉ đa thức bậc n-1 trên [0,1])
fprintf('\n%s\nPHẦN B: XẤP XỈ ĐA THỨC TRÊN [0,1] - hệ chuẩn tắc H.c = b, H là ma trận Hilbert\n', repmat('=', 1, 76));
fprintf('  Nghiệm đúng c = (1,...,1); b được làm nhiễu với ||db||_2/||b||_2 = 1e-8\n');
eps_b = 1e-8;
ds_n = 2:10;
ds_ss = zeros(size(ds_n));  ds_can = ds_ss;
fprintf('  %3s%16s%14s%22s%12s  Đánh giá\n', 'n', 'k_2 (tự tính)', 'k_2 cond()', 'Sai số ||dc||/||c||', 'Cận k.eps');
for t = 1:numel(ds_n)
    n = ds_n(t);
    H = hilb(n);
    c = ones(n, 1);
    b = H * c;
    huong = randn(n, 1);
    db = eps_b * norm(b) * huong / norm(huong);
    ss = norm(H \ (b + db) - c) / norm(c);
    k_tu_tinh = chuan_2(H) * chuan_2(inv(H));
    k_ml = cond(H, 2);
    ds_ss(t) = ss;  ds_can(t) = k_ml * eps_b;
    dg = danh_gia_on_dinh(k_ml);
    fprintf('  %3d%16.3e%14.3e%22.3e%12.3e  %s\n', n, k_tu_tinh, k_ml, ss, k_ml * eps_b, extractBefore(dg, ','));
end

fig = figure('Position', [100 100 750 450]);
semilogy(ds_n, ds_can, '-o', 'Color', mau_do_thi('cam'), 'LineWidth', 2, 'MarkerFaceColor', mau_do_thi('cam'));
hold on;
semilogy(ds_n, ds_ss, '-o', 'Color', mau_do_thi('xanh'), 'LineWidth', 2, 'MarkerFaceColor', mau_do_thi('xanh'));
yline(1, '--', 'sai số 100%', 'Color', mau_do_thi('truc'));
xlabel('Cấp n của ma trận Hilbert (bậc đa thức + 1)');
ylabel('Sai số tương đối (thang log)');
title('Nhiễu 10^{-8} ở vế phải được khuếch đại theo số điều kiện');
legend('Cận trên lý thuyết  k_2(H)·||δb||/||b||', 'Sai số tương đối thực tế ||δc||/||c||', 'Location', 'northwest');
dinh_dang_truc(gca);
luu_hinh(fig, 'bai2_hilbert_sai_so_matlab.png');

function s = danh_gia_on_dinh(k)
mat = max(0, log10(k));
if k < 1e2
    muc = 'TỐT (ổn định)';
elseif k < 1e6
    muc = 'TRUNG BÌNH';
else
    muc = 'XẤU (kém ổn định)';
end
s = sprintf('%s, có thể mất ~%.1f chữ số chính xác', muc, mat);
end
