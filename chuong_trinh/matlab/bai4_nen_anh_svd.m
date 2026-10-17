%% BÀI TOÁN 4 - Đo độ chênh lệch dữ liệu ảnh bằng chuẩn Frobenius (mục 4.3)
% Ảnh xám m x n = ma trận A (độ sáng 0..255).
% (a) Nén ảnh bằng SVD, giữ k thành phần -> A_k. Định lý Eckart-Young:
%       ||A - A_k||_F = sqrt(sigma_(k+1)^2 + ... + sigma_r^2),   ||A - A_k||_2 = sigma_(k+1)
% (b) Đo độ chênh lệch giữa ảnh gốc và ảnh nhiễu; cắt SVD để khử nhiễu.
clc; clear;
NGUONG = 0.05;                                   % yêu cầu sai số tương đối <= 5%

% Đổi đường dẫn này để dùng ảnh của bạn (ảnh màu sẽ được chuyển sang ảnh xám)
duong_dan_anh = fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'anh_mau.png');
if exist(duong_dan_anh, 'file')
    I = double(imread(duong_dan_anh));
    if size(I, 3) >= 3
        I = 0.299 * I(:, :, 1) + 0.587 * I(:, :, 2) + 0.114 * I(:, :, 3);
    end
    A = I;
else
    A = tao_anh_mau(256);
end
[m, n] = size(A);
fprintf('BÀI TOÁN 4: NÉN VÀ SO SÁNH ẢNH BẰNG CHUẨN MA TRẬN\n');

%% (a) NÉN ẢNH BẰNG SVD
[U, S, V] = svd(A, 'econ');
s = diag(S);
nA = chuan_F(A);
Ak_fun = @(k) U(:, 1:k) * S(1:k, 1:k) * V(:, 1:k)';
fprintf('%s\n(a) NÉN ẢNH BẰNG SVD   ảnh %dx%d, hạng = %d, ||A||_F = %.1f, ||A||_2 = sigma_1 = %.1f\n', ...
    repmat('=', 1, 92), m, n, sum(s > 1e-8 * s(1)), nA, s(1));
fprintf('  %4s%13s%15s%13s%11s%10s%9s%14s\n', 'k', '||A-Ak||_F', 'sqrt(sum s^2)', '||A-Ak||_2', ...
    's_(k+1)', 'Sai số %', 'PSNR dB', 'Dung lượng %');
for k = [1 2 5 10 20 30 50 80 120]
    E = A - Ak_fun(k);
    eF = chuan_F(E);
    psnr_k = 20 * log10(255 / sqrt(eF^2 / (m * n)));
    fprintf('  %4d%13.2f%15.2f%13.2f%11.2f%10.2f%9.2f%14.1f\n', k, eF, sqrt(sum(s(k+1:end).^2)), ...
        chuan_2(E), s(k+1), eF / nA * 100, psnr_k, k * (m + n + 1) / (m * n) * 100);
end
% Sai số tương đối e_k = sqrt(sum_{i>k} sigma_i^2)/||A||_F, k = 1..r (chỉ cần các sigma)
duoi = flipud(cumsum(flipud(s.^2)));             % duoi(i) = sum_{j>=i} sigma_j^2
e_k = sqrt([duoi(2:end); 0] / sum(s.^2));
k_chon = find(e_k <= NGUONG, 1);
fprintf('  => k nhỏ nhất để sai số <= %.0f%%: k = %d (sai số %.2f%%, lưu %.1f%% dung lượng ảnh gốc)\n', ...
    NGUONG * 100, k_chon, e_k(k_chon) * 100, k_chon * (m + n + 1) / (m * n) * 100);

%% (b) ẢNH NHIỄU
fprintf('%s\n(b) ĐỘ CHÊNH LỆCH GIỮA ẢNH GỐC VÀ ẢNH NHIỄU (nhiễu Gauss, độ lệch chuẩn 25)\n', repmat('=', 1, 92));
rng(7);
B = A + 25 * randn(size(A));
fprintf('  ||B - A||_F / ||A||_F = %.2f%%\n', chuan_F(B - A) / nA * 100);
[Ub, Sb, Vb] = svd(B, 'econ');
k_tot = 0; e_tot = inf;
for k = [5 10 20 30 40 60 100]
    Bk = Ub(:, 1:k) * Sb(1:k, 1:k) * Vb(:, 1:k)';
    e = chuan_F(Bk - A) / nA;
    fprintf('  k = %3d: ||B_k - A||_F / ||A||_F = %.2f%%\n', k, e * 100);
    if e < e_tot, e_tot = e; k_tot = k; B_tot = Bk; end
end
fprintf('  => Cắt SVD ở k = %d đưa ảnh nhiễu về gần ảnh gốc nhất (chênh lệch %.2f%%).\n', k_tot, e_tot * 100);

%% HÌNH
ds_k = unique([5 k_chon 50]);
fig = figure('Position', [50 100 1300 380], 'Color', 'w');
subplot(1, 4, 1); hien_anh(A, 'Ảnh gốc');
for t = 1:3
    Ak = Ak_fun(ds_k(t));
    subplot(1, 4, t + 1); hien_anh(Ak, sprintf('k = %d   sai số %.1f%%', ds_k(t), chuan_F(A - Ak) / nA * 100));
end
luu_hinh(fig, 'bai4_anh_nen_matlab.png');

r = min(numel(s) - 1, 150);
fig = figure('Position', [100 100 750 450]);
plot(1:r, e_k(1:r) * 100, '-', 'Color', mau_do_thi('xanh'), 'LineWidth', 2); hold on;
plot(1:r, s(2:r+1) / s(1) * 100, '-', 'Color', mau_do_thi('cam'), 'LineWidth', 2);
plot(k_chon, e_k(k_chon) * 100, 'o', 'Color', mau_do_thi('xanh'), 'MarkerFaceColor', mau_do_thi('xanh'), ...
    'HandleVisibility', 'off');
yline(NGUONG * 100, '--', 'Color', mau_do_thi('truc'), 'HandleVisibility', 'off');
text(k_chon + 3, e_k(k_chon) * 100 + 1.5, sprintf('k = %d đạt sai số ≤ %.0f%%', k_chon, NGUONG * 100));
xlabel('Số thành phần SVD giữ lại k'); ylabel('Sai số tương đối (%)');
title('Sai số nén ảnh theo chuẩn Frobenius và chuẩn 2');
legend('||A - A_k||_F / ||A||_F', '||A - A_k||_2 / ||A||_2 = \sigma_{k+1}/\sigma_1');
ylim([0 inf]);
dinh_dang_truc(gca);
luu_hinh(fig, 'bai4_sai_so_theo_k_matlab.png');

fig = figure('Position', [50 100 1000 380], 'Color', 'w');
subplot(1, 3, 1); hien_anh(A, 'Ảnh gốc A');
subplot(1, 3, 2); hien_anh(B, sprintf('Ảnh nhiễu B (%.1f%%)', chuan_F(B - A) / nA * 100));
subplot(1, 3, 3); hien_anh(B_tot, sprintf('B cắt SVD, k = %d (%.1f%%)', k_tot, e_tot * 100));
luu_hinh(fig, 'bai4_khu_nhieu_matlab.png');

function hien_anh(X, tieu_de)
image(min(max(X, 0), 255)); colormap(gca, gray(256));
axis image off; title(tieu_de);
end
