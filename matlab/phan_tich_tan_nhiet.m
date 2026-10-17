%% PHÂN TÍCH NHIỆT ĐỘ TRÊN TẤM TẢN NHIỆT CPU BẰNG CHUẨN MA TRẬN
%   Bước 1 (mục 3.3, 4.1)     : Lập hệ A*theta = b, tính 4 chuẩn của A bằng hàm tự cài đặt.
%   Bước 2 (mục 3.4, 3.5, 4.2): Sai số cảm biến lan sang nghiệm; số điều kiện và lưới.
%   Bước 3 (mục 3.6)          : Giải bằng Jacobi / Gauss-Seidel, dừng bằng đánh giá theo chuẩn.
%   Bước 4 (mục 4.3)          : So sánh với "ảnh camera nhiệt", nén và khử nhiễu bằng SVD.
% Toàn bộ chương trình nằm trong file này: phần chạy chính ở trên, các hàm ở cuối file
%   (Phần 1: chuẩn tự cài đặt, Phần 2: mô hình tấm tản nhiệt, Phần 3: vẽ hình, Phần 4: hàm phụ).
% Chạy: chọn thư mục matlab làm Current Folder rồi gõ  phan_tich_tan_nhiet
clc; clear; close all;
EPS = 1e-3;                 % sai số cho phép khi giải lặp (°C)
SAI_SO_CAM_BIEN = 0.5;      % sai số mỗi điểm đo dọc cạnh CPU (°C)
NHIEU_CAMERA = 2.0;         % độ lệch chuẩn nhiễu của camera nhiệt (°C)
DO_PHAN_GIAI = 128;         % ảnh camera 128 x 128 điểm ảnh
rng(2026);

fprintf('PHÂN TÍCH NHIỆT ĐỘ TRÊN TẤM TẢN NHIỆT CPU BẰNG CHUẨN MA TRẬN\n');
ts = thong_so();
if ~batchStartupOptionUsed
    fprintf('Nhập thông số (Enter để giữ giá trị mặc định):\n');
    v = input(sprintf('  Số điểm lưới mỗi chiều [%d]: ', ts.N));               if ~isempty(v), ts.N = v; end
    v = input(sprintf('  Nhiệt độ CPU (°C) [%g]: ', ts.T_CPU));                if ~isempty(v), ts.T_CPU = v; end
    v = input(sprintf('  Hệ số toả nhiệt của quạt (W/m^2.K) [%g]: ', ts.h));   if ~isempty(v), ts.h = v; end
end
N = ts.N;

%% ============================ BƯỚC 1 ============================
tieu_de('BƯỚC 1 - LẬP MÔ HÌNH VÀ TÍNH CHUẨN CỦA MA TRẬN A  (mục 3.3, 4.1)');
[A, b, c] = lap_he(ts);
fprintf('  Tấm nhôm %g x %g cm, dày %g mm, T_CPU = %g°C, T_phòng = %g°C, h = %g W/m^2.K\n', ...
    ts.L*100, ts.L*100, ts.t*1000, ts.T_CPU, ts.T_inf, ts.h);
fprintf('  m^2 = 2h/(k.t) = %g 1/m^2, lưới %dx%d -> A cỡ %dx%d, c = m^2.dx^2 = %.5f\n', ...
    he_so_m2(ts), N, N, N^2, N^2, c);

fprintf('\n  a) Kiểm tra hàm tự cài đặt trên các ví dụ của phần lý thuyết:\n');
ten_vd = {'[2 5 4; 6 4 3; 5 3 5]', '[-3 5 -2; 1 -7 4; 2 0 9]', '[5 2; 2 5]', '[1 3; 0 2]'};
vd = {[2 5 4; 6 4 3; 5 3 5], [-3 5 -2; 1 -7 4; 2 0 9], [5 2; 2 5], [1 3; 0 2]};
ghi_chu = {'chuẩn 1 = 13', 'chuẩn ∞ = 12', 'chuẩn 2 = 7', 'chuẩn 2 = 3.702'};
fprintf('     %-26s%9s%11s%10s%10s   Báo cáo\n', 'Ma trận', '||.||_1', '||.||_inf', '||.||_F', '||.||_2');
for i = 1:4
    M = vd{i};
    fprintf('     %-26s%9.4f%11.4f%10.4f%10.4f   %s\n', ten_vd{i}, chuan_1(M), chuan_vo_cung(M), ...
        chuan_frobenius(M), chuan_2(M), ghi_chu{i});
end

fprintf('\n  b) Các chuẩn của ma trận A (%dx%d):\n', N^2, N^2);
[lam_min, lam_max] = tri_rieng_A(ts, N);
[n2, so_lan] = chuan_2(A);
nF_ct = sqrt(N^2 * (4 + c)^2 + 4 * N * (N - 1));
ten_chuan = {'Chuẩn 1', 'Chuẩn vô cùng', 'Frobenius (đ.nghĩa)', 'Frobenius (vết)', 'Chuẩn 2 (lũy thừa)'};
tu_tinh = [chuan_1(A), chuan_vo_cung(A), chuan_frobenius(A), chuan_frobenius_vet(A), n2];
co_san  = [norm(A, 1), norm(A, inf), norm(A, 'fro'), norm(A, 'fro'), norm(A, 2)];
cong_thuc = [8 + c, 8 + c, nF_ct, nF_ct, lam_max];
fprintf('     %-22s%14s%14s%14s\n', 'Chuẩn', 'Tự cài đặt', 'norm()', 'Công thức');
for i = 1:5
    fprintf('     %-22s%14.6f%14.6f%14.6f\n', ten_chuan{i}, tu_tinh(i), co_san(i), cong_thuc(i));
end
fprintf('     (lũy thừa hội tụ sau %d vòng; công thức: ||A||_1 = ||A||_inf = 8 + c, ||A||_2 = lambda_max)\n', so_lan);
fprintf('     Kiểm chứng: ||A||_2 = %.4f <= ||A||_F = %.4f <= sqrt(r)||A||_2 = %.4f;  ||A||_2 <= sqrt(||A||_1.||A||_inf) = %.4f\n', ...
    n2, chuan_frobenius(A), N * n2, sqrt(chuan_1(A) * chuan_vo_cung(A)));

%% ============================ BƯỚC 2 ============================
tieu_de('BƯỚC 2 - SAI SỐ CẢM BIẾN VÀ ĐỘ ỔN ĐỊNH CỦA HỆ  (mục 3.4, 3.5, 4.2)');
theta = A \ b;
fprintf('  a) Số điều kiện của A:\n');
fprintf('     kappa_1   = %10.4f\n     kappa_2   = %10.4f\n     kappa_inf = %10.4f\n', ...
    so_dieu_kien(A, 1), so_dieu_kien(A, 2), so_dieu_kien(A, inf));
k2 = so_dieu_kien(A, 2);
fprintf('     Công thức: kappa_2 = lambda_max/lambda_min = %.4f -> mất khoảng %.1f chữ số chính xác: hệ ỔN ĐỊNH\n', ...
    lam_max / lam_min, log10(k2));

fprintf('\n  b) Mỗi điểm đo dọc cạnh CPU sai ngẫu nhiên tới +-%g°C (nhiễu vế phải b):\n', SAI_SO_CAM_BIEN);
db = zeros(size(b));
db(1:N) = SAI_SO_CAM_BIEN * (2 * rand(N, 1) - 1);
dtheta = A \ db;
fprintf('     %-14s%14s%14s%12s%12s  Đúng?\n', 'Chuẩn', '||db||/||b||', '||dθ||/||θ||', 'Cận dưới', 'Cận trên');
ds_p = {1, 2, inf};  ten_p = {'Chuẩn 1', 'Chuẩn 2', 'Chuẩn vô cùng'};
for t = 1:3
    p = ds_p{t};
    k = so_dieu_kien(A, p);
    eb = norm(db, p) / norm(b, p);
    ex = norm(dtheta, p) / norm(theta, p);
    kq = 'SAI'; if eb / k <= ex && ex <= k * eb, kq = 'ĐÚNG'; end
    fprintf('     %-14s%14.3e%14.3e%12.3e%12.3e  %s\n', ten_p{t}, eb, ex, eb / k, k * eb, kq);
end
fprintf('     Theo °C: sai số lớn nhất thực tế = %.4f°C <= ||A^-1||_inf.||db||_inf = %.4f°C\n', ...
    norm(dtheta, inf), chuan_vo_cung(inv(A)) * norm(db, inf));
[U, S, ~] = svd(A);
b_xau = U(:, 1);  db_xau = 1e-3 * U(:, end);
he_so = (norm(A \ db_xau) / norm(A \ b_xau)) / (norm(db_xau) / norm(b_xau));
fprintf('     Trường hợp xấu nhất (b theo u_1, nhiễu theo u_n): hệ số khuếch đại = %.4f = kappa_2 -> cận trên đạt được\n', he_so);

fprintf('\n  c) Chia lưới mịn hơn: mô hình chính xác hơn nhưng kappa tăng và phép lặp chậm hơn\n');
fprintf('     %4s%8s%19s%21s%11s%16s\n', 'N', 'Số ẩn', 'kappa_2 (tự tính)', 'kappa_2 (công thức)', 'q Jacobi', 'Sai số mô hình');
ds_N_bang = [5 10 20 40];
k_bang = zeros(size(ds_N_bang));
for t = 1:numel(ds_N_bang)
    Nn = ds_N_bang(t);
    [An, bn, cn] = lap_he(ts, Nn);
    [X, Y] = toa_do_nut(ts, Nn);
    chinh_xac = reshape(nghiem_giai_tich(ts, X, Y)', [], 1);    % đúng thứ tự ẩn (j-1)*N + i
    sai_so = chuan_frobenius(An \ bn - chinh_xac) / chuan_frobenius(chinh_xac);
    [lmin, lmax] = tri_rieng_A(ts, Nn);
    k_bang(t) = chuan_2(An) * chuan_2(inv(An));
    fprintf('     %4d%8d%19.4f%21.4f%11.5f%15.3f%%\n', Nn, Nn^2, k_bang(t), lmax / lmin, 4 / (4 + cn), sai_so * 100);
end
fprintf('     (Sai số mô hình = ||θ_lưới - θ_chính xác||_F / ||θ_chính xác||_F, θ_chính xác tính bằng chuỗi Fourier)\n');
fprintf('     Đối chứng: ma trận Hilbert 10x10 có kappa_2 = %.2e -> dù có %d ẩn, hệ tấm tản nhiệt vẫn ổn định hơn rất nhiều.\n', ...
    cond(hilb(10)), N^2);

%% ============================ BƯỚC 3 ============================
tieu_de(sprintf('BƯỚC 3 - GIẢI LẶP VÀ THUẬT TOÁN DỪNG  (mục 3.6),  eps = %g°C', EPS));
fprintf('  Dừng khi đánh giá hậu nghiệm  q/(1-q).||x_k - x_(k-1)||_inf <= eps,  q = hệ số co (chuẩn vô cùng)\n');
ds_pp = {'Jacobi', 'Gauss-Seidel'};
kq_lap = struct();
for t = 1:2
    pp = ds_pp{t};
    [x, q, k_tien, ls] = giai_lap(A, b, pp, EPS, theta);
    k_nt = find(ls(:, 2) <= EPS, 1);
    kq_lap(t).x = x;  kq_lap(t).q = q;  kq_lap(t).ls = ls;
    fprintf('\n  %s:  q = %.5f < 1 => hội tụ;  q/(1-q) = %.1f\n', upper(pp), q, q / (1 - q));
    fprintf('     %6s%18s%16s%14s\n', 'k', '||x_k - x_k-1||', 'q/(1-q).||.||', 'Sai số thật');
    Nk = size(ls, 1);
    for i = [1 2 3 0 Nk-2 Nk-1 Nk]
        if i == 0, fprintf('     %6s\n', '...'); continue; end
        fprintf('     %6d%18.4e%16.4e%14.4e\n', ls(i, :));
    end
    fprintf('     * Dừng ở k = %d: sai số thật %.2e <= eps  (tiên nghiệm dự báo cần %d vòng)\n', Nk, ls(end, 4), k_tien);
    fprintf('     * Nếu dừng "ngây thơ" khi ||x_k - x_k-1|| <= eps: dừng ở k = %d, sai số thật %.2e = %.0f lần eps  => KHÔNG đạt\n', ...
        k_nt, ls(k_nt, 4), ls(k_nt, 4) / EPS);
end
T_luoi = reshape(kq_lap(2).x, N, N)' + ts.T_inf;      % hàng = y, cột = x
fprintf('\n  Kết quả: nhiệt độ cao nhất trong tấm %.2f°C (sát CPU), tại tâm tấm %.2f°C\n', ...
    max(T_luoi(:)), T_luoi(floor(N/2) + 1, floor(N/2) + 1));

%% ============================ BƯỚC 4 ============================
tieu_de('BƯỚC 4 - SO SÁNH VỚI ẢNH CAMERA NHIỆT, NÉN VÀ KHỬ NHIỄU BẰNG SVD  (mục 4.3)');
n = DO_PHAN_GIAI;
[X, Y] = toa_do_nut(ts, n);
that = nghiem_giai_tich(ts, X, Y);                   % trường nhiệt thật (độ tăng nhiệt θ)
camera = that + NHIEU_CAMERA * randn(size(that));    % ảnh camera có nhiễu
nT = chuan_frobenius(that);
fprintf('  Ảnh nhiệt %dx%d điểm ảnh = ma trận θ (độ tăng nhiệt so với phòng, °C).\n', n, n);
fprintf('  Camera có nhiễu %g°C:  ||θ_camera - θ_thật||_F / ||θ_thật||_F = %.2f%%\n', NHIEU_CAMERA, chuan_frobenius(camera - that) / nT * 100);
[Uc, Sc, Vc] = svd(camera);
s = diag(Sc);
nC = chuan_frobenius(camera);
fprintf('\n  SVD ảnh camera, giữ k thành phần (Eckart-Young: ||C - C_k||_F = sqrt(sum s_i^2, i>k), ||C - C_k||_2 = s_(k+1))\n');
fprintf('     %4s%12s%14s%12s%10s%17s%18s%12s\n', 'k', '||C-Ck||_F', 'sqrt(sum s^2)', '||C-Ck||_2', 's_(k+1)', ...
    'Lệch ảnh camera', 'Lệch trường thật', 'Dung lượng');
ds_k = 1:60;
lech_camera = zeros(size(ds_k));  lech_that = lech_camera;
for k = ds_k
    Ck = Uc(:, 1:k) * Sc(1:k, 1:k) * Vc(:, 1:k)';
    lech_camera(k) = chuan_frobenius(camera - Ck) / nC;
    lech_that(k) = chuan_frobenius(Ck - that) / nT;
    if ismember(k, [1 2 3 5 10 20 50])
        E = camera - Ck;
        fprintf('     %4d%12.2f%14.2f%12.2f%10.2f%16.2f%%%17.2f%%%11.1f%%\n', k, chuan_frobenius(E), sqrt(sum(s(k+1:end).^2)), ...
            chuan_2(E), s(k+1), lech_camera(k) * 100, lech_that(k) * 100, k * (2*n + 1) / n^2 * 100);
    end
end
[lech_min, k_tot] = min(lech_that);
Ck_tot = Uc(:, 1:k_tot) * Sc(1:k_tot, 1:k_tot) * Vc(:, 1:k_tot)';
fprintf('  => Giữ k = %d thành phần: ảnh chỉ còn lệch trường thật %.2f%% (so với %.2f%% của ảnh gốc) và chỉ cần lưu %.1f%% dữ liệu.\n', ...
    k_tot, lech_min * 100, chuan_frobenius(camera - that) / nT * 100, k_tot * (2*n + 1) / n^2 * 100);
fprintf('     Trường nhiệt trơn nên có hạng thấp; nhiễu trải đều trên mọi thành phần nên bị cắt bỏ.\n\n');

%% ============================ HÌNH VẼ ============================
Lcm = ts.L * 100;
% Hình 1: bản đồ nhiệt
T_day_du = ts.T_inf * ones(N + 2);
T_day_du(2:end-1, 2:end-1) = T_luoi;
T_day_du(1, 2:end-1) = ts.T_CPU;
fig = figure('Position', [100 100 620 520]);
imagesc([0 Lcm], [0 Lcm], T_day_du); set(gca, 'YDir', 'normal'); axis image;
colormap(gca, mau_nhiet()); cb = colorbar; cb.Label.String = 'Nhiệt độ (°C)';
xlabel(sprintf('x (cm)  —  cạnh dưới áp vào CPU (%g°C)', ts.T_CPU)); ylabel('y (cm)');
title(sprintf('Phân bố nhiệt độ trên tấm (lưới %dx%d, Gauss-Seidel)', N, N));
luu_hinh(fig, 'hinh1_ban_do_nhiet_matlab.png');

% Hình 2: số điều kiện theo lưới
ds_N = [4 5 7 10 14 20 28 40 57 80 113 160];
k_ct = zeros(size(ds_N));
for t = 1:numel(ds_N), [a1, a2] = tri_rieng_A(ts, ds_N(t)); k_ct(t) = a2 / a1; end
fig = figure('Position', [100 100 750 450]);
loglog(ds_N, k_ct, '-', 'Color', mau_ve('xanh'), 'LineWidth', 2); hold on;
loglog(ds_N_bang, k_bang, 'o', 'Color', mau_ve('cam'), 'MarkerFaceColor', mau_ve('cam'), 'MarkerSize', 8);
xlabel('Số điểm lưới mỗi chiều N  (số ẩn = N^2)'); ylabel('Số điều kiện \kappa_2(A)');
title('Lưới càng mịn, hệ càng nhạy với sai số (\kappa_2 tăng cỡ N^2)');
legend('Công thức \lambda_{max}/\lambda_{min}', 'Tự tính ||A||_2·||A^{-1}||_2 (lũy thừa)', 'Location', 'northwest');
dinh_dang_truc(gca);
luu_hinh(fig, 'hinh2_so_dieu_kien_matlab.png');

% Hình 3: hội tụ
fig = figure('Position', [100 100 750 450]);
mau_duong = {mau_ve('xanh'), mau_ve('cam')};
for t = 1:2
    ls = kq_lap(t).ls;
    semilogy(ls(:, 1), ls(:, 4), '-', 'Color', mau_duong{t}, 'LineWidth', 2); hold on;
end
for t = 1:2
    ls = kq_lap(t).ls;
    semilogy(ls(end, 1), ls(end, 4), 'o', 'Color', mau_duong{t}, 'MarkerFaceColor', mau_duong{t}, 'HandleVisibility', 'off');
    text(ls(end, 1), ls(end, 4) * 2, sprintf('dừng k = %d', ls(end, 1)), 'HorizontalAlignment', 'right');
end
yline(EPS, '--', sprintf('\\epsilon = %g °C', EPS), 'Color', mau_ve('truc'), 'HandleVisibility', 'off');
xlabel('Số vòng lặp k'); ylabel('Sai số thật ||\theta^{(k)} - \theta^*||_\infty  (°C)');
title('Hội tụ của phép lặp - dừng theo đánh giá hậu nghiệm');
legend(sprintf('Jacobi (q = %.4f)', kq_lap(1).q), sprintf('Gauss-Seidel (q = %.4f)', kq_lap(2).q));
dinh_dang_truc(gca);
luu_hinh(fig, 'hinh3_hoi_tu_matlab.png');

% Hình 4: ảnh camera
fig = figure('Position', [50 100 1200 400], 'Color', 'w');
anh = {that, camera, Ck_tot};
td = {'Trường nhiệt thật \theta', sprintf('Ảnh camera (lệch %.1f%%)', chuan_frobenius(camera - that) / nT * 100), ...
      sprintf('Camera cắt SVD, k = %d (lệch %.1f%%)', k_tot, lech_min * 100)};
for t = 1:3
    subplot(1, 3, t);
    imagesc(anh{t}, [min(that(:)) max(that(:))]); set(gca, 'YDir', 'normal');
    colormap(gca, mau_nhiet()); axis image off; title(td{t});
end
luu_hinh(fig, 'hinh4_anh_camera_matlab.png');

% Hình 5: độ lệch theo k
fig = figure('Position', [100 100 750 450]);
plot(ds_k, lech_that * 100, '-', 'Color', mau_ve('cam'), 'LineWidth', 2); hold on;
plot(ds_k, lech_camera * 100, '-', 'Color', mau_ve('xanh'), 'LineWidth', 2);
plot(k_tot, lech_min * 100, 'o', 'Color', mau_ve('cam'), 'MarkerFaceColor', mau_ve('cam'), 'HandleVisibility', 'off');
text(k_tot + 2, lech_min * 100 - 0.6, sprintf('k = %d: gần trường thật nhất', k_tot));
xlabel('Số thành phần SVD giữ lại k'); ylabel('Độ lệch tương đối (%)');
title('Cắt SVD: giữ đủ thông tin, bỏ phần nhiễu');
legend('||C_k - \theta_{thật}||_F / ||\theta_{thật}||_F', '||C - C_k||_F / ||C||_F  (so với ảnh camera)');
ylim([0 inf]);
dinh_dang_truc(gca);
luu_hinh(fig, 'hinh5_do_lech_svd_matlab.png');


%% ====================== PHẦN 1. CÁC CHUẨN MA TRẬN TỰ CÀI ĐẶT ======================
% Tự cài đặt theo đúng công thức lý thuyết (không dùng hàm norm có sẵn):
%     ||A||_1   = max_j  sum_i |a_ij|          (chuẩn cột)
%     ||A||_inf = max_i  sum_j |a_ij|          (chuẩn hàng)
%     ||A||_F   = sqrt( sum_ij |a_ij|^2 ) = sqrt( tr(A' A) )
%     ||A||_2   = sqrt( lambda_max(A' A) ) = sigma_max(A)

function n = chuan_1(A)
    % Chuẩn 1: tổng trị tuyệt đối lớn nhất theo CỘT.
    n = max(sum(abs(A), 1));
end

function n = chuan_vo_cung(A)
    % Chuẩn vô cùng: tổng trị tuyệt đối lớn nhất theo HÀNG.
    n = max(sum(abs(A), 2));
end

function n = chuan_frobenius(A)
    % Chuẩn Frobenius theo công thức định nghĩa.
    n = sqrt(sum(abs(A(:)).^2));
end

function n = chuan_frobenius_vet(A)
    % Chuẩn Frobenius qua vết: sqrt(trace(A'*A)).
    n = sqrt(trace(A' * A));
end

function [n, k] = chuan_2(A, tol, max_lap)
    % Chuẩn 2 (chuẩn phổ) = sqrt(lambda_max(A'*A)).
    % lambda_max tìm bằng PHƯƠNG PHÁP LŨY THỪA:
    %     y = M*x,  x = y/||y||,  lambda = x'*M*x
    % Dừng khi |lambda_k - lambda_(k-1)| <= tol*lambda_k.  [n, k] trả thêm số vòng lặp.
    if nargin < 2, tol = 1e-12; end
    if nargin < 3, max_lap = 10000; end
    M = A' * A;                      % đối xứng, xác định không âm
    m = size(M, 1);
    % Vector đầu ngẫu nhiên (cố định seed): vector (1,...,1) có thể vuông góc với
    % vector riêng cần tìm (xảy ra với ma trận lưới đối xứng) và cho kết quả sai.
    x = randn(RandStream('mt19937ar', 'Seed', 0), m, 1);
    x = x / norm(x);
    lam_cu = 0;
    for k = 1:max_lap
        y = M * x;
        ny = sqrt(y' * y);
        if ny == 0, n = 0; return; end      % A là ma trận không
        x = y / ny;
        lam = x' * M * x;
        if abs(lam - lam_cu) <= tol * abs(lam), break; end
        lam_cu = lam;
    end
    n = sqrt(lam);
end

function k = so_dieu_kien(A, p)
    % Số điều kiện kappa_p(A) = ||A||_p * ||A^-1||_p,  p = 1, 2, inf hoặc 'fro'.
    if nargin < 2, p = 2; end
    Ai = inv(A);
    if ischar(p) || isstring(p)
        k = chuan_frobenius(A) * chuan_frobenius(Ai);
    elseif p == 1
        k = chuan_1(A) * chuan_1(Ai);
    elseif p == 2
        k = chuan_2(A) * chuan_2(Ai);
    else
        k = chuan_vo_cung(A) * chuan_vo_cung(Ai);
    end
end

%% ====================== PHẦN 2. MÔ HÌNH TẤM TẢN NHIỆT CPU ======================
% Tấm nhôm vuông cạnh L, dày t. Cạnh dưới (y = 0) áp vào CPU ở nhiệt độ T_CPU, ba cạnh
% còn lại ở nhiệt độ phòng T_inf, hai mặt tấm toả nhiệt (quạt, hệ số h). theta = T - T_inf:
%     theta_xx + theta_yy - m^2 theta = 0,     m^2 = 2h / (k t)
% Sai phân 5 điểm, lưới N x N điểm trong, dx = L/(N+1), c = m^2 dx^2:
%     (4 + c) theta_ij - (4 điểm lân cận) = 0
% => A*theta = b,  A = kron(I, K) + kron(K, I) + c I,  K = tridiag(-1, 2, -1).

function ts = thong_so()
    % Thông số mặc định (có thể nhập lại khi chạy chương trình).
    ts.L     = 0.10;     % cạnh tấm (m)
    ts.t     = 0.001;    % độ dày tấm (m)
    ts.k     = 200;      % hệ số dẫn nhiệt của nhôm (W/m.K)
    ts.h     = 100;      % hệ số toả nhiệt khi có quạt (W/m^2.K)
    ts.T_CPU = 80;       % nhiệt độ cạnh áp vào CPU (°C)
    ts.T_inf = 25;       % nhiệt độ phòng (°C)
    ts.N     = 20;       % số điểm lưới trong mỗi chiều
end

function m2 = he_so_m2(ts)
    m2 = 2 * ts.h / (ts.k * ts.t);
end

function [A, b, c] = lap_he(ts, N)
    % Lập hệ A*theta = b. Ẩn đánh số theo hàng: chỉ số = (j-1)*N + i (i theo x, j theo y).
    if nargin < 2, N = ts.N; end
    dx = ts.L / (N + 1);
    c = he_so_m2(ts) * dx^2;
    K = 2 * eye(N) - diag(ones(N-1, 1), 1) - diag(ones(N-1, 1), -1);
    A = kron(eye(N), K) + kron(K, eye(N)) + c * eye(N^2);
    b = zeros(N^2, 1);
    b(1:N) = ts.T_CPU - ts.T_inf;    % hàng j = 1 kề cạnh CPU (theta = T_CPU - T_inf)
end                                  % ba cạnh còn lại có theta = 0 nên không đóng góp vào b

function [lam_min, lam_max] = tri_rieng_A(ts, N)
    % Trị riêng nhỏ nhất / lớn nhất của A theo công thức giải tích (để kiểm chứng):
    %   lambda_pq = c + 4 sin^2(p pi / (2(N+1))) + 4 sin^2(q pi / (2(N+1)))
    c = he_so_m2(ts) * (ts.L / (N + 1))^2;
    s = @(p) sin(p * pi / (2 * (N + 1)))^2;
    lam_min = c + 8 * s(1);
    lam_max = c + 8 * s(N);
end

function theta = nghiem_giai_tich(ts, X, Y, so_hang)
    % Nghiệm chính xác theta(x, y) bằng chuỗi Fourier:
    %   theta = sum_n a_n sin(n pi x/L) sinh(beta_n (L-y)) / sinh(beta_n L)
    %   a_n = 2 theta0 (1 - (-1)^n)/(n pi),  beta_n = sqrt((n pi/L)^2 + m^2)
    if nargin < 4, so_hang = 400; end
    L = ts.L;  theta0 = ts.T_CPU - ts.T_inf;  m2 = he_so_m2(ts);
    theta = zeros(size(X));
    for n = 1:2:so_hang               % chỉ n lẻ có a_n khác 0
        beta = sqrt((n * pi / L)^2 + m2);
        % sinh(beta(L-y))/sinh(beta L) viết dạng mũ để tránh tràn số
        ti_so = exp(-beta * Y) .* (1 - exp(-2 * beta * (L - Y))) / (1 - exp(-2 * beta * L));
        theta = theta + 4 * theta0 / (n * pi) * sin(n * pi * X / L) .* ti_so;
    end
end

function [B, q] = ma_tran_lap(A, pp)
    % Ma trận lặp B = -D^-1 (L + U) và hệ số co q tương ứng (chuẩn vô cùng).
    D = diag(A);
    B = -(A - diag(D)) ./ D;
    if strcmp(pp, 'Jacobi')
        q = chuan_vo_cung(B);      % q = ||B||_inf
    else                                        % Gauss-Seidel
        q = 0;
        for i = 1:numel(D)
            q = max(q, sum(abs(B(i, i+1:end))) / (1 - sum(abs(B(i, 1:i-1)))));
        end
    end
end

function [x, q, k_tien, ls] = giai_lap(A, b, pp, EPS, x_dung)
    % Giải A x = b bằng 'Jacobi' hoặc 'Gauss-Seidel'.
    % Dừng khi đánh giá HẬU NGHIỆM  q/(1-q) ||x_k - x_(k-1)||_inf <= EPS.
    % ls: mỗi hàng [k, ||x_k - x_(k-1)||, đánh giá hậu nghiệm, sai số thật].
    n = numel(b);
    [~, q] = ma_tran_lap(A, pp);
    D = diag(A);
    R = A - diag(D);                  % phần ngoài đường chéo
    if nargin < 5, x_dung = nan(n, 1); end
    x = zeros(n, 1);
    ls = zeros(0, 4);  k_tien = NaN;
    for k = 1:20000
        if strcmp(pp, 'Jacobi')
            x_moi = (b - R * x) ./ D;
        else
            x_moi = x;
            for i = 1:n               % dùng ngay các thành phần vừa tính
                x_moi(i) = (b(i) - R(i, :) * x_moi) / D(i);
            end
        end
        buoc = norm(x_moi - x, inf);
        if k == 1
            k_tien = ceil(log(EPS * (1 - q) / buoc) / log(q));
        end
        hau_nghiem = q / (1 - q) * buoc;
        ls(end+1, :) = [k, buoc, hau_nghiem, norm(x_moi - x_dung, inf)]; %#ok<AGROW>
        x = x_moi;
        if hau_nghiem <= EPS, break; end
    end
end

%% ====================== PHẦN 3. CẤU HÌNH ĐỒ THỊ ======================
function c = mau_ve(ten)
    % Bảng màu dùng chung.
    switch ten
        case 'xanh', c = [42 120 214]/255;     % chuỗi số liệu 1
        case 'cam',  c = [235 104 52]/255;     % chuỗi số liệu 2
        case 'truc', c = [195 194 183]/255;
        otherwise,   c = [82 81 78]/255;       % chữ phụ
    end
end

function cmap = mau_nhiet(n)
    % Thang màu tuần tự một sắc độ cho bản đồ nhiệt (nhạt = nguội, đậm = nóng).
    if nargin < 1, n = 256; end
    moc = [255 245 235; 253 141 60; 127 39 4] / 255;
    cmap = interp1([0 0.5 1], moc, linspace(0, 1, n));
end

function dinh_dang_truc(ax)
    % Định dạng trục đồ thị đồng bộ với hình bên Python.
    set(ax, 'Color', [252 252 251]/255, 'XColor', [82 81 78]/255, 'YColor', [82 81 78]/255, ...
        'GridColor', [225 224 217]/255, 'GridAlpha', 1, 'Box', 'off', 'FontSize', 11, ...
        'XGrid', 'on', 'YGrid', 'on', 'TickDir', 'out');
    set(ancestor(ax, 'figure'), 'Color', [252 252 251]/255);
    lg = get(ax, 'Legend');
    if ~isempty(lg), lg.Box = 'off'; end
end

function luu_hinh(fig, ten_file)
    % Lưu hình vào thư mục ket_qua (cạnh thư mục matlab).
    thu_muc = fullfile(fileparts(mfilename('fullpath')), '..', 'ket_qua');
    if ~exist(thu_muc, 'dir'), mkdir(thu_muc); end
    duong_dan = fullfile(thu_muc, ten_file);
    try
        theme(fig, 'light');   % MATLAB R2025a+ mặc định giao diện tối -> chữ bị nhạt trên nền sáng
    catch
        % bản MATLAB cũ không có hàm theme, vốn đã dùng giao diện sáng
    end
    set(findall(fig, 'Type', 'text'), 'Color', [11 11 11]/255);
    exportgraphics(fig, duong_dan, 'Resolution', 150);
    fprintf('  -> Đã lưu hình: %s\n', duong_dan);
end

%% ====================== PHẦN 4. HÀM PHỤ ======================
function tieu_de(chu)
fprintf('\n%s\n%s\n%s\n', repmat('=', 1, 90), chu, repmat('=', 1, 90));
end

function [X, Y] = toa_do_nut(ts, N)
% Toạ độ các điểm lưới trong; X, Y cỡ N x N với hàng = y, cột = x
g = (1:N) * ts.L / (N + 1);
[X, Y] = meshgrid(g, g);
end
