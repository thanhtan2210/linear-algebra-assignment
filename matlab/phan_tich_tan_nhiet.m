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
ts = nhap_thong_so(thong_so());   % mặc định / tình huống mẫu / trả lời từng câu (Phần 5)
N = ts.N;

%% ============================ BƯỚC 1 ============================
tieu_de('BƯỚC 1 - LẬP MÔ HÌNH VÀ TÍNH CHUẨN CỦA MA TRẬN A  (mục 3.3, 4.1)');
[A, b, c] = lap_he(ts);
ds_cau = cau_hoi();
vat_lieu = lower(ten_lua_chon(ds_cau(1), ts.k));
if isempty(vat_lieu), vat_lieu = 'vật liệu'; end
fprintf(['  Tấm %s %g x %g cm, dày %g mm, k = %g W/m.K, T_CPU = %g°C, T_phòng = %g°C,' ...
    ' h = %g W/m^2.K\n'], vat_lieu, ts.L*100, ts.L*100, ts.t*1000, ts.k, ts.T_CPU, ts.T_inf, ts.h);
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
% Hình 0: sơ đồ đề bài
ve_de_bai(ts);

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

function ve_de_bai(ts)
    % Hình 0: sơ đồ đề bài - tấm nhìn từ trên xuống và nhìn ngang.
    L = ts.L * 100;
    N = ts.N;
    ds = cau_hoi();
    vat_lieu = lower(ten_lua_chon(ds(1), ts.k));
    if isempty(vat_lieu), vat_lieu = 'kim loại'; end
    nong = mau_ve('cam'); lanh = mau_ve('xanh'); phu = mau_ve('chu_phu');
    chu = [11 11 11]/255; tam = [217 213 204]/255; cpu = [58 57 54]/255;
    giua = {'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle'};
    fig = figure('Position', [50 100 1300 540], 'Color', [252 252 251]/255);

    % (a) Nhìn từ trên xuống: cạnh dưới áp vào CPU, ba cạnh kia ở nhiệt độ phòng
    ax = axes(fig, 'Position', [0.02 0.03 0.40 0.86]); hold(ax, 'on');
    rectangle(ax, 'Position', [0 0 L L], 'FaceColor', tam, 'EdgeColor', 'none');
    g = linspace(0, L, 10); g = g(2:end-1);
    [X, Y] = meshgrid(g, g);
    plot(ax, X(:), Y(:), 'o', 'MarkerSize', 3, 'MarkerFaceColor', phu, 'MarkerEdgeColor', phu);
    text(ax, L/2, L/2, {'Nhiệt độ tại mỗi', 'điểm lưới = ?'}, giua{:}, 'FontSize', 12, 'Color', chu, ...
        'BackgroundColor', 'w', 'EdgeColor', mau_ve('truc'), 'Margin', 6);
    plot(ax, [0 L], [0 0], 'Color', nong, 'LineWidth', 6);
    plot(ax, [0 0 L L], [0 L L 0], 'Color', lanh, 'LineWidth', 4);
    rectangle(ax, 'Position', [0 -0.16*L L 0.16*L], 'FaceColor', cpu, 'EdgeColor', 'none');
    text(ax, L/2, -0.08*L, sprintf('CPU  %g°C', ts.T_CPU), giua{:}, 'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold');
    text(ax, L/2, 0.05*L, sprintf('cạnh dưới: %g°C', ts.T_CPU), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'Color', nong, 'FontSize', 10);
    text(ax, L/2, 0.95*L, sprintf('%g°C (nhiệt độ phòng)', ts.T_inf), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', 'Color', lanh, 'FontSize', 10);
    text(ax, 0.05*L, L/2, sprintf('%g°C', ts.T_inf), giua{:}, 'Rotation', 90, 'Color', lanh, 'FontSize', 10);
    text(ax, 0.95*L, L/2, sprintf('%g°C', ts.T_inf), giua{:}, 'Rotation', -90, 'Color', lanh, 'FontSize', 10);
    mui_ten(ax, [0 L], [1.07 1.07]*L, phu, 1, 0.03*L, true);
    text(ax, L/2, 1.09*L, sprintf('%g cm', L), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'Color', phu, 'FontSize', 11);
    mui_ten(ax, [-0.07 -0.07]*L, [0 L], phu, 1, 0.03*L, true);
    text(ax, -0.09*L, L/2, sprintf('%g cm', L), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'Rotation', 90, 'Color', phu, 'FontSize', 11);
    axis(ax, 'equal', 'off');
    xlim(ax, [-0.2 1.05]*L); ylim(ax, [-0.2 1.2]*L);
    title(ax, sprintf('(a) Nhìn từ trên xuống: tấm %s vuông', vat_lieu), 'Color', chu, ...
        'FontWeight', 'normal', 'FontSize', 13);

    % (b) Nhìn ngang: nhiệt lan dọc tấm (k) và thoát ra hai mặt nhờ quạt (h); độ dày vẽ phóng to
    ax = axes(fig, 'Position', [0.45 0.03 0.54 0.86]); hold(ax, 'on');
    rectangle(ax, 'Position', [-1.6 -1.1 1.6 2.8], 'FaceColor', cpu, 'EdgeColor', 'none');
    text(ax, -0.8, 0.3, 'CPU', giua{:}, 'Rotation', 90, 'Color', 'w', 'FontWeight', 'bold', 'FontSize', 11);
    rectangle(ax, 'Position', [0 0 10 0.6], 'FaceColor', tam, 'EdgeColor', mau_ve('truc'));
    plot(ax, [0 0], [0 0.6], 'Color', nong, 'LineWidth', 5);
    plot(ax, [10 10], [0 0.6], 'Color', lanh, 'LineWidth', 5);
    mui_ten(ax, [0.4 9.3], [0.3 0.3], nong, 2, 0.45, false);
    for x = [2 4 6 8]
        mui_ten(ax, [x x], [0.75 1.9], lanh, 1.8, 0.3, false);
        mui_ten(ax, [x x], [-0.15 -1.3], lanh, 1.8, 0.3, false);
    end
    text(ax, 5, 2.1, {sprintf('Quạt thổi hai mặt tấm: nhiệt thoát ra không khí (%g°C)', ts.T_inf), ...
        sprintf('hệ số toả nhiệt h = %g W/m^2.K', ts.h)}, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'Color', lanh, 'FontSize', 10);
    text(ax, 5, -1.5, {sprintf('Nhiệt từ CPU (%g°C) lan dọc tấm', ts.T_CPU), ...
        sprintf('hệ số dẫn nhiệt k = %g W/m.K', ts.k)}, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', 'Color', nong, 'FontSize', 10);
    text(ax, 10.3, 0.3, {sprintf('dày %g mm', ts.t*1000), '(vẽ phóng to)', sprintf('đầu xa: %g°C', ts.T_inf)}, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'Color', phu, 'FontSize', 9);
    text(ax, -1.6, -3.6, { ...
        sprintf('Biết:  tấm %s %g x %g cm, dày %g mm, k = %g W/m.K, h = %g W/m^2.K,', ...
                vat_lieu, L, L, ts.t*1000, ts.k, ts.h), ...
        sprintf('          CPU %g°C, phòng %g°C', ts.T_CPU, ts.T_inf), ...
        sprintf('Tìm:  nhiệt độ tại %d x %d = %d điểm lưới trên tấm', N, N, N^2), ...
        sprintf('  =>   hệ %d phương trình  A·θ = b  (A cỡ %d x %d)', N^2, N^2, N^2)}, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'Color', chu, 'FontSize', 10.5, ...
        'BackgroundColor', 'w', 'EdgeColor', mau_ve('truc'), 'Margin', 8);
    axis(ax, 'equal', 'off');
    xlim(ax, [-2 13.2]); ylim(ax, [-7.4 3.2]);
    title(ax, '(b) Nhìn ngang: nhiệt vào từ CPU, thoát ra hai mặt', 'Color', chu, ...
        'FontWeight', 'normal', 'FontSize', 13);
    luu_hinh(fig, 'hinh0_de_bai_matlab.png', true);
end

function mui_ten(ax, x, y, mau, lw, dau, hai_dau)
    % Mũi tên từ (x(1), y(1)) tới (x(2), y(2)); 'dau' là chiều dài đầu mũi (đơn vị dữ liệu).
    plot(ax, x, y, 'Color', mau, 'LineWidth', lw);
    dau_mui_ten(ax, [x(1) y(1)], [x(2) y(2)], mau, dau);
    if hai_dau, dau_mui_ten(ax, [x(2) y(2)], [x(1) y(1)], mau, dau); end
end

function dau_mui_ten(ax, goc, ngon, mau, dau)
    u = (ngon - goc) / norm(ngon - goc);       % hướng mũi tên
    v = [-u(2) u(1)];                          % hướng vuông góc
    day = ngon - dau * u;
    p = [ngon; day + 0.45 * dau * v; day - 0.45 * dau * v];
    patch(ax, p(:, 1), p(:, 2), mau, 'EdgeColor', 'none');
end

function luu_hinh(fig, ten_file, giu_mau_chu)
    % Lưu hình vào thư mục ket_qua (cạnh thư mục matlab).
    % giu_mau_chu = true: giữ màu chữ đã đặt riêng (hình sơ đồ), không đổi hết sang màu đen.
    thu_muc = fullfile(fileparts(mfilename('fullpath')), '..', 'ket_qua');
    if ~exist(thu_muc, 'dir'), mkdir(thu_muc); end
    duong_dan = fullfile(thu_muc, ten_file);
    try
        theme(fig, 'light');   % MATLAB R2025a+ mặc định giao diện tối -> chữ bị nhạt trên nền sáng
    catch
        % bản MATLAB cũ không có hàm theme, vốn đã dùng giao diện sáng
    end
    if nargin < 3 || ~giu_mau_chu
        set(findall(fig, 'Type', 'text'), 'Color', [11 11 11]/255);
    end
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

%% ====================== PHẦN 5. NHẬP THÔNG SỐ ======================
% Người dùng không cần biết số liệu, mỗi câu hỏi nhận nhiều kiểu trả lời:
%   Enter              -> giữ giá trị mặc định
%   a / b / c          -> chọn một lựa chọn có sẵn
%   tên lựa chọn       -> giống chọn chữ cái, gõ có dấu hay không dấu đều được (nhôm, nhom)
%   một con số         -> tự nhập giá trị, được kèm đơn vị (10 cm, 1mm, 80 độ) và dùng dấu phẩy (0,5)
%   ?                  -> xem giải thích
% Giá trị lựa chọn và khoảng hợp lý ghi theo đơn vị hiển thị; doi_so(1) đổi về đơn vị của thong_so().

function ds = cau_hoi()
    do_c = {{'', 'c', '°c', 'do', 'doc'}, [1 1 1 1 1]};
    ds = [ ...
        tao_cau('k', 'Vật liệu', 'Tấm tản nhiệt làm bằng gì?', 'hệ số dẫn nhiệt k', 'W/m.K', ...
            {{'', 'w/m.k', 'w/mk'}, [1 1 1]}, [1 500], ...
            {'Nhôm', 'Đồng', 'Thép'}, [200 400 50], {{'al'}, {'cu'}, {'sat'}}, ...
            ['Vật liệu quyết định nhiệt lan trong tấm nhanh hay chậm. Đồng dẫn nhiệt tốt nhất,' ...
             ' thép kém nhất. Tản nhiệt CPU thường làm bằng nhôm.']), ...
        tao_cau('h', 'Làm mát', 'Tản nhiệt được làm mát thế nào?', 'hệ số toả nhiệt h', 'W/m^2.K', ...
            {{'', 'w/m2.k', 'w/m^2.k', 'w/m².k'}, [1 1 1 1]}, [1 2000], ...
            {'Không quạt', 'Quạt thường', 'Quạt mạnh'}, [10 100 250], {{'khong'}, {'quat', 'thuong'}, {'manh'}}, ...
            ['Quạt thổi càng mạnh thì tấm mất nhiệt ra không khí càng nhanh. Không có quạt,' ...
             ' nhiệt chỉ thoát nhờ không khí tự lưu thông nên rất chậm.']), ...
        tao_cau('T_inf', 'Nhiệt độ phòng', 'Phòng đặt máy nóng hay lạnh?', 'nhiệt độ phòng', '°C', ...
            do_c, [0 50], ...
            {'Phòng máy lạnh', 'Bình thường', 'Trời nóng'}, [20 25 35], {{'lanh', 'may lanh'}, {'thuong'}, {'nong'}}, ...
            'Nhiệt độ không khí quanh máy. Ba cạnh tấm xa CPU được coi là bằng nhiệt độ này.'), ...
        tao_cau('T_CPU', 'Nhiệt độ CPU', 'CPU đang chạy nặng cỡ nào?', 'nhiệt độ CPU', '°C', ...
            do_c, [1 150], ...
            {'Máy rảnh', 'Tải nặng', 'Quá tải'}, [45 80 95], {{'ranh'}, {'nang', 'game'}, {'qua tai'}}, ...
            ['Nhiệt độ của cạnh tấm áp vào CPU. Máy để không khoảng 40-50°C, chơi game hay' ...
             ' tính toán nặng khoảng 80°C, gần 100°C là CPU sắp tự giảm tốc.']), ...
        tao_cau('L', 'Cạnh tấm', 'Tấm tản nhiệt to cỡ nào?', 'cạnh tấm', 'cm', ...
            {{'', 'cm', 'mm', 'm'}, [0.01 0.01 0.001 1]}, [1 50], ...
            {'Nhỏ', 'Vừa', 'Lớn'}, [5 10 15], {{'be'}, {'trung binh'}, {'to'}}, ...
            'Tấm hình vuông, đây là độ dài một cạnh. Tản nhiệt CPU máy bàn thường 8-12 cm.'), ...
        tao_cau('t', 'Độ dày', 'Tấm dày hay mỏng?', 'độ dày', 'mm', ...
            {{'', 'mm', 'cm', 'm'}, [0.001 0.001 0.01 1]}, [0.1 10], ...
            {'Mỏng', 'Vừa', 'Dày'}, [0.5 1 2], {{}, {'trung binh'}, {}}, ...
            'Tấm càng dày thì nhiệt càng dễ lan xa khỏi CPU trước khi thoát ra không khí.'), ...
        tao_cau('N', 'Lưới tính', 'Muốn tính chi tiết đến mức nào?', 'số điểm lưới mỗi chiều N', 'điểm', ...
            {{'', 'diem'}, [1 1]}, [5 40], ...
            {'Thô', 'Vừa', 'Mịn'}, [10 20 40], {{'nhanh'}, {'trung binh'}, {'chi tiet', 'cham'}}, ...
            ['Chương trình chia tấm thành lưới N x N điểm và tính nhiệt độ tại mỗi điểm.' ...
             ' Lưới càng mịn thì kết quả càng sát thực tế nhưng chạy càng lâu (mịn mất cỡ 1 phút).'])];
end

function c = tao_cau(khoa, nhan, hoi, ten_so, don_vi, doi, khoang, ten_lc, so_lc, bi_danh, giai_thich)
    c = struct('khoa', khoa, 'nhan', nhan, 'hoi', hoi, 'ten_so', ten_so, 'don_vi', don_vi, ...
        'doi_ten', {doi{1}}, 'doi_so', doi{2}, 'khoang', khoang, 'ten_lc', {ten_lc}, ...
        'so_lc', so_lc, 'bi_danh', {bi_danh}, 'giai_thich', giai_thich);
end

function ds = tinh_huong()
    ds = struct('mo_ta', {'Máy văn phòng: tản nhôm không quạt, CPU 45°C', ...
                          'Máy chơi game: tản nhôm có quạt, CPU 80°C (mặc định)', ...
                          'Máy trạm: tản đồng quạt mạnh, CPU 95°C, phòng 35°C'}, ...
                'bi_danh', {{'van phong'}, {'game'}, {'may tram'}}, ...
                'thay_doi', {struct('k', 200, 'h', 10, 'T_CPU', 45), struct(), ...
                             struct('k', 400, 'h', 250, 'T_CPU', 95, 'T_inf', 35)});
end

function s = hoi(loi_nhac)
    % Đọc một dòng; khi chạy -batch (không có bàn phím) thì coi như bấm Enter.
    if batchStartupOptionUsed
        s = '';
    else
        s = strtrim(input(loi_nhac, 's'));
    end
end

function s = bo_dau(s)
    % Chữ thường, bỏ dấu tiếng Việt, gộp khoảng trắng: ' Đồng ' -> 'dong'.
    s = lower(strtrim(s));
    nhom = {'aàáạảãâầấậẩẫăằắặẳẵ', 'eèéẹẻẽêềếệểễ', 'iìíịỉĩ', 'oòóọỏõôồốộổỗơờớợởỡ', ...
            'uùúụủũưừứựửữ', 'yỳýỵỷỹ', 'dđĐ'};
    for i = 1:numel(nhom)
        for ch = nhom{i}(2:end)
            s(s == ch) = nhom{i}(1);
        end
    end
    s = regexprep(s, '\s+', ' ');
end

function ten = ten_lua_chon(cau, gia_tri)
    % Tên lựa chọn có sẵn ứng với giá trị (đơn vị thong_so()), '' nếu là giá trị tự nhập.
    i = find(abs(cau.so_lc - gia_tri / cau.doi_so(1)) < 1e-9 * max(1, abs(cau.so_lc)), 1);
    if isempty(i), ten = ''; else, ten = cau.ten_lc{i}; end
end

function gia_tri = hoi_mot_cau(cau, gia_tri)
    he_so = cau.doi_so(1);
    chu = 'abc';
    fprintf('\n  %s\n', cau.hoi);
    for i = 1:3
        fprintf('     %s) %-16s(%s = %g %s)\n', chu(i), cau.ten_lc{i}, cau.ten_so, cau.so_lc(i), cau.don_vi);
    end
    fprintf('     Hoặc gõ thẳng %s theo %s.  Enter = giữ %g %s,  ? = giải thích\n', ...
        cau.ten_so, cau.don_vi, gia_tri / he_so, cau.don_vi);
    while true
        tra_loi = bo_dau(hoi('     > '));
        if isempty(tra_loi), return; end
        if strcmp(tra_loi, '?')
            fprintf('     %s\n', cau.giai_thich);
            continue
        end
        for i = 1:3
            if any(strcmp(tra_loi, [{chu(i), bo_dau(cau.ten_lc{i})}, cau.bi_danh{i}]))
                fprintf('     -> %s: %s = %g %s\n', cau.ten_lc{i}, cau.ten_so, cau.so_lc(i), cau.don_vi);
                gia_tri = cau.so_lc(i) * he_so;
                return
            end
        end
        khop = regexp(tra_loi, '^([+-]?\d+(?:[.,]\d+)?)\s*(.*)$', 'tokens', 'once');
        if isempty(khop)
            fprintf('     Chưa hiểu câu trả lời. Gõ a, b, c, tên lựa chọn, một con số, hoặc ? để xem giải thích.\n');
            continue
        end
        don_vi = strrep(khop{2}, ' ', '');
        if isempty(don_vi), don_vi = ''; end
        j = find(strcmp(don_vi, cau.doi_ten), 1);
        if isempty(j)
            fprintf('     Không nhận ra đơn vị ''%s''. Dùng: %s (hoặc bỏ trống = %s)\n', ...
                khop{2}, strjoin(cau.doi_ten(2:end), ', '), cau.don_vi);
            continue
        end
        so = str2double(strrep(khop{1}, ',', '.')) * cau.doi_so(j) / he_so;   % theo đơn vị hiển thị
        if so < cau.khoang(1) || so > cau.khoang(2)
            fprintf('     %g %s nằm ngoài khoảng hợp lý %g - %g %s. Nhập lại.\n', ...
                so, cau.don_vi, cau.khoang(1), cau.khoang(2), cau.don_vi);
            continue
        end
        if strcmp(cau.khoa, 'N') && so ~= round(so)
            fprintf('     Cần một số nguyên. Nhập lại.\n');
            continue
        end
        fprintf('     -> %s = %g %s\n', cau.ten_so, so, cau.don_vi);
        gia_tri = so * he_so;
        return
    end
end

function ts = hoi_tung_cau(ts)
    ds = cau_hoi();
    for cau = ds
        ts.(cau.khoa) = hoi_mot_cau(cau, ts.(cau.khoa));
        while strcmp(cau.khoa, 'T_CPU') && ts.T_CPU <= ts.T_inf
            fprintf('     CPU phải nóng hơn phòng (%g°C) thì mới có nhiệt cần tản. Chọn lại.\n', ts.T_inf);
            ts.T_CPU = hoi_mot_cau(cau, ts.T_CPU);
        end
    end
end

function ts = chon_tinh_huong(ts)
    ds = tinh_huong();
    chu = 'abc';
    fprintf('\n  Chọn tình huống (gõ a, b, c hoặc tên; Enter = b):\n');
    for i = 1:3, fprintf('     %s) %s\n', chu(i), ds(i).mo_ta); end
    while true
        tra_loi = bo_dau(hoi('     > '));
        if isempty(tra_loi), tra_loi = 'b'; end
        for i = 1:3
            if strcmp(tra_loi, chu(i)) || any(strcmp(tra_loi, ds(i).bi_danh))
                fprintf('     -> %s\n', ds(i).mo_ta);
                for ten = fieldnames(ds(i).thay_doi)'
                    ts.(ten{1}) = ds(i).thay_doi.(ten{1});
                end
                return
            end
        end
        fprintf('     Chưa hiểu câu trả lời. Gõ a, b hoặc c.\n');
    end
end

function in_tom_tat(ts)
    fprintf('\n  Thông số sẽ dùng:\n');
    for cau = cau_hoi()
        ten = ten_lua_chon(cau, ts.(cau.khoa));
        if ~isempty(ten), ten = sprintf('  (%s)', ten); end
        fprintf('     %-16s: %g %s%s\n', cau.nhan, ts.(cau.khoa) / cau.doi_so(1), cau.don_vi, ten);
    end
end

function ts = nhap_thong_so(ts)
    fprintf('Chọn cách nhập thông số:\n');
    fprintf('  1) Chạy ngay với thông số mặc định      (Enter)\n');
    fprintf('  2) Chọn nhanh một tình huống mẫu        (không cần biết số liệu)\n');
    fprintf('  3) Trả lời từng câu hỏi                 (chọn a/b/c, gõ tên, hoặc gõ số kèm đơn vị)\n');
    while true
        cach = hoi('  > ');
        if isempty(cach), cach = '1'; end
        if any(strcmp(cach, {'1', '2', '3'})), break; end
        fprintf('  Gõ 1, 2 hoặc 3 (hoặc Enter để chạy với thông số mặc định).\n');
    end
    if strcmp(cach, '1'), return; end
    while true
        if strcmp(cach, '2'), ts = chon_tinh_huong(ts); else, ts = hoi_tung_cau(ts); end
        in_tom_tat(ts);
        fprintf('\n');
        if ~strcmp(bo_dau(hoi('  Enter = bắt đầu tính,  s = sửa lại từng câu: ')), 's'), return; end
        cach = '3';
    end
end
