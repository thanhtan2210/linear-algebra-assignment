%% BÀI TOÁN 3 - Thuật toán dừng trong phương pháp lặp (mục 3.6)
% Phân bố nhiệt độ trên thanh kim loại tản nhiệt - Jacobi và Gauss-Seidel.
% Thanh nhôm L = 0.5 m, hai đầu T_A = 100°C, T_B = 40°C, môi trường T_inf = 25°C:
%     T'' - m^2 (T - T_inf) = 0,  m^2 = 100
% Sai phân n điểm trong, h = L/(n+1), c = m^2 h^2:
%     -T_(i-1) + (2+c) T_i - T_(i+1) = c.T_inf
% Lặp x^(k) = B x^(k-1) + g hội tụ nếu q = ||B|| < 1:
%     tiên nghiệm: ||x^(k) - x*|| <= q^k/(1-q) ||x^(1) - x^(0)||
%     hậu nghiệm : ||x^(k) - x*|| <= q/(1-q)  ||x^(k) - x^(k-1)||   <= eps  -> DỪNG
clc; clear;
L = 0.5; n = 9;
T_A = 100; T_B = 40; T_inf = 25; m2 = 100;
EPS = 1e-3;

h = L / (n + 1);
c = m2 * h^2;
A = (2 + c) * eye(n) - diag(ones(n-1, 1), 1) - diag(ones(n-1, 1), -1);
b = c * T_inf * ones(n, 1);
b(1) = b(1) + T_A;  b(end) = b(end) + T_B;
x_dung = A \ b;                          % nghiệm đúng của hệ rời rạc

fprintf('BÀI TOÁN 3: n = %d nút trong, h = %g m, c = m^2h^2 = %.4f, eps = %g °C\n', n, h, c, EPS);
fprintf('Ma trận A trội đường chéo chặt: |a_ii| = %.2f > tổng |a_ij| = 2\n', 2 + c);

kq_j  = giai_lap(A, b, x_dung, T_inf, EPS, 'Jacobi');
kq_gs = giai_lap(A, b, x_dung, T_inf, EPS, 'Gauss-Seidel');
in_ket_qua('Jacobi', kq_j, EPS);
in_ket_qua('Gauss-Seidel', kq_gs, EPS);
fprintf('%s\n  Nhiệt độ tại các nút (°C): %s\n', repmat('=', 1, 78), mat2str(round(kq_gs.x', 3)));

%% Hình 1: sai số thật theo vòng lặp
fig = figure('Position', [100 100 750 450]);
semilogy(kq_j.ls(:, 1), kq_j.ls(:, 4), '-', 'Color', mau_do_thi('xanh'), 'LineWidth', 2); hold on;
semilogy(kq_gs.ls(:, 1), kq_gs.ls(:, 4), '-', 'Color', mau_do_thi('cam'), 'LineWidth', 2);
semilogy(kq_j.ls(end, 1), kq_j.ls(end, 4), 'o', 'Color', mau_do_thi('xanh'), ...
    'MarkerFaceColor', mau_do_thi('xanh'), 'HandleVisibility', 'off');
semilogy(kq_gs.ls(end, 1), kq_gs.ls(end, 4), 'o', 'Color', mau_do_thi('cam'), ...
    'MarkerFaceColor', mau_do_thi('cam'), 'HandleVisibility', 'off');
yline(EPS, '--', sprintf('\\epsilon = %g °C', EPS), 'Color', mau_do_thi('truc'), 'HandleVisibility', 'off');
text(kq_j.ls(end, 1), kq_j.ls(end, 4) * 1.8, sprintf('dừng k = %d', kq_j.ls(end, 1)));
text(kq_gs.ls(end, 1), kq_gs.ls(end, 4) * 1.8, sprintf('dừng k = %d', kq_gs.ls(end, 1)));
xlabel('Số vòng lặp k'); ylabel('||x^{(k)} - x*||_\infty  (°C, thang log)');
title('Sai số theo vòng lặp - dừng bằng đánh giá hậu nghiệm');
legend(sprintf('Jacobi (q = %.3f)', kq_j.q), sprintf('Gauss-Seidel (q = %.3f)', kq_gs.q));
dinh_dang_truc(gca);
luu_hinh(fig, 'bai3_hoi_tu_matlab.png');

%% Hình 2: phân bố nhiệt độ
m = sqrt(m2);
T_gt = @(xx) T_inf + ((T_A - T_inf) * sinh(m * (L - xx)) + (T_B - T_inf) * sinh(m * xx)) / sinh(m * L);
xs = linspace(0, L, 200);
x_nut = linspace(0, L, n + 2);
T_nut = [T_A; kq_gs.x; T_B];
fig = figure('Position', [100 100 750 450]);
plot(xs, T_gt(xs), '-', 'Color', mau_do_thi('xanh'), 'LineWidth', 2); hold on;
plot(x_nut, T_nut, 'o', 'Color', mau_do_thi('cam'), 'MarkerFaceColor', mau_do_thi('cam'), 'MarkerSize', 8);
xlabel('Vị trí trên thanh x (m)'); ylabel('Nhiệt độ T (°C)');
title('Phân bố nhiệt độ trên thanh tản nhiệt');
legend('Nghiệm giải tích', 'Nghiệm lặp Gauss-Seidel (các nút)');
dinh_dang_truc(gca);
luu_hinh(fig, 'bai3_nhiet_do_matlab.png');
fprintf('  Sai số rời rạc hoá tại các nút: %.3f °C\n', max(abs(T_nut - T_gt(x_nut)')));

%% ---------------- Hàm cục bộ ----------------
function q = he_so_co(A, pp)
% q của Jacobi = ||B_J||_inf ; q của Gauss-Seidel theo công thức hàng
D = diag(A);
B = -A ./ D;
B(logical(eye(size(A)))) = 0;
if strcmp(pp, 'Jacobi')
    q = chuan_vc(B);
else
    q = 0;
    for i = 1:size(A, 1)
        p = sum(abs(B(i, 1:i-1)));       % phần dưới đường chéo
        r = sum(abs(B(i, i+1:end)));     % phần trên đường chéo
        q = max(q, r / (1 - p));
    end
end
end

function kq = giai_lap(A, b, x_dung, x0, EPS, pp)
nA = numel(b);
q = he_so_co(A, pp);
x = x0 * ones(nA, 1);                    % x^(0): cả thanh bằng nhiệt độ môi trường
ls = [];  k_ngay_tho = [];
for k = 1:1000
    x_moi = x;
    for i = 1:nA
        if strcmp(pp, 'Gauss-Seidel'), nguon = x_moi; else, nguon = x; end   % GS dùng ngay giá trị mới
        tong = A(i, :) * nguon - A(i, i) * nguon(i);
        x_moi(i) = (b(i) - tong) / A(i, i);
    end
    buoc = norm(x_moi - x, inf);
    if k == 1
        k_tien = ceil(log(EPS * (1 - q) / buoc) / log(q));
    end
    hau_nghiem = q / (1 - q) * buoc;
    that = norm(x_moi - x_dung, inf);
    ls(end+1, :) = [k, buoc, hau_nghiem, that]; %#ok<AGROW>
    x = x_moi;
    if isempty(k_ngay_tho) && buoc <= EPS, k_ngay_tho = k; end
    if hau_nghiem <= EPS, break; end
end
kq = struct('q', q, 'x', x, 'ls', ls, 'k_tien', k_tien, 'k_ngay_tho', k_ngay_tho);
end

function in_ket_qua(ten, kq, EPS)
ls = kq.ls;
N = size(ls, 1);
fprintf('%s\nPHƯƠNG PHÁP %s   (q = %.4f < 1 => hội tụ)\n', repmat('=', 1, 78), upper(ten), kq.q);
fprintf('  %4s%18s%17s%15s\n', 'k', '||x_k - x_k-1||', 'q/(1-q)*||.||', 'Sai số thật');
for i = unique([1:min(5, N), max(1, N-2):N])
    if i == N-2 && N > 8, fprintf('  %4s\n', '...'); end
    fprintf('  %4d%18.6e%17.6e%15.6e\n', ls(i, :));
end
kn = kq.k_ngay_tho;
fprintf('  * Dừng theo đánh giá hậu nghiệm sau k = %d vòng; sai số thật = %.2e <= %g\n', N, ls(end, 4), EPS);
fprintf('  * Ước lượng tiên nghiệm (trước khi lặp) cần k >= %d vòng\n', kq.k_tien);
if ls(kn, 4) > EPS, nx = '> eps  => KHÔNG đạt yêu cầu!'; else, nx = '<= eps'; end
fprintf('  * Nếu dừng "ngây thơ" (||x_k - x_k-1|| <= eps): dừng ở k = %d, sai số thật = %.2e %s\n', ...
    kn, ls(kn, 4), nx);
end
