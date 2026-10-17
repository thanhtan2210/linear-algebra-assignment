%% BÀI TOÁN 5 - Dự báo thị phần bằng chuỗi Markov
% pi_(k+1) = pi_k * P, P ngẫu nhiên hàng => ||P||_inf = 1.
% Hệ số co Dobrushin delta(P) = 1/2 max_{i,j} ||P_i - P_j||_1 < 1:
%     ||pi_k - pi*||_1 <= 2 delta^k                               (tiên nghiệm)
%     ||pi_k - pi*||_1 <= delta/(1-delta) ||pi_k - pi_(k-1)||_1   (hậu nghiệm)
clc; clear;
P = [0.80 0.10 0.10;
     0.20 0.70 0.10;
     0.10 0.20 0.70];
pi0 = [0.50 0.30 0.20];                 % thị phần hiện tại
EPS = 1e-6;

fprintf('BÀI TOÁN 5: CHUỖI MARKOV - THỊ PHẦN 3 NHÀ MẠNG\n  P =\n'); disp(P);
fprintf('  ||P||_inf = %.4f   (= 1 vì tổng mỗi hàng bằng 1)\n', chuan_vc(P));
fprintf('  ||P^T||_1 = %.4f   (chuẩn 1 của ma trận chuyển vị cũng = 1)\n', chuan_1(P'));
d = 0;
for i = 1:3
    for j = 1:3
        d = max(d, 0.5 * sum(abs(P(i, :) - P(j, :))));
    end
end
fprintf('  Hệ số co delta(P) = %.4f < 1 => hội tụ về trạng thái dừng duy nhất\n', d);
k_tien = ceil(log(EPS / 2) / log(d));
fprintf('  Đánh giá tiên nghiệm: cần k >= %d tháng để sai số <= %g\n', k_tien, EPS);

M = (P - eye(3))';  M(end, :) = 1;      % pi(P - I) = 0 và tổng pi = 1
pi_sao = (M \ [0; 0; 1])';

fprintf('\n  %5s%10s%10s%10s%20s%16s%13s\n', 'Tháng', 'Mạng A', 'Mạng B', 'Mạng C', ...
    '||pi_k-pi_k-1||_1', 'Cận hậu nghiệm', 'Sai số thật');
fprintf('  %5d%10.4f%10.4f%10.4f\n', 0, pi0);
p_k = pi0;
for k = 1:1000
    p_moi = p_k * P;
    buoc = sum(abs(p_moi - p_k));
    can = d / (1 - d) * buoc;
    that = sum(abs(p_moi - pi_sao));
    p_k = p_moi;
    if k <= 6 || can <= EPS
        fprintf('  %5d%10.4f%10.4f%10.4f%20.3e%16.3e%13.3e\n', k, p_k, buoc, can, that);
    elseif k == 7
        fprintf('  %5s\n', '...');
    end
    if can <= EPS, break; end
end
fprintf('\n  Dừng ở tháng k = %d (hậu nghiệm), sớm hơn ước lượng tiên nghiệm (%d tháng).\n', k, k_tien);
fprintf('  Thị phần dài hạn: A: %.2f%%  B: %.2f%%  C: %.2f%%\n', 100 * pi_sao);
