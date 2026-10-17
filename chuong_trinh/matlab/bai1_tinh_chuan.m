%% BÀI TOÁN 1 - Chương trình tính chuẩn ma trận từ dữ liệu nhập vào (mục 3.3, 4.1)
% - Nhập ma trận A (cú pháp MATLAB, vd: [1 2 3; 4 5 6]).
% - Tính 4 chuẩn bằng hàm TỰ CÀI ĐẶT rồi đối chiếu với hàm norm có sẵn.
% - Kiểm chứng bằng số các tính chất của chuẩn ma trận.
% Nhấn Enter (không nhập gì) để chạy các ví dụ trong báo cáo.
clc; clear;
fprintf('BÀI TOÁN 1: TÍNH CHUẨN MA TRẬN\n');

ten_vd = {'VD chuẩn 1 (mục 2.3.1)', 'VD chuẩn vô cùng (2.3.2)', ...
          'VD1 chuẩn 2 - đối xứng', 'VD2 chuẩn 2 - không đ.xứng'};
vd = {[2 5 4; 6 4 3; 5 3 5], [-3 5 -2; 1 -7 4; 2 0 9], [5 2; 2 5], [1 3; 0 2]};

A = [];
if ~batchStartupOptionUsed
    A = input('Nhập ma trận A (vd [1 2 3; 4 5 6]), Enter để chạy ví dụ: ');
end
if ~isempty(A)
    xu_ly('Ma trận nhập vào', A);
else
    for i = 1:numel(vd)
        xu_ly(ten_vd{i}, vd{i});
    end
end

function xu_ly(ten, A)
fprintf('%s\n%s:  A =\n', repmat('=', 1, 72), ten);
disp(A);
[n2, so_lan] = chuan_2(A);
ten_chuan = {'Chuẩn 1', 'Chuẩn vô cùng', 'Frobenius (đ.nghĩa)', 'Frobenius (vết)', 'Chuẩn 2 (lũy thừa)'};
tu_tinh = [chuan_1(A), chuan_vc(A), chuan_F(A), chuan_F_vet(A), n2];
co_san  = [norm(A, 1), norm(A, inf), norm(A, 'fro'), norm(A, 'fro'), norm(A, 2)];
fprintf('  %-22s%14s%14s%13s\n', 'Chuẩn', 'Tự cài đặt', 'norm()', 'Chênh lệch');
for i = 1:5
    fprintf('  %-22s%14.6f%14.6f%13.1e\n', ten_chuan{i}, tu_tinh(i), co_san(i), abs(tu_tinh(i) - co_san(i)));
end
fprintf('  (Phương pháp lũy thừa hội tụ sau %d vòng lặp)\n', so_lan);
kiem_tra_tinh_chat(A);
end

function kiem_tra_tinh_chat(A)
rng(1);
[m, n] = size(A);
B = randi([-5 5], m, n);         % cùng cỡ để cộng
C = randi([-5 5], n, n);         % để nhân A*C
x = randn(n, 1);
a = -2.5;
fprintf('  Kiểm chứng tính chất (dùng chuẩn 2; B, C, x sinh ngẫu nhiên):\n');
ten = {'||aA|| = |a|.||A||', '||A+B|| <= ||A||+||B||', '||AC|| <= ||A||.||C||', ...
       '||Ax|| <= ||A||.||x||', '||A||_2 <= ||A||_F', '||A||_F <= sqrt(r)||A||_2', ...
       '||A||_2 <= sqrt(||A||_1||A||_inf)'};
trai = [chuan_2(a*A), chuan_2(A+B), chuan_2(A*C), norm(A*x), chuan_2(A), chuan_F(A), chuan_2(A)];
phai = [abs(a)*chuan_2(A), chuan_2(A)+chuan_2(B), chuan_2(A)*chuan_2(C), chuan_2(A)*norm(x), ...
        chuan_F(A), sqrt(rank(A))*chuan_2(A), sqrt(chuan_1(A)*chuan_vc(A))];
for i = 1:numel(ten)
    if i == 1
        dau = ' ='; dung = abs(trai(i) - phai(i)) <= 1e-9 * phai(i);
    else
        dau = '<='; dung = trai(i) <= phai(i) * (1 + 1e-12);
    end
    kq = 'SAI'; if dung, kq = 'ĐÚNG'; end
    fprintf('    %-36s %10.4f %s %10.4f   %s\n', ten{i}, trai(i), dau, phai(i), kq);
end
end
