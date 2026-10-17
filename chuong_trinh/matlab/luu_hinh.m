function luu_hinh(fig, ten_file)
%LUU_HINH  Lưu hình vào thư mục ket_qua (cạnh thư mục matlab).
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
