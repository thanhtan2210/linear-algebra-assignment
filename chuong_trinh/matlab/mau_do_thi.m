function c = mau_do_thi(ten)
%MAU_DO_THI  Bảng màu dùng chung cho các đồ thị.
switch ten
    case 'xanh', c = [42 120 214]/255;
    case 'cam',  c = [235 104 52]/255;
    case 'truc', c = [195 194 183]/255;
    otherwise,   c = [82 81 78]/255;
end
end
