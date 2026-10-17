function dinh_dang_truc(ax)
%DINH_DANG_TRUC  Định dạng trục đồ thị đồng bộ với hình bên Python.
set(ax, 'Color', [252 252 251]/255, 'XColor', [82 81 78]/255, 'YColor', [82 81 78]/255, ...
    'GridColor', [225 224 217]/255, 'GridAlpha', 1, 'Box', 'off', 'FontSize', 11, ...
    'XGrid', 'on', 'YGrid', 'on', 'TickDir', 'out');
set(ancestor(ax, 'figure'), 'Color', [252 252 251]/255);
lg = get(ax, 'Legend');
if ~isempty(lg), lg.Box = 'off'; end
end
