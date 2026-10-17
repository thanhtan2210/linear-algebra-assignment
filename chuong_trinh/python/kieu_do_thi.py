"""Cấu hình chung cho đồ thị (màu, font, lưới) để các hình trong báo cáo đồng bộ."""
import os
import matplotlib

matplotlib.use("Agg")                  # lưu file, không cần cửa sổ
import matplotlib.pyplot as plt

MAU = {
    "xanh": "#2a78d6",     # chuỗi số liệu 1
    "cam": "#eb6834",      # chuỗi số liệu 2
    "ngoc": "#1baf7a",     # chuỗi số liệu 3
    "chu": "#0b0b0b",
    "chu_phu": "#52514e",
    "luoi": "#e1e0d9",
    "truc": "#c3c2b7",
    "nen": "#fcfcfb",
}

THU_MUC_KET_QUA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ket_qua")


def cai_dat():
    plt.rcParams.update({
        "font.family": ["Segoe UI", "DejaVu Sans"],   # cả hai đều có dấu tiếng Việt
        "font.size": 11,
        "axes.facecolor": MAU["nen"],
        "figure.facecolor": MAU["nen"],
        "axes.edgecolor": MAU["truc"],
        "axes.labelcolor": MAU["chu_phu"],
        "axes.titlecolor": MAU["chu"],
        "axes.titlesize": 13,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": MAU["luoi"],
        "grid.linewidth": 0.8,
        "xtick.color": MAU["chu_phu"],
        "ytick.color": MAU["chu_phu"],
        "lines.linewidth": 2,
        "lines.markersize": 6,
        "legend.frameon": False,
        "savefig.dpi": 150,
        "savefig.bbox": "tight",
    })
    os.makedirs(THU_MUC_KET_QUA, exist_ok=True)


def luu_hinh(fig, ten_file):
    duong_dan = os.path.normpath(os.path.join(THU_MUC_KET_QUA, ten_file))
    fig.savefig(duong_dan)
    plt.close(fig)
    print(f"  -> Đã lưu hình: {duong_dan}")
