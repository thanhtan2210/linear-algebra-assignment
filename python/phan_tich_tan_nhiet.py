"""
BÀI TOÁN: PHÂN TÍCH NHIỆT ĐỘ TRÊN TẤM TẢN NHIỆT CPU BẰNG CHUẨN MA TRẬN
(Toàn bộ chương trình trong một file - chỉ cần chạy file này.)

  Bước 1 (mục 3.3, 4.1)     : Lập hệ A.theta = b, tính 4 chuẩn của A bằng code tự cài đặt.
  Bước 2 (mục 3.4, 3.5, 4.2): Sai số cảm biến lan sang nghiệm thế nào? Số điều kiện và lưới.
  Bước 3 (mục 3.6)          : Giải bằng Jacobi / Gauss-Seidel, dừng bằng đánh giá theo chuẩn.
  Bước 4 (mục 4.3)          : So sánh với "ảnh camera nhiệt", nén và khử nhiễu bằng SVD.

Cấu trúc file:
  Phần 1. Các chuẩn ma trận tự cài đặt
  Phần 2. Mô hình tấm tản nhiệt
  Phần 3. Cấu hình đồ thị
  Phần 4. Chương trình chính (4 bước + vẽ hình)

Chạy:  python phan_tich_tan_nhiet.py      (Enter để dùng thông số mặc định; có thể chọn
       tình huống mẫu hoặc trả lời từng câu hỏi bằng a/b/c, tên, hay con số kèm đơn vị)
Cần:   pip install numpy matplotlib
"""
import os
import re
import sys
import unicodedata
import numpy as np
import matplotlib

matplotlib.use("Agg")                  # lưu hình ra file, không cần cửa sổ
import matplotlib.pyplot as plt

# In tiếng Việt có dấu trên cửa sổ dòng lệnh Windows
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


# ======================================================================================
# PHẦN 1. CÁC CHUẨN MA TRẬN TỰ CÀI ĐẶT
# ======================================================================================
# Tự cài đặt theo đúng công thức lý thuyết (không dùng hàm norm có sẵn):
#     ||A||_1   = max_j  sum_i |a_ij|          (chuẩn cột)
#     ||A||_inf = max_i  sum_j |a_ij|          (chuẩn hàng)
#     ||A||_F   = sqrt( sum_ij |a_ij|^2 ) = sqrt( tr(A^T A) )
#     ||A||_2   = sqrt( lambda_max(A^T A) ) = sigma_max(A)


def chuan_1(A):
    """Chuẩn 1: tổng trị tuyệt đối lớn nhất theo CỘT."""
    A = np.asarray(A, dtype=float)
    return np.abs(A).sum(axis=0).max()


def chuan_vo_cung(A):
    """Chuẩn vô cùng: tổng trị tuyệt đối lớn nhất theo HÀNG."""
    A = np.asarray(A, dtype=float)
    return np.abs(A).sum(axis=1).max()


def chuan_frobenius(A):
    """Chuẩn Frobenius theo công thức định nghĩa."""
    A = np.asarray(A, dtype=float)
    return np.sqrt((np.abs(A) ** 2).sum())


def chuan_frobenius_vet(A):
    """Chuẩn Frobenius qua vết: sqrt(tr(A^T A))."""
    A = np.asarray(A, dtype=float)
    return np.sqrt(np.trace(A.T @ A))


def chuan_2(A, tol=1e-12, max_lap=10000, tra_ve_so_lan_lap=False):
    """
    Chuẩn 2 (chuẩn phổ) = sqrt(lambda_max(A^T A)).

    lambda_max được tìm bằng PHƯƠNG PHÁP LŨY THỪA (power method):
        y_k = M x_{k-1},  x_k = y_k / ||y_k||,  lambda_k = x_k^T M x_k
    Điều kiện dừng: |lambda_k - lambda_{k-1}| <= tol * lambda_k.
    """
    A = np.asarray(A, dtype=float)
    M = A.T @ A                       # ma trận đối xứng, xác định không âm
    n = M.shape[0]
    # Vector đầu ngẫu nhiên (cố định seed): vector (1,...,1) có thể vuông góc với
    # vector riêng cần tìm (xảy ra với ma trận lưới đối xứng) và cho kết quả sai.
    x = np.random.default_rng(0).standard_normal(n)
    x /= np.linalg.norm(x)
    lam_cu = 0.0
    for k in range(1, max_lap + 1):
        y = M @ x
        ny = np.sqrt(y @ y)
        if ny == 0:                   # A là ma trận không
            return (0.0, k) if tra_ve_so_lan_lap else 0.0
        x = y / ny
        lam = x @ M @ x
        if abs(lam - lam_cu) <= tol * abs(lam):
            break
        lam_cu = lam
    gia_tri = np.sqrt(lam)
    return (gia_tri, k) if tra_ve_so_lan_lap else gia_tri


CAC_CHUAN = {
    "1": chuan_1,
    "2": chuan_2,
    "inf": chuan_vo_cung,
    "fro": chuan_frobenius,
}


def so_dieu_kien(A, p="2"):
    """Số điều kiện kappa_p(A) = ||A||_p * ||A^{-1}||_p (A vuông khả nghịch)."""
    A = np.asarray(A, dtype=float)
    chuan = CAC_CHUAN[p]
    return chuan(A) * chuan(np.linalg.inv(A))


# ======================================================================================
# PHẦN 2. MÔ HÌNH TẤM TẢN NHIỆT CPU
# ======================================================================================
# Tấm nhôm vuông cạnh L, dày t. Cạnh dưới (y = 0) áp vào CPU ở nhiệt độ T_CPU,
# ba cạnh còn lại gắn vào khung ở nhiệt độ phòng T_inf, hai mặt tấm toả nhiệt
# ra không khí (quạt thổi, hệ số toả nhiệt h). Đặt theta = T - T_inf:
#
#     d²theta/dx² + d²theta/dy² - m² theta = 0,     m² = 2h / (k t)
#
# Sai phân 5 điểm trên lưới N x N điểm trong, bước dx = L/(N+1), c = m² dx²:
#
#     (4 + c) theta_ij - theta_(i-1)j - theta_(i+1)j - theta_i(j-1) - theta_i(j+1) = 0
#
# => hệ A.theta = b với A cỡ N² x N²,  A = kron(I, K) + kron(K, I) + c.I,
#    K = tridiag(-1, 2, -1).

# Thông số mặc định (có thể nhập lại khi chạy chương trình)
THONG_SO = {
    "L": 0.10,        # cạnh tấm (m)
    "t": 0.001,       # độ dày tấm (m)
    "k": 200.0,       # hệ số dẫn nhiệt của nhôm (W/m.K)
    "h": 100.0,       # hệ số toả nhiệt khi có quạt (W/m².K)
    "T_CPU": 80.0,    # nhiệt độ cạnh áp vào CPU (°C)
    "T_inf": 25.0,    # nhiệt độ phòng (°C)
    "N": 20,          # số điểm lưới trong mỗi chiều
}


def he_so_m2(ts):
    return 2 * ts["h"] / (ts["k"] * ts["t"])


def lap_he(ts, N=None):
    """Lập hệ A.theta = b. Ẩn được đánh số theo hàng: chỉ số = j*N + i (i theo x, j theo y)."""
    N = ts["N"] if N is None else N
    dx = ts["L"] / (N + 1)
    c = he_so_m2(ts) * dx ** 2
    K = 2 * np.eye(N) - np.eye(N, k=1) - np.eye(N, k=-1)
    I = np.eye(N)
    A = np.kron(I, K) + np.kron(K, I) + c * np.eye(N * N)
    b = np.zeros(N * N)
    b[:N] += ts["T_CPU"] - ts["T_inf"]        # hàng j = 0 kề cạnh CPU (theta = T_CPU - T_inf)
    return A, b, c                             # ba cạnh còn lại có theta = 0 nên không đóng góp vào b


def tri_rieng_A(ts, N):
    """Trị riêng nhỏ nhất / lớn nhất của A theo công thức giải tích (để kiểm chứng)."""
    c = he_so_m2(ts) * (ts["L"] / (N + 1)) ** 2
    s = lambda p: np.sin(p * np.pi / (2 * (N + 1))) ** 2
    return c + 8 * s(1), c + 8 * s(N)


def nghiem_giai_tich(ts, x, y, so_hang=400):
    """
    Nghiệm chính xác theta(x, y) bằng chuỗi Fourier:
      theta = sum_n a_n sin(n pi x / L) sinh(beta_n (L - y)) / sinh(beta_n L),
      a_n = 2 theta0 (1 - (-1)^n) / (n pi),   beta_n = sqrt((n pi / L)² + m²)
    """
    L, theta0, m2 = ts["L"], ts["T_CPU"] - ts["T_inf"], he_so_m2(ts)
    x, y = np.asarray(x, float), np.asarray(y, float)
    kq = np.zeros(np.broadcast(x, y).shape)
    for n in range(1, so_hang + 1, 2):        # chỉ n lẻ có a_n khác 0
        beta = np.sqrt((n * np.pi / L) ** 2 + m2)
        # sinh(beta(L-y))/sinh(beta L) viết dạng mũ để tránh tràn số
        ti_so = np.exp(-beta * y) * (1 - np.exp(-2 * beta * (L - y))) / (1 - np.exp(-2 * beta * L))
        kq += 4 * theta0 / (n * np.pi) * np.sin(n * np.pi * x / L) * ti_so
    return kq


def ma_tran_lap(A, phuong_phap):
    """Ma trận lặp B = -D^-1 (L + U) và hệ số co q tương ứng (chuẩn vô cùng)."""
    D = np.diag(A)
    B = -A / D[:, None]
    np.fill_diagonal(B, 0)
    if phuong_phap == "Jacobi":
        q = np.abs(B).sum(axis=1).max()                                   # q = ||B||_inf
    else:                                                                 # Gauss-Seidel
        n = len(D)
        duoi = np.array([np.abs(B[i, :i]).sum() for i in range(n)])
        tren = np.array([np.abs(B[i, i + 1:]).sum() for i in range(n)])
        q = (tren / (1 - duoi)).max()
    return B, q


def giai_lap(A, b, phuong_phap, eps, x_dung=None, x0=None, max_lap=20000):
    """
    Giải A x = b bằng Jacobi hoặc Gauss-Seidel.
    Dừng khi đánh giá HẬU NGHIỆM  q/(1-q) ||x_k - x_(k-1)||_inf <= eps.
    Trả về nghiệm, q, số vòng tiên nghiệm, và lịch sử từng vòng:
      (k, ||x_k - x_(k-1)||, đánh giá hậu nghiệm, sai số thật nếu biết x_dung).
    """
    n = len(b)
    _, q = ma_tran_lap(A, phuong_phap)
    D = np.diag(A)
    R = A - np.diag(D)                       # phần ngoài đường chéo
    x = np.zeros(n) if x0 is None else x0.astype(float).copy()
    lich_su, k_tien = [], None
    for k in range(1, max_lap + 1):
        if phuong_phap == "Jacobi":
            x_moi = (b - R @ x) / D
        else:
            x_moi = x.copy()
            for i in range(n):               # dùng ngay các thành phần vừa tính
                x_moi[i] = (b[i] - R[i] @ x_moi) / D[i]
        buoc = np.abs(x_moi - x).max()
        if k == 1:
            k_tien = int(np.ceil(np.log(eps * (1 - q) / buoc) / np.log(q)))
        hau_nghiem = q / (1 - q) * buoc
        that = np.abs(x_moi - x_dung).max() if x_dung is not None else np.nan
        lich_su.append((k, buoc, hau_nghiem, that))
        x = x_moi
        if hau_nghiem <= eps:
            break
    return x, q, k_tien, np.array(lich_su)


# ======================================================================================
# PHẦN 3. CẤU HÌNH ĐỒ THỊ (màu, font, lưới, lưu hình)
# ======================================================================================
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

# Thang màu tuần tự một sắc độ (nhạt = nguội, đậm = nóng) cho bản đồ nhiệt
BANG_MAU_NHIET = matplotlib.colors.LinearSegmentedColormap.from_list(
    "nhiet", ["#fff5eb", "#fd8d3c", "#7f2704"])

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


# ======================================================================================
# PHẦN 4. CHƯƠNG TRÌNH CHÍNH
# ======================================================================================
EPS = 1e-3                  # sai số cho phép khi giải lặp (°C)
SAI_SO_CAM_BIEN = 0.5       # sai số mỗi điểm đo nhiệt độ dọc cạnh CPU (°C)
NHIEU_CAMERA = 2.0          # độ lệch chuẩn nhiễu của camera nhiệt (°C)
DO_PHAN_GIAI_CAMERA = 128   # ảnh camera 128 x 128 điểm ảnh
rng = np.random.default_rng(2026)


def tieu_de(chu):
    print("\n" + "=" * 90 + f"\n{chu}\n" + "=" * 90)


# ---------------------------------------------------------------------------------------
# Nhập thông số. Người dùng không cần biết số liệu, mỗi câu hỏi nhận nhiều kiểu trả lời:
#   Enter              -> giữ giá trị mặc định
#   a / b / c          -> chọn một lựa chọn có sẵn
#   tên lựa chọn       -> giống chọn chữ cái, gõ có dấu hay không dấu đều được (nhôm, nhom)
#   một con số         -> tự nhập giá trị, được kèm đơn vị (10 cm, 1mm, 80 độ) và dùng dấu phẩy (0,5)
#   ?                  -> xem giải thích
# Giá trị lựa chọn và khoảng hợp lý ghi theo đơn vị hiển thị; doi[""] đổi về đơn vị trong THONG_SO.
# ---------------------------------------------------------------------------------------
DO_C = {"": 1, "c": 1, "°c": 1, "do": 1, "doc": 1}
CAU_HOI = [
    {"khoa": "k", "nhan": "Vật liệu", "hoi": "Tấm tản nhiệt làm bằng gì?",
     "ten_so": "hệ số dẫn nhiệt k", "don_vi": "W/m.K", "doi": {"": 1, "w/m.k": 1, "w/mk": 1},
     "khoang": (1, 500),
     "lua_chon": [("Nhôm", 200, ["al"]), ("Đồng", 400, ["cu"]), ("Thép", 50, ["sat"])],
     "giai_thich": "Vật liệu quyết định nhiệt lan trong tấm nhanh hay chậm. Đồng dẫn nhiệt tốt nhất,"
                   " thép kém nhất. Tản nhiệt CPU thường làm bằng nhôm."},
    {"khoa": "h", "nhan": "Làm mát", "hoi": "Tản nhiệt được làm mát thế nào?",
     "ten_so": "hệ số toả nhiệt h", "don_vi": "W/m².K", "doi": {"": 1, "w/m2.k": 1, "w/m².k": 1},
     "khoang": (1, 2000),
     "lua_chon": [("Không quạt", 10, ["khong"]), ("Quạt thường", 100, ["quat", "thuong"]),
                  ("Quạt mạnh", 250, ["manh"])],
     "giai_thich": "Quạt thổi càng mạnh thì tấm mất nhiệt ra không khí càng nhanh. Không có quạt,"
                   " nhiệt chỉ thoát nhờ không khí tự lưu thông nên rất chậm."},
    {"khoa": "T_inf", "nhan": "Nhiệt độ phòng", "hoi": "Phòng đặt máy nóng hay lạnh?",
     "ten_so": "nhiệt độ phòng", "don_vi": "°C", "doi": DO_C, "khoang": (0, 50),
     "lua_chon": [("Phòng máy lạnh", 20, ["lanh", "may lanh"]), ("Bình thường", 25, ["thuong"]),
                  ("Trời nóng", 35, ["nong"])],
     "giai_thich": "Nhiệt độ không khí quanh máy. Ba cạnh tấm xa CPU được coi là bằng nhiệt độ này."},
    {"khoa": "T_CPU", "nhan": "Nhiệt độ CPU", "hoi": "CPU đang chạy nặng cỡ nào?",
     "ten_so": "nhiệt độ CPU", "don_vi": "°C", "doi": DO_C, "khoang": (1, 150),
     "lua_chon": [("Máy rảnh", 45, ["ranh"]), ("Tải nặng", 80, ["nang", "game"]),
                  ("Quá tải", 95, ["qua tai"])],
     "giai_thich": "Nhiệt độ của cạnh tấm áp vào CPU. Máy để không khoảng 40-50°C, chơi game hay"
                   " tính toán nặng khoảng 80°C, gần 100°C là CPU sắp tự giảm tốc."},
    {"khoa": "L", "nhan": "Cạnh tấm", "hoi": "Tấm tản nhiệt to cỡ nào?",
     "ten_so": "cạnh tấm", "don_vi": "cm", "doi": {"": 0.01, "cm": 0.01, "mm": 0.001, "m": 1},
     "khoang": (1, 50),
     "lua_chon": [("Nhỏ", 5, ["be"]), ("Vừa", 10, ["trung binh"]), ("Lớn", 15, ["to"])],
     "giai_thich": "Tấm hình vuông, đây là độ dài một cạnh. Tản nhiệt CPU máy bàn thường 8-12 cm."},
    {"khoa": "t", "nhan": "Độ dày", "hoi": "Tấm dày hay mỏng?",
     "ten_so": "độ dày", "don_vi": "mm", "doi": {"": 0.001, "mm": 0.001, "cm": 0.01, "m": 1},
     "khoang": (0.1, 10),
     "lua_chon": [("Mỏng", 0.5, []), ("Vừa", 1, ["trung binh"]), ("Dày", 2, [])],
     "giai_thich": "Tấm càng dày thì nhiệt càng dễ lan xa khỏi CPU trước khi thoát ra không khí."},
    {"khoa": "N", "nhan": "Lưới tính", "hoi": "Muốn tính chi tiết đến mức nào?",
     "ten_so": "số điểm lưới mỗi chiều N", "don_vi": "điểm", "doi": {"": 1, "diem": 1},
     "khoang": (5, 40), "so_nguyen": True,
     "lua_chon": [("Thô", 10, ["nhanh"]), ("Vừa", 20, ["trung binh"]), ("Mịn", 40, ["chi tiet", "cham"])],
     "giai_thich": "Chương trình chia tấm thành lưới N x N điểm và tính nhiệt độ tại mỗi điểm."
                   " Lưới càng mịn thì kết quả càng sát thực tế nhưng chạy càng lâu (mịn mất cỡ 1 phút)."},
]

TINH_HUONG = [
    ("Máy văn phòng: tản nhôm không quạt, CPU 45°C", ["van phong"], {"k": 200, "h": 10, "T_CPU": 45}),
    ("Máy chơi game: tản nhôm có quạt, CPU 80°C (mặc định)", ["game"], {}),
    ("Máy trạm: tản đồng quạt mạnh, CPU 95°C, phòng 35°C", ["may tram"],
     {"k": 400, "h": 250, "T_CPU": 95, "T_inf": 35}),
]


def hoi(loi_nhac):
    """Đọc một dòng; khi không có bàn phím (chạy tự động) thì coi như bấm Enter."""
    try:
        return input(loi_nhac).strip()
    except EOFError:
        return ""


def bo_dau(chuoi):
    """Chữ thường, bỏ dấu tiếng Việt, gộp khoảng trắng: ' Đồng ' -> 'dong'."""
    chuoi = unicodedata.normalize("NFD", chuoi.lower()).replace("đ", "d")
    chuoi = "".join(ch for ch in chuoi if unicodedata.category(ch) != "Mn")
    return " ".join(chuoi.split())


def ten_lua_chon(cau, gia_tri):
    """Tên lựa chọn có sẵn ứng với giá trị (đơn vị THONG_SO), None nếu là giá trị tự nhập."""
    so = gia_tri / cau["doi"][""]
    return next((ten for ten, v, _ in cau["lua_chon"] if np.isclose(v, so)), None)


def hoi_mot_cau(cau, gia_tri):
    he_so = cau["doi"][""]
    print(f"\n  {cau['hoi']}")
    for chu, (ten, so, _) in zip("abc", cau["lua_chon"]):
        print(f"     {chu}) {ten:<16}({cau['ten_so']} = {so:g} {cau['don_vi']})")
    print(f"     Hoặc gõ thẳng {cau['ten_so']} theo {cau['don_vi']}."
          f"  Enter = giữ {gia_tri / he_so:g} {cau['don_vi']},  ? = giải thích")
    while True:
        tra_loi = bo_dau(hoi("     > "))
        if tra_loi == "":
            return gia_tri
        if tra_loi == "?":
            print(f"     {cau['giai_thich']}")
            continue
        for chu, (ten, so, bi_danh) in zip("abc", cau["lua_chon"]):
            if tra_loi in [chu, bo_dau(ten)] + bi_danh:
                print(f"     -> {ten}: {cau['ten_so']} = {so:g} {cau['don_vi']}")
                return so * he_so
        khop = re.fullmatch(r"([+-]?\d+(?:[.,]\d+)?)\s*(.*)", tra_loi)
        if khop is None:
            print("     Chưa hiểu câu trả lời. Gõ a, b, c, tên lựa chọn, một con số, hoặc ? để xem giải thích.")
            continue
        don_vi = khop.group(2).replace(" ", "")
        if don_vi not in cau["doi"]:
            print(f"     Không nhận ra đơn vị '{khop.group(2)}'. Dùng: "
                  + ", ".join(dv for dv in cau["doi"] if dv) + " (hoặc bỏ trống = " + cau["don_vi"] + ")")
            continue
        so = float(khop.group(1).replace(",", ".")) * cau["doi"][don_vi] / he_so   # theo đơn vị hiển thị
        thap, cao = cau["khoang"]
        if not thap <= so <= cao:
            print(f"     {so:g} {cau['don_vi']} nằm ngoài khoảng hợp lý {thap:g} - {cao:g} {cau['don_vi']}."
                  " Nhập lại.")
            continue
        if cau.get("so_nguyen"):
            if so != int(so):
                print("     Cần một số nguyên. Nhập lại.")
                continue
            so = int(so)
        print(f"     -> {cau['ten_so']} = {so:g} {cau['don_vi']}")
        return so * he_so


def hoi_tung_cau(ts):
    for cau in CAU_HOI:
        ts[cau["khoa"]] = hoi_mot_cau(cau, ts[cau["khoa"]])
        while cau["khoa"] == "T_CPU" and ts["T_CPU"] <= ts["T_inf"]:
            print(f"     CPU phải nóng hơn phòng ({ts['T_inf']:g}°C) thì mới có nhiệt cần tản. Chọn lại.")
            ts["T_CPU"] = hoi_mot_cau(cau, ts["T_CPU"])
    return ts


def chon_tinh_huong(ts):
    print("\n  Chọn tình huống (gõ a, b, c hoặc tên; Enter = b):")
    for chu, (mo_ta, _, _) in zip("abc", TINH_HUONG):
        print(f"     {chu}) {mo_ta}")
    while True:
        tra_loi = bo_dau(hoi("     > ")) or "b"
        for chu, (mo_ta, bi_danh, thay_doi) in zip("abc", TINH_HUONG):
            if tra_loi == chu or tra_loi in bi_danh:
                print(f"     -> {mo_ta}")
                ts.update(thay_doi)
                return ts
        print("     Chưa hiểu câu trả lời. Gõ a, b hoặc c.")


def in_tom_tat(ts):
    print("\n  Thông số sẽ dùng:")
    for cau in CAU_HOI:
        so = ts[cau["khoa"]] / cau["doi"][""]
        ten = ten_lua_chon(cau, ts[cau["khoa"]])
        print(f"     {cau['nhan']:<16}: {so:g} {cau['don_vi']}" + (f"  ({ten})" if ten else ""))


def nhap_thong_so():
    ts = dict(THONG_SO)
    print("Chọn cách nhập thông số:")
    print("  1) Chạy ngay với thông số mặc định      (Enter)")
    print("  2) Chọn nhanh một tình huống mẫu        (không cần biết số liệu)")
    print("  3) Trả lời từng câu hỏi                 (chọn a/b/c, gõ tên, hoặc gõ số kèm đơn vị)")
    while True:
        cach = hoi("  > ") or "1"
        if cach in ("1", "2", "3"):
            break
        print("  Gõ 1, 2 hoặc 3 (hoặc Enter để chạy với thông số mặc định).")
    if cach == "1":
        return ts
    while True:
        ts = chon_tinh_huong(ts) if cach == "2" else hoi_tung_cau(ts)
        in_tom_tat(ts)
        if bo_dau(hoi("\n  Enter = bắt đầu tính,  s = sửa lại từng câu: ")) != "s":
            return ts
        cach = "3"


def toa_do_nut(ts, N):
    """Toạ độ các điểm lưới trong, xếp đúng thứ tự ẩn (j*N + i)."""
    g = np.arange(1, N + 1) * ts["L"] / (N + 1)
    X, Y = np.meshgrid(g, g)                   # hàng = y, cột = x
    return X, Y


# =====================================================================================
def buoc_1(ts):
    tieu_de("BƯỚC 1 - LẬP MÔ HÌNH VÀ TÍNH CHUẨN CỦA MA TRẬN A  (mục 3.3, 4.1)")
    N = ts["N"]
    A, b, c = lap_he(ts)
    vat_lieu = ten_lua_chon(CAU_HOI[0], ts["k"])
    print(f"  Tấm {vat_lieu.lower() if vat_lieu else 'vật liệu'} {ts['L'] * 100:g} x {ts['L'] * 100:g} cm,"
          f" dày {ts['t'] * 1000:g} mm, k = {ts['k']:g} W/m.K, T_CPU = {ts['T_CPU']:g}°C,"
          f" T_phòng = {ts['T_inf']:g}°C, h = {ts['h']:g} W/m².K")
    print(f"  m² = 2h/(k.t) = {he_so_m2(ts):g} 1/m²,  lưới {N}x{N} -> A cỡ {N * N}x{N * N},"
          f"  c = m².dx² = {c:.5f}")

    print("\n  a) Kiểm tra thư viện tự cài đặt trên các ví dụ của phần lý thuyết:")
    vi_du = [("[2 5 4; 6 4 3; 5 3 5]", [[2, 5, 4], [6, 4, 3], [5, 3, 5]], "chuẩn 1 = 13"),
             ("[-3 5 -2; 1 -7 4; 2 0 9]", [[-3, 5, -2], [1, -7, 4], [2, 0, 9]], "chuẩn ∞ = 12"),
             ("[5 2; 2 5]", [[5, 2], [2, 5]], "chuẩn 2 = 7"),
             ("[1 3; 0 2]", [[1, 3], [0, 2]], "chuẩn 2 = 3.702")]
    print(f"     {'Ma trận':<26}{'||.||_1':>9}{'||.||_inf':>11}{'||.||_F':>10}{'||.||_2':>10}   Báo cáo")
    for ten, M, ghi_chu in vi_du:
        M = np.array(M, float)
        print(f"     {ten:<26}{chuan_1(M):9.4f}{chuan_vo_cung(M):11.4f}{chuan_frobenius(M):10.4f}"
              f"{chuan_2(M):10.4f}   {ghi_chu}")

    print(f"\n  b) Các chuẩn của ma trận A ({N * N}x{N * N}):")
    lam_min, lam_max = tri_rieng_A(ts, N)
    n2, so_lan = chuan_2(A, tra_ve_so_lan_lap=True)
    dong = [
        ("Chuẩn 1",            chuan_1(A),             np.linalg.norm(A, 1),      8 + c),
        ("Chuẩn vô cùng",      chuan_vo_cung(A),       np.linalg.norm(A, np.inf), 8 + c),
        ("Frobenius (đ.nghĩa)", chuan_frobenius(A),    np.linalg.norm(A, "fro"),
         np.sqrt(N * N * (4 + c) ** 2 + 4 * N * (N - 1))),
        ("Frobenius (vết)",    chuan_frobenius_vet(A), np.linalg.norm(A, "fro"),
         np.sqrt(N * N * (4 + c) ** 2 + 4 * N * (N - 1))),
        ("Chuẩn 2 (lũy thừa)", n2,                     np.linalg.norm(A, 2),      lam_max),
    ]
    print(f"     {'Chuẩn':<22}{'Tự cài đặt':>14}{'NumPy':>14}{'Công thức':>14}")
    for ten, ta, np_, ct in dong:
        print(f"     {ten:<22}{ta:14.6f}{np_:14.6f}{ct:14.6f}")
    print(f"     (lũy thừa hội tụ sau {so_lan} vòng; công thức: ||A||_1 = ||A||_inf = 8 + c,"
          f" ||A||_2 = lambda_max)")
    nF, n1, ni = chuan_frobenius(A), chuan_1(A), chuan_vo_cung(A)
    print(f"     Kiểm chứng: ||A||_2 = {n2:.4f} <= ||A||_F = {nF:.4f} <= sqrt(r)||A||_2 = {N * n2:.4f};"
          f"  ||A||_2 <= sqrt(||A||_1.||A||_inf) = {np.sqrt(n1 * ni):.4f}")
    return A, b, c


# =====================================================================================
def buoc_2(ts, A, b):
    tieu_de("BƯỚC 2 - SAI SỐ CẢM BIẾN VÀ ĐỘ ỔN ĐỊNH CỦA HỆ  (mục 3.4, 3.5, 4.2)")
    N = ts["N"]
    theta = np.linalg.solve(A, b)
    A_nghich = np.linalg.inv(A)
    lam_min, lam_max = tri_rieng_A(ts, N)

    print("  a) Số điều kiện của A:")
    for p, ten in [("1", "kappa_1  "), ("2", "kappa_2  "), ("inf", "kappa_inf")]:
        print(f"     {ten} = {so_dieu_kien(A, p):10.4f}")
    k2 = so_dieu_kien(A, "2")
    print(f"     Công thức: kappa_2 = lambda_max/lambda_min = {lam_max / lam_min:.4f}"
          f"  -> mất khoảng {np.log10(k2):.1f} chữ số chính xác: hệ ỔN ĐỊNH")

    print(f"\n  b) Mỗi điểm đo dọc cạnh CPU sai ngẫu nhiên tới +-{SAI_SO_CAM_BIEN}°C (nhiễu vế phải b):")
    db = np.zeros_like(b)
    db[:N] = rng.uniform(-SAI_SO_CAM_BIEN, SAI_SO_CAM_BIEN, N)
    dtheta = np.linalg.solve(A, db)
    print(f"     {'Chuẩn':<14}{'||db||/||b||':>14}{'||dθ||/||θ||':>14}{'Cận dưới':>12}{'Cận trên':>12}  Đúng?")
    for p, chuan, ten in [("1", chuan_1, "Chuẩn 1"), ("2", chuan_2, "Chuẩn 2"),
                          ("inf", chuan_vo_cung, "Chuẩn vô cùng")]:
        k = so_dieu_kien(A, p)
        eb = chuan(db[:, None]) / chuan(b[:, None])
        ex = chuan(dtheta[:, None]) / chuan(theta[:, None])
        print(f"     {ten:<14}{eb:14.3e}{ex:14.3e}{eb / k:12.3e}{k * eb:12.3e}"
              f"  {'ĐÚNG' if eb / k <= ex <= k * eb else 'SAI'}")
    can_tuyet_doi = chuan_vo_cung(A_nghich) * np.abs(db).max()
    print(f"     Theo °C: sai số lớn nhất thực tế = {np.abs(dtheta).max():.4f}°C"
          f" <= ||A^-1||_inf.||db||_inf = {can_tuyet_doi:.4f}°C")

    U, s, Vt = np.linalg.svd(A)
    b_xau, db_xau = U[:, 0], 1e-3 * U[:, -1]
    he_so = (np.linalg.norm(np.linalg.solve(A, db_xau)) / np.linalg.norm(np.linalg.solve(A, b_xau))) \
        / (np.linalg.norm(db_xau) / np.linalg.norm(b_xau))
    print(f"     Trường hợp xấu nhất (b theo u_1, nhiễu theo u_n): hệ số khuếch đại = {he_so:.4f}"
          f" = kappa_2 -> cận trên đạt được")

    print("\n  c) Chia lưới mịn hơn: mô hình chính xác hơn nhưng kappa tăng và phép lặp chậm hơn")
    print(f"     {'N':>4}{'Số ẩn':>8}{'kappa_2 (tự tính)':>19}{'kappa_2 (công thức)':>21}"
          f"{'q Jacobi':>11}{'Sai số mô hình':>16}")
    bang = []
    for Nn in [5, 10, 20, 40]:
        An, bn, cn = lap_he(ts, Nn)
        X, Y = toa_do_nut(ts, Nn)
        chinh_xac = nghiem_giai_tich(ts, X, Y).ravel()
        th = np.linalg.solve(An, bn)
        sai_so = chuan_frobenius(th - chinh_xac) / chuan_frobenius(chinh_xac)
        lmin, lmax = tri_rieng_A(ts, Nn)
        k_tu_tinh = chuan_2(An) * chuan_2(np.linalg.inv(An))
        bang.append((Nn, k_tu_tinh))
        print(f"     {Nn:>4}{Nn * Nn:>8}{k_tu_tinh:19.4f}{lmax / lmin:21.4f}{4 / (4 + cn):11.5f}{sai_so:16.3%}")
    print("     (Sai số mô hình = ||θ_lưới - θ_chính xác||_F / ||θ_chính xác||_F, θ_chính xác tính bằng chuỗi Fourier)")

    i = np.arange(1, 11)
    H = 1.0 / (i[:, None] + i[None, :] - 1)
    print(f"     Đối chứng: ma trận Hilbert 10x10 có kappa_2 = {so_dieu_kien(H, '2'):.2e}"
          f" -> dù có {N * N} ẩn, hệ tấm tản nhiệt vẫn ổn định hơn rất nhiều.")
    return theta, bang


# =====================================================================================
def buoc_3(ts, A, b, theta):
    tieu_de(f"BƯỚC 3 - GIẢI LẶP VÀ THUẬT TOÁN DỪNG  (mục 3.6),  eps = {EPS}°C")
    print("  Dừng khi đánh giá hậu nghiệm  q/(1-q).||x_k - x_(k-1)||_inf <= eps,  q = hệ số co (chuẩn vô cùng)")
    ket_qua = {}
    for pp in ["Jacobi", "Gauss-Seidel"]:
        x, q, k_tien, ls = giai_lap(A, b, pp, EPS, x_dung=theta)
        k_ngay_tho = int(np.argmax(ls[:, 1] <= EPS)) + 1
        ket_qua[pp] = (x, q, ls)
        print(f"\n  {pp.upper()}:  q = {q:.5f} < 1 => hội tụ;  q/(1-q) = {q / (1 - q):.1f}")
        print(f"     {'k':>6}{'||x_k - x_k-1||':>18}{'q/(1-q).||.||':>16}{'Sai số thật':>14}")
        chon = list(range(3)) + [None] + list(range(len(ls) - 3, len(ls)))
        for i in chon:
            if i is None:
                print(f"     {'...':>6}")
            else:
                print(f"     {int(ls[i, 0]):>6}{ls[i, 1]:18.4e}{ls[i, 2]:16.4e}{ls[i, 3]:14.4e}")
        print(f"     * Dừng ở k = {len(ls)}: sai số thật {ls[-1, 3]:.2e} <= eps  (tiên nghiệm dự báo cần {k_tien} vòng)")
        sai = ls[k_ngay_tho - 1, 3]
        print(f"     * Nếu dừng 'ngây thơ' khi ||x_k - x_k-1|| <= eps: dừng ở k = {k_ngay_tho},"
              f" sai số thật {sai:.2e} = {sai / EPS:.0f} lần eps  => KHÔNG đạt")
    T = theta + ts["T_inf"]
    print(f"\n  Kết quả: nhiệt độ cao nhất trong tấm {T.max():.2f}°C (sát CPU),"
          f" tại tâm tấm {T.reshape(ts['N'], ts['N'])[ts['N'] // 2, ts['N'] // 2]:.2f}°C")
    return ket_qua


# =====================================================================================
def buoc_4(ts):
    tieu_de("BƯỚC 4 - SO SÁNH VỚI ẢNH CAMERA NHIỆT, NÉN VÀ KHỬ NHIỄU BẰNG SVD  (mục 4.3)")
    n = DO_PHAN_GIAI_CAMERA
    X, Y = toa_do_nut(ts, n)
    that = nghiem_giai_tich(ts, X, Y)                       # trường nhiệt thật (độ tăng nhiệt θ)
    camera = that + rng.normal(0, NHIEU_CAMERA, that.shape)  # ảnh camera có nhiễu
    nT = chuan_frobenius(that)
    print(f"  Ảnh nhiệt {n}x{n} điểm ảnh = ma trận θ (độ tăng nhiệt so với phòng, °C).")
    print(f"  Camera có nhiễu {NHIEU_CAMERA}°C:  ||θ_camera - θ_thật||_F / ||θ_thật||_F ="
          f" {chuan_frobenius(camera - that) / nT:.2%}")

    U, s, Vt = np.linalg.svd(camera)
    nC = chuan_frobenius(camera)
    print("\n  SVD ảnh camera, giữ k thành phần (Eckart-Young: ||C - C_k||_F = sqrt(sum s_i², i>k), ||C - C_k||_2 = s_(k+1))")
    print(f"     {'k':>4}{'||C-Ck||_F':>12}{'sqrt(sum s²)':>14}{'||C-Ck||_2':>12}{'s_(k+1)':>10}"
          f"{'Lệch ảnh camera':>17}{'Lệch trường thật':>18}{'Dung lượng':>12}")
    ds_k = np.arange(1, 61)
    lech_camera, lech_that = [], []
    for k in ds_k:
        Ck = (U[:, :k] * s[:k]) @ Vt[:k]
        lech_camera.append(chuan_frobenius(camera - Ck) / nC)
        lech_that.append(chuan_frobenius(Ck - that) / nT)
        if k in (1, 2, 3, 5, 10, 20, 50):
            E = camera - Ck
            print(f"     {k:>4}{chuan_frobenius(E):12.2f}{np.sqrt((s[k:] ** 2).sum()):14.2f}{chuan_2(E):12.2f}"
                  f"{s[k]:10.2f}{lech_camera[-1]:17.2%}{lech_that[-1]:18.2%}{k * (2 * n + 1) / n ** 2:12.1%}")
    k_tot = int(ds_k[np.argmin(lech_that)])
    Ck_tot = (U[:, :k_tot] * s[:k_tot]) @ Vt[:k_tot]
    print(f"  => Giữ k = {k_tot} thành phần: ảnh chỉ còn lệch trường thật {min(lech_that):.2%}"
          f" (so với {chuan_frobenius(camera - that) / nT:.2%} của ảnh gốc) và chỉ cần lưu"
          f" {k_tot * (2 * n + 1) / n ** 2:.1%} dữ liệu.")
    print("     Trường nhiệt trơn nên có hạng thấp; nhiễu trải đều trên mọi thành phần nên bị cắt bỏ.")
    return that, camera, Ck_tot, k_tot, np.array(lech_camera), np.array(lech_that)


# =====================================================================================
def ve_de_bai(ts):
    """Hình 0: sơ đồ đề bài - tấm nhìn từ trên xuống và nhìn ngang."""
    L = ts["L"] * 100
    vat_lieu = ten_lua_chon(CAU_HOI[0], ts["k"]) or "kim loại"
    nong, lanh, tam = MAU["cam"], MAU["xanh"], "#d9d5cc"
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.4), gridspec_kw={"width_ratios": [1, 1.3]})
    for ax in (ax1, ax2):
        ax.set_aspect("equal")
        ax.axis("off")

    # (a) Nhìn từ trên xuống: cạnh dưới áp vào CPU, ba cạnh kia ở nhiệt độ phòng
    ax1.add_patch(plt.Rectangle((0, 0), L, L, facecolor=tam, edgecolor="none"))
    g = np.linspace(0, L, 10)[1:-1]
    X, Y = np.meshgrid(g, g)
    ax1.plot(X.ravel(), Y.ravel(), "o", color=MAU["chu_phu"], ms=3)
    ax1.text(L / 2, L / 2, "Nhiệt độ tại mỗi\nđiểm lưới = ?", ha="center", va="center", fontsize=12,
             color=MAU["chu"], bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=MAU["truc"]))
    ax1.plot([0, L], [0, 0], color=nong, lw=6, solid_capstyle="butt")
    ax1.plot([0, 0, L, L], [0, L, L, 0], color=lanh, lw=4)
    ax1.add_patch(plt.Rectangle((0, -0.16 * L), L, 0.16 * L, facecolor="#3a3936", edgecolor="none"))
    ax1.text(L / 2, -0.08 * L, f"CPU  {ts['T_CPU']:g}°C", ha="center", va="center", color="white",
             fontsize=12, fontweight="bold")
    ax1.text(L / 2, 0.05 * L, f"cạnh dưới: {ts['T_CPU']:g}°C", ha="center", va="bottom", color=nong, fontsize=10)
    ax1.text(L / 2, 0.95 * L, f"{ts['T_inf']:g}°C (nhiệt độ phòng)", ha="center", va="top", color=lanh, fontsize=10)
    for x, goc in [(0.05 * L, 90), (0.95 * L, -90)]:
        ax1.text(x, L / 2, f"{ts['T_inf']:g}°C", ha="center", va="center", rotation=goc, color=lanh, fontsize=10)
    ax1.annotate("", (0, 1.07 * L), (L, 1.07 * L), arrowprops=dict(arrowstyle="<->", color=MAU["chu_phu"]))
    ax1.text(L / 2, 1.09 * L, f"{L:g} cm", ha="center", va="bottom", color=MAU["chu_phu"])
    ax1.annotate("", (-0.07 * L, 0), (-0.07 * L, L), arrowprops=dict(arrowstyle="<->", color=MAU["chu_phu"]))
    ax1.text(-0.09 * L, L / 2, f"{L:g} cm", ha="right", va="center", rotation=90, color=MAU["chu_phu"])
    ax1.set_xlim(-0.2 * L, 1.05 * L)
    ax1.set_ylim(-0.2 * L, 1.2 * L)
    ax1.set_title(f"(a) Nhìn từ trên xuống: tấm {vat_lieu.lower()} vuông")

    # (b) Nhìn ngang: nhiệt lan dọc tấm (k) và thoát ra hai mặt nhờ quạt (h); độ dày vẽ phóng to
    ax2.add_patch(plt.Rectangle((-1.6, -1.1), 1.6, 2.8, facecolor="#3a3936", edgecolor="none"))
    ax2.text(-0.8, 0.3, "CPU", ha="center", va="center", color="white", fontweight="bold", rotation=90)
    ax2.add_patch(plt.Rectangle((0, 0), 10, 0.6, facecolor=tam, edgecolor=MAU["truc"]))
    ax2.plot([0, 0], [0, 0.6], color=nong, lw=5)
    ax2.plot([10, 10], [0, 0.6], color=lanh, lw=5)
    ax2.annotate("", (9.3, 0.3), (0.4, 0.3), arrowprops=dict(arrowstyle="-|>", color=nong, lw=2))
    for x in (2, 4, 6, 8):
        ax2.annotate("", (x, 1.9), (x, 0.75), arrowprops=dict(arrowstyle="-|>", color=lanh, lw=1.8))
        ax2.annotate("", (x, -1.3), (x, -0.15), arrowprops=dict(arrowstyle="-|>", color=lanh, lw=1.8))
    ax2.text(5, 2.1, f"Quạt thổi hai mặt tấm: nhiệt thoát ra không khí ({ts['T_inf']:g}°C)\n"
             f"hệ số toả nhiệt h = {ts['h']:g} W/m².K", ha="center", va="bottom", color=lanh, fontsize=10)
    ax2.text(5, -1.5, f"Nhiệt từ CPU ({ts['T_CPU']:g}°C) lan dọc tấm\n"
             f"hệ số dẫn nhiệt k = {ts['k']:g} W/m.K", ha="center", va="top", color=nong, fontsize=10)
    ax2.text(10.3, 0.3, f"dày {ts['t'] * 1000:g} mm\n(vẽ phóng to)\nđầu xa: {ts['T_inf']:g}°C",
             ha="left", va="center", color=MAU["chu_phu"], fontsize=9)
    N = ts["N"]
    ax2.text(-1.6, -3.6,
             f"Biết:  tấm {vat_lieu.lower()} {L:g} x {L:g} cm, dày {ts['t'] * 1000:g} mm,"
             f" k = {ts['k']:g} W/m.K, h = {ts['h']:g} W/m².K,\n"
             f"         CPU {ts['T_CPU']:g}°C, phòng {ts['T_inf']:g}°C\n"
             f"Tìm:  nhiệt độ tại {N} x {N} = {N * N} điểm lưới trên tấm\n"
             f"  =>   hệ {N * N} phương trình  A·θ = b  (A cỡ {N * N} x {N * N})",
             ha="left", va="top", color=MAU["chu"], fontsize=10.5, linespacing=1.6,
             bbox=dict(boxstyle="round,pad=0.6", fc="white", ec=MAU["truc"]))
    ax2.set_xlim(-2, 13.2)
    ax2.set_ylim(-7.4, 3.2)
    ax2.set_anchor("N")                        # đặt sát mép trên để hai tiêu đề thẳng hàng
    ax2.set_title("(b) Nhìn ngang: nhiệt vào từ CPU, thoát ra hai mặt")
    luu_hinh(fig, "hinh0_de_bai.png")


def ve_hinh(ts, theta, bang_kappa, ket_qua_lap, kq4):
    cai_dat()
    ve_de_bai(ts)
    N, L = ts["N"], ts["L"] * 100

    # Hình 1: bản đồ nhiệt tính được
    T = np.full((N + 2, N + 2), ts["T_inf"])
    T[1:-1, 1:-1] = theta.reshape(N, N) + ts["T_inf"]
    T[0, 1:-1] = ts["T_CPU"]
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    anh = ax.imshow(T, origin="lower", extent=[0, L, 0, L], cmap=BANG_MAU_NHIET, interpolation="bilinear")
    fig.colorbar(anh, ax=ax, label="Nhiệt độ (°C)")
    ax.set_xlabel(f"x (cm)  —  cạnh dưới áp vào CPU ({ts['T_CPU']:g}°C)")
    ax.set_ylabel("y (cm)")
    ax.set_title(f"Phân bố nhiệt độ trên tấm (lưới {N}x{N}, Gauss-Seidel)")
    ax.grid(False)
    luu_hinh(fig, "hinh1_ban_do_nhiet.png")

    # Hình 2: số điều kiện theo kích thước lưới
    ds_N = np.array([4, 5, 7, 10, 14, 20, 28, 40, 57, 80, 113, 160])
    k_ct = [lmax / lmin for lmin, lmax in (tri_rieng_A(ts, n) for n in ds_N)]
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.loglog(ds_N, k_ct, "-", color=MAU["xanh"], label="Công thức λ_max / λ_min")
    ax.loglog([r[0] for r in bang_kappa], [r[1] for r in bang_kappa], "o", color=MAU["cam"], ms=8,
              label="Tự tính ‖A‖₂·‖A⁻¹‖₂ (lũy thừa)")
    ax.set_xlabel("Số điểm lưới mỗi chiều N  (số ẩn = N²)")
    ax.set_ylabel("Số điều kiện κ₂(A)")
    ax.set_title("Lưới càng mịn, hệ càng nhạy với sai số (κ₂ tăng cỡ N²)")
    ax.legend()
    luu_hinh(fig, "hinh2_so_dieu_kien.png")

    # Hình 3: hội tụ của phép lặp
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    for pp, mau in [("Jacobi", MAU["xanh"]), ("Gauss-Seidel", MAU["cam"])]:
        _, q, ls = ket_qua_lap[pp]
        ax.semilogy(ls[:, 0], ls[:, 3], "-", color=mau, label=f"{pp} (q = {q:.4f})")
        ax.plot(ls[-1, 0], ls[-1, 3], "o", color=mau, ms=8)
        ax.annotate(f"dừng k = {int(ls[-1, 0])}", (ls[-1, 0], ls[-1, 3]), textcoords="offset points",
                    xytext=(-10, 10), ha="right", color=MAU["chu_phu"], fontsize=9)
    ax.axhline(EPS, color=MAU["truc"], lw=1, ls="--")
    ax.text(5, EPS * 1.4, f"ε = {EPS:g} °C", color=MAU["chu_phu"], fontsize=9)
    ax.set_xlabel("Số vòng lặp k")
    ax.set_ylabel("Sai số thật ‖θ⁽ᵏ⁾ − θ*‖∞  (°C)")
    ax.set_title("Hội tụ của phép lặp - dừng theo đánh giá hậu nghiệm")
    ax.legend()
    luu_hinh(fig, "hinh3_hoi_tu.png")

    # Hình 4: ảnh camera nhiệt
    that, camera, Ck, k_tot, lech_camera, lech_that = kq4
    nT = chuan_frobenius(that)
    vmin, vmax = that.min(), that.max()
    fig, axs = plt.subplots(1, 3, figsize=(12, 4))
    for ax, M, td in [(axs[0], that, "Trường nhiệt thật θ"),
                      (axs[1], camera, f"Ảnh camera (lệch {chuan_frobenius(camera - that) / nT:.1%})"),
                      (axs[2], Ck, f"Camera cắt SVD, k = {k_tot} (lệch {chuan_frobenius(Ck - that) / nT:.1%})")]:
        ax.imshow(M, origin="lower", cmap=BANG_MAU_NHIET, vmin=vmin, vmax=vmax, extent=[0, L, 0, L])
        ax.set_title(td)
        ax.axis("off")
    luu_hinh(fig, "hinh4_anh_camera.png")

    # Hình 5: độ lệch theo k
    k = np.arange(1, len(lech_that) + 1)
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.plot(k, lech_that * 100, "-", color=MAU["cam"], label="‖Cₖ − θ_thật‖F / ‖θ_thật‖F")
    ax.plot(k, lech_camera * 100, "-", color=MAU["xanh"], label="‖C − Cₖ‖F / ‖C‖F  (so với ảnh camera)")
    ax.plot(k_tot, lech_that[k_tot - 1] * 100, "o", color=MAU["cam"], ms=8)
    ax.annotate(f"k = {k_tot}: gần trường thật nhất", (k_tot, lech_that[k_tot - 1] * 100),
                textcoords="offset points", xytext=(12, -14), color=MAU["chu_phu"], fontsize=9)
    ax.set_xlabel("Số thành phần SVD giữ lại k")
    ax.set_ylabel("Độ lệch tương đối (%)")
    ax.set_title("Cắt SVD: giữ đủ thông tin, bỏ phần nhiễu")
    ax.set_ylim(0, None)
    ax.legend()
    luu_hinh(fig, "hinh5_do_lech_svd.png")


if __name__ == "__main__":
    print("PHÂN TÍCH NHIỆT ĐỘ TRÊN TẤM TẢN NHIỆT CPU BẰNG CHUẨN MA TRẬN")
    ts = nhap_thong_so()
    A, b, c = buoc_1(ts)
    theta, bang_kappa = buoc_2(ts, A, b)
    ket_qua_lap = buoc_3(ts, A, b, theta)
    kq4 = buoc_4(ts)
    print()
    ve_hinh(ts, ket_qua_lap["Gauss-Seidel"][0], bang_kappa, ket_qua_lap, kq4)
