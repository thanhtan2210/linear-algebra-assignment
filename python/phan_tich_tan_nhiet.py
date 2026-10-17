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

Chạy:  python phan_tich_tan_nhiet.py      (Enter để dùng thông số mặc định)
Cần:   pip install numpy matplotlib
"""
import os
import sys
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


def nhap_thong_so():
    ts = dict(THONG_SO)
    print("Nhập thông số (Enter để giữ giá trị mặc định):")
    for ten, mo_ta, kieu in [("N", "Số điểm lưới mỗi chiều", int),
                             ("T_CPU", "Nhiệt độ CPU (°C)", float),
                             ("h", "Hệ số toả nhiệt của quạt (W/m².K)", float)]:
        try:
            chuoi = input(f"  {mo_ta} [{ts[ten]}]: ").strip()
        except EOFError:
            chuoi = ""
        if chuoi:
            ts[ten] = kieu(chuoi)
    return ts


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
    print(f"  Tấm nhôm {ts['L'] * 100:g} x {ts['L'] * 100:g} cm, dày {ts['t'] * 1000:g} mm,"
          f" T_CPU = {ts['T_CPU']}°C, T_phòng = {ts['T_inf']}°C, h = {ts['h']} W/m².K")
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
def ve_hinh(ts, theta, bang_kappa, ket_qua_lap, kq4):
    cai_dat()
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
