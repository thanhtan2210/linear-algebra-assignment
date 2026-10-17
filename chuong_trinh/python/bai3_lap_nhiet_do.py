"""
BÀI TOÁN 3 - Thuật toán dừng trong phương pháp lặp  (mục 3.6)
Tính phân bố nhiệt độ trên thanh kim loại tản nhiệt bằng phương pháp Jacobi và Gauss-Seidel.

Mô hình: thanh nhôm dài L = 0.5 m, hai đầu giữ ở T_A = 100°C và T_B = 40°C,
toả nhiệt ra không khí T_inf = 25°C:      T'' - m^2 (T - T_inf) = 0,   m^2 = 100 (1/m^2)
Sai phân với n điểm trong, bước h = L/(n+1), c = m^2 h^2:
      -T_{i-1} + (2 + c) T_i - T_{i+1} = c.T_inf          (ma trận 3 đường chéo A)

Phương pháp lặp  x^(k) = B x^(k-1) + g  hội tụ nếu q = ||B|| < 1, và:
  - Đánh giá tiên nghiệm : ||x^(k) - x*|| <= q^k/(1-q) . ||x^(1) - x^(0)||
  - Đánh giá hậu nghiệm  : ||x^(k) - x*|| <= q/(1-q)  . ||x^(k) - x^(k-1)||
Điều kiện dừng đúng: q/(1-q).||x^(k) - x^(k-1)|| <= eps   (đảm bảo sai số thật <= eps)
Điều kiện dừng "ngây thơ": ||x^(k) - x^(k-1)|| <= eps     (KHÔNG đảm bảo khi q gần 1)
"""
import numpy as np
from chuan_ma_tran import chuan_vo_cung
import kieu_do_thi
from kieu_do_thi import MAU, plt

# ---------------- Dữ liệu bài toán ----------------
L, n = 0.5, 9
T_A, T_B, T_inf, m2 = 100.0, 40.0, 25.0, 100.0
EPS = 1e-3                                    # sai số cho phép (°C)

h = L / (n + 1)
c = m2 * h ** 2
A = (2 + c) * np.eye(n) - np.eye(n, k=1) - np.eye(n, k=-1)
b = np.full(n, c * T_inf)
b[0] += T_A
b[-1] += T_B
x_dung = np.linalg.solve(A, b)                # nghiệm đúng của hệ rời rạc (để so sánh)


def chuan_vec(v):
    return np.abs(v).max()                    # chuẩn vô cùng của vector


def he_so_co(A, phuong_phap):
    """q của Jacobi = ||B_J||_inf ; q của Gauss-Seidel theo công thức hàng."""
    D = np.diag(A)
    B = -A / D[:, None]
    np.fill_diagonal(B, 0)
    if phuong_phap == "Jacobi":
        return chuan_vo_cung(B)
    p = np.array([np.abs(B[i, :i]).sum() for i in range(len(D))])        # phần dưới
    r = np.array([np.abs(B[i, i + 1:]).sum() for i in range(len(D))])    # phần trên
    return (r / (1 - p)).max()


def buoc_lap(A, b, x, phuong_phap):
    x_moi = x.copy()
    for i in range(len(b)):
        nguon = x_moi if phuong_phap == "Gauss-Seidel" else x     # GS dùng ngay giá trị mới
        tong = A[i] @ nguon - A[i, i] * nguon[i]
        x_moi[i] = (b[i] - tong) / A[i, i]
    return x_moi


def giai_lap(phuong_phap, max_lap=1000):
    q = he_so_co(A, phuong_phap)
    x = np.full(n, T_inf)                     # x^(0): cả thanh bằng nhiệt độ môi trường
    lich_su = []                              # (k, ||x_k - x_{k-1}||, đánh giá hậu nghiệm, sai số thật)
    k_ngay_tho = None
    for k in range(1, max_lap + 1):
        x_moi = buoc_lap(A, b, x, phuong_phap)
        buoc = chuan_vec(x_moi - x)
        if k == 1:
            k_tien_nghiem = int(np.ceil(np.log(EPS * (1 - q) / buoc) / np.log(q)))
        hau_nghiem = q / (1 - q) * buoc
        that = chuan_vec(x_moi - x_dung)
        lich_su.append((k, buoc, hau_nghiem, that))
        x = x_moi
        if k_ngay_tho is None and buoc <= EPS:
            k_ngay_tho = k
        if hau_nghiem <= EPS:
            break
    return dict(q=q, x=x, lich_su=lich_su, k_tien=k_tien_nghiem, k_ngay_tho=k_ngay_tho)


def in_ket_qua(ten, kq):
    ls = kq["lich_su"]
    print("=" * 78)
    print(f"PHƯƠNG PHÁP {ten.upper()}   (q = {kq['q']:.4f} < 1 => hội tụ)")
    print(f"  {'k':>4}{'||x_k - x_k-1||':>18}{'q/(1-q)*||.||':>17}{'Sai số thật':>15}")
    hien = ls[:5] + [None] + ls[-3:] if len(ls) > 9 else ls
    for dong in hien:
        if dong is None:
            print(f"  {'...':>4}")
            continue
        print(f"  {dong[0]:>4}{dong[1]:18.6e}{dong[2]:17.6e}{dong[3]:15.6e}")
    k_ht = ls[-1][0]
    kn = kq["k_ngay_tho"]
    print(f"  * Dừng theo đánh giá hậu nghiệm sau k = {k_ht} vòng; sai số thật = {ls[-1][3]:.2e} <= {EPS}")
    print(f"  * Ước lượng tiên nghiệm (trước khi lặp) cần k >= {kq['k_tien']} vòng")
    print(f"  * Nếu dừng 'ngây thơ' (||x_k - x_k-1|| <= eps): dừng ở k = {kn},"
          f" sai số thật = {ls[kn - 1][3]:.2e} "
          f"{'> eps  => KHÔNG đạt yêu cầu!' if ls[kn - 1][3] > EPS else '<= eps'}")


def ve_hinh(kq_j, kq_gs):
    kieu_do_thi.cai_dat()

    # Hình 1: sai số thật theo số vòng lặp
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    for ten, kq, mau in [("Jacobi", kq_j, MAU["xanh"]), ("Gauss-Seidel", kq_gs, MAU["cam"])]:
        ls = np.array(kq["lich_su"])
        ax.semilogy(ls[:, 0], ls[:, 3], "-", color=mau, label=f"{ten} (q = {kq['q']:.3f})")
        ax.plot(ls[-1, 0], ls[-1, 3], "o", color=mau, ms=8)
        ax.annotate(f"dừng k = {int(ls[-1, 0])}", (ls[-1, 0], ls[-1, 3]), textcoords="offset points",
                    xytext=(6, 8), color=MAU["chu_phu"], fontsize=9)
    ax.axhline(EPS, color=MAU["truc"], lw=1, ls="--")
    ax.text(1, EPS * 1.4, f"ε = {EPS:g} °C", color=MAU["chu_phu"], fontsize=9)
    ax.set_xlabel("Số vòng lặp k")
    ax.set_ylabel("‖x⁽ᵏ⁾ − x*‖∞  (°C, thang log)")
    ax.set_title("Sai số theo vòng lặp - dừng bằng đánh giá hậu nghiệm")
    ax.legend()
    kieu_do_thi.luu_hinh(fig, "bai3_hoi_tu.png")

    # Hình 2: phân bố nhiệt độ trên thanh
    m = np.sqrt(m2)
    xs = np.linspace(0, L, 200)
    T_giai_tich = T_inf + ((T_A - T_inf) * np.sinh(m * (L - xs)) + (T_B - T_inf) * np.sinh(m * xs)) / np.sinh(m * L)
    x_nut = np.linspace(0, L, n + 2)
    T_nut = np.concatenate([[T_A], kq_gs["x"], [T_B]])
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.plot(xs, T_giai_tich, "-", color=MAU["xanh"], label="Nghiệm giải tích")
    ax.plot(x_nut, T_nut, "o", color=MAU["cam"], ms=8, label="Nghiệm lặp Gauss-Seidel (các nút)")
    ax.set_xlabel("Vị trí trên thanh x (m)")
    ax.set_ylabel("Nhiệt độ T (°C)")
    ax.set_title("Phân bố nhiệt độ trên thanh tản nhiệt")
    ax.legend()
    kieu_do_thi.luu_hinh(fig, "bai3_nhiet_do.png")
    sai_so_roi_rac = np.abs(T_nut - (T_inf + ((T_A - T_inf) * np.sinh(m * (L - x_nut))
                                              + (T_B - T_inf) * np.sinh(m * x_nut)) / np.sinh(m * L))).max()
    print(f"  Sai số rời rạc hoá (so với nghiệm giải tích) tại các nút: {sai_so_roi_rac:.3f} °C")


if __name__ == "__main__":
    print(f"BÀI TOÁN 3: n = {n} nút trong, h = {h} m, c = m²h² = {c:.4f}, eps = {EPS} °C")
    print("Ma trận A có trội đường chéo chặt: |a_ii| = %.2f > tổng |a_ij| = 2" % (2 + c))
    kq_j = giai_lap("Jacobi")
    kq_gs = giai_lap("Gauss-Seidel")
    in_ket_qua("Jacobi", kq_j)
    in_ket_qua("Gauss-Seidel", kq_gs)
    print("=" * 78)
    print("  Nhiệt độ tại các nút (°C):", np.round(kq_gs["x"], 3))
    ve_hinh(kq_j, kq_gs)
