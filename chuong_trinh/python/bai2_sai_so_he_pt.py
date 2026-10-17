"""
BÀI TOÁN 2 - Đánh giá sai số và độ ổn định của hệ phương trình Ax = b  (mục 3.4, 3.5, 4.2)

Lý thuyết dùng:
  Nếu b bị nhiễu thành b + db thì nghiệm x đổi thành x + dx và
      (1/k) * ||db||/||b||  <=  ||dx||/||x||  <=  k * ||db||/||b||,   k = ||A||.||A^-1||
  Nếu A bị nhiễu dA (với k.||dA||/||A|| < 1):
      ||dx||/||x||  <=  k/(1 - k.||dA||/||A||) * ( ||dA||/||A|| + ||db||/||b|| )

Phần A: Mạch điện 3 vòng (định luật Kirchhoff) - ma trận điện trở, điều kiện TỐT.
Phần B: Xấp xỉ đa thức bình phương tối thiểu trên [0,1] - ma trận Gram chính là
        ma trận Hilbert H_ij = 1/(i+j-1), điều kiện RẤT XẤU.
"""
import numpy as np
from chuan_ma_tran import CAC_CHUAN, TEN_CHUAN, so_dieu_kien, chuan_2
import kieu_do_thi
from kieu_do_thi import MAU, plt

rng = np.random.default_rng(2026)


def danh_gia_on_dinh(k):
    """Phân loại mức ổn định theo số điều kiện (quy tắc kinh nghiệm)."""
    mat = max(0.0, np.log10(k))
    if k < 1e2:
        muc = "TỐT (ổn định)"
    elif k < 1e6:
        muc = "TRUNG BÌNH"
    else:
        muc = "XẤU (kém ổn định)"
    return f"{muc}, có thể mất ~{mat:.1f} chữ số chính xác"


def nhieu_ve_phai(A, b, p, ten_he):
    """Thêm nhiễu vào b, so sánh sai số thực tế với cận lý thuyết theo chuẩn p."""
    chuan = CAC_CHUAN[p]
    x = np.linalg.solve(A, b)
    db = b * rng.uniform(-0.01, 0.01, size=b.shape)          # nhiễu đo +-1%
    dx = np.linalg.solve(A, b + db) - x
    k = so_dieu_kien(A, p)
    sai_so_b = chuan(db.reshape(-1, 1)) / chuan(b.reshape(-1, 1))
    sai_so_x = chuan(dx.reshape(-1, 1)) / chuan(x.reshape(-1, 1))
    return k, sai_so_b, sai_so_x


def phan_A_mach_dien():
    print("=" * 76)
    print("PHẦN A: MẠCH ĐIỆN 3 VÒNG - phương pháp dòng mắt lưới  R.I = E")
    # R1=10, R2=5, R3=15 (riêng từng vòng); R12=5 (chung vòng 1-2), R23=10 (chung vòng 2-3)
    R = np.array([[15.0, -5.0, 0.0],
                  [-5.0, 20.0, -10.0],
                  [0.0, -10.0, 25.0]])
    E = np.array([12.0, 0.0, 6.0])         # nguồn 12V ở vòng 1, 6V ở vòng 3
    I = np.linalg.solve(R, E)
    print("  R =\n" + "\n".join("    " + d for d in str(R).splitlines()))
    print(f"  E = {E}  (V)")
    print(f"  Dòng điện mắt lưới I = {np.round(I, 5)}  (A)")

    print("\n  A1) Nguồn đo có sai số +-1%  (nhiễu vế phải b)")
    print(f"  {'Chuẩn':<16}{'k(R)':>9}{'||dE||/||E||':>15}{'||dI||/||I||':>15}"
          f"{'Cận dưới':>12}{'Cận trên':>12}  Đúng?")
    for p in ["1", "2", "inf"]:
        k, eb, ex = nhieu_ve_phai(R, E, p, "mạch")
        print(f"  {TEN_CHUAN[p]:<16}{k:9.4f}{eb:15.3e}{ex:15.3e}{eb / k:12.3e}{k * eb:12.3e}"
              f"  {'ĐÚNG' if eb / k <= ex <= k * eb else 'SAI'}")

    print("\n  A2) Điện trở có dung sai +-5%  (nhiễu ma trận A), dùng chuẩn vô cùng")
    dR = R * rng.uniform(-0.05, 0.05, size=R.shape)
    dI = np.linalg.solve(R + dR, E) - I
    k = so_dieu_kien(R, "inf")
    eA = CAC_CHUAN["inf"](dR) / CAC_CHUAN["inf"](R)
    ex = np.abs(dI).max() / np.abs(I).max()
    can = k / (1 - k * eA) * eA
    print(f"    ||dR||/||R|| = {eA:.4f},  sai số thực tế ||dI||/||I|| = {ex:.4f},"
          f"  cận trên = {can:.4f}  -> {'ĐÚNG' if ex <= can else 'SAI'}")

    print("\n  A3) Trường hợp XẤU NHẤT (chuẩn 2): E theo hướng u_1, nhiễu theo hướng u_n (SVD)")
    U, s, Vt = np.linalg.svd(R)
    E_xau, dE_xau = U[:, 0], 1e-3 * U[:, -1]
    x_xau = np.linalg.solve(R, E_xau)
    dx_xau = np.linalg.solve(R, dE_xau)
    tile = (np.linalg.norm(dx_xau) / np.linalg.norm(x_xau)) / (np.linalg.norm(dE_xau) / np.linalg.norm(E_xau))
    print(f"    Hệ số khuếch đại sai số = {tile:.4f} = k_2(R) = {s[0] / s[-1]:.4f}"
          "  -> cận trên đạt được, không thể làm tốt hơn")
    print(f"\n  Đánh giá: k_2(R) = {so_dieu_kien(R, '2'):.3f} -> {danh_gia_on_dinh(so_dieu_kien(R, '2'))}")


def hilbert(n):
    i = np.arange(1, n + 1)
    return 1.0 / (i[:, None] + i[None, :] - 1)


def phan_B_hilbert():
    print("\n" + "=" * 76)
    print("PHẦN B: XẤP XỈ ĐA THỨC BẬC n-1 TRÊN [0,1] - hệ chuẩn tắc H.c = b, H là ma trận Hilbert")
    print("  Nghiệm đúng c = (1,1,...,1); b được làm nhiễu với ||db||_2/||b||_2 = 1e-8")
    eps = 1e-8
    print(f"  {'n':>3}{'k_2 (tự tính)':>16}{'k_2 (NumPy)':>14}{'Sai số ||dc||/||c||':>22}"
          f"{'Cận k.eps':>12}  Đánh giá")
    ds_n, ds_ss, ds_can = [], [], []
    for n in range(2, 11):
        H = hilbert(n)
        c = np.ones(n)
        b = H @ c
        huong = rng.standard_normal(n)
        db = eps * np.linalg.norm(b) * huong / np.linalg.norm(huong)
        c_moi = np.linalg.solve(H, b + db)
        ss = np.linalg.norm(c_moi - c) / np.linalg.norm(c)
        k_tu_tinh = chuan_2(H) * chuan_2(np.linalg.inv(H))
        k_np = np.linalg.cond(H, 2)
        ds_n.append(n); ds_ss.append(ss); ds_can.append(k_np * eps)
        print(f"  {n:>3}{k_tu_tinh:16.3e}{k_np:14.3e}{ss:22.3e}{k_np * eps:12.3e}"
              f"  {danh_gia_on_dinh(k_np).split(',')[0]}")
    return ds_n, ds_ss, ds_can


def ve_hinh(ds_n, ds_ss, ds_can):
    kieu_do_thi.cai_dat()
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.semilogy(ds_n, ds_can, "-o", color=MAU["cam"], label="Cận trên lý thuyết  k₂(H)·‖δb‖/‖b‖")
    ax.semilogy(ds_n, ds_ss, "-o", color=MAU["xanh"], label="Sai số tương đối thực tế  ‖δc‖/‖c‖")
    ax.axhline(1, color=MAU["truc"], lw=1, ls="--")
    ax.text(ds_n[0], 1.6, "sai số 100% - nghiệm vô nghĩa", color=MAU["chu_phu"], fontsize=9)
    ax.set_xlabel("Cấp n của ma trận Hilbert (bậc đa thức + 1)")
    ax.set_ylabel("Sai số tương đối (thang log)")
    ax.set_title("Nhiễu 10⁻⁸ ở vế phải được khuếch đại theo số điều kiện")
    ax.legend(loc="upper left")
    kieu_do_thi.luu_hinh(fig, "bai2_hilbert_sai_so.png")


if __name__ == "__main__":
    phan_A_mach_dien()
    ve_hinh(*phan_B_hilbert())
