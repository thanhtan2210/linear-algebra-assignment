"""
BÀI TOÁN 4 - Đo độ chênh lệch dữ liệu ảnh bằng chuẩn Frobenius  (mục 4.3)

Ảnh xám kích thước m x n chính là một ma trận A (mỗi phần tử = độ sáng 0..255).
  (a) NÉN ẢNH bằng SVD:  A = sum sigma_i u_i v_i^T,  giữ k thành phần đầu -> A_k.
      Định lý Eckart-Young:
          ||A - A_k||_F = sqrt(sigma_{k+1}^2 + ... + sigma_r^2)
          ||A - A_k||_2 = sigma_{k+1}
      Sai số tương đối  e_k = ||A - A_k||_F / ||A||_F ;  MSE = ||A - A_k||_F^2 / (m n)
  (b) ĐO ĐỘ CHÊNH LỆCH giữa ảnh gốc và ảnh bị nhiễu, và kiểm tra việc
      cắt SVD có làm ảnh nhiễu "gần" ảnh gốc hơn không (khử nhiễu).

Chạy:  python bai4_nen_anh_svd.py [đường_dẫn_ảnh]   (bỏ trống -> dùng data/anh_mau.png)
"""
import os
import sys
import numpy as np
from chuan_ma_tran import chuan_frobenius, chuan_2
import kieu_do_thi
from kieu_do_thi import MAU, plt

THU_MUC_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
NGUONG = 0.05                                  # yêu cầu sai số tương đối <= 5%


def tao_anh_mau(N=256):
    """Ảnh thử nghiệm tổng hợp (cùng công thức với tao_anh_mau.m bên MATLAB)."""
    y, x = np.mgrid[0:N, 0:N] / (N - 1)
    A = 90 + 80 * x                                            # nền chuyển sắc
    A[(x - 0.30) ** 2 + (y - 0.35) ** 2 < 0.18 ** 2] = 235     # hình tròn sáng
    A[(x - 0.72) ** 2 + (y - 0.70) ** 2 < 0.14 ** 2] = 30      # hình tròn tối
    vung = (x > 0.55) & (x < 0.92) & (y > 0.08) & (y < 0.40)   # vùng sọc chéo
    A[vung] = 128 + 100 * np.sin(2 * np.pi * 10 * (x[vung] + y[vung]))
    A[y > x + 0.45] = 190                                      # tam giác góc dưới trái
    return np.round(A)


def doc_anh(duong_dan):
    anh = plt.imread(duong_dan).astype(float)
    if anh.max() <= 1.0:                                       # PNG được đọc ở dạng [0,1]
        anh = anh * 255
    if anh.ndim == 3:                                          # ảnh màu -> ảnh xám
        anh = anh[..., :3] @ np.array([0.299, 0.587, 0.114])
    return anh


def xap_xi_hang_k(U, s, Vt, k):
    return (U[:, :k] * s[:k]) @ Vt[:k, :]


def phan_a_nen_anh(A, U, s, Vt):
    m, n = A.shape
    nA = chuan_frobenius(A)
    print("=" * 92)
    print(f"(a) NÉN ẢNH BẰNG SVD   ảnh {m}x{n}, hạng = {np.sum(s > 1e-8 * s[0])},"
          f" ||A||_F = {nA:.1f}, ||A||_2 = sigma_1 = {s[0]:.1f}")
    print(f"  {'k':>4}{'||A-Ak||_F':>13}{'sqrt(sum s^2)':>15}{'||A-Ak||_2':>13}{'s_(k+1)':>11}"
          f"{'Sai số %':>10}{'PSNR dB':>9}{'Dung lượng %':>14}")
    ds_k = [1, 2, 5, 10, 20, 30, 50, 80, 120]
    for k in ds_k:
        Ak = xap_xi_hang_k(U, s, Vt, k)
        E = A - Ak
        eF = chuan_frobenius(E)
        eF_ly_thuyet = np.sqrt((s[k:] ** 2).sum())
        e2 = chuan_2(E)
        psnr = 20 * np.log10(255 / np.sqrt(eF ** 2 / (m * n)))
        dung_luong = k * (m + n + 1) / (m * n) * 100
        print(f"  {k:>4}{eF:13.2f}{eF_ly_thuyet:15.2f}{e2:13.2f}{s[k]:11.2f}"
              f"{eF / nA * 100:10.2f}{psnr:9.2f}{dung_luong:14.1f}")

    # Chọn k nhỏ nhất đạt yêu cầu -> chỉ cần các giá trị kỳ dị, không phải dựng lại ảnh
    sai_so_tuong_doi = np.sqrt(np.cumsum((s ** 2)[::-1])[::-1] / (s ** 2).sum())   # e_k với k = 0..r-1
    sai_so_tuong_doi = np.append(sai_so_tuong_doi[1:], 0.0)                          # dịch thành k = 1..r
    k_chon = int(np.argmax(sai_so_tuong_doi <= NGUONG)) + 1
    print(f"  => k nhỏ nhất để sai số <= {NGUONG:.0%}: k = {k_chon}  (sai số {sai_so_tuong_doi[k_chon - 1]:.2%},"
          f" lưu {k_chon * (m + n + 1) / (m * n):.1%} dung lượng ảnh gốc)")
    return sai_so_tuong_doi, k_chon


def phan_b_nhieu(A, rng):
    print("=" * 92)
    print("(b) ĐỘ CHÊNH LỆCH GIỮA ẢNH GỐC VÀ ẢNH NHIỄU (nhiễu Gauss, độ lệch chuẩn 25)")
    B = A + rng.normal(0, 25, size=A.shape)
    nA = chuan_frobenius(A)
    print(f"  ||B - A||_F / ||A||_F = {chuan_frobenius(B - A) / nA:.2%}   (ảnh nhiễu khác ảnh gốc)")
    U, s, Vt = np.linalg.svd(B, full_matrices=False)
    tot_nhat = None
    for k in [5, 10, 20, 30, 40, 60, 100]:
        e = chuan_frobenius(xap_xi_hang_k(U, s, Vt, k) - A) / nA
        print(f"  k = {k:>3}: ||B_k - A||_F / ||A||_F = {e:.2%}")
        if tot_nhat is None or e < tot_nhat[1]:
            tot_nhat = (k, e)
    print(f"  => Cắt SVD ở k = {tot_nhat[0]} đưa ảnh nhiễu về gần ảnh gốc nhất "
          f"(chênh lệch {tot_nhat[1]:.2%}), tức là đã khử bớt nhiễu.")
    return B, xap_xi_hang_k(U, s, Vt, tot_nhat[0]), tot_nhat


def ve_hinh(A, U, s, Vt, sai_so_tuong_doi, k_chon, B, B_k, tot_nhat):
    kieu_do_thi.cai_dat()
    nA = chuan_frobenius(A)

    # Hình 1: ảnh gốc và các ảnh nén
    ds_k = sorted({5, k_chon, 50})
    fig, axs = plt.subplots(1, 4, figsize=(13, 3.8))
    axs[0].imshow(A, cmap="gray", vmin=0, vmax=255)
    axs[0].set_title("Ảnh gốc")
    for ax, k in zip(axs[1:], ds_k):
        Ak = xap_xi_hang_k(U, s, Vt, k)
        ax.imshow(np.clip(Ak, 0, 255), cmap="gray", vmin=0, vmax=255)
        ax.set_title(f"k = {k}   sai số {chuan_frobenius(A - Ak) / nA:.1%}")
    for ax in axs:
        ax.axis("off")
    kieu_do_thi.luu_hinh(fig, "bai4_anh_nen.png")

    # Hình 2: sai số tương đối theo k
    r = min(len(s), 150)
    k = np.arange(1, r + 1)
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.plot(k, sai_so_tuong_doi[:r] * 100, "-", color=MAU["xanh"], label="‖A − Aₖ‖F / ‖A‖F")
    ax.plot(k, s[1:r + 1] / s[0] * 100, "-", color=MAU["cam"], label="‖A − Aₖ‖₂ / ‖A‖₂ = σₖ₊₁/σ₁")
    ax.axhline(NGUONG * 100, color=MAU["truc"], lw=1, ls="--")
    ax.plot(k_chon, sai_so_tuong_doi[k_chon - 1] * 100, "o", color=MAU["xanh"], ms=8)
    ax.annotate(f"k = {k_chon} đạt sai số ≤ {NGUONG:.0%}", (k_chon, sai_so_tuong_doi[k_chon - 1] * 100),
                textcoords="offset points", xytext=(8, 10), color=MAU["chu_phu"], fontsize=9)
    ax.set_xlabel("Số thành phần SVD giữ lại k")
    ax.set_ylabel("Sai số tương đối (%)")
    ax.set_title("Sai số nén ảnh theo chuẩn Frobenius và chuẩn 2")
    ax.set_ylim(0, None)
    ax.legend()
    kieu_do_thi.luu_hinh(fig, "bai4_sai_so_theo_k.png")

    # Hình 3: khử nhiễu
    fig, axs = plt.subplots(1, 3, figsize=(10, 3.8))
    for ax, anh, tieu_de in [(axs[0], A, "Ảnh gốc A"),
                             (axs[1], B, f"Ảnh nhiễu B  ({chuan_frobenius(B - A) / nA:.1%})"),
                             (axs[2], B_k, f"B cắt SVD, k = {tot_nhat[0]}  ({tot_nhat[1]:.1%})")]:
        ax.imshow(np.clip(anh, 0, 255), cmap="gray", vmin=0, vmax=255)
        ax.set_title(tieu_de)
        ax.axis("off")
    kieu_do_thi.luu_hinh(fig, "bai4_khu_nhieu.png")


def main():
    if len(sys.argv) > 1:
        A = doc_anh(sys.argv[1])
    else:
        duong_dan = os.path.normpath(os.path.join(THU_MUC_DATA, "anh_mau.png"))
        if not os.path.exists(duong_dan):
            os.makedirs(THU_MUC_DATA, exist_ok=True)
            plt.imsave(duong_dan, tao_anh_mau(), cmap="gray", vmin=0, vmax=255)
            print(f"Đã tạo ảnh mẫu: {duong_dan}")
        A = doc_anh(duong_dan)
    print("BÀI TOÁN 4: NÉN VÀ SO SÁNH ẢNH BẰNG CHUẨN MA TRẬN")
    U, s, Vt = np.linalg.svd(A, full_matrices=False)
    sai_so_tuong_doi, k_chon = phan_a_nen_anh(A, U, s, Vt)
    B, B_k, tot_nhat = phan_b_nhieu(A, np.random.default_rng(7))
    ve_hinh(A, U, s, Vt, sai_so_tuong_doi, k_chon, B, B_k, tot_nhat)


if __name__ == "__main__":
    main()
