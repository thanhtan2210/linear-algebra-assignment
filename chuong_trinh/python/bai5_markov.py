"""
BÀI TOÁN 5 - Dự báo thị phần bằng chuỗi Markov, dùng chuẩn để chứng minh hội tụ
và quyết định khi nào dừng  (minh hoạ ý nghĩa chuẩn vô cùng / chuẩn 1)

Ba nhà mạng A, B, C. Mỗi tháng khách hàng chuyển mạng theo ma trận chuyển P
(hàng i = xác suất khách của mạng i chuyển sang các mạng, tổng mỗi hàng = 1):
      pi_(k+1) = pi_k . P
  - ||P||_inf = 1 (ma trận ngẫu nhiên hàng) -> thị phần không bao giờ "bùng nổ".
  - Hệ số co Dobrushin  delta(P) = 1/2 . max_{i,j} ||P_i - P_j||_1  < 1  thì
        ||pi_k - pi*||_1 <= delta^k . ||pi_0 - pi*||_1 <= 2 delta^k     (tiên nghiệm)
        ||pi_k - pi*||_1 <= delta/(1-delta) . ||pi_k - pi_(k-1)||_1     (hậu nghiệm)
"""
import numpy as np
from chuan_ma_tran import chuan_vo_cung, chuan_1

P = np.array([[0.80, 0.10, 0.10],
              [0.20, 0.70, 0.10],
              [0.10, 0.20, 0.70]])
pi0 = np.array([0.50, 0.30, 0.20])          # thị phần hiện tại
EPS = 1e-6
TEN = ["A", "B", "C"]


def he_so_dobrushin(P):
    n = P.shape[0]
    return 0.5 * max(np.abs(P[i] - P[j]).sum() for i in range(n) for j in range(n))


def phan_phoi_dung(P):
    """Giải pi (P - I) = 0, sum(pi) = 1 bằng cách thay 1 phương trình bởi điều kiện tổng."""
    n = P.shape[0]
    M = (P - np.eye(n)).T
    M[-1, :] = 1
    rhs = np.zeros(n)
    rhs[-1] = 1
    return np.linalg.solve(M, rhs)


if __name__ == "__main__":
    print("BÀI TOÁN 5: CHUỖI MARKOV - THỊ PHẦN 3 NHÀ MẠNG")
    print("  P =\n" + "\n".join("    " + d for d in str(P).splitlines()))
    print(f"  ||P||_inf = {chuan_vo_cung(P):.4f}   (= 1 vì tổng mỗi hàng bằng 1)")
    print(f"  ||P^T||_1 = {chuan_1(P.T):.4f}   (chuẩn 1 của ma trận chuyển vị cũng = 1)")
    d = he_so_dobrushin(P)
    print(f"  Hệ số co delta(P) = {d:.4f} < 1  => thị phần hội tụ về trạng thái dừng duy nhất")
    k_tien = int(np.ceil(np.log(EPS / 2) / np.log(d)))
    print(f"  Đánh giá tiên nghiệm: cần k >= ln(eps/2)/ln(delta) = {k_tien} tháng để sai số <= {EPS:g}")

    pi_sao = phan_phoi_dung(P)
    pi = pi0.copy()
    print(f"\n  {'Tháng':>5}" + "".join(f"{'Mạng ' + t:>10}" for t in TEN)
          + f"{'||pi_k-pi_k-1||_1':>20}{'Cận hậu nghiệm':>16}{'Sai số thật':>13}")
    print(f"  {0:>5}" + "".join(f"{v:10.4f}" for v in pi))
    for k in range(1, 1000):
        pi_moi = pi @ P
        buoc = np.abs(pi_moi - pi).sum()
        can = d / (1 - d) * buoc
        that = np.abs(pi_moi - pi_sao).sum()
        pi = pi_moi
        if k <= 6 or can <= EPS:
            print(f"  {k:>5}" + "".join(f"{v:10.4f}" for v in pi) + f"{buoc:20.3e}{can:16.3e}{that:13.3e}")
        elif k == 7:
            print(f"  {'...':>5}")
        if can <= EPS:
            break
    print(f"\n  Dừng ở tháng k = {k} (hậu nghiệm), sớm hơn ước lượng tiên nghiệm ({k_tien} tháng).")
    print("  Thị phần dài hạn (giải trực tiếp):", "  ".join(f"{t}: {v:.2%}" for t, v in zip(TEN, pi_sao)))
