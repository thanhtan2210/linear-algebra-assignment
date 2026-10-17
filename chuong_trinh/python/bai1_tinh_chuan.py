"""
BÀI TOÁN 1 - Chương trình tính chuẩn ma trận từ dữ liệu nhập vào  (mục 3.3, 4.1)

- Nhập ma trận A từ bàn phím (cú pháp giống MATLAB: 1 2 3; 4 5 6).
- Tính 4 chuẩn bằng code TỰ CÀI ĐẶT theo công thức, rồi đối chiếu với numpy.linalg.norm.
- Kiểm chứng bằng số các tính chất của chuẩn ma trận.
Nhấn Enter (không nhập gì) để chạy các ví dụ trong báo cáo.
"""
import numpy as np
from chuan_ma_tran import (chuan_1, chuan_vo_cung, chuan_frobenius, chuan_frobenius_vet,
                           chuan_2, doc_ma_tran)

# Các ma trận ví dụ lấy từ phần lý thuyết của báo cáo
VI_DU = {
    "VD chuẩn 1 (mục 2.3.1)":     np.array([[2, 5, 4], [6, 4, 3], [5, 3, 5]]),
    "VD chuẩn vô cùng (2.3.2)":   np.array([[-3, 5, -2], [1, -7, 4], [2, 0, 9]]),
    "VD1 chuẩn 2 - đối xứng":     np.array([[5, 2], [2, 5]]),
    "VD2 chuẩn 2 - không đ.xứng": np.array([[1, 3], [0, 2]]),
}


def bang_chuan(A):
    """In bảng so sánh chuẩn tự cài đặt với NumPy."""
    n2, so_lan = chuan_2(A, tra_ve_so_lan_lap=True)
    dong = [
        ("Chuẩn 1",              chuan_1(A),         np.linalg.norm(A, 1)),
        ("Chuẩn vô cùng",        chuan_vo_cung(A),   np.linalg.norm(A, np.inf)),
        ("Frobenius (đ.nghĩa)",  chuan_frobenius(A), np.linalg.norm(A, "fro")),
        ("Frobenius (vết)",      chuan_frobenius_vet(A), np.linalg.norm(A, "fro")),
        ("Chuẩn 2 (lũy thừa)",   n2,                 np.linalg.norm(A, 2)),
    ]
    print(f"  {'Chuẩn':<22}{'Tự cài đặt':>14}{'NumPy':>14}{'Chênh lệch':>13}")
    for ten, ta, np_ in dong:
        print(f"  {ten:<22}{ta:>14.6f}{np_:>14.6f}{abs(ta - np_):>13.1e}")
    print(f"  (Phương pháp lũy thừa hội tụ sau {so_lan} vòng lặp)")


def kiem_tra_tinh_chat(A):
    """Kiểm chứng số học các tính chất/bất đẳng thức của chuẩn."""
    rng = np.random.default_rng(1)
    m, n = A.shape
    B = rng.integers(-5, 6, size=(m, n)).astype(float)   # cùng cỡ để cộng
    C = rng.integers(-5, 6, size=(n, n)).astype(float)   # để nhân A @ C
    x = rng.standard_normal(n)
    alpha = -2.5
    print("  Kiểm chứng tính chất (dùng chuẩn 2; B, C, x sinh ngẫu nhiên):")
    ds = [
        ("||aA|| = |a|.||A||",       chuan_2(alpha * A),  abs(alpha) * chuan_2(A), "="),
        ("||A+B|| <= ||A||+||B||",   chuan_2(A + B),      chuan_2(A) + chuan_2(B), "<="),
        ("||AC|| <= ||A||.||C||",    chuan_2(A @ C),      chuan_2(A) * chuan_2(C), "<="),
        ("||Ax|| <= ||A||.||x||",    np.linalg.norm(A @ x), chuan_2(A) * np.linalg.norm(x), "<="),
        ("||A||_2 <= ||A||_F",       chuan_2(A),          chuan_frobenius(A), "<="),
        ("||A||_F <= sqrt(r)||A||_2", chuan_frobenius(A),
         np.sqrt(np.linalg.matrix_rank(A)) * chuan_2(A), "<="),
        ("||A||_2 <= sqrt(||A||_1||A||_inf)", chuan_2(A),
         np.sqrt(chuan_1(A) * chuan_vo_cung(A)), "<="),
    ]
    for ten, trai, phai, dau in ds:
        dung = np.isclose(trai, phai) if dau == "=" else trai <= phai * (1 + 1e-12)
        print(f"    {ten:<36} {trai:10.4f} {dau:>2} {phai:10.4f}   {'ĐÚNG' if dung else 'SAI'}")


def xu_ly(ten, A):
    print("=" * 72)
    print(f"{ten}:  A =")
    print("\n".join("  " + dong for dong in str(A).splitlines()))
    bang_chuan(A)
    kiem_tra_tinh_chat(A)


def main():
    print("BÀI TOÁN 1: TÍNH CHUẨN MA TRẬN")
    try:
        chuoi = input("Nhập ma trận A (vd: 1 2 3; 4 5 6), Enter để chạy ví dụ: ")
    except EOFError:
        chuoi = ""
        print()
    if chuoi.strip():
        xu_ly("Ma trận nhập vào", doc_ma_tran(chuoi))
    else:
        for ten, A in VI_DU.items():
            xu_ly(ten, A.astype(float))


if __name__ == "__main__":
    main()
