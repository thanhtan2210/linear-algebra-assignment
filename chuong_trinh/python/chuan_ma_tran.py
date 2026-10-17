"""
Thư viện tự cài đặt các chuẩn ma trận theo đúng công thức lý thuyết
(không dùng hàm norm có sẵn), dùng chung cho các bài toán 1-5.

    ||A||_1   = max_j  sum_i |a_ij|          (chuẩn cột)
    ||A||_inf = max_i  sum_j |a_ij|          (chuẩn hàng)
    ||A||_F   = sqrt( sum_ij |a_ij|^2 ) = sqrt( tr(A^T A) )
    ||A||_2   = sqrt( lambda_max(A^T A) ) = sigma_max(A)
"""
import sys
import numpy as np

# In tiếng Việt có dấu trên cửa sổ dòng lệnh Windows
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


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
    x = np.ones(n) / np.sqrt(n)
    if np.allclose(M @ x, 0):         # vector đầu trùng nhân của M -> đổi vector
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

TEN_CHUAN = {"1": "Chuẩn 1", "2": "Chuẩn 2", "inf": "Chuẩn vô cùng", "fro": "Chuẩn Frobenius"}


def so_dieu_kien(A, p="2"):
    """Số điều kiện kappa_p(A) = ||A||_p * ||A^{-1}||_p (A vuông khả nghịch)."""
    A = np.asarray(A, dtype=float)
    chuan = CAC_CHUAN[p]
    return chuan(A) * chuan(np.linalg.inv(A))


def doc_ma_tran(chuoi):
    """Đọc ma trận kiểu MATLAB: '1 2 3; 4 5 6' hoặc '1,2;3,4'."""
    hang = [h.replace(",", " ").split() for h in chuoi.strip().strip("[]").split(";")]
    hang = [h for h in hang if h]
    if len({len(h) for h in hang}) != 1:
        raise ValueError("Các hàng phải có cùng số phần tử")
    return np.array([[float(v) for v in h] for h in hang])
