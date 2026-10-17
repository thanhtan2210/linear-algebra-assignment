# BTL Đại số tuyến tính – Chuẩn của ma trận: phần chương trình ứng dụng

Lớp L01 – Nhóm 22 – HK261. Thư mục này gồm **5 bài toán ứng dụng chuẩn ma trận**, mỗi bài có cả bản **Python** và **MATLAB** cho ra cùng kết quả. Các bài được chọn để khớp với dàn ý báo cáo (mục 3.3 → 4.3).

| Bài | Bài toán thực tế | Ma trận dùng | Chuẩn dùng | Mục báo cáo |
|---|---|---|---|---|
| 1 | Chương trình tính chuẩn từ ma trận nhập vào | Ma trận bất kỳ (ví dụ trong báo cáo) | 1, ∞, F, 2 | 3.3, 4.1 |
| 2 | Sai số khi giải mạch điện và khi xấp xỉ đa thức | Ma trận điện trở 3×3; ma trận Hilbert | 1, 2, ∞ (số điều kiện) | 3.4, 3.5, 4.2 |
| 3 | Phân bố nhiệt độ trên thanh tản nhiệt (phương pháp lặp) | Ma trận 3 đường chéo 9×9 | ∞ (hệ số co q) | 3.6 |
| 4 | Nén ảnh và đo độ chênh lệch ảnh | Ma trận ảnh xám 256×256 | Frobenius, 2 | 4.3 |
| 5 | Dự báo thị phần nhà mạng (chuỗi Markov) | Ma trận chuyển trạng thái 3×3 | ∞, 1 | 2.3.2 (ý nghĩa), 5 |

---

## Cấu trúc thư mục

```
python/
  chuan_ma_tran.py      # thư viện TỰ CÀI ĐẶT 4 chuẩn + số điều kiện (dùng chung)
  kieu_do_thi.py        # cấu hình màu/font cho hình vẽ
  bai1_tinh_chuan.py  …  bai5_markov.py
  requirements.txt
matlab/
  chuan_1.m  chuan_vc.m  chuan_F.m  chuan_F_vet.m  chuan_2.m  so_dieu_kien.m   # hàm tự cài đặt
  tao_anh_mau.m  dinh_dang_truc.m  mau_do_thi.m  luu_hinh.m                     # hàm phụ
  bai1_tinh_chuan.m  …  bai5_markov.m
data/anh_mau.png        # ảnh thử nghiệm cho bài 4 (có thể thay bằng ảnh của nhóm)
ket_qua/                # hình xuất ra (*.png của Python, *_matlab.png của MATLAB)
```

## Cách chạy

**Python** (cần Python 3.9+), chạy từ thư mục gốc của repo:

```bash
pip install -r chuong_trinh/python/requirements.txt
```

```bash
cd chuong_trinh/python
python bai1_tinh_chuan.py
```

Tương tự cho `bai2_sai_so_he_pt.py`, `bai3_lap_nhiet_do.py`, `bai4_nen_anh_svd.py [đường_dẫn_ảnh]`, `bai5_markov.py`.

**MATLAB** (R2020b trở lên, không cần toolbox): mở thư mục `chuong_trinh/matlab` làm Current Folder rồi gõ `bai1_tinh_chuan`, `bai2_sai_so_he_pt`, … trong Command Window. Bài 4 muốn dùng ảnh khác thì sửa biến `duong_dan_anh` ở đầu file.

Nguyên tắc chung: **các chuẩn đều được tự cài đặt theo công thức lý thuyết** (chuẩn 2 tính bằng phương pháp lũy thừa trên AᵀA), sau đó đối chiếu với `numpy.linalg.norm` / `norm` của MATLAB để chứng minh là đúng.

---

## Bài 1 – Chương trình tính chuẩn ma trận (mục 3.3, 4.1)

**Đầu vào:** ma trận A nhập từ bàn phím theo cú pháp MATLAB (`1 2 3; 4 5 6`). Nếu nhấn Enter, chương trình chạy 4 ví dụ có sẵn trong phần lý thuyết.

**Thuật toán:**
- ‖A‖₁ = max cột Σ|aᵢⱼ| ; ‖A‖∞ = max hàng Σ|aᵢⱼ| ; ‖A‖_F = √(Σ aᵢⱼ²) = √tr(AᵀA)
- ‖A‖₂ = √λ_max(AᵀA), trong đó λ_max tìm bằng **phương pháp lũy thừa**: yₖ = M xₖ₋₁, xₖ = yₖ/‖yₖ‖, λₖ = xₖᵀ M xₖ, dừng khi |λₖ − λₖ₋₁| ≤ 10⁻¹² λₖ.
- Kiểm chứng bằng số các tính chất: thuần nhất, bất đẳng thức tam giác, ‖AB‖ ≤ ‖A‖‖B‖, ‖Ax‖ ≤ ‖A‖‖x‖, ‖A‖₂ ≤ ‖A‖_F ≤ √r‖A‖₂, ‖A‖₂ ≤ √(‖A‖₁‖A‖∞).

**Kết quả (khớp với các ví dụ trong báo cáo):**

| Ma trận | ‖A‖₁ | ‖A‖∞ | ‖A‖_F | ‖A‖₂ |
|---|---|---|---|---|
| VD chuẩn 1: [2 5 4; 6 4 3; 5 3 5] | **13** | 13 | 12.8452 | 12.3897 |
| VD chuẩn ∞: [−3 5 −2; 1 −7 4; 2 0 9] | 15 | **12** | 13.7477 | 11.7485 |
| VD1 chuẩn 2: [5 2; 2 5] | 7 | 7 | 7.6158 | **7** |
| VD2 chuẩn 2: [1 3; 0 2] | 5 | 4 | 3.7417 | **3.7025** (= √(7+3√5)) |

Chênh lệch giữa code tự viết và thư viện ≤ 3·10⁻¹³; mọi tính chất đều được kiểm chứng ĐÚNG.

---

## Bài 2 – Đánh giá sai số và độ ổn định của hệ phương trình (mục 3.4, 3.5, 4.2)

**Cơ sở lý thuyết:** với số điều kiện κ(A) = ‖A‖·‖A⁻¹‖,

- Nhiễu vế phải: (1/κ)·‖δb‖/‖b‖ ≤ ‖δx‖/‖x‖ ≤ κ·‖δb‖/‖b‖
- Nhiễu ma trận: ‖δx‖/‖x‖ ≤ κ/(1 − κ‖δA‖/‖A‖) · ‖δA‖/‖A‖
- Quy tắc kinh nghiệm: giải hệ có thể **mất khoảng log₁₀κ chữ số chính xác**.

### Phần A – Mạch điện 3 vòng (hệ có điều kiện tốt)

Phương pháp dòng mắt lưới (Kirchhoff) cho hệ **R·I = E**:

```
R = [ 15  -5   0 ;      E = [12; 0; 6] (V)     ->   I = [0.94884; 0.44651; 0.41860] (A)
      -5  20 -10 ;
       0 -10  25 ]
```

- **A1.** Nguồn đo sai ±1% → sai số tương đối của dòng điện luôn nằm giữa cận dưới và cận trên, theo cả chuẩn 1, 2, ∞ (κ₁ = κ∞ = 4.23; κ₂ = 3.60).
- **A2.** Điện trở có dung sai ±5% (nhiễu ma trận) → sai số thực tế 2.9% ≤ cận trên 7.7%.
- **A3.** Chọn E theo hướng véc-tơ kỳ dị u₁ và nhiễu theo hướng uₙ → hệ số khuếch đại sai số **đúng bằng κ₂ = 3.5954**, tức là cận trên đạt được, không thể làm tốt hơn.
- Kết luận: κ₂(R) ≈ 3.6 → hệ **ổn định**, chỉ mất khoảng 0.6 chữ số.

### Phần B – Xấp xỉ đa thức bình phương tối thiểu (hệ có điều kiện rất xấu)

Tìm đa thức bậc n−1 xấp xỉ hàm trên [0,1] dẫn tới hệ chuẩn tắc **H·c = b**, với H là **ma trận Hilbert** Hᵢⱼ = 1/(i+j−1). Cho nghiệm đúng c = (1,…,1) và làm nhiễu b với ‖δb‖/‖b‖ = 10⁻⁸:

| n | κ₂(H) | Sai số thực tế ‖δc‖/‖c‖ | Cận trên κ·10⁻⁸ | Đánh giá |
|---|---|---|---|---|
| 2 | 1.93·10¹ | 1.5·10⁻⁷ | 1.9·10⁻⁷ | Tốt |
| 4 | 1.55·10⁴ | 1.1·10⁻⁴ | 1.6·10⁻⁴ | Trung bình |
| 6 | 1.50·10⁷ | 6.3·10⁻² | 1.5·10⁻¹ | Xấu |
| 8 | 1.53·10¹⁰ | 1.1·10¹ | 1.5·10² | Xấu |
| 10 | 1.60·10¹³ | 6.5·10⁴ | 1.6·10⁵ | Xấu |

Nhiễu chỉ 10⁻⁸ nhưng từ n = 7 trở đi sai số đã vượt 100%, tức là nghiệm vô nghĩa. Đây là lý do thực tế người ta không giải hệ chuẩn tắc trực tiếp mà dùng QR/SVD. Hình: `ket_qua/bai2_hilbert_sai_so.png`.

---

## Bài 3 – Thuật toán dừng trong phương pháp lặp (mục 3.6)

**Bài toán:** thanh nhôm dài 0.5 m, hai đầu giữ ở 100°C và 40°C, toả nhiệt ra không khí 25°C. Phương trình T'' − m²(T − T∞) = 0 (m² = 100), sai phân với 9 nút trong:

```
-T(i-1) + 2.25·T(i) - T(i+1) = 0.25·25        (ma trận A 9×9 ba đường chéo, trội đường chéo chặt)
```

**Vai trò của chuẩn:** đưa về dạng lặp x⁽ᵏ⁾ = Bx⁽ᵏ⁻¹⁾ + g. Nếu **q = ‖B‖∞ < 1** thì phương pháp hội tụ, và:

- Tiên nghiệm: ‖x⁽ᵏ⁾ − x*‖ ≤ qᵏ/(1−q) · ‖x⁽¹⁾ − x⁽⁰⁾‖ → biết trước cần bao nhiêu vòng.
- Hậu nghiệm: ‖x⁽ᵏ⁾ − x*‖ ≤ q/(1−q) · ‖x⁽ᵏ⁾ − x⁽ᵏ⁻¹⁾‖ → **điều kiện dừng đúng**: q/(1−q)·‖x⁽ᵏ⁾ − x⁽ᵏ⁻¹⁾‖ ≤ ε.

**Kết quả (ε = 10⁻³ °C):**

| | Jacobi | Gauss-Seidel |
|---|---|---|
| Hệ số co q | 0.8889 | 0.8000 |
| Ước lượng tiên nghiệm | ≥ 108 vòng | ≥ 54 vòng |
| Dừng theo hậu nghiệm | **64 vòng**, sai số thật 3.7·10⁻⁴ ✔ | **30 vòng**, sai số thật 6.2·10⁻⁴ ✔ |
| Dừng "ngây thơ" ‖x⁽ᵏ⁾ − x⁽ᵏ⁻¹⁾‖ ≤ ε | 52 vòng, sai số thật **2.8·10⁻³ > ε** ✘ | 26 vòng, sai số thật **2.4·10⁻³ > ε** ✘ |

Điểm nhấn cho báo cáo: nếu chỉ dừng khi hai lần lặp liên tiếp gần nhau thì **không bảo đảm** sai số, vì khi q gần 1 hệ số q/(1−q) lớn (= 8 với Jacobi). Phải dùng chuẩn của B để đánh giá. Nghiệm lặp khớp nghiệm giải tích (sai số rời rạc hoá 0.29°C). Hình: `ket_qua/bai3_hoi_tu.png`, `ket_qua/bai3_nhiet_do.png`.

---

## Bài 4 – Đo độ chênh lệch dữ liệu ảnh (mục 4.3)

**Mô hình:** ảnh xám 256×256 là ma trận A. Phân tích SVD A = Σσᵢuᵢvᵢᵀ và giữ k thành phần đầu → Aₖ. Theo **định lý Eckart–Young**:

- ‖A − Aₖ‖_F = √(σ²ₖ₊₁ + … + σ²ᵣ)  và  ‖A − Aₖ‖₂ = σₖ₊₁ (Aₖ là xấp xỉ hạng k tốt nhất)
- Sai số tương đối eₖ = ‖A − Aₖ‖_F / ‖A‖_F; MSE = ‖A − Aₖ‖²_F/(mn); PSNR = 20·log₁₀(255/√MSE)

**(a) Nén ảnh** – kết quả với ảnh mẫu `data/anh_mau.png`:

| k | ‖A−Aₖ‖_F | √Σσ² (lý thuyết) | ‖A−Aₖ‖₂ = σₖ₊₁ | Sai số | PSNR | Dung lượng |
|---|---|---|---|---|---|---|
| 5 | 4782.4 | 4782.4 | 2406.1 | 12.0% | 22.7 dB | 3.9% |
| 20 | 1843.4 | 1843.4 | 454.5 | 4.6% | 31.0 dB | 15.7% |
| 50 | 925.6 | 925.6 | 184.1 | 2.3% | 37.0 dB | 39.1% |

Cột tính trực tiếp và cột lý thuyết trùng nhau, kiểm chứng định lý Eckart–Young. **k nhỏ nhất để sai số ≤ 5% là k = 18**, khi đó chỉ cần lưu 14.1% dung lượng ảnh gốc.

**(b) Ảnh nhiễu** – thêm nhiễu Gauss (σ = 25): ‖B − A‖_F/‖A‖_F = 16.1%. Cắt SVD ảnh nhiễu ở k = 20 giảm độ chênh lệch còn **8.6%**, tức là khử được khoảng một nửa nhiễu. Hình: `ket_qua/bai4_anh_nen.png`, `bai4_sai_so_theo_k.png`, `bai4_khu_nhieu.png`.

---

## Bài 5 – Dự báo thị phần bằng chuỗi Markov

**Bài toán:** 3 nhà mạng A, B, C; mỗi tháng khách hàng chuyển mạng theo ma trận P (tổng mỗi hàng = 1), thị phần hiện tại π₀ = (50%, 30%, 20%), πₖ₊₁ = πₖP.

```
P = [0.8 0.1 0.1 ;
     0.2 0.7 0.1 ;
     0.1 0.2 0.7 ]
```

- ‖P‖∞ = 1 (ma trận ngẫu nhiên hàng): thị phần không bao giờ "bùng nổ", đúng như ý nghĩa của chuẩn vô cùng nêu ở mục 2.3.2.
- Hệ số co δ(P) = ½·maxᵢ,ⱼ‖Pᵢ − Pⱼ‖₁ = 0.7 < 1 → hội tụ về trạng thái dừng duy nhất, với ‖πₖ − π*‖₁ ≤ 2δᵏ.
- Tiên nghiệm cần 41 tháng; dừng theo hậu nghiệm δ/(1−δ)·‖πₖ − πₖ₋₁‖₁ ≤ 10⁻⁶ ở **tháng 27**.
- Thị phần dài hạn: **A 43.75%, B 31.25%, C 25.00%**.

---

## Ghi chú

- Phần Python và MATLAB có cùng thuật toán và cùng dữ liệu. Các phép thử có yếu tố ngẫu nhiên (nhiễu ở bài 2 phần A1/A2 và bài 4b, ma trận B, C ở bài 1) dùng bộ sinh số ngẫu nhiên khác nhau, nên các con số đó lệch nhau chút ít giữa hai ngôn ngữ, nhưng kết luận vẫn giữ nguyên.
- Các tham số vật lý (điện trở, nhiệt độ, ma trận chuyển P) là số liệu giả định hợp lý để minh hoạ; nhóm có thể thay bằng số liệu khác trong code.
