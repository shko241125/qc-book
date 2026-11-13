# R12 — independent technical review of chapter 20 (`ch20/ch20_linear_systems.*`)

Reviewer: R12 · 2026-11-11 · read-only on the repository. All scratch code and outputs are in
`/tmp/claude-0/-home-claude/fd63146d-511d-5ae9-8dab-3100fb85f0db/scratchpad/mid/R12/` (`mine/*.py` are my own checks, `run/` is a scratch copy of `ch20` + `tools`).

Line numbers refer to the md as it stands now (one line, the `![…](ch20_results.png)` image at md:341, was added while I was reviewing; nothing else changed; the py, stdout, json and png are byte-identical to what I started from).

## Verdict

The core of the chapter is sound. I re-derived the polynomial, the phases, the convention conversion, the main experiment (all 15 rows), the pre-registration exponents, the post hoc table and the gate-level instance with my own code and they reproduce. "Both criteria pass" is justified and robust (N = 256 total-call exponent 1.64; 1.62–1.65 under every alternative reading I tried; 1.59–1.69 over 60 other random S).

Two things are wrong and must be fixed before the chapter is used:

1. **A1 — the 2D Poisson CG cost at N = 16 is wrong (1,120 instead of 448)** because `A.nnz` counts 96 explicitly stored zeros. With the true value the CG exponent for 2D Poisson is 1.58, not 1.25, and the sentence "Poisson 행렬에서는 양자 쪽의 호출 수가 CG의 곱셈 수보다 N에 대해 더 가파르게 자란다" is no longer supported for 2D (asymptotically both are N^1.5).
2. **A2 — the PennyLane snippet will very probably fail at `qp.poly_to_angles(..., angle_solver="root-finding")`.** I could not run PennyLane; this is a prediction from a numpy emulation of the 0.45.1 solver that reproduces the documented example to 16 digits. The chapter's polynomial has ‖c‖₁ = 1.56 > 1, which that solver cannot handle.

Everything else is B/C: honesty of the pre-registration wording (ch19's own caveat is not repeated), sample-noise-dominated exponents in 5절 (7), an off-by-one Newton count, a citation attribution, unsupported side remarks, and small rounding/wording items.

Counts: **A 2 · B 8 · C 19**.

---

## A — wrong, must fix

### A1. 2D Poisson, N = 16: "CG nnz × 반복" = 1,120 is wrong; the 2D exponent 1.25 and the "steeper than CG" claim depend on it

- Location: md:315 (`| 2D Poisson | 16 | 9.5 | 125 | 0.0424 | 607 | 7.0 | 1,120 |`), md:325 (`| 2D Poisson | 0.91 | 1.52 | 1.25 |`), md:328 ("Poisson 행렬에서는 양자 쪽의 호출 수가 CG의 곱셈 수보다 N에 대해 더 가파르게 자란다. 호출 하나의 비용을 넣기도 전에 진다."), md:397 ("κ가 N^1.96(1차원), N^0.91(2차원)으로 자라 양자 호출 수가 CG보다 가파르다"), figure panel (5) (grey square at N = 16), stdout lines 106 and 116. Script: `ch20_linear_systems.py:493` (`poisson2d`) and `:508` (`A.nnz * np.mean(its)`).
- What is wrong: for L = 4, `sp.kron(I_, T)` returns BSR (scipy switches to BSR when `2*B.nnz >= B.shape[0]*B.shape[1]`: 20 ≥ 16), and the BSR blocks carry explicit zeros through the sum. `A.nnz` is the number of *stored* entries: 160. The matrix has 64 nonzeros. L = 8 and 16 are unaffected (288, 1,216 are correct).
- Evidence (`mine/c_nnz.py`, `mine/h_rep78.py`, `mine/g_poisson.py`):
  - `L 4 N 16 A.nnz (stored) = 160  true nonzeros = 64  formats: bsr coo  explicit zeros: 96`
  - Same random stream as the chapter, my own CG: 7.0 iterations × 64 = **448** (chapter: 1,120).
  - Log-log fit with the chapter's own other two points: (448, 4,349, 35,507) → **N^1.58**; the printed 1.25 comes from (1,120, 4,349, 35,507).
  - Quantum side for 2D: 1.52 from the chapter's 10-sample means, 1.41 with the exact E[p_s] (see B2). Either way it is *not* steeper than 1.58.
  - Asymptotically both are N^{3/2}: κ ∝ N, CG = N·√κ; quantum = d/√p_s with d ∝ κ ∝ N and E[p_s] = (1/4N)·Σ(λ_min/λ_kl)² → C/N. With N = 64, 256, 1,024 I get quantum N^1.45 and CG N^1.51.
- Proposed fix:
  - Script: count true nonzeros, e.g. `A.eliminate_zeros()` after building, or `nnz = int((Ad != 0).sum())` (as the script already does for `nnz6`), and add `assert A.nnz == np.count_nonzero(Ad)` in the family loop.
  - md:315 → `… | 7.0 | 448 |`; md:325 → CG exponent `1.58` (re-run to get the exact figure).
  - md:328 replacement: "1차원 Poisson에서는 양자 쪽의 호출 수(N^2.3\~2.5)가 CG의 곱셈 수(N^1.99)보다 가파르게 자란다. 2차원에서는 두 지수가 거의 같다(양자 약 1.4\~1.5, CG 약 1.6; 큰 N에서는 둘 다 N^1.5). 같은 지수라면 호출 하나가 2큐비트 게이트 수백 개인 양자 쪽이 지고, 고전 쪽에는 O(N)의 멀티그리드가 있다."
  - md:397: restrict "양자 호출 수가 CG보다 가파르다" to 1D and say "2차원에서는 지수가 같다" (the multigrid sentence that follows already carries the 2D argument).
  - Summary item 6 (md:18) only quotes the 1D numbers and survives, but see B2 for 2.33.
  - Regenerate the figure (panel 5).

### A2. PennyLane snippet: `poly_to_angles(poly, "QSVT", angle_solver="root-finding")` is expected to raise for this polynomial

- Location: md:372 and the note at md:383 ("차수 25에서는 문제가 없지만 … 그래서 스니펫은 κ = 2 사례만 쓴다").
- What is wrong: PennyLane's root-finding solver completes the target with its *analytic* partner. It puts the odd Chebyshev coefficients into P(z) = Σ c_{2j+1} z^j and needs a complementary polynomial with |P|² + |Q|² = 1 on the whole unit circle, i.e. it needs max_θ |Σ c_k e^{ikθ}| ≤ 1. That is stronger than |P(x)| ≤ 1 on [−1, 1]. The CKS coefficients alternate in sign (+−+−…), so at θ = π/2 the sum equals ‖c‖₁ = 1.564 > 1. Then 1 − |P|² changes sign on the circle, R(z) has roots on |z| = 1, the inside/outside split is not 25/25, and `np.array([P, complementary])` raises `ValueError: … inhomogeneous shape`. (The chapter itself notes that ‖c‖₁ = 1.56 is outside the fixed-point range 0.861; the same number breaks this solver.)
- Evidence:
  - `mine/i_pl.py`, `i_pl2.py`, `i_pl3.py`: numpy emulation of `_complementary_poly`, `_compute_qsp_angle`, the QSP→QSVT shift (`update_vals[0] = 3π/4 − (3 + n % 4)π/2`, inner `π/2`, last `−π/4`; confirmed against the 0.45.1 source page) and the `QSVT`/`PCPhase` matrix. On the documented example `poly = [0, 1, 0, −1/2, 0, 1/3]` it gives angles `[-5.498 1.571 1.571 0.583 1.611 0.748]` and `output.real = 0.19610666666647042` (docs: `0.19610666666647059`), so the emulation is faithful.
  - Chapter polynomial: `max over unit circle of |sum_k c_k e^{ik theta}| = 1.5639`; roots of R: 13 at 0, 11 inside, **2 on the circle**, 11 outside → `ValueError: setting an array element with a sequence. The requested array has an inhomogeneous shape`.
  - Scaling the polynomial by ≤ 1/1.564 = 0.639 makes it work: `scale 0.6: max|Re block/scale − P(A)| = 3.0e-11; p_s = 0.1249; sol err = 5.970e-04` (scale 0.64 and above fail).
  - The same error message from `poly_to_angles` (root-finding, degree-25 polynomial, PennyLane 0.41.1) is reported on the PennyLane forum, with staff answering that root-finding "in very degenerate polynomials … is an unstable method" and recommending `angle_solver="iterative"`: https://discuss.pennylane.ai/t/error-when-using-qml-qsvt/8612
  - The rest of the snippet is right: with exact angles pushed through PennyLane's conversion (my phases with both ends shifted by −π/4) the real part of the top-left block equals P(A) to 4e-15, so `PCPhase(a, dim=16, wires)`, the `QubitUnitary` wire order (wire 0 = ancilla = MSB) and "real part" are all correct for this Hermitian dilation.
- Confidence: high but not certain; PennyLane is not installed here and I did not execute the real function.
- Proposed fix (either; (a) is the one I could check numerically):
  - (a) scale and undo:
    ```python
    scale = 0.6                                                               # ‖c‖₁ = 1.56 > 1 이라 root-finding 해법은 그대로는 실패한다. 0.6 × 1.56 < 1
    angles = qp.poly_to_angles(scale * poly, "QSVT", angle_solver="root-finding")
    …
    PA = M[:16, :16].real / scale
    ```
    Expected outputs stay `True`, `0.1249`, `약 6.0e-04` (emulation: 3e-11, 0.1249, 5.970e-04).
  - (b) `angle_solver="iterative"` (SciPy L-BFGS-B, symmetric phases). Not checked by me; the `atol=1e-6` in `np.allclose` may then be too tight.
  - md:383 note (1): replace "차수 25에서는 문제가 없지만 차수가 높아지면 이 변환부터 정밀도가 무너진다. 그래서 스니펫은 κ = 2 사례만 쓴다." with something like "단항식 변환(차 6.5×10⁻¹¹)보다 먼저 걸리는 것은 위상 해법이다. `root-finding` 해법은 Σ c\_k e^{ikθ}의 크기가 단위원 전체에서 1 이하여야 하는데(‖c‖₁ ≤ 1이면 충분) 이 장의 다항식은 ‖c‖₁ = 1.56이다. 그래서 다항식을 0.6배로 줄여 넣고 결과를 0.6으로 나눈다."
  - Add `scale`-related value(s) to the script output if the text quotes them (house rule: numbers only from scripts).

---

## B — should fix (misleading, unsupported, fragile)

### B1. The chapter says the pre-registration was run "그대로" but omits ch19's own caveat that the ch20 thresholds were tuned after a preliminary calculation

- Location: md:9 ("19장의 사전 등록 표(20장 행)를 그대로 실행했다 … **두 통과 기준은 모두 통과했다.**"), md:257–263, learning note md:433.
- What is wrong: `ch19/ch19_candidates.md:345` states "20장과 23a장의 기준은 엄밀한 뜻의 사전 등록이 아니다. 독립 검토가 … 작은 행렬로 예비 계산을 해 본 뒤, 판정이 … 맞춤 구간에 따라 뒤집히지 않도록 구간과 문턱을 고쳤기 때문이다 … 20·23a장이 다시 잰다." ch20 never mentions this (no "예비", no "엄밀한 뜻"). A reader of ch20 alone takes "1.3 이상 2.3 이하" as a blind prediction.
- Also not flagged as analysis choices made after registration (all harmless for the verdict, but they are choices): (i) the degree is the lowest one whose *polynomial relative error on a 4,201-point grid* is ≤ 10⁻³, not the lowest one that reaches solution error 10⁻³ on the instances (actual solution errors are 5.1–7.5×10⁻⁴); (ii) b = ⌈κ² ln(2/ε)⌉ rather than CKS's κ² log(κ/ε); (iii) total calls = d/√(mean p_s), not mean(d/√p_s); (iv) p_s is the QSVT value ‖P(A)b‖², not the registered formula; (v) CG stops on the true solution error; (vi) the fit uses all five κ.
- Evidence that the verdict does not depend on them (`mine/d_main.py`): degree exponent 1.021 (min degree), 1.052 (rigorous tail-sum degree), 1.045 (CKS j₀ with this b), 1.104 (CKS b and j₀). Total-call exponent at N = 256: 1.643 (chapter), 1.650 (mean of calls), 1.650 (median), 1.623 (exact (2m+1)d counts), 1.636 (exact E[p_s]); κ = 2–16 only 1.673, κ = 4–32 only 1.643; bootstrap over the 10 b: 1.643 ± 0.012 (range 1.61–1.68); 60 other random S at N = 256: 1.636 ± 0.023 (range 1.59–1.69).
- Proposed fix: add after md:263, e.g. "19장 7절이 밝힌 대로 이 행의 문턱(1.3\~2.3)은 독립 검토의 예비 계산을 본 뒤에 정한 것이어서 엄밀한 뜻의 사전 등록은 아니다. 이 장이 정한 세부(차수는 격자에서 상대 오차가 ε 이하인 최저 차수, 총 호출은 d/√(p\_s의 평균), 맞춤은 κ 다섯 점 전부)를 바꿔도 지수는 1.62\~1.65에 머문다." (print those alternative exponents from the script first), and soften md:9 to "… 표(20장 행)의 설정대로 실행했다".

### B2. 5절 (7): the quantum N-exponents (2.33, 1.52, 0.15) are fitted to 10-sample means of p_s whose scatter is 20–40%; the chapter's own exact formula gives different values

- Location: md:308–326, md:18 ("양자 총 호출이 N^2.33으로 자란다").
- What is wrong: for Poisson matrices p_s is dominated by the b-component on the lowest mode (χ²₁-distributed), so a mean of 10 has a relative s.d. of 0.3–0.4 (1D) and 0.2 (2D). The chapter's 1D N = 16 value 0.0073 is less than half of the exact expectation 0.0170. The exponent is then quoted to three figures. The exact route (E[p_s] = (1/N)Σ1/(2κλ_i)², no sampling) is already used in the 사후 확인 of (4).
- Evidence (`mine/g_poisson.py`; chapter sample → exact E[p_s]): 1D: 0.0073/0.0058/0.0009 → 0.01696/0.00423/0.00106; 2D: 0.0424/0.0073/0.0014 → 0.02869/0.00706/0.00171. Exponents with exact E[p_s]: 1D **2.46** (chapter 2.33; asymptote 2.5), 2D **1.41** (chapter 1.52; asymptote 1.5), I + γS 0.15 (unchanged). "sd of 10-sample means/E": 0.29–0.40 (1D), 0.18–0.24 (2D).
- Proposed fix: compute the p_s column and the quantum exponent in (7) from E[p_s] (keep the samples only for CG), or state the sampling error next to the exponents. The 1D conclusion gets stronger (2.46 vs 1.99); the 2D one changes as in A1.

### B3. "끝의 회귀 검사가 본문이 인용하는 값을 고정한다" overstates what the asserts pin

- Location: md:209; script lines 611–628.
- Pinned: degrees, b, CKS column, ε-row degrees, degree slope, verdict, the N = 256 slopes and two N = 16/64 slopes, one row of (3) (N = 256, κ = 32), κ_eff degrees, the N = 16 α values, gate counts/depth/T/R, p_s and κ_eff of (6), family degrees, exercise 3, snippet p_s, monomial max.
- Quoted in the text but not pinned: all other rows of table (3) (p_s, solution errors, CG iterations, nnz); table (2) (Newton counts, residuals, 1.3×10⁻¹², 0.062, 4.93×10⁻⁵); max|P| and ‖c‖₁; the whole 사후 확인 table (`EXACT` is computed at py:311–322 and never asserted); 1.9–16.5배; SDD row counts 256/238/167; s, max|a_ij|, α, κ_eff of (5); γ, α, 4.58×10⁻⁴, rG values of (6); every p_s, call count, CG value and all nine N-exponents of (7); the readout table (8); exercises 2 and 4; 2.4×10⁻¹⁵, 6.5×10⁻¹¹, 1.4897.
- Why it matters: A1 would have been caught by one assert; a scipy change in `sp.random`'s stream would silently change most tables.
- Proposed fix: either add asserts for the values above (at least the 사후 table, the (7) exponents, the SDD counts and `A.nnz == count_nonzero`) or change md:209 to "끝의 회귀 검사가 판정에 쓰는 값과 표의 대표 값을 고정한다".

### B4. Newton iteration count is off by one ("5\~6회")

- Location: md:14, md:225–231 (column "Newton 반복": 6, 5, 5, 5, 5), md:302 ("Newton 5회"), md:432 ("Newton 법이 5\~6회에 끝난다"); script py:112–127, 137.
- What is wrong: `len(hist)` counts residual evaluations, including the last one that only detects convergence (the loop breaks before the update). The number of Newton updates is `len(hist) − 1`.
- Evidence (`mine/d_main.py`, own Newton with the same start and tolerance): update steps 5, 4, 4, 4, 4 for κ = 2…32, residuals 2.2e-16 … 4.0e-15.
- Proposed fix: print `len(hist) - 1` and write "Newton 법 4\~5회" (md:14, table, md:302 → 4회, md:432), or rename the column to "잔차 평가 횟수".

### B5. max|P| = "1.00" at ε = 10⁻⁴ hides that it is above 1

- Location: md:221 ("max|P|가 0.74, 0.88, 1.00, 1.22로 커진다. 1을 넘으면 그대로는 QSVT에 쓸 수 없어"), md:399 ("ε = 10⁻⁴이면 max|P| = 1.00").
- What is wrong: the value is 1.0048 (`mine/a_poly2.py`: `kappa=8 eps=0.0001 … max|P|=1.0048`). Printed as 1.00, the reader cannot tell whether 10⁻⁴ is still usable; the paragraph title says "1을 넘는다".
- Extra, useful for the text: max|P| ≈ 0.638·√(ln(2/ε))/2 independent of κ (0.880 at 10⁻³, 1.004 at 10⁻⁴, 1.215 at 10⁻⁶), so the 1/(2κ) scaling fails for every κ once ε ≲ 1.1×10⁻⁴.
- Proposed fix: print three decimals in the script and write "0.740, 0.881, 1.005, 1.216" and "ε = 10⁻⁴이면 max|P| = 1.005로 1을 막 넘고".

### B6. Existence of symmetric phases is attributed to the wrong paper

- Location: md:77 ("실수 계수에 패리티가 정해진 |P| ≤ 1인 다항식에는 이런 위상이 있다(Dong et al. 2021)").
- What is wrong: Dong–Meng–Whaley–Lin (2021) proposed the symmetric ansatz and the optimisation method. Existence (and uniqueness in the admissible class) is Theorem 1 of Wang, Dong, Lin, "On the energy landscape of symmetric quantum signal processing", *Quantum* 6, 850 (2022), arXiv:2110.04993, which itself says Ref. [5] (DMWL) "suggests that the solution can be restricted to be symmetric". Its standing assumption is the strict bound ‖f‖∞ < 1.
- Evidence: https://arxiv.org/pdf/2110.04993 (fetched; introduction and Theorem 1).
- Proposed fix: "… ‖P‖∞ < 1인 다항식에는 이런 위상이 있다(Wang–Dong–Lin 2022, 정리 1; 대칭 위상을 처음 제안한 것은 Dong et al. 2021)" and add the reference to 10절. The chapter's polynomials satisfy the strict bound (max 0.898).

### B7. Learning note asserts "16장의 최소제곱 맞춤이 d = 105에서 멈춰서" — a specific result that no script produces, and it contradicts 5절 (2)

- Location: md:432 vs md:233 ("초안 단계에서 이 장의 차수에 써 보았으나 수렴이 느려 버렸다. 그 시험은 스크립트에 넣지 않았으므로 수치는 적지 않는다").
- What is wrong: (2) says "slow, no numbers"; the note says "stalled at d = 105". "d = 105" is a number with no source in `ch20_stdout.txt` (house rule). For what it is worth, a plain `scipy.optimize.least_squares` on the reduced symmetric phases from a zero start solves d = 25 to 2e-16 immediately (`mine/b_qsp.py`), so the reader cannot reproduce the claim.
- Proposed fix: either add the test to the script and quote its output, or make md:432 match md:233: "이 장의 분기점: 16장의 최소제곱 맞춤 대신 Newton 법을 썼다(이유는 5절 (2))."

### B8. Summary item 5 mixes two sources under "부록 A에서 잰 CZ 오류율"

- Location: md:17 ("부록 A에서 잰 CZ 오류율이면 rG = 28\~43이다").
- What is wrong: only 28.4 (r = 2.0×10⁻³, the middle of 1.58–2.46×10⁻³) comes from appendix A; 42.5 uses the ch14 public median 3.0×10⁻³ (md:305 says so correctly). And those appendix-A rates are the backend's calibration values at job time (appA: "이 작업이 쓴 연결의 CZ 오류는 …"), not something the book measured; ch00/CLAUDE.md reserve "잰/실측" language for QPU results.
- Proposed fix: "14장의 공개 중앙값(3.0×10⁻³)이면 rG = 42.5, 부록 A 실행 당일의 보정값(가운데 2.0×10⁻³)이면 28.4다."

---

## C — minor / style

- **C1. md:248** p_s (N = 64, κ = 32) is printed 0.0071; the value is 0.0070499 → **0.0070** (double rounding from stdout's 0.00705).
- **C2. md:235** "0.3% 안에서 같다" is not in stdout; the true maximum over the 15 rows is 0.18% (0.19% per individual b). Print it and write "0.2% 안".
- **C3. md:305** "약 60\~140배" is not in stdout (28.4/0.5 = 56.8, 42.5/0.3 = 141.7). Print it or write "약 57\~142배".
- **C4. md:164, 306** "(π/2)/√p\_s = 7.3배 … 약 102,953개": six significant figures on an estimate, and the chapter's own exercise 3 gives the exact count m = 3 → 2m + 1 = 7 → 7 × 14,178 = 99,246. Suggest "약 7배(정확히는 2m + 1 = 7, 99,246개; 어림 (π/2)/√p\_s = 7.3)".
- **C5. md:421** "어림은 1\~4% 안에서 맞는다": the three errors are −1.1%, −0.4%, +3.7% → "4% 안에서".
- **C6. md:15** "무작위 b는 N = 256에서 최악보다 1.9\~16.5배 낫지만": this is the ratio of success probabilities (calls differ by its square root, 1.4–4.1×). Write "성공 확률이 … 1.9\~16.5배 크지만". The 1.71 in the same sentence is the post hoc number and should carry "(사후 확인)".
- **C7. md:79** "이 장의 두 가지 블록 인코딩이 모두 그렇다": the 9-qubit gate-list U of (6) is Hermitian only on the work = 00 sector (‖U − U†‖ = 1.2 on the full space; 1.1e-15 and U² = I to 9e-15 on work = 00; `mine/e_gates.py`). The circuit never leaves that sector, so the result is right, but the script has no assert for the property the formula relies on. Add "(작업 큐비트 = 00인 부분공간에서)" and an assert.
- **C8. md:141, 387** "만들려면 A를 대각화해야 한다": strictly one needs √(I − A²) (or any factor B with BB† = I − A², e.g. Cholesky), which can be computed without diagonalisation but costs at least as much as solving A x = b; and a dense 2N × 2N unitary needs O(N²) gates anyway. Suggest "√(I − A²)를 고전적으로 계산해야 하고(대각화나 그에 맞먹는 O(N³) 계산) 그 비용은 선형계를 직접 푸는 것 이상이다".
- **C9. md:391** "네 단서 — b의 준비, A의 접근 모형, κ, 출력 — 가 이 장의 5절 (5)·(7)·(8)에 하나씩 대응한다": four caveats, three sections; b-preparation is the one the chapter does not cost. → "뒤의 셋이 5절 (5)·(7)·(8)에 대응하고, b의 준비는 이 장이 세지 않았다".
- **C10. md:401** "횟수를 늘려 가며 시도하거나 고정점 진폭 증폭을 쓴다(9장)": ch09 covers the randomised/BBHT schedule and quantum counting; fixed-point amplitude amplification does not appear in `ch09_grover.md`. Move "(9장)" after "시도하거나" or cite Yoder–Low–Chuang 2014.
- **C11. md:273** "CG 반복 수의 지수 0.54는 이론값 1/2과 맞는다": true for N = 256 only (local slopes 0.74, 0.51, 0.49, 0.47); N = 16 and 64 give 0.31 and 0.37 because CG terminates in ≤ N steps on a discrete spectrum. Add half a sentence.
- **C12. md:211–219, figure (2)** "CKS 상한": with this chapter's b the CKS cut-off only bounds the truncation error by ε in absolute terms, so the total relative error bound is 1.5ε, not ε; it is not an upper bound on the degree for the chapter's criterion. Call it "CKS 끊는 자리" (the script's rigorous column 29/59/123/257/533 is the real upper bound and could be shown instead).
- **C13. md:303** "i^(1−d)는 회로에서 뺀 전역 위상이다": i^(−d) is the omitted global phase; the remaining i comes from the 1/(2i) of the Im extraction (the script comment at py:470 says this correctly).
- **C14. md:157, 163** T estimate: with ε = 10⁻⁶ per rotation and R = 1,892 the synthesis error can accumulate to Rε ≈ 1.9×10⁻³, larger than the 10⁻³ target (ch00 3절 makes exactly this point). Either note it or use ε ≈ 10⁻⁷ (≈ 70 T per rotation). In the same paragraph ε means the synthesis error while everywhere else it is the solution error.
- **C15. Notation clashes inside the chapter** (not covered by ch00): b = right-hand side and b = CKS exponent (table column md:213, 3절, exercise 2 next to exercise 1's "b = (1, 1)/√2"); d = degree and d = row sparsity in the AKP sentence (md:395, flagged in parentheses); s = sparsity (md:287) and s = ±1 (md:154); δ = π/4-or-π/2 shift (md:154) and readout error (md:330). For a first-time reader the b/b clash is the worst; consider writing the exponent as β or "지수 b".
- **C16. Figure** (`ch20_results.png`): panel (1) x-label "x (A의 고유값)" should be "A/α의 고유값"; panel (5) y-label "횟수 (로그)" mixes block-encoding calls and multiplications without the "단위가 서로 다르다" note that panel (4) has, and the grey CG markers have no legend entry (only the title explains them); the grey square at N = 16 is the A1 error. Titles, legends and the caption at md:339 otherwise match the data (2: κ^1.02 / κ^0.54; 4: κ^1.64 / κ^2.02; 6: 8.0 / 38.6 / 22.8).
- **C17. md:383 note (2)** hedges about `qp.BlockEncode` ("어느 노름인지는 적혀 있지 않다 … 행 합 노름으로 판단하면"). The 0.45.1 source settles it: `normalization = maximum(norm(A @ A†, ord=pnp.inf), norm(A† @ A, ord=pnp.inf))` and `A / maximum(normalization, 1)` (https://docs.pennylane.ai/en/stable/_modules/pennylane/ops/qubit/matrix_ops.html). So A *is* divided by 1.4897 here; the note can say so and cite the source page.
- **C18. md escaping.** No unescaped `~ _ [ *` outside code; table cells use `\|`; all 15 tables have consistent column counts. Bare `|` remains in running text at md:42, 63, 69, 77, 93, 130, 221, 287, 330, 391, 399, 415, 429. This complies with CLAUDE.md as it reads now (escape in table cells, no paragraph starting with `|`) and renders correctly, but not with the stricter wording CLAUDE.md had when the review started.
- **C19. Script docstring** py:4 says "약 1분"; the run takes 25 s (md:452 and stdout say 25 s).

Unverified nuance (not graded): md:353 "`qiskit.algorithms`가 1.0에서 빠졌다" — the module was removed in 1.0, but I believe `linear_solvers`/HHL was deprecated and removed earlier (Terra 0.22 → 0.25); I could not confirm the version on the web. The conclusion ("지금의 Qiskit 본체에 없다") is right.

---

## Verified as correct

Numbers (md vs `ch20_stdout.txt`): I checked every number in all ten sections, all tables, the exercises and the learning notes. Apart from A1, C1–C5 and B4/B5 they match stdout or follow by trivial arithmetic (e.g. 15,567 = 2·7,784 − 1; 14,178 = 61·206 + 62·26; 363× = 14,178/39; 875 ≈ "약 900배"). Re-running the script in a scratch copy reproduces `ch20_stdout.txt` (only the two timing fields differ), `ch20_case.json` and `ch20_results.png` byte for byte. Cited hardware numbers match ch14/appA (3.0×10⁻³, 68 ns, budget 0.3–0.5, 1.58–2.46×10⁻³, CZ 3–39, 0.81–0.92).

Independent re-derivation (my code, `mine/`):
- Chebyshev formula for g_b checked against the function itself (≤ 1e-12 for b = 3, 7, 31), coefficients from exact integer binomials. b = ⌈κ² ln(2/ε)⌉ gives untruncated relative error (1 − 1/κ²)^b = 1.3e-4 … 4.98e-4 ≤ ε/2. Minimal degrees 25, 51, 105, 211, 423 are the global minima (full scan; the acceptable set is contiguous) and are unchanged on a 10× finer grid; relative errors, max|P| (0.898 … 0.880), ‖c‖₁ = 1.56, slope 1.021, ε row 67/105/139/211, CKS column 41/81/169/353/735, and the extra degrees 61, 251, 125, 425, 33 all reproduce.
- QSP: with my own 2×2 products, Im⟨0|U_Φ|0⟩ = P to 4e-15 for the case-file phases; ⟨0|U_{−Φ}|0⟩ is the exact conjugate; W = i·e^{−iπ/4 Z} R e^{−iπ/4 Z}; the reflection-form formula with i^d equals the W_x product for random non-symmetric phases (d = 1, 2, 5, 8; ≤ 4e-16); (U_Φ − U_{−Φ})/(2i) has Π-block P(A) for the dilation (2.7e-15) and for a generic 3-ancilla complex Hermitian reflection with indefinite spectrum (8e-15). Jacobian at 0 = 2I (1.3e-12). An independent least-squares solve gives the same phases to 1.4e-16. Consistent with ch16 6절.
- Main experiment: all 15 rows (d, nnz, solution error, p_s mean, worst p_s, both call counts, CG iterations and cost), SDD counts, s and max|a_ij|, the four exponent columns for the three sizes, the 사후 table for N = 16 … 4,096 (1.38/1.58/1.64/1.69/1.71, E[p_s], floor shares) reproduce with my own QSVT, Newton and CG. p_s = λ_min²‖A⁻¹b‖²/4 and its independence of α hold (exercise 4: 0.04679 vs 0.04681 vs 0.04683); E[p_s] ≥ 1/(4N) is exact algebra; (2m + 1)d and the exercise-3 values are right.
- Gate-level instance: with my own expansion rules (written from the md paragraph) and my own state-vector simulator, the gate list gives ⟨000,00|U|000,00⟩ = A/α to 2.0e-15, no leakage to work ≠ 00 (3e-16), final block = i^(1−d)·P(A/α)b to 4.2e-14, p_s = 0.046794, error 4.575e-4. Counts: U 206 two-qubit / depth 327 / T 196 / R 30; projector phase 26 / 44 / 28 / R 1 (the second Rz is Clifford); full circuit 35,208 gates, 14,178 two-qubit, depth 22,675, T 13,692, R 1,892, 10 qubits. Nothing is undercounted under the stated convention: X gates of `and3` are one-qubit gates (in the depth, not in G), CZ and CX count 1 each, Rz(−π) and Rz(−π/2) are correctly Clifford. Z_Π = −Z_{w₂} and the rzz/rz angle signs are right. The construction is the plain (non-optimised) one, as the text says.
- κ_eff table, Pauli LCU α = 2.844 with 88 terms, sparse α = 4.824; Poisson κ values (cot² formula); readout arithmetic (⟨M⟩ = −0.0361, var 0.9987, 99.9 / 9,987 shots, 1,110 calls).

Pre-registration: family, κ values, sizes (16, 64, 256 within "16\~256"), ε = 10⁻³, 10 random unit b, metric definitions, baseline and both thresholds match the ch19 row; "차수가 N과 무관" is reported as a check only, as registered; the post hoc analysis is clearly labelled "사후 확인(사전 등록 밖)" in 5절 (but see C6 for the summary).

Concepts and citations checked on the web: CKS Lemmas 17–19 and Theorems 4, 5 (arXiv:1511.02306); Gilyén et al. Lemma 40, Theorem 41, Lemma 48 (√(s_r s_c), |a_ij| ≤ 1), Corollary 69 (arXiv:1806.01838); Dong–Lin–Ni–Wang Lemma 16 (DF(0) = 2I), Theorem 5 (0.861), Im-part target and W(x) convention (arXiv:2209.10162); Andoni–Krauthgamer–Pogrow Theorems 1.1, 1.2 (κ ≤ 3, 1/5·‖x‖∞, n^{Ω(1/d²)}), 5.1 (arXiv:1809.02995); PennyLane 0.45.1 pages for `poly_to_angles`, `QSVT`, `BlockEncode`. The complexity table (HHL κ²/ε, Ambainis κ log³κ/ε³, CKS, QSVT, Costa et al., CG), the κ^{1−δ} ⇒ BQP = PSPACE statement, the Orsucci–Dunjko remark, the dequantisation and multigrid remarks, and "N&C/KLM have no linear-systems algorithm" are consistent with what I know; the max|P| > 1 discussion is right (and is quantified in B5).

SDK snippets: the Qiskit snippet's semantics (RZZ(θ) = exp(−iθ/2 Z⊗Z), RZ(θ) = exp(−iθ/2 Z), `UnitaryGate` first-listed qubit = LSB, `Statevector` indexing, `1j ** (1 − d)`) reproduce 0.1249 and 5.97e-4 in a little-endian numpy emulation (`mine/j_qiskit_emul.py`); I see nothing that should fail on Qiskit 2.x. PennyLane: signatures, wire order, `PCPhase(phi, dim, wires)`, "real part of the top-left block" and the 26 angles are right; the failure is confined to the angle solver (A2).

House rules: ten sections present and in order; one other SDK (PennyLane) with a one-sentence reason; both snippets start with the status comment; "numpy 구현, 실행 완료" label present; cross-references (5절 (2)/(6), 16장 6절·7절, 17장 1·2단계, 19장 학습 노트, 10장 gate-cost rule, 부록 A) point to things that exist and say what is claimed, except C10. The html is byte-identical to a fresh build from the current md (15 tables, figure embedded, no stray escapes).

## Not checked / could not verify

- Qiskit and PennyLane are not installed: neither snippet was executed. A2 rests on a numpy emulation plus the docs/source page and a forum report, not on a real run. Option (b) of A2 (`"iterative"`) is untested.
- Citations taken from memory only (not re-opened): HHL, Ambainis, Costa et al., Orsucci–Dunjko (abstract-level), Aaronson, Montanaro–Pallister, Clader et al., Hestenes–Stiefel, Briggs et al., Babbush et al., Ross–Selinger, Dong–Meng–Whaley–Lin page numbers, the arXiv:2307.12468 abstract, N&C/KLM section numbers, the ITCS venue of AKP, and the Qiskit version in which HHL was removed.
- I did not test the Levenberg–Marquardt claim of B7 in ch16's convention, did not re-run ch14/ch16/ch19/appA scripts, and did not view the html in a browser (structure only).
- The "약 5초" / "0.3초" timings are machine-dependent and were only checked for plausibility (4.5–4.7 s, 0.3 s here).
