# R5 findings — ch18 "탐색 도구와 반례" (independent review)

Files reviewed: `/home/claude/qc-book/ch18/ch18_search_tools.md` (451 lines), `ch18_search_tools.py` (709 lines),
`ch18_stdout.txt`, `ch18_results.png`. Line numbers are of those files as of 2026-11-10.
Scratch code and outputs: `/tmp/claude-0/-home-claude/fd63146d-511d-5ae9-8dab-3100fb85f0db/scratchpad/mid/R5/`
(`c1_synth.py`, `c2_ccz.py`+`c2_out.txt`, `c3_adapt.py`+`c3_out.txt`, `c4_ff.py`, `c5_fkv.py`, `c6_misc.py`, `c7_lin7.py`+`c7_out.txt`,
`c8_misc.py`, `c9_hva_low.py`, `c10_warm.py`, `c11_hva_hi.py`+`c11_out.txt`, `patched/patched_B.py`+`patched_out.txt`, `run/rerun_stdout.txt`).

Re-run: a copy of the script was run once (126 s). Its stdout is identical to `ch18_stdout.txt` except the four timing fields
(D3 9.9 s, E3 SVD 0.14 s vs 0.15 s, per-step timings, total). The png is byte-identical.

Counts: **A 4 · B 9 · C 17** (C1, C8 and C17 partly overlap with A1, B3 and B1).

---

## A — factually or numerically wrong

### [A1] ch18 §5(C) "어디서 깨지는가" and §1 item 5 — md:298–305, md:17; py:455, 469–471
Current text:
> `| 계수 40, 평평한 스펙트럼 | 40 | 0.966 | 1.6×10⁵ | 읽는 것보다 쌌 |`
> `오차는 안정 계수에 비례해 커진다. 세 번째 열의 r은 (안정 계수/0.1)²로 잡은 상한 추정이다.`
> (md:17) `안정 계수가 300을 넘는 행렬에서는 같은 방법이 전체를 읽는 것보다 비싸다.`

What is wrong.
1. The verdict column compares two different units. `need = (srank/0.1)**2` (py:469) is a number of sampled **rows**; the script
   compares it with `1e6`, the number of **entries** (py:471). The matrix has only m = 1000 rows. For rank 40 the estimate is
   1.6×10⁵ rows = the whole matrix 160 times over (r·n = 1.6×10⁸ entries; with the chapter's own "읽은 성분은 r²개" accounting 2.6×10¹⁰).
   "읽는 것보다 쌈" is therefore false for rank 40, and §1's threshold "300을 넘는" is wrong: the method already loses at stable rank 40.
2. This is not an artefact of the loose bound. Measured Gram error ‖AᵀA − SᵀS‖₂/‖A‖² on the same four matrix types (`c5_fkv.py`):

   | matrix (1000×1000) | r = 200 | r = 1 000 | r = 4 000 | r = 16 000 |
   | --- | --- | --- | --- | --- |
   | rank 5 | 0.138 | 0.059 | 0.020 | 0.016 |
   | rank 40 flat | 0.906 | 0.372 | 0.189 | 0.094 |
   | rank 300 flat | 4.38 | 1.57 | 0.75 | 0.33 |
   | tridiagonal | 5.75 | 2.03 | 0.81 | 0.40 |

   Rank 40 needs about 16 000 sampled rows for 10 % — 16 times the number of rows the matrix has. Only the rank-5 matrix
   (about 400 rows, r² ≈ 1.6×10⁵ entries) is cheaper than reading everything.
3. "오차는 안정 계수에 비례해 커진다" is not what the table shows: stable rank ×125 (3 → 375) gives error ×43 (0.136 → 5.84);
   error/(안정 계수/√r) = 0.66, 0.32, 0.21, 0.22 for the four rows. (The spectral error behaves like 2√(s/r) + s/r, s = stable rank; the
   s/√r form is the Frobenius upper bound of §3.)
4. "세 번째 열" is the fourth column. "읽는 것보다 쌌" is a truncated word (stdout: "읽는 것보다 쌈").

Proposed replacement (table and the first two sentences after it, md:298–305):
```
| 행렬 | 안정 계수 | Gram 오차(r = 200) | 10% 오차에 필요한 표본 행 수 r(상한 추정) | 행렬의 행 수 1,000과 비교 |
| --- | --- | --- | --- | --- |
| 계수 5(추천 행렬형) | 3.0 | 0.136 | 9×10² | 행렬보다 적다 — 표본이 이득 |
| 계수 40, 평평한 스펙트럼 | 40 | 0.966 | 1.6×10⁵ | 행렬 전체보다 많다 |
| 계수 300, 평평한 스펙트럼 | 300 | 4.41 | 9×10⁶ | 행렬 전체보다 많다 |
| 3중 대각(희소, 최대 계수) | 375 | 5.84 | 1.4×10⁷ | 행렬 전체보다 많다 |

오차는 안정 계수와 함께 커진다(안정 계수 3 → 375에서 0.14 → 5.8). 넷째 열의 r은 3절의 상한 (안정 계수)/√r ≤ 0.1에서 나온 (안정 계수/0.1)²이다. 견줄 상대는 성분 수 10⁶이 아니라 행 수 1,000이다. 표본 행이 1,000개를 넘으면 행렬을 통째로 읽은 것과 같기 때문이다. 이 기준으로는 계수 5인 행렬만 이득이고(상한 추정으로는 881행이라 아슬아슬하지만, 실제 오차는 r = 200에서 이미 0.14다), 안정 계수가 40만 돼도 표본 방법이 전체를 읽는 것보다 비싸다.
```
md:17, last sentence:
```
1000×1000 행렬에서는 안정 계수가 40만 돼도 필요한 표본 행이 행렬의 행 수를 넘어, 같은 방법이 전체를 읽는 것보다 비싸다.
```
Script (py:471 and the header string at py:456): compare `need` with the row count, e.g.
`'행렬의 행 수보다 적음' if need < N else '행렬 전체(1000행)보다 많음'`, header `"행 수 m = 1000 과 비교"`. If the author wants the
measured crossover in the book, add the r = 1000/4000/16000 loop above to C3 and quote it (rank 40: 0.37, 0.19, 0.09).
Also change the related phrase in md:86 if desired: "크면 표본이 행렬 전체보다 많이 필요하다" is correct as it stands.

### [A2] ch18 §8 문제 3 해답 (a) — md:420 (and the problem statement md:418)
Current text:
> `(a) r ≈ (0.65 × 3/0.05)² ≈ 1,500이고 읽는 성분은 r² ≈ 2.3×10⁶개다. 행렬 전체가 1.92×10⁶개이므로 표본 방법조차 이득이 없다.`

What is wrong. The matrix is 10,000 × 192. "r²" assumes c = r = 1,500 sampled columns, but there are only 192 columns; reading the
1,500 sampled rows in full costs 1,500 × 192 = 288,000 entries = 15 % of the matrix. Checked (`c5_fkv.py`): a 10,000×192 matrix of
stable rank 2.97, r = 1,500 full rows → Gram error 0.055. So the sampling method does read about 7 times less; the claim
"표본 방법조차 이득이 없다" is false. (The conclusion of (b) still holds for the other reason given there: the full SVD takes 0.15 s.)
Also, "5절에서 … c = 0.65" (md:418): §5 never states c; it is printed only by stdout E3.

Proposed replacement for the first three sentences of (a):
```
(a) r ≈ (0.65 × 3/0.05)² ≈ 1,500이다. 열이 192개뿐이므로 열을 다시 뽑을 필요 없이 표본 행을 통째로 읽으면 되고, 읽는 성분은 1,500 × 192 ≈ 2.9×10⁵개로 전체 1.92×10⁶개의 15%다. 표본 방법의 이득은 7배쯤이다(r² ≈ 2.3×10⁶으로 세면 전체보다 많아 보이지만, 그것은 열도 1,500개 뽑는다고 가정한 셈이라 이 행렬에는 맞지 않는다).
```
and in (b), first sentences:
```
(b) 1단계(문제 정의)다. 7배의 이득은 의미가 없다. 이 크기의 행렬은 전체 SVD가 0.15초면 끝난다(스크립트 E3, 이 컨테이너 기준).
```
Problem statement (md:418), parenthesis:
```
(스크립트 E3가 5절의 계수 5 행렬, r = 200에서 맞춘 식: 오차 ≈ c·안정 계수/√r, c = 0.65)
```

### [A3] ch18 §3(나), §5(B), §8, §9 — "대칭 장벽" is attributed to the wrong mechanism — md:76, 282, 414, 426 (and md:139)
Current text:
> (md:76) `둘째, 풀의 모든 P가 어떤 대칭 Q와 교환하면 회로는 Q의 구역별 가중치를 바꿀 수 없다. 시작 상태가 답과 다른 구역에 걸쳐 있으면 거기서 막힌다(대칭 장벽, Shkolnikov et al. 2023).`
> (md:282) `셋째 줄이 대칭 장벽이다.`

What is wrong. Shkolnikov et al. (arXiv:2109.05340, abstract and §4, fetched today from arxiv.org/pdf/2109.05340 and
ar5iv.arxiv.org/html/2109.05340) describe the opposite situation: the reference state lies in one symmetry sector, H conserves
the symmetry, and a pool operator that **changes** the sector has zero gradient at every iteration — "every Pauli string that
changes at least one of these parities will always have zero gradient and will never be picked" — so a complete pool effectively
shrinks; "complete pools can fail to yield convergent results, unless the pool is chosen to obey certain symmetry rules"
(abstract). Their remedy is a pool that **preserves** the symmetry. The scenario ch18 calls 대칭 장벽 (symmetry-preserving pool,
start state spread over two sectors) does not appear in the part of the paper I could read (the HTML fetch is truncated after the
induction argument of §4).
In ch18's own table the Shkolnikov mechanism is the **first** row, not the third: start |+…+⟩ (sector ∏X = +1), pool {Y_i},
{Y_i, ∏X} = 0, all gradients exactly 0 from the start. Checked (`c3_adapt.py`): {Y_i,Q} = 0, [Y_iZ_j,Q] = 0, {Y_iX_j,Q} = 0,
max|g_{Y_i}| at |+…+⟩ = 0.0 for n = 4, 6, 8. The text gives no reason for that row ("처음부터 기울기가 전부 0").

Proposed replacements.
md:76, from "둘째," to the end of the paragraph:
```
둘째, 대칭이 두 가지 방식으로 탐색을 막는다. (i) H가 대칭 Q와 교환하고 지금 상태가 Q의 한 구역 안에 있으면, Q와 반교환하는(구역을 바꾸는) P의 기울기는 언제나 0이다. P|ψ⟩는 다른 구역에 있고 H는 구역을 섞지 않기 때문이다. 그런 연산자는 풀에 있어도 끝내 뽑히지 않는다. Shkolnikov et al.(2023)이 대칭 장벽(symmetry roadblock)이라 부른 것이 이것이고, 처방은 풀을 대칭을 지키는 연산자로 짜는 것이다. (ii) 반대로 풀의 모든 P가 Q와 교환하면 회로는 Q의 구역별 가중치를 바꿀 수 없다. 시작 상태가 답과 다른 구역에 걸쳐 있으면 거기서 막힌다(구역 가중치 보존). 5절의 함정 표에서 첫째 줄이 (i), 셋째 줄이 (ii)다.
```
md:276 (table row 1, last cell):
```
시작 상태의 에너지 오차 그대로. Y\_i가 ∏X와 반교환해 처음부터 기울기가 전부 0(대칭 장벽)
```
md:278 (table row 3, last cell):
```
(E₁ − E₀)/2. 풀이 패리티를 보존한다(구역 가중치 보존)
```
md:282, replace "셋째 줄이 대칭 장벽이다." by:
```
첫째 줄이 Shkolnikov et al.의 대칭 장벽이다. |+…+⟩는 패리티 ∏X = +1 구역에 있고 Y\_i는 ∏X와 반교환하므로 기울기가 처음부터 전부 0이다. 셋째 줄은 그 반대 경우다.
```
md:414: `**문제 2 (자작 — 구역 가중치 보존).**`
md:426, the "깨지는 곳" sentence:
```
깨지는 곳: 기울기가 모두 0이어도 답이 아닐 수 있다. 대칭을 깨는 연산자는 기울기가 0이라 뽑히지 않고(대칭 장벽), 대칭을 지키는 풀은 시작 상태의 구역 가중치를 바꾸지 못해 오차가 (E₁ − E₀)/2에서 멈춘다.
```
md:139 (3′ 탈락 신호) may stay; it describes case (ii) without naming it.

### [A4] ch18 §5(B) — md:272 (small)
Current text: `모든 기울기가 10⁻⁸ 이하인데 오차는 큰 상태에서 멈춘다.`
stdout B2, n = 6: max gradient 1.2e-08 ({Y_i} from |0…0⟩) and 1.5e-08 ({Y_iZ_j} from |0…0⟩). Replacement:
```
모든 기울기가 2×10⁻⁸ 이하인데(멈춤 문턱값은 10⁻⁶) 오차는 큰 상태에서 멈춘다.
```

---

## B — misleading, incomplete, or contradicting another chapter

### [B1] ch18 §5(B) HVA baseline — md:260–268, figure panel (3); py:298–308, 349–354
Current text:
> `| 6 | 6 | 12 | 60 | 2.3×10⁻³ | 0.19 |` · `| 6 | 8 | 16 | 80 | 5.9×10⁻⁵ | 0.05 |`
> `HVA p = 8에서는 시작점 21개 중 1개만 최선에 갔고`

What I found (own HVA implementation with its own adjoint gradient, checked against central differences to 1e-9; L-BFGS-B,
random starts in [−π, π]; `c3_adapt.py`, `c9_hva_low.py`, `c11_hva_hi.py`):

| n | p | book "최선 오차" (21 starts) | my best (random starts) | starts | fraction within 2× of my best |
| --- | --- | --- | --- | --- | --- |
| 4 | 3 | 1.56e-4 | 1.56e-4 | 1500 | 0.31 |
| 4 | 4 | 2.3e-14 | < 1e-13 | 300 | 0.58 |
| 6 | 3 | 2.68e-2 | 2.68e-2 | 300 | 0.44 |
| 6 | 4 | 1.77e-2 | **8.4e-3** | 300 | 0.003 |
| 6 | 5 | 6.79e-3 | 6.79e-3 | 400 | 0.26 |
| 6 | 6 | 2.26e-3 | **1.87e-4** | 400 | 0.008 |
| 6 | 7 | 7.35e-4 | **1.82e-5** | 600 | 0.003 |
| 6 | 8 | 5.90e-5 | **3.7e-6** | 350 | 0.011 |
| 6 | 9 | 1.1e-13 | < 1e-13 | 120 | 0.23 |

- The thresholds are genuine, not optimiser luck: n = 4 exact at p = 4 (58 % of starts) and not at p = 3 (floor 1.56e-4, 1500
  starts); n = 6 exact at p = 9 (23 % of starts) and not reached at p = 8 in 350 starts (best 3.7e-6). They coincide with
  2p = n²/2 (8 and 18 parameters). The headline "HVA는 n = 6에서 파라미터 18개(CNOT 90)" and "62 대 90" therefore stand.
- The sub-threshold entries for n = 6 are optimiser luck. p = 4, 6, 7, 8 are 2×, 12×, 40×, 16× above values that more starts
  find. "최선의 2배 이내로 간 비율" is measured against a "best" that is not the optimum, and the HVA triangles of figure (3) at
  12–16 parameters sit one to one and a half decades too high.
- The "warm start" (py:302, previous optimum + [0, 0]) is a stationary point of the deeper circuit, so it never moves. Reason:
  at θ_new = 0 both new derivatives are ±⟨i[ΣZZ, ΣX]⟩, which is the derivative with respect to the last β and vanishes at the
  previous optimum. Check (`c10_warm.py`): new-angle gradients (6.8e-9, −6.8e-9) at n = 4, p = 3 → 4. In an instrumented copy of
  the book script (`patched/patched_out.txt`) the warm start's final error equals the previous p's best in all 11 cases. So "시작점
  21개" is 20 random starts plus one dead one; at p = 8 it was 1 of 20 random starts.

Proposed replacement for the third "읽는 법" (md:268, from "셋째,"), keeping the book's stdout numbers:
```
셋째, 고정된 틀을 무작위로 시작하면 국소 최소가 많다. 표의 "최선 오차"는 시작점 21개 가운데 최선일 뿐 그 층 수에서의 최적값이 아니다. 믿을 수 있는 것은 문턱이다. 오차 10⁻⁸ 아래로 가는 것은 n = 4에서 p = 4, n = 6에서 p = 9부터이고, 이는 파라미터 수 2p가 n²/2에 닿는 층 수와 같다. 그 아래 층 수의 값은 시작점을 늘리면 더 내려간다. 그리고 n = 6의 p = 1\~9 전체에 에너지·기울기 평가가 60,622번 들었다. ADAPT는 새 각도를 0에서 시작해 에너지가 단조로 내려가므로 2,091번으로 끝났다. 그림 (3).
```
and change the column header (md:260) to `최선 오차(시작점 21개 중)`.
Recommended script change (then regenerate the table and figure): for n = 6 use `starts=300` (about 6 more minutes), and either drop
the warm start or perturb it (`warm + 0.05·rng.normal`). With that the table rows would read approximately p = 6: 1.9×10⁻⁴,
p = 8: 3.7×10⁻⁶, p = 9: 10⁻¹³ (my numbers; they must come from the book's script before they are printed).

### [B2] ch18 §5(A) CCZ search — md:242, 244, 248
Current text:
> `| 6 | 70 | 0 (4묶음이 닿음) | 20 | 0.1464 |`
> `여기서는 바닥값이 구조와 무관하게 같은 상수(0.4375, 0.375, 0.1464)로 나오고 알려진 정리와 일치한다는 점이 보강 근거다.`

What I found (`c2_ccz.py`, `c6_misc.py`, `c7_lin7.py`; the book's `Template` was reused after checking its unitary and gradient against
my own builder: 5e-16 and 5e-12).
- With 12 starts instead of 4, **six** all-to-all classes reach CCZ at k = 6, not four. The two the book misses are
  `01 02 12 01 12 02` and `01 02 12 12 01 02`. Per-start success probability of the six classes (60 single starts each):
  0.63, 0.62, 0.12, 0.38, 0.68, 0.10 — a solvable structure is missed by 4 starts with probability up to 0.66. For the four linear
  8-CNOT classes: 0.20, 0.12, 0.35, 0.55. "4묶음" is a property of the run, not of CCZ.
- "구조와 무관하게 같은 상수" is not true per structure: at k = 4–5 the classes end at 0.3750, 0.2714 or 0.1464; at k = 6 also at
  0.1323 and 0.1461. Only the minimum over structures takes the three quoted values.
- The symmetry reduction itself is valid: CCZ is invariant under qubit permutations and CCZᵀ = CCZ; the template has exactly one
  U3 on every wire segment, so the reversed/permuted template reaches the same set. Nothing is lost in coverage: all 243 raw
  k = 5 structures, 3 starts each, give min 0.14645. What is lost is starts: a class gets 4 starts in total.
- Strengthening the linear claim: all 36 linear k = 7 classes with **40** starts each reach nothing below 0.146447 (106 s). With a
  per-start success of even 0.10, 40 starts miss with probability 1.5 %, against 66 % for 4 starts.

Proposed replacements.
md:242 cell: `0 (시작점 4개로 4묶음이 닿음)` · md:244 cell: `0 (시작점 4개로 4묶음이 닿음)`
md:248:
```
주의: "닿지 못했다"는 수치 탐색의 결과이지 증명이 아니다. 닿을 수 있는 구조에서도 시작점 하나가 성공하는 비율은 10\~70%여서, 시작점 4개로는 닿는 구조를 놓칠 수 있다(전결합 6개에서 시작점을 12개로 늘리면 닿는 묶음이 4개에서 6개로 는다). 여기서는 CNOT 수마다 구조 전체의 최솟값이 0.4375, 0.375, 0.1464라는 닫힌 꼴의 상수로 나오고(구조 하나하나의 바닥은 0.375, 0.2714 등으로 다르다) 알려진 정리와 일치한다는 점이 보강 근거다. 선형 연결 7개는 근거가 탐색뿐이므로 시작점을 40개로 늘려 다시 확인하는 편이 안전하다.
```
(The 12-start and 40-start numbers are mine; if quoted, add the runs to the script: `starts=12` for the two last k, about 2 minutes.)

### [B3] ch18 §4 "고전 경계 검사 한 벌" — md:122–131 (contradicts ch17:150–160, 192)
Current text: `후보는 여섯 줄을 모두 빠져나가야 한다. 하나라도 걸리면 5단계 탈락이다.`
What is wrong. The rows do not all apply to every kind of candidate.
- stoquastic row: absence of a sign problem helps ground-state and thermal (imaginary-time) problems only. Real-time dynamics
  under a stoquastic H is not classically easy for that reason. ch17's surviving candidate "스핀 격자 급랭 동역학" uses the Kim–Huse
  Ising model −ΣZZ − 0.9045ΣX − 0.8090ΣZ (ch17:150), whose off-diagonal elements are all −0.9045 or 0 (checked at n = 6): by ch18's D5 test it is stoquastic, and
  by "하나라도 걸리면 탈락" it would be dropped — contradicting ch17's verdict "가능".
- Even for ground states "stoquastic → QMC 가능" means "no sign problem", not a polynomial-time guarantee; ch15:58 says
  "몬테카를로로 시뮬레이션 가능한 경우가 많다".
- free-fermion row: needs a Gaussian (or few-Gaussian) initial state as well as quadratic terms (see B4).
- dequantization row: applies only to QRAM-input linear-algebra candidates.

Proposed replacement for md:122:
```
**고전 경계 검사 한 벌.** 15장 3절·7절의 경계마다 검사를 하나씩 붙였다. 후보는 자기 문제 종류에 해당하는 줄을 모두 빠져나가야 하고, 하나라도 걸리면 5단계 탈락이다. 해당 범위는 줄마다 다르다. stoquastic 줄은 바닥상태·열평형 문제에만 쓴다. 실시간 동역학 e^{−iHt}는 H가 stoquastic이어도 몬테카를로에 위상 문제가 생기므로 이 줄로 탈락시키지 않는다(17장의 비적분 Ising은 비대각 원소가 모두 음수인 stoquastic 모형이지만 그 급랭 동역학은 "가능" 판정이다). 부호 문제가 없다는 것도 다항 시간 수렴까지 보장하지는 않는다. 자유 페르미온 줄은 초기 상태가 가우스 상태일 때, 입력 모형 줄은 QRAM 입력을 쓰는 선형대수 후보에만 해당한다.
```
and md:129 (table row), second and third cells:
```
| 부호 문제 없음(stoquastic) — 바닥상태·열평형만 | 비대각 원소의 부호(기저 회전 포함) | 비대각 원소가 모두 0 이하 | 이 장 5절 (D) |
```

### [B4] ch18 §1 item 6, §5(D), §6 vs ch17 §5(A) — md:18, 318–327, 337, 384; ch17:150
Current text:
> (ch17:150) `이 급랭은 2n×2n 행렬로 정확히 풀린다(18장 5절: n = 1000, t = 40까지 약 10초).`
> (md:18) `(여기서는 |+…+⟩에서 시작)` · (md:337) `차이가 0이 아니면 상태가 가우스가 아니라는 뜻이고 2n×2n 계산은 쓸 수 없다.`

What is wrong. ch17's quench starts from |0…0⟩; ch18 solves a different one (|+…+⟩) and does not say why. The reason matters:
with the X-string Jordan–Wigner of §3(라), |0…0⟩ is an equal superposition of the two ∏X sectors and is **not** a Gaussian state,
so the Γ(t) calculation of ch18 does not reproduce ch17's numbers, and ch18's own Wick test rejects ch17's quench although it is
classically easy. Check (`c4_ff.py`, `c8_misc.py`; n = 12, critical TFIM, half chain):

| t | S true, from \|0…0⟩ | S "Gaussian" from Γ | max ⟨c_a⟩ | S true = S Gaussian for (\|0…0⟩ ± \|1…1⟩)/√2 |
| --- | --- | --- | --- | --- |
| 1 | 0.6013 | 1.6013 | 0.645 | 1.6013 / 1.6013 |
| 2 | 1.4980 | 2.4979 | 0.561 | 2.4977 / 2.4979 |
| 3 | 2.4364 | 3.3507 | 0.501 | 3.2282 / 3.3189 |

(The true values agree with ch17's n = 16 row 0.601, 1.498.) The two parity components are each Gaussian and pass the Wick test,
so the |0…0⟩ quench is "two Gaussian calculations plus a cross term", still polynomial — but the sentence at md:337 tells the
reader that a non-zero difference means the 2n×2n method cannot be used. The CUDA-Q snippet of §6 also starts from |0…0⟩ and
md:384 says "이 회로는 자유 페르미온이라 … 고전으로 쉽다" without this qualification.
What is delivered exactly as ch17 says: n = 1000, t = 40, "약 10초" (9.9 s), 2000×2000. I reproduced S(2), S(10), S(40) =
2.0015, 9.3443, 36.8963 with an independent implementation (scipy `expm` + repeated squaring instead of `eigh(iA)`).

Proposed insertion after md:318 (end of the n = 12 check paragraph):
```
17장의 급랭은 |0…0⟩에서 시작했는데 여기서는 |+…+⟩로 바꿨다. 끈을 X로 잡은 Jordan–Wigner에서 |0…0⟩은 패리티 ∏X의 짝·홀 구역에 절반씩 걸친 상태 (|GHZ₊⟩ + |GHZ₋⟩)/√2, |GHZ±⟩ = (|0…0⟩ ± |1…1⟩)/√2 이고, ⟨c₀⟩ = ⟨Z₀⟩ = 1 ≠ 0이라 가우스 상태 하나가 아니기 때문이다. 그래도 고전으로 쉽다. |GHZ±⟩ 각각은 가우스 상태이고 H가 패리티를 보존하므로, 가우스 계산 두 번과 두 구역 사이의 교차항으로 풀린다. |+…+⟩는 처음부터 가우스 상태 하나라서 식이 가장 간단하다.
```
md:337, last sentence:
```
차이가 0이 아니면 그 상태는 가우스 상태 하나가 아니고, 3절의 Γ(t) 계산을 그대로 쓸 수 없다. 그렇다고 곧 어렵다는 뜻은 아니다. 가우스 상태 몇 개의 중첩(위의 |0…0⟩ 급랭)은 여전히 다항 시간에 풀린다.
```
Suggested wording for the pointer in ch17:150 (not my file to edit):
```
(18장 5절: 같은 모형을 |+…+⟩에서 급랭해 n = 1000, t = 40까지 약 10초)
```
If the author wants a number in the inserted paragraph, add the n = 12 |0…0⟩ row to D4 of the script (true 1.4980, Gaussian 2.4979 at t = 2).
Note on my own proposal: that the cross term between the two sectors is polynomial-time is the standard result (matrix elements of
Gaussian operators between Gaussian states are Pfaffians); I verified only that each GHZ± component is Gaussian, not an
implementation of the cross term. If the author prefers not to assert it, end the inserted paragraph after "가우스 상태 하나가 아니기 때문이다." and
add "짝·홀 구역의 성분 |GHZ±⟩ 각각은 가우스 상태다(n = 12에서 Wick 검사 통과)."

### [B5] ch18 §1 message and §7 first paragraph — md:20, 394, 429
Current text:
> (md:20) `TFIM의 변분 탐색이 잘 되는 까닭과 그 동역학이 고전으로 쉬운 까닭은 같은 구조(자유 페르미온)다.`
> (md:394) `**잘 훈련되는 이유가 고전으로 쉬운 이유다.** 이 장의 변분 탐색은 TFIM에서 잘 돌았다. … 우연이 아니다.`

What is overstated (see also the closing judgement).
- The variational search of this chapter is ADAPT with the pool {Y_iZ_j}, and that tool is not free-fermionic. Majorana degree of
  Y_iZ_j (X string) is 2 only for |i − j| = 1; it is 4 for |i − j| = 2 and 6 for |i − j| = 3 (`c4_ff.py`, using the book's
  `majorana_degree`, itself verified against explicit Majorana operators on all 256 Pauli strings of n = 4 for all three axes).
  Four of the seven gates chosen at n = 4 (Y₃Z₁, Y₀Z₂, Y₀Z₃, Y₁Z₃) are such gates. The dynamical Lie algebra of the pool has
  dimension 12 (n = 3) and 56 (n = 4) = 2·dim so(2ⁿ⁻¹) — exponential if that pattern continues (I computed n = 3, 4 only) — against 9 and 16 = n² for the HVA generators
  (`c8_misc.py`). ADAPT's success here is not explained by free-fermion structure of its circuit, and it is outside the
  "provably no barren plateau" class the Cerezo argument is about.
- The free-fermion argument of md:394 is about HVA, and HVA is the method that did **not** train well in this chapter: below
  p = n²/4, about 1 % of random starts get near the best value (B1). "No barren plateau" is not "trains well".
- ADAPT itself did not "잘 돌았다" at n = 8 (60 gates, 4.8×10⁻⁵).
- What the chapter does show is weaker and still worth saying: the *target problem* (TFIM ground state and dynamics) is
  classically easy, so success on it is no evidence of advantage.

Proposed replacements.
md:20, second and third sentences:
```
이 장의 변분 탐색이 연습 문제로 삼은 TFIM은 자유 페르미온 모형이라 2n×2n 행렬로 풀린다. 구조를 아는 틀 HVA는 바로 그 구조(Majorana 2차식) 덕분에 불모 고원이 없고, 같은 구조 때문에 고전으로 시뮬레이션된다. 탐색이 성공했다는 사실만으로는 양자 이득의 증거가 되지 않는다. 그래서 탐색 도구와 반례 검사는 한 쌍으로 쓴다.
```
md:394:
```
**불모 고원이 없는 이유가 고전으로 쉬운 이유일 수 있다.** 이 장의 변분 탐색은 TFIM 바닥상태를 표적으로 삼았다. 그런데 5절 (D)는 같은 모형이 자유 페르미온이라 2n×2n 행렬로 풀린다고 말한다. HVA에서는 둘이 같은 사실이다. HVA의 게이트(ZZ 회전과 X 회전)는 모두 Majorana 2차식의 지수라서 회로 전체가 자유 페르미온 회로다. 동역학 리 대수가 다항 크기(n = 4에서 16차원)여서 불모 고원이 없고(13장), 같은 이유로 고전 시뮬레이션이 된다. 단 불모 고원이 없다는 것과 잘 훈련된다는 것은 다르다. 5절의 HVA는 층이 모자랄 때 무작위 시작점 대부분이 국소 최소에 빠졌다. ADAPT는 사정이 다르다. 풀 {Y\_iZ\_j}의 먼 쌍은 Majorana 4차 이상이라 자유 페르미온 회로가 아니고, 풀의 리 대수도 HVA보다 훨씬 크다(n = 4에서 56차원 대 16차원). ADAPT의 성공은 "회로가 고전으로 쉽다"로는 설명되지 않는다. 쉬운 것은 표적 문제(TFIM 바닥상태) 쪽이다. Cerezo et al.(2023)은 앞의 경우가 일반적 경향이라고 본다. 불모 고원이 없다고 증명된 회로족은 대개 그 증명에 쓴 구조로 손실함수를 고전에서 계산할 수 있다는 것이다(경우에 따라 양자 장치에서 초기 데이터를 한 번 모아야 한다). 어느 쪽이든 변분 탐색이 성공했다는 사실은 양자 이득의 증거가 아니다. 휴리스틱 트랙의 5′ 단계가 그래서 있다.
```
(The DLA dimensions 16 and 56 are mine; add a `dla_dim` call to the script if they are to be printed. ch13 already has `dla_dim`.)
md:429 may stay; its "깨지는 곳" already says "경향이지 정리가 아니다".

### [B6] ch18 §5(C) 추천 질의 — md:294; py:433–446
Current text:
> `추천 질의도 표본만으로 답했다. … 질의 하나가 160,400개를 읽었다. 사용자 20명의 상대 오차는 평균 0.17 …, 상위 10개 추천의 일치율은 평균 0.81이다. 양자 추천 알고리즘이 내놓는 것(저계수 근사 행에서의 표본)과 같은 종류의 출력이다.`

What is misleading. The 160,400 counted reads produce only the five coefficients λ_t. The quantities reported (relative error of
the whole row, top-10 overlap) are computed from `xhat = A[idx].T @ (coef @ lam)` (py:440), which reads 400 × 4000 = 1.6×10⁶
further entries (10 % of the matrix) outside the counter; the code comment says "검증용". The inner-product estimator itself is
correct and unbiased (checked: mean −0.00902 vs exact −0.00899, std as predicted). But no sample from x̂ is ever drawn, so the
output is not "the same kind" as the quantum algorithm's, and a top-10 list is not obtainable from 160,400 reads.

Proposed replacement for md:294:
```
추천 질의의 계수도 표본만으로 구했다. 4000×4000, r = 400에서 부분공간 단계가 성분 160,000개(전체의 1%)를, 질의 하나가 160,400개를 읽어 사용자 i의 계수-5 근사 행 x̂ = Σ\_t λ\_t v̂\_t의 계수 λ\_t 다섯 개를 얻는다. x̂의 성분 하나는 여기에 400개를 더 읽으면 나온다. 품질은 검증용으로 x̂ 전체를 만들어 쟀다(성분 1.6×10⁶개를 더 읽으므로 위 개수에는 넣지 않았다). 사용자 20명의 상대 오차는 평균 0.17(내적을 정확히 계산하면 0.14), 상위 10개 추천의 일치율은 평균 0.81이다. 양자 추천 알고리즘과 Tang의 알고리즘은 x̂ 전체가 아니라 x̂에서 뽑은 표본 하나를 내놓는다. 그 마지막 단계(기각 표본 추출)는 여기서 구현하지 않았다.
```

### [B7] ch18 §5(C) — md:305 "조건수도 정해져 있다"
Current text: `3중 대각 행렬은 행마다 성분이 3개뿐인 희소 행렬이고 조건수도 정해져 있다. … "희소 접근 + 높은 계수"(HHL의 원래 설정, 20장)는 남는다.`
The matrix (py:458, tridiag(−1, 2, −1), N = 1000) has condition number 4.06×10⁵ ≈ (2(N+1)/π)² (`np.linalg.cond`), growing as N².
HHL's cost is linear in κ at best (ch17:195), so this particular matrix is an example of where sampling fails, not of where a
quantum speedup survives; §7 (md:400) correctly says "잘 조건화된 행렬". Proposed replacement of the first sentence:
```
3중 대각 행렬은 행마다 성분이 3개뿐인 희소 행렬이다. 다만 조건수가 N²에 비례해(N = 1000에서 4×10⁵) 양자 선형계 알고리즘에도 유리한 예는 아니다. 여기서는 "표본이 통하지 않는 희소·최대 계수 입력"의 예로만 쓴다.
```

### [B8] ch18 §7 and §4 — Cerezo et al. caveat omitted — md:394, 141
The paper (arXiv:2312.09121 abstract; Nature Communications 16, 7907 (2025)) says such models can be classically simulated
"provided that one can collect some classical data from quantum devices during an initial data acquisition phase", and lists
caveats (smart initialisation, models outside its assumptions). ch18 states only "고전에서 계산할 수 있다". The parenthesis added in
the B5 replacement covers md:394. For md:141 (5′ 탈락 신호):
```
같은 손실함수를 고전으로(또는 양자 장치에서 초기 데이터만 한 번 모은 뒤 고전으로) 계산할 수 있음(Cerezo et al. 2023)
```

### [B9] ch18 §4 휴리스틱 트랙 vs ch15:43 and ch17:239 — md:135–142
- ch15:43 sends the reader to "커널 집중 — 18장 휴리스틱 트랙". ch18 only names "커널 집중" (md:141) and defers it to ch22; the term is
  defined nowhere in ch01–18 and has no English gloss.
- ch17:239 describes the track as three steps "(훈련 가능성 → 고전 대리 계산 검사 → 같은 예산의 고전 휴리스틱과 실측 비교)" applied by
  "21장"; ch18 has four steps (4′ 측정 예산 is new) and says "19장이 … 적용한다".
Proposed addition after the track table (md:142):
```
커널 집중(kernel concentration)은 큐비트 수가 늘수록 서로 다른 두 입력의 커널 값이 한 상수로 몰려, 그 차이를 읽는 데 필요한 샷 수가 지수로 느는 현상이다. 불모 고원의 커널 판에 해당한다. 검사는 n을 키우며 커널 값의 분산을 재는 것이고, 22장에서 한다.
```
(I did not verify a source for this definition today; it is the standard usage of Thanasilp et al., "Exponential concentration in
quantum kernel methods" — unverified citation, add only after checking.) And ask the author to align ch17:239 with ch18 (four steps;
"19장이 평가하고 21·22·23장이 실행").

---

## C — minor

- **[C1] md:384** `이 검사가 놀치는 것이 5절 (D)다.` → `이 검사가 놓치는 것이 5절 (D)다.` (garbled syllable). md:300–301 `쌌` → see A1.
- **[C2] md:14** `CCZ(Toffoli와 H 하나 차이)` — two H gates (target conjugation), as md:233 says. → `CCZ(표적 큐비트 앞뒤에 H를 붙이면 Toffoli)`.
- **[C3] md:45, 246, 437** "QSearch" is not the name used in arXiv:1912.02727 (it calls it "our synthesis algorithm … inspired by A*";
  fetched today). The name and the BQSKit pass cite the published version: Davis et al., "Towards Optimal Topology Aware Quantum
  Circuit Synthesis", IEEE QCE 2020 (BQSKit `QSearchSynthesisPass` docs). The Toffoli numbers are right: "the optimal CNOT count for
  Toffoli is six, which our algorithm matches", "Our linear topology Toffoli contains only eight CNOTs" (found, not proven optimal).
  §10 entry: append `출판판: "Towards Optimal Topology Aware Quantum Circuit Synthesis", IEEE QCE 2020(여기서 QSearch라는 이름을 쓴다).`
- **[C4] §10** (a) Valiant 2002 is cited at md:107 but is not in §10 — add `Valiant, "Quantum circuits that can be simulated classically in polynomial time", *SIAM J. Comput.* 31, 1229 (2002).` (reference from memory, not re-fetched). (b) Bravyi–DiVincenzo–Oliveira–Terhal 2008 is listed (md:448) but never cited in the text; cite it where "stoquastic" first appears (md:129 or 339). (c) md:441: give the full reference `*Nature Communications* 16, 7907 (2025)`. (d) md:435: the arXiv title is longer ("… and Timing Hamiltonians to Compute Controlled-Not Gates"); the PRA title as printed is fine; add the link `https://arxiv.org/abs/quant-ph/0308045` like the neighbours.
- **[C5] md:404** `Shor 알고리즘의 게이트 수를 O(n^1.5)로 줄였고` — Regev's abstract: a circuit of Õ(n^{3/2}) gates that is run √n + 4 times. → `회로 한 번의 게이트 수를 Õ(n^1.5)로 줄였고(대신 회로를 √n + 4번 돌린다; 11장)`. (ch11:227 has the same shorthand.)
- **[C6] md:88** `현재 가장 좋은 고전 표본 알고리즘` is time-dependent and I did not check for improvements after 2023 → `2023년 기준 가장 좋은`. The formulas themselves match the abstract of arXiv:2303.01492 exactly (classical Õ(d¹¹‖A‖_F⁴/(ε²‖A‖⁴)), previous d²²‖A‖_F⁶/(ε⁶‖A‖⁶), QSVT Õ(d‖A‖_F/(ε‖A‖)) for estimating one entry of p(A)b after QRAM loading).
- **[C7] md:33** `계산 기저뿐 아니라 모든 상태를 바꾸므로 유니터리로서 SWAP이다.` — one superposition input with the lower wire at |0⟩ does not show that. → `계산 기저 상태뿐 아니라 중첩 상태도 옮긴다. 네 기저 상태에서 모두 맞고 선형이므로 유니터리로서 SWAP이다.` Both Quirk links decode to what the text says (cols `["X"],["•","X"],["X","•"],["•","X"],["Chance2"]` and `["H"],["Z^¼"],["•","X"],["X","•"],["•","X"],["Bloch","Bloch"]`).
- **[C8] md:122, 348** "15장 3절의 경계마다" — the dequantization row comes from ch15 §7 and the 가속 분류 table, not from the §3 boundary list (which says "이 네 가지"); fixed in the B3 text. md:348 `첫째·둘째 줄` → `둘째 줄` (ch17 drops 1D TFIM by DMRG, not by QMC).
- **[C9] terms (rule: 한국어(English) on first use)** no English gloss for 대칭 장벽 (symmetry roadblock), 연산자 풀 (operator pool), 평균장 (mean field), 공분산 행렬 (covariance matrix), 가우스 상태 (Gaussian state), 전진 선택법 (forward selection), 공정 비충실도 (process infidelity), 좌절 (frustration), 커널 집중 (kernel concentration). md:394 "동역학적 Lie 대수" vs ch13 "동역학 리 대수(DLA)" — use ch13's form.
- **[C10] figure** panel (3): `ax.set_ylim(1e-11, 5)` (py:678) hides exactly the HVA points the text is about (n = 4, p = 4 at 2e-14; n = 6, p = 9 at 1e-13); clip them to the floor as the ADAPT curves effectively are, e.g. plot `max(err, 2e-11)`. The x-axis "파라미터(게이트) 수" mixes HVA parameters with ADAPT gates. Panel (6): the label "MPS라면 χ ≥ 2^S ≈ 1e+11" is drawn across the curve near t = 40.
- **[C11] numbers not printed by the script** md:220 `0.0000은 10⁻¹⁵ 이하다` (only visible in figure (1); stdout prints 4 decimals); md:418 `c = 0.65` attributed to §5 (see A2); md:146/149 "약 1초" — the cell took 2.1 s here including imports (output 0.7500, 0.7500, 0.5000, 0.0000 confirmed).
- **[C12] md:214** `drain_pool=True는 게이트 이름과 와이어만 비교해 풀에서 빼므로, 이름이 모두 PauliRot인 풀에서는 False로 둔다.` The first half is right (source: `g.name != op.name or g.wires != op.wires`). But in this pool every element has a distinct ordered wire pair, so name + wires does identify it; the actual reason to keep `False` is different. → `drain_pool=True는 한 번 고른 게이트를 풀에서 빼는데, ADAPT는 같은 연산자를 다시 고를 수 있어야 하므로 False로 둔다. 또 비교가 이름과 와이어만으로 이뤄져서, 같은 와이어에 파울리 문자열이 다른 PauliRot(예: YZ와 YX)을 함께 둔 풀에서는 하나를 고르면 둘 다 빠진다.`
- **[C13] md:384** the CUDA-Q MPS target also truncates by `CUDAQ_MPS_ABS_CUTOFF` and `CUDAQ_MPS_RELATIVE_CUTOFF` (both default 1e-5), so "χ 상한을 늘려도 변하지 않으면" should mention them: add `(특이값 절단 CUDAQ_MPS_ABS_CUTOFF·CUDAQ_MPS_RELATIVE_CUTOFF의 기본값 10⁻⁵도 함께 작용한다)`. The docs do not say the variable must be set before `set_target`; setting it before `import cudaq` as the snippet does is the safe choice — the comment "타깃 설정 전에 둔다" is an assumption.
- **[C14] md:109, 410** N&C 4.5.4 is a covering/counting argument (ch04:201 calls it "세는 논법"), not parameter (dimension) counting → `세는 논법으로 보인다(N&C 4.5.4)`; exercise 1's title likewise. Exercise mix is N&C + 자작 + 자작 (no KLM-derived exercise; same deviation R4-C4 recorded for other chapters).
- **[C15] md:107** `세로장 Z_j나 다음 이웃 결합 Z_jZ_{j+2}는 Majorana 차수가 2를 넘어` — Z_j has degree 2j + 1 (Z₀ has degree 1, an odd term); → `세로장 Z_j(차수 2j + 1, 홀수)나 다음 이웃 결합 Z_jZ_{j+2}(차수 4)는 2차식이 아니어서`.
- **[C16] md:16, 276** `시작 상태의 에너지` — 0.759 is the energy error E(|+…+⟩) − E₀ = −4 + 4.7588, not the energy. → `시작 상태의 에너지 오차`.
- **[C17] md:268** the evaluation-count comparison (60,622 vs 2,091) leaves out ADAPT's 960 pool-gradient measurements listed two tables earlier and sets nine depths × 21 starts against one ADAPT run; say so or drop the comparison. (No replacement proposed; it is folded into the B1 text if that is adopted.)

---

## Overall judgement on the central message

"도구가 '된다'고 답한 이유가 곧 고전으로도 되는 이유일 수 있다" — the hedged sentence is sound. Two stronger statements built on it
are not supported by what the chapter ran. (1) §1's "TFIM의 변분 탐색이 잘 되는 까닭 … 같은 구조(자유 페르미온)" and §7's unhedged
header: the search tool actually demonstrated, ADAPT with {Y_iZ_j}, uses non-free-fermion gates and has a much larger Lie
algebra (56 vs 16 at n = 4), so its success is not explained by, and does not illustrate, the free-fermion/small-DLA argument; that
argument applies to HVA, which in this chapter trained badly from random starts and whose sub-threshold numbers are optimiser-limited.
(2) Cerezo et al. is quoted without its data-acquisition caveat. What the chapter does establish is the weaker, correct point: the
benchmark problem is classically easy (verified: n = 1000 by 2000×2000 matrices), so success on it proves nothing about advantage.
Separately, the free-fermion counterexample is demonstrated for a different initial state than the one ch17 points to (B4), and the
"stoquastic" row of the screening table is applied more broadly than it is valid (B3).

## Summary (most important first)

1. A1 — the stable-rank table compares sampled rows with matrix entries; the sampling method already loses at stable rank 40, not 300.
2. A3 — "대칭 장벽" is attached to the wrong trap; Shkolnikov et al.'s roadblock is the chapter's first trap (symmetry-breaking pool operators have zero gradient).
3. B1 — HVA: thresholds p = 4 and p = 9 are real; the n = 6 sub-threshold errors are 2–40× too high (too few starts; the warm start is a saddle point).
4. B4/B5 — ch18 solves a different quench than ch17 cites (|0…0⟩ is not Gaussian in this JW and fails ch18's own Wick test), and the central claim is overstated for ADAPT.
5. A2/B3 — exercise 3 counts r² reads on a 192-column matrix; "하나라도 걸리면 탈락" contradicts ch17's surviving stoquastic dynamics candidate.

## Checked and found correct

- **Numbers**: every number in §1, all §5 tables and sentences, and the §8 answers against `ch18_stdout.txt` (exceptions: A4, C11, C16). 43 = 1+1+2+4+10+25 classes; (2−√2)/4 = 0.1464; 4·sin(π/18) = 0.6946; pool-measurement counts 96/960/3416 = (gates+1) × pool size; 1/21, 5/21, 15/21, 8/21 fractions; 2^S column; r²/mn = 10.24, 0.64, 0.04; exercise 1(b),(c), exercise 2.
- **Analytic gradient** of `Template.cost_grad` (own builder, 5e-12) and of `energy_grad`/the ADAPT formula g_P = (i/2)⟨[P,H]⟩ (−0.52538437 both ways).
- **γ(U) criterion** as implemented and as tabulated in §3, on 21 gates with known cost (I, local, CNOT, CZ, exp(−iπ/4 XX), iSWAP, √iSWAP, DCNOT, CP(φ), B gate, generic 2-CNOT product, SWAP, √SWAP, generic canonical gate, Haar, QFT₂): criterion = known cost = numerically found minimum in all 21. It agrees with Props. III.1–III.3 of quant-ph/0308045 (χ[γ] = (x±1)⁴; (x+i)²(x−i)²; tr γ real); γ is unitary, so the book's "γ = ±I", "tr γ = 0 and γ² = −I" are equivalent forms; the det^{1/4} branch only flips the sign of γ.
- **CNOT direction** is absorbed (H⊗H conjugation; three direction choices give identical minima on 4 random targets). **Parameter count**: Jacobian rank of the 2-qubit template is 6, 10, 14, 15, 15 for k = 0…4 = min(15, 6 + 4k); bound ⌈(4ⁿ−3n−1)/4⌉ = 3, 14, 61, 252 and the Rz-control/Rx-target argument match Prop. III.1 of quant-ph/0308033 ("almost all n-qubit operators").
- **Symmetry reduction** in the CCZ search is valid (see B2); floors 0.4375/0.375/0.1464 reproduced with 16 starts; 01·01·02·12·02·12 is the reversed N&C Fig. 4.9 order.
- **Traps**: E(|+…+⟩) − E₀ = 0.7588, 1.2962; best product state (general complex product state, 30 starts) 0.3551, 0.3937; sector minima by projection give (E₁−E₀)/2 = 0.3473, 0.2411, 0.1845 (n = 8) and gap = 4·sin(π/(4n+2)) for n = 4, 6, 8 (and the lowest mode at n = 200); sector weights of |0…0⟩ are 1/2.
- **FKV**: E[SᵀS] = AᵀA (0.3 % Monte-Carlo error); E‖AᵀA − SᵀS‖_F² = (‖A‖_F⁴ − ‖AᵀA‖_F²)/r ≤ ‖A‖_F⁴/r (0.8314 vs 0.8305 vs bound 1.048); the column mixture equals norm² sampling of the columns of S because all rows of S have norm² ‖A‖_F²/r; the rescaling in `fkv` gives E[WWᵀ] = SSᵀ; the inner-product estimator is unbiased; read counts r·c and 160,400.
- **Free fermions**: Majorana algebra and the identities of §3(라); `majorana_degree` = brute-force expansion on all 256 strings of n = 4 for axes X, Y, Z; H = (i/4)ΣA c c equals the TFIM matrix; i[H, c_a] = ΣA_ab c_b; Γ(t) and the entropy formula against a state vector at n = 12 (6 digits, independent `expm`); n = 1000 values; the D1 table; Wick-test numbers; stoquastic verdicts.
- **Literature** (fetched today): quant-ph/0308033, quant-ph/0308045, 0803.2316 (first proof that six CNOTs are necessary), 1912.02727 (authors, 6/8), 1812.11173 and 1911.10205 (titles, authors; "minimal pool size … linear"), 2109.05340 (title, authors, 2n−2; Quantum 7, 1040 via the journal URL), 2312.09121 (12 authors, DOI, Nat. Commun. 16, 7907), 2303.01492, 1807.04271, 1910.06151 (six authors), Arrazola et al. Quantum 4, 307, Chapman–Flammia Quantum 4, 278, Marvian–Lidar–Hen (Nat. Commun., April 2019), Bravyi et al. QIC 8, 361–385, FKV 2004, Regev 2308.06572.
- **SDK** (official docs, PennyLane 0.45.1 pages; CUDA-Q "latest"): `two_qubit_decomposition(U, wires)` returns a list of operations and chooses 0–3 CNOTs by the same γ criterion; `AdaptiveOptimizer(param_steps=10, stepsize=0.5)`; `step_and_cost(circuit, operator_pool, drain_pool=False, params_zero=True)` returns (circuit, cost before the step, largest gradient); only the newly selected gate's parameter is optimised (source code); `lie_closure` exists; `tensornet-mps`, `CUDAQ_MPS_MAX_BOND` default 64, NVIDIA GPU required; Stim backend listed; BQSKit `QSearchSynthesisPass`; Qiskit `UnitarySynthesis` (current API page exists). Signs in the CUDA-Q kernel: exp(+i dt ZZ) = CX·Rz(−2dt)·CX and exp(+i dt X) = Rx(−2dt) are right for H = −ΣZZ − ΣX, consistent with ch13:82; `int(t/dt)` gives 20, 40, 60.
- **Cross-chapter**: TFIM convention and E₀ (−4.758770, −7.296230) as in ch12/13/17; HVA n = 4 (p = 3: 1.56×10⁻⁴; CNOT 24, 8 parameters) as in ch13; ch14:43's 8-CX linear Toffoli; ch17's χ = 24 at t = 3 and rate 0.91; ch17's "n = 1000, t = 40까지 약 10초" (9.9 s); promises at ch06:211, ch11:227, ch15:60, ch15:111, ch16:315, ch17:13, ch17:52, ch17:192 are kept (ch15:43 and ch17:239 only partly — B9; ch17:150 — B4).
- **Template**: 10 sections; 교과서 비교 in §3; three analogies in §2 and six notes in §9, each with a "깨지는 곳"; three exercises with answers; one other SDK with a stated reason; no "…" elisions in code cells; no table with a wrong cell count; versions in §10 match this container (Python 3.13.17, numpy 2.5.3, scipy 1.18.1, matplotlib 3.11.2).

## Could not check

- Running the PennyLane and CUDA-Q snippets (no PyPI); whether `pennylane==0.45.x`/`cudaq==0.16.x` are the right pins beyond the docs pages seen.
- N&C §4.5/Fig. 4.9 and KLM ch. 4 (no access to the books); the plan document ("작업계획서 후보 2번", PennyLane as this chapter's planned tool).
- Journal details not opened today (consistent with what I know, not re-fetched): Grimsley et al. Nat. Commun. 10, 3007; PRX Quantum 2, 020310; Shende–Markov QIC 9, 461; HHL PRL 103, 150502; Bittel–Kliesch PRL 127, 120502; Gidney 2505.15917; Kerenidis–Prakash ITCS 2017; STOC venues of Tang and Chia et al.; Lieb–Schultz–Mattis; Pfeuty; Terhal–DiVincenzo; Valiant. The body of Arrazola et al. (only the abstract: good only for low rank, low condition number, very large dimension — the "직접 해법보다 느리다" wording at md:400 goes slightly beyond the abstract).
- The remainder of Shkolnikov et al. beyond §4 (fetch truncated), and whether a better classical bound than Bakshi–Tang has appeared since 2023.
- Whether HVA at n = 6, p = 8 can reach zero with still more starts (350 starts: best 3.7×10⁻⁶; parameter counting says no), and whether a 7-CNOT linear CCZ exists (40 starts per class found none; no proof known to me).
- One of my own runs (`c3_adapt.py`) was stopped by the 20-minute background limit before its n = 6, p = 8/9 part; those two cases were rerun separately with fewer starts (`c11_out.txt`).
