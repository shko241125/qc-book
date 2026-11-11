# R2 findings — exercises ch10–ch17 (24) + second review of ch17

Reviewer R2. Scratch code and outputs: `/tmp/claude-0/-home-claude/fd63146d-511d-5ae9-8dab-3100fb85f0db/scratchpad/mid/R2/`
(`ex11.py`, `ex13.py`, `ex14.py`, `ex16.py`, `ch17_indep.py`, `ch17.out` = fresh stdout of a copy of `ch17_discovery.py`).
Nothing under `/home/claude/qc-book/` was modified.

Counts: **A = 1, B = 3, C = 10.**

---

## A — factually / numerically wrong

### A1. ch11 §8 문제 2 (and §3, §1, script) — the "≥ 3/4 for N = pq" bound is false; the chapter's own N = 21 violates it
- `ch11/ch11_shor.md:241` (answer to 문제 2), `ch11/ch11_shor.md:29` (§3), `ch11/ch11_shor.md:13` (§1), `ch11/ch11_shor.py:60` (print string).
- Current text (241): `N&C 정리 5.3은 "r이 짝수이고 a^{r/2} ≢ −1"의 확률을 ≥ 1 − 1/2^m(N = pq에서 3/4)으로 주며, 이것이 Shor 원논문의 "성공 ≥ 1/2"보다 강한 근거.`
  (29): `확률은 ≥ 1 − 1/2^m (m = N의 서로 다른 소인수 개수; N&C 정리 5.3 = A4.13), N = pq면 ≥ 3/4 — Shor 원논문의 느슨한 판은 ≥ 1/2.`
  (13): `하한 3/4`
- What is wrong: 1 − 1/2^m is a known misprint in N&C. The correct bound is 1 − 1/2^{m−1} (Shor), i.e. 1/2 for N = pq, and it is tight. Shor's bound is not the "loose version"; it is the right one.
- Evidence:
  - Enumeration (`R2/ex11.py`), fraction of a coprime to N with r even and a^{r/2} ≢ −1: N = 15: 6/8 = 0.75; **N = 21: 6/12 = 0.50**; N = 33: 10/20; N = 57: 18/36; N = 77: 30/60 (all 0.50 < 0.75). For N = 21 the failing a are 1, 4, 16 (odd r) and 5, 17, 20 (a^{r/2} ≡ −1). N = 21 is the chapter's own second example.
  - Analytic reason: failure ⇔ the 2-adic valuations of ord_p(a) and ord_q(a) are equal; for p = 3, q = 7 each valuation is 0 or 1 with probability 1/2 → failure 1/2.
  - Official N&C errata: "Theorem 5.3 should be amended so that the lower bound in Eq. (5.60) is 1 − 1/2^{m−1}", and the same for Theorem A4.13 / Eq. (A4.32)–(A4.33). https://michaelnielsen.org/qcqi/errata/errata/errata.html
- Proposed replacement (241, from "−1일 확률이 ≤ 1/2 —" up to "하한은 하한일 뿐."):
  > 정확히 세면, a mod p와 a mod q의 위수에 들어 있는 2의 지수가 같을 때만 실패한다(둘 다 0이면 r이 홀수, 같은 양수이면 a^{r/2} ≡ −1). 그 확률은 ≤ 1/2이므로 "r이 짝수이고 a^{r/2} ≢ −1"의 확률은 ≥ 1 − 1/2^{m−1}(N = pq에서 1/2)이다 — Shor 원논문의 하한이다. N&C 정리 5.3의 1 − 1/2^m은 인쇄 오류이고 공식 정오표가 1 − 1/2^{m−1}로 고쳤다. 이 하한은 실제로 닿는다: N = 21, 33, 57, 77에서는 서로소 a의 정확히 절반이 실패한다(N = 21: 12개 중 a = 1, 4, 5, 16, 17, 20). N = 15는 p − 1 = 2와 q − 1 = 4의 2의 지수가 달라 8개 중 2개(a = 1, 14)만 실패한다 — a = 1을 빼면 7개 중 1개.
- Proposed replacement (29):
  > 난수 a에 대해 이 조건이 성립할 확률은 ≥ 1 − 1/2^{m−1} (m = N의 서로 다른 소인수 개수; Shor 1997. N&C 정리 5.3 = A4.13에는 1 − 1/2^m로 인쇄되어 있으나 공식 정오표가 1 − 1/2^{m−1}로 고쳤다), N = pq면 ≥ 1/2. 이 하한은 N = 21에서 실제로 닿는다(서로소 a 12개 중 6개만 성공). 연습문제 2가 실패 예를 다룬다.
- Proposed replacement (13): `하한 3/4` → `일반 하한은 1/2`
- Proposed replacement (py:60 string): `(일반 하한 1/2 = 1 − 1/2^(m−1): Shor 1997, N&C 정리 5.3 정오표)`

---

## B — misleading / incomplete / contradicts another chapter

### B1. ch16 §8 문제 3 — the closed form for the MSE has the wrong constant; the exact result is simpler and reproduces the script
- `ch16/ch16_primitives.md:341`
- Current text: `따라서 MSE ≈ Σ_{1≤|k|≤M/2} [sin²(πf)/(π²k²)]·(π k sin 2θ/M)² = sin²(πf)·sin²2θ/M … (k가 M에 가까우면 선형 근사가 깨지지만 차수는 같다).` followed by `a를 (0, 1)에 고르게 평균하면 RMSE·√M = 0.50~0.52`.
- What is wrong: the linearised formula does not reproduce the numbers quoted next to it. Exact identity (no approximation): MSE = sin²(Mθ)/(2M) = sin²(πf)/(2M), independent of θ.
  Derivation: sin²A − sin²θ = sin(A+θ)sin(A−θ) with A = πy/M; one peak has P(y) = sin²(Mθ)/(2M² sin²(A−θ)), so sin²(A−θ) cancels exactly and Σ_y sin²(A±θ) = M/2.
- Evidence (`R2/ex16.py`, exact QAE distribution, a = 0.3):

  | m | exact RMSE | √(sin²πf/(2M)) | book formula √(sin²πf·sin²2θ/M) |
  | --- | --- | --- | --- |
  | 4 | 2.651e-2 | 2.651e-2 | 3.436e-2 |
  | 8 | 2.957e-2 | 2.957e-2 | 3.832e-2 |
  | 12 | 8.225e-3 | 8.225e-3 | 1.066e-2 |

  The exact column equals the chapter's own stdout (2.65e-02, 2.96e-02, 8.22e-03). Averaged over a ∈ (0,1): exact RMSE·√M = 0.501 / 0.500 / 0.494 (m = 4/8/12) = √(1/4); the book formula predicts √(1/2 · 2/3) = 0.577, which does not match the 0.50–0.52 the answer cites as confirmation.
- Proposed replacement (from "k칸 떨어진 출력의 추정 오차는" through "RMSE·√M = 0.50\~0.52(m = 4\~12)로 일정하다."):
  > k칸 떨어진 출력의 추정 오차는 선형 근사로 Δa ≈ sin 2θ·πk/M이다. 꼬리 확률 ∝ 1/k²과 곱하면 각 칸이 k²를 상쇄해 칸 수(≈ M)만큼 쌓이므로 MSE ∝ 1/M이다. 선형 근사 없이도 닫힌 식이 나온다. A = πy/M으로 두면 sin²A − sin²θ = sin(A+θ)·sin(A−θ)이고, 봉우리 하나의 분포는 sin²(Mθ)/(2M² sin²(A−θ))이므로 분모가 정확히 상쇄된다. 거울상 봉우리도 같은 꼴이고 Σ\_y sin²(A±θ) = M/2이므로 MSE = sin²(Mθ)/(2M) = sin²(πf)/(2M) — θ에 무관하다. RMSE = |sin πf|/√(2M) ∝ 1/√M ∝ 1/√N, 몬테카를로와 같은 차수다. 5절의 a = 0.3에 넣으면 m = 4에서 2.65×10⁻², m = 12에서 8.22×10⁻³으로 표의 값과 일치한다. a를 (0, 1)에 고르게 평균하면 sin²(πf)의 평균이 1/2이라 RMSE·√M = 1/2이다(스크립트 0.50\~0.52, m = 4\~12).

### B2. Overlap symbol γ means a probability in ch10/ch12 and an amplitude in ch17
- `ch10/ch10_qpe.md:255` (`겹침이 γ = |c_0|²인 시험 상태 … 확률 γ … 1/γ번 반복 … 진폭 증폭으로 1/√γ번`), `ch12/ch12_hamiltonian.md:231` (`겹침 γ = |⟨ψ₀|φ⟩|²`), ch12 stdout (`겹침 γ=0.813 … 1/γ ≈ 1.2 회`)
  versus `ch17/ch17_discovery.md:59` (`QPE가 바닥상태를 뽑을 확률이 γ²이므로 반복 횟수는 1/γ²(진폭 증폭으로 1/γ)`), lines 101, 179 (`1/γ²`), 225.
- What is wrong: same symbol, same quantity name ("겹침"), two definitions; a reader carrying ch12's γ = 0.81 into ch17's 1/γ² gets the wrong repetition count. Both are internally correct.
- Proposed fix (smallest edit — one sentence appended in ch17:59 after "…(진폭 증폭으로 1/γ)."):
  > 주의: 10·12장에서는 겹침을 확률 |⟨ψ₀|φ⟩|²로 정의해 γ라 썼다. 이 장의 γ는 진폭 겹침 |⟨ψ₀|φ⟩|이고, 10·12장의 γ는 이 장의 γ²에 해당한다.

### B3. ch17 — "구조 없음" is both a stop signal and a row with a primitive; the SAT row fills steps 3–5 after failing at step 2
- `ch17/ch17_discovery.md:32` (step 2 탈락 신호 `구조 없음(→ 최선이 √N)`), `:27` (`앞 단계에서 떨어지면 뒤 단계는 하지 않는다`), `:50` (structure table row `구조 없는 탐색 | Grover / 진폭 증폭 | 제곱 | SAT, 무차별 대입`), `:182` (candidate row), `:185`.
- What is wrong (internal consistency, introduced by the post-review edit that moved SAT to "탈락(2단계)"):
  1. The SAT row still carries step 3 (Grover), step 4 (√N 호출), step 5 (SAT 솔버), and its verdict text gives a step-5 reason (`고전 솔버는 실제 입력의 구조를 씀`) together with a step-2 label.
  2. Derivative pricing has the same quadratic speedup, passes step 2 and is rejected at step 5. The text never says why "quadratic via amplitude estimation" passes step 2 while "quadratic via Grover" stops there.
  3. The structure→primitive table lists the stop signal as an ordinary structure with a primitive and an expected speedup.
- Proposed replacement (row 182):
  > | 구조 없는 SAT 탐색 | 구조 없음 | —(참고: Grover) | —(참고: √N 호출) | —(참고: 최신 SAT 솔버는 실제 입력의 구조를 씀) | 탈락(2단계) — 구조 없음, 최선이 √N. 2단계에서 멈추므로 뒤 칸은 참고용 |
- Proposed sentence to add after the five-step table (after line 35):
  > 제곱 가속이라는 사실만으로는 2단계에서 떨어지지 않는다. 확률·기댓값 구조(진폭 추정)는 2단계를 통과하고 5단계의 손익분기로 판정한다. "구조 없음"은 양자 쪽에 쓸 구조가 없는데 고전 쪽은 실제 입력의 구조를 쓸 수 있는 경우이므로 2단계에서 멈춘다.
- Proposed replacement (row 50, first cell): `구조 없는 탐색` → `구조 없는 탐색(2단계 탈락 신호 — 참고)`

---

## C — minor

### C1. ch10 §8 문제 3 — question asks 4 bits and "how many shots", answer analyses 6 bits; "≈ 0.69" should be 0.68
- `ch10/ch10_qpe.md:257–259`. Current: `그래서 6비트 전부 맞을 확률은 ≈ 0.75 × 0.93 × 0.98 × … ≈ 0.69이고 … 더하면 실측 0.845가 된다.`
- Evidence: 0.75 × 0.9330 × 0.9830 × 0.9957 × 0.9989 × 0.9997 = 0.6840, identical to the §5 table value P(최빈) = 0.6840 at t = 6 (single-shot IPE samples the QPE distribution exactly). Exact P(오차 ≤ 1/64) = P(k=21) + P(k=22) = 0.6840 + 0.1710 = 0.855; the quoted 0.845 is the 200-trial sample (SE 0.025), fine but it is not "what you get by adding". Binomial values 0.10 / 0.006 and Hoeffding 0.535 / 0.072 are correct.
- Proposed replacement of that sentence:
  > 그래서 6비트 전부 맞을 확률은 0.75 × 0.933 × 0.983 × … = 0.684다 — 5절 표의 t = 6 최빈 확률과 같은 값이다(단일 샷 IPE는 QPE와 같은 분포에서 뽑는다). 최하위 비트를 틀려도 상위 보정이 흡수해 최종 오차 ≤ 1/64가 되는 경우(k = 22, 확률 0.171)를 더하면 0.855이고, 200회 표본의 실측은 0.845다.
- Proposed replacement of the question's last sentence: `φ = 1/3일 때 비트당 몇 샷이 필요한가?` → `φ = 1/3일 때(5절처럼 6비트) 비트당 샷 수에 따라 성공 확률이 어떻게 달라지는가?`

### C2. ch10 §8 문제 1 — grouping in the last inequality; attribution
- `ch10/ch10_qpe.md:251`. `P ≥ (2·2^tδ/π ÷ 2^tδ)²` gives the right value 4/π² but does not show the two bounds used (numerator ≥ 2·2^t|δ|, denominator ≤ 2^t·π|δ|).
- Proposed: `P ≥ (2·2^t|δ| ÷ (2^t·π|δ|))² = 4/π² ≈ 0.405`
- Attribution `(N&C 5.2.1)`: the geometric-sum derivation is N&C 5.2.1, but to my recollection the 4/π² statement itself is from Cleve–Ekert–Macchiavello–Mosca 1998 / KLM 7.1, not N&C. **Unverified suspicion** (no access to the book text). If confirmed: `(N&C 5.2.1의 유도, 4/π² 하한은 KLM 7.1)`.

### C3. ch11 §8 문제 3 — the reason given for the t = 5 failure is not the actual mechanism
- `ch11/ch11_shor.md:245`. Current: `t = 5(2^t = 32 < r² = 36)는 연분수가 r̂ = 6을 한 번도 못 줌 — 1/6과 1/5, 1/7을 32등분으로 구분할 수 없다.`
- Evidence: `Fraction(5, 32).limit_denominator(21) = 3/19`, `Fraction(27, 32).limit_denominator(21) = 16/19` (script uses `limit_denominator(N)`). The peak k = 5 is closer to 3/19 (0.0016) than to 1/6 (0.0104), so r̂ = 19. At t = 7, k = 21 → 1/6 and k = 107 → 5/6 (2 × 0.684/6 = 0.228 ✓ 23%).
- Proposed: `— 1/6과 1/5, 1/7을 32등분으로 구분할 수 없다.` → `— 봉우리 k = 5의 5/32는 1/6(차이 0.0104)보다 3/19(차이 0.0016)에 가까워 연분수가 r̂ = 19를 낸다. 격자 간격 1/32이 분모 ≤ 21인 분수들의 간격보다 넓기 때문이다.`

### C4. ch12 §8 문제 2 — operator order opposite to the script's circuit
- `ch12/ch12_hamiltonian.md:243`. `2항 분해 e^{−iH_{ZZ}dt}e^{−iH_X dt}`. The script (`ch12_hamiltonian.py:27–28, 47`) applies the ZZ terms first, so the step operator is e^{−iH_X dt}·e^{−iH_ZZ dt}. Verified the answer's numbers with this order: ‖[H_ZZ, H_X]‖ = 8.944, bound 0.559 at r = 32, actual error 0.0836.
- Proposed: `e^{−iH\_X dt}e^{−iH\_{ZZ}dt}(ZZ 항을 먼저 적용)`

### C5. ch13 §8 문제 1 — "유한차분은 쓸 수 없다" is an overstatement
- `ch13/ch13_variational.md:286`. Current: `샷이 있는 한 유한차분은 쓸 수 없다`. The algebra (a + b cos θ + c sin θ, shift rule, σ/(√2ε) vs σ/√2) is correct. With a larger step the finite-difference estimator trades bias for variance and is usable.
- Proposed: `작은 ε의 유한차분은 샷 잡음에 묻힌다. ε를 키우면 분산은 줄지만 편향 O(ε²)이 생긴다 — 이동 규칙은 편향 없이 σ/√2를 준다.`

### C6. ch13 §8 문제 3 — "거리 2p+1 이내"
- `ch13/ch13_variational.md:294` (same wording in `:48`, `:202`). Numbers verified with an independent implementation (`R2/ex13.py`): n = 12 ring, p = 1, 2, 3 → ⟨C⟩ = 9.000, 10.000, 10.500. The light cone of one edge is the set of vertices within distance p of its endpoints: 2p+2 vertices, 2p+1 edges.
- Proposed: `각 변에서 거리 2p+1 이내의 그래프만 "보고"` → `각 변의 양 끝에서 거리 p 이내(고리에서는 정점 2p+2개·변 2p+1개의 사슬)만 "보고"`

### C7. ch14 §8 문제 1 — the budget omits the 1-qubit-gate term that §5's budget column includes
- `ch14/ch14_hardware.md:190`. `오류 예산 1 − 0.997¹⁴⁴·0.987⁴ = 0.38`. Arithmetic is right (0.384), as are θ = 0.2527, k = 3, P = 0.961, CZ 144, n = 5 → 320 / 0.64, and 231. But §5's budget (e.g. Grover n = 3: 0.147 = 1 − 0.997³²·0.99975⁹⁵·0.987³) includes 1q gates. Scaling 95 1q gates per 4 CCZ to 18 Toffoli-class gates gives roughly 430–480 1q gates → budget ≈ 0.45 (my estimate, not a script run).
- Proposed: `오류 예산 1 − 0.997¹⁴⁴·0.987⁴ = 0.38` → `오류 예산(CZ와 판독만) 1 − 0.997¹⁴⁴·0.987⁴ = 0.38, 1큐비트 게이트 수백 개(×2.5×10⁻⁴)까지 넣으면 0.45 안팎`

### C8. ch17 §4 code excerpt is stale against the script and the text
- `ch17/ch17_discovery.md:84–87`. Excerpt: `def qpe_qubitization_cost(L, alpha, n_sys, eps, b_rot=10)` and `qubits=n_sys + 2*m + b_rot + 1`. Script (`ch17_discovery.py:85`) and §3 text use `mu` (keep 비트 μ).
- Proposed: `b_rot=10` → `mu=10`, `b_rot + 1` → `mu + 1`.

### C9. ch17 §3 — the formula attributed to Babbush et al. 2018
- `ch17/ch17_discovery.md:59`. `최적 창을 쓰는 문헌 추정은 N_W ≈ πα/(2ε)(Babbush et al. 2018)`.
- Evidence: the paper gives m = ⌈log(√2πλ/(2ΔE))⌉ (Eq. 24) and 2^m < √2πλ/ΔE (Eq. 26); πλ/(2ΔE) does not appear (ar5iv full text, https://ar5iv.arxiv.org/html/1805.03662). The √2 comes from splitting the error budget. "상수 2배 안팎 차이" remains true (book: πα/ε ≤ 2^b < 2πα/ε; paper: ≈ 2.22α/ε).
- Proposed: `N\_W ≈ √2πα/(2ε)(Babbush et al. 2018, 식 24; QPE 몫의 오차만 ε\_QPE로 두면 πα/(2ε\_QPE))`

### C10. ch17 §5 candidate table — linear-systems cost does not belong to the listed assembly
- `ch17/ch17_discovery.md:183`. `블록 인코딩 + QSVT(1/x) | κ·log(1/ε)에 비례`. The 1/x polynomial has degree O(κ log(κ/ε)) per run; with plain amplitude amplification of the ~1/κ success amplitude the total is O(κ² log(κ/ε)). O(κ log(1/ε)) is the optimal query complexity of a different construction (Costa et al. 2022). Stated from my knowledge of the literature; no source fetched in this session.
- Proposed cell: `QSVT 차수 κ·log(κ/ε), 진폭 증폭 포함 총 κ²·log(κ/ε)(최적 알고리즘은 κ·log(1/ε))`

---

## Summary
- One real error: ch11 repeats the N&C Theorem 5.3 misprint (1 − 1/2^m) and calls Shor's correct bound the loose one; the chapter's own N = 21 succeeds for exactly 1/2 of a, and the official errata confirm 1 − 1/2^{m−1}.
- ch16 exercise 3 gives a linearised MSE formula that is ~30% off the script values it cites; the exact identity MSE = sin²(πf)/(2M) reproduces them to four digits.
- ch17 reproduces completely: every number in §1, §5 and §8 matches stdout and my independent recomputation; all formulas in §3–§4 check out. Remaining ch17 issues are consistency ones (γ definition vs ch10/12, the SAT row vs the step-2 rule, a stale code excerpt, one citation constant).
- The other 22 exercise answers agree with my own solutions except for the C-level items (one rounded product in ch10 P3, an omitted 1-qubit term in ch14 P1, wording elsewhere).
- Literature and PennyLane API claims in ch17 match primary sources; the Qualtran snippet could be confirmed only in part.

## Checked and found correct
**Exercises (solved independently, then compared):**
- ch10: P1 Dirichlet-kernel bound 4/π², 0.684 limit; P2 linearity/sampling; P3 IPE steps for T (0.0010₂), residual phases 1/12, 1/24, 1/48 and errors 0.067, 0.017, 0.0043, 0.0011, Hoeffding and binomial tails.
- ch11: P1 eigenstate and Σ_s|u_s⟩ = |1⟩ (N&C Exercise 5.13 exists and has no erratum); P2's N = 15 specifics (a = 14, a = 4); P3 t ≥ 2n argument and the 0 / 0.23 / 1/3 pattern.
- ch12: P1 ZZ identity and Y basis change (S† then H); P2 commutator bound; P3 N ≈ 1.69×10⁷, M ≈ 2.58×10⁴, ratio 654.
- ch13: P1 full derivation and noise ratio 1/ε; P2 2 groups (TFIM) / 3 groups (Heisenberg); P3 0.7500 / 0.8333 / 0.8750 by my own simulator; arXiv:2606.29562 (Marwaha, "The QAOA on the ring of disagrees") exists and says what the chapter says.
- ch14: P1 all numbers (and consistency with ch9's "보조 n−3, Toffoli 2n−5"); P2 distribution {001: .427, 010: .213, 011: .213, 101: .125, 100: .021} reproduced by my own simulation, correct rule gives 001 with probability 1; P3 M, 0.9245 → 0.9800 (script stdout), residual ≈ 0.017 + thermal.
- ch15: P1 conjugation rules, TXT† = (X+Y)/√2; P2 stabilizer entropy = |A| − k; P3 amplitude = Φ, E|Φ| = √(2/π)/√N (0.0249 at n = 10), √(2/π) for the forrelated pair.
- ch16: P1 α = 2, probability 1/2, phases ±π/4, ±3π/4, α/‖H‖ = √2 (numerically); P2 P(x) = i^d T_d(x) for d = 3, 4, 5 at x = 0.3, −0.7 (0.792i ✓); P3's f values, the 3× vs 16× statement and the three remedies.
- ch17: P1 (N log₂N = 4,608); P2 (N* = 10²², T* = 10¹¹ s ≈ 3,170 yr, 10⁷ s ≈ 3.8 months, 1 s / 5.4×10⁷ = 18.5 ns); P3 (b = 16, 65,536 × 603 = 3.95×10⁷, 8×, n² scaling).

**ch17 second review:**
- Re-ran a copy of `ch17_discovery.py` (8.4 s): all of §1, §5 (A)(B)(C) tables, the §4 example row, §8 and §9 numbers match stdout, including hand conversions (49 s, 8.2 min, 2.3 min, 8.9 min, 37 min; 1.2 d, 12 d, 32 yr, 2.8 h, 3.6 min; 2^S ≈ 1.5 and ≈ 14).
- Independent recomputation of (A) with a differently built Hamiltonian (`R2/ch17_indep.py`): E₀, S, χ for n = 6…18, slope 0.0966, quench S/χ at t = 0…6, rate 0.910 — identical.
- Formulas: b = ⌈log₂(πα/ε)⌉ from rounding error π/2^b × α|sin θ|; Toffoli/W = (L−1) + 2L + (m−1) and qubits n + 2m + μ + 1 recomputed for all five models; L and α for each lattice (280, 540, 2280; J₂ adds 486 → 1,026); N* = (Cτ/t_c)^{p/(p−1)}, T* = (Cτ)²/t_c at p = 2, the 10⁷ sensitivity factor; equivalence to the Babbush 2021 break-even form; χ ≥ 2^S stated as a lower bound; c/6 = 1/12 in bits per log₂n for the open-chain half cut.
- Verdict table: counts 2 / 2 / 4 and "3 at step 5, 1 at step 2" agree in §3, §5, §9; checklist figure agrees with the text; cross-references to ch15 (Tang, stoquastic, 0.64 / 0.50) and ch16 (QPE bound, cost table, input-model caveat) exist and agree.
- Literature (primary sources fetched): Chakrabarti et al. published abstract = 8k logical qubits, T-depth 54 million, "order of a second"; arXiv v2 = 7.5k, 46 million, 10 MHz. Gidney 2025 = < 1 million noisy qubits, < 1 week, 0.1% gate error; 2019 estimate 20 million / 8 hours. Babbush et al. 2021 = quadratic insufficient, quartic far more practical, PRX Quantum 2, 010103. Lee et al. 2023 = Nat. Commun. 14, 1952; evidence "has yet to be found"; overlap / initial-state argument as described. Stoquastic statements about TFIM are correct.
- PennyLane 0.45.1 docs: `PauliHamiltonian(num_qubits, pauli_terms, one_norm=None)`, `UniformStatePrep(num_states)`, `SelectPauli(pauli_ham)`, `Qubitization(prep_op, select_op)`, `QPE(base, num_estimation_wires, adj_qft_op=None)`, `TrotterPauli(pauli_ham, num_steps, order)`, `estimate(workflow, gate_set=None, …)`, `Resources.gate_counts`, default gate set includes Toffoli and T — the snippet's names and argument order are all valid. Both cited doc URLs and the demo URL resolve.
- Qualtran: `GateCounts` fields (`t`, `toffoli`, `cswap`, `and_bloq`, `clifford`, `rotation`, `measurement`), `QECGatesCost` and `get_cost_value` in `qualtran.resource_counting`, module path `bloqs/multiplexers/select_pauli_lcu`, and "unary iteration uses L−1 AND operations" confirmed.

## Could not check
- Qualtran `SelectPauliLCU` constructor argument names (`selection_bitsize`, `target_bitsize`, `select_unitaries`) and the claimed output `and_bloq ≈ L − 2 = 197`: the class page was not reachable from this session; they match my recollection of the source only.
- No SDK snippet was executed (PyPI unreachable); the PennyLane check is against documentation, not a run. Whether `Qubitization(UniformStatePrep(199), SelectPauli(H))` accepts the two operators' wire counts without an explicit `wires` argument is untested.
- N&C book text for the "문제 1 (N&C 5.2.1)" attribution of 4/π² (C2) and section titles 4.5.5 / 6.6 / 10.5.2 — from memory only.
- Page / issue details of Hoefler–Häner–Troyer (CACM 66(5), 82), Aaronson (Nat. Phys. 11, 291), Calabrese–Cardy, Fattal et al. — from memory only; Costa et al. 2022 complexity in C10 likewise.
- Quirk links were checked by algebra, not opened.
- ch14 P3's 0.9245 → 0.9800 and the C7 1-qubit-gate count rely on the chapter script's circuits; I did not rebuild the noise model.
- "작업계획서 후보 1번" / "22장 반례 박스" references in ch17 P1 point to a plan document outside the chapter directories.
