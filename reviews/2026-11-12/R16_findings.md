# R16 — verification of the post-R15 revision of chapter 21 (`ch21/ch21_hp_folding.*`)

Reviewer: R16 · 2026-11-12 (UTC) · read-only on the repository. Own code and outputs: `/tmp/claude-0/-home-claude/fd63146d-511d-5ae9-8dab-3100fb85f0db/scratchpad/r16/` (`lib.py` own tables + own QAOA state, `check.py`, `replay.py`, `replay2.py`, `author_kernel.py` = py:248–276 copied for testing, `*.out`). The chapter script was not run.
Line numbers: `md:N` = `ch21_hp_folding.md` (421 lines, mtime 01:57:53), `py:N` = `ch21_hp_folding.py` (696 lines, 01:37:02), `out:N` = `ch21_stdout.txt` (144 lines, 01:57:36).

**Counts: A 0 · B 4 · C 15.**

## Verdict

The R15 delta is applied correctly. Every changed or added number in the md equals `ch21_stdout.txt`; every number quoted from an earlier draft is attributed and matches R13 / R15 / `ch21_stdout.before_R15.txt` (one exception, C7). The new kernel is right (checked against my own per-qubit implementation at 10, 14, 18 and 22 qubits), and I reproduced the new headline numbers with my own state-vector code: all 20 N = 12 SPSA seeds at p = 2 and the 5 seeds at p = 3 (p = 2: 7.24× (4.07–20.09), 169,302 (62,647–299,062), first five 13.75× / 90,328, three seeds shorter than that; p = 3: 13.29× (11.04–19.78), 93,403 (63,591–111,906) — identical to the stdout). The registered verdict (방식 1: 불통과) is stated once in §1, once in §5 (6), and nothing in §1, 4, 5, 7, 8, 9 contradicts it; 25 of 25 seeds are longer than SA, the best one 4.0×.

Nothing found is wrong enough to be an A. Four things should still be fixed: a comparison in the new p = 3 paragraph that leans on the subsample the same paragraph has just called unrepresentative (B1); the "거의 다 내려왔다" sentence omits that two of the three restarts were cut at the cap, which the script prints (B2); a failed regression run now leaves a complete-looking stdout (B3); the two lists of "what changed after the reviews" disagree (B4).

---

## B — should fix

### B1. md:252 compares the p = 3 run with "p = 2의 처음 다섯 시드(90,328)" — the subsample the paragraph has just shown to be favourable — and concludes "차이를 말할 수 없다"

- Location: md:252 ("p = 2의 처음 다섯 시드(90,328)와 비슷하고 역시 모두 담금질보다 길다. 시드 다섯 개로는 p = 2와 p = 3의 차이를 말할 수 없다."); out:85–87.
- Evidence: the fair reference for five fresh p = 3 seeds is the 20-seed p = 2 sample. From the printed numbers alone: the p = 3 maximum (111,906) is below the p = 2 median (169,302), i.e. 5 of 5 below the median (chance 1/32 if the two had the same distribution). My replay (`replay.out`, `replay2.out`): only 6 of the 20 p = 2 seeds are shorter than the p = 3 maximum and 4 shorter than the p = 3 median, and three of those six are among the first five (seeds 1, 2, 4); the median of p = 2 seeds 5–19 is 175,621. Mann–Whitney, two-sided: p = 0.024 for the five p = 3 seeds against all twenty p = 2 seeds, 0.008 against seeds 5–19. So the chapter's own data do distinguish p = 3 from p = 2 at about the 2 % level; they "look alike" only against the lucky five. R13's three p = 3 seeds were 133k, 70k, 157k (R13:125) against 125k–275k for its four p = 2 seeds — same direction.
- Why it matters: the registered hypothesis is "QAOA(p ≤ 3)". The better member of the registered family is probably p = 3 (6.0× SA, best seed 4.1×), not the p = 2 headline (10.9×). The verdict is the same either way; the sentence just should not flatten the difference by picking the lucky five.
- Proposed replacement for the two sentences:
  > 다섯 시드 모두 담금질보다 길다(가장 좋은 시드 63,591도 담금질의 4.1배다). 다섯 시드의 최대(111,906)가 p = 2 스무 시드의 중앙값(169,302)보다 짧으니 p = 3 쪽이 나아 보인다. 다만 시드가 다섯 개뿐이라 얼마나 나은지는 말할 수 없다.
  (all four numbers are in out:85, out:87, out:104.)

### B2. md:221 "거의 다 내려왔다" leaves out the flag the script prints: two of the three second searches stopped at the evaluation cap

- Location: md:221; out:61 ("p = 1: 3.3e-07 (횟수 상한에서 끊김), p = 2: 0.0e+00 (횟수 상한에서 끊김), p = 3: 8.3e-08 (수렴 기준 충족)"); py:305–309.
- `nm2` is correct: gain = f(end of first half) − f(end of second half), and the flag is `success` of the second Nelder–Mead only (the evaluation counts 62 / 162 / 362 at out:56–58 confirm 30+30, 80+80, 200+160). But a gain of exactly 0 from a restart that was cut after 80 evaluations (p = 2, four parameters) is weaker evidence than a converged restart: the fresh 5 % simplex had not yet shrunk to the scale of the remaining error. The values themselves are fine — R13's independent searches gave ⟨C⟩ = 2.6604 / 1.1746 / 0.3863 (R13:244) and 35.5× / 33,973 at p = 2 (R13:48; chapter 35.53× / 33,978).
- Proposed replacement for "둘째 탐색에서 ⟨C⟩가 더 내려간 폭은 … 전역 최솟값의 증명은 아니다.":
  > 둘째 탐색에서 ⟨C⟩가 더 내려간 폭은 p = 1에서 3.3×10⁻⁷, p = 2에서 0, p = 3에서 8.3×10⁻⁸이다. 수렴 기준까지 간 것은 p = 3뿐이고, p = 1·2의 둘째 탐색은 횟수 상한에서 끊겼다. 세 값은 검토 R13이 별도 코드로 찾은 값(2.6604, 1.1746, 0.3863)과 넷째 자리까지 같다. 거의 다 내려왔다는 뜻이고 전역 최솟값의 증명은 아니다.

### B3. A failed regression no longer shows in `ch21_stdout.txt`

- Location: py:692–696. Before R15 the asserts ran before the two closing prints, so a complete stdout implied a passed regression. Now `chk`/`near` only collect, the closing lines are printed, and `assert not _bad` comes last.
- Evidence: `scratchpad/c21/full15/ch21/out.txt` is the author's first full run of this revision with the old pins: `rc=1`, yet 144 complete lines ending in "(실행 시간: 총 1228 s)", numerically identical to the final stdout. Only `rc.txt` / stderr tell the two apart. `assert not _bad` itself is reached in every run (also under `CH21_QUICK`, where only the `if not QUICK:` block at py:664 is skipped).
- Fix (needs one re-run because stdout changes): before py:696 add
  `print(f"    회귀 검사: {'모두 통과' if not _bad else '불일치 %d건' % len(_bad)}" + (" (CH21_QUICK: 고정값 검사는 건너뜀)" if QUICK else ""))`.

### B4. The two lists of "what changed after the reviews" disagree, and md:277 points to §4 for something §4 does not say

- md:129: "검토 뒤에 바뀐 것은 넷이다" = 출발점의 부호 / 시드 수와 통계량 / 담금질의 실행 수와 시드 / p = 3 훈련. md:277: "출발점의 부호, 시드 수, 담금질의 실행 수, 탐색 예산은 검토 뒤에 바뀌었다(4절)" — 탐색 예산 is not in §4 (it is in md:221), and the p = 3 run is missing here.
- Proposed replacement for that sentence in md:277:
  > 출발점의 부호, 시드 수와 적는 통계량, 담금질의 실행 수와 시드는 검토 뒤에 바뀌었고 잔기 12개의 p = 3 훈련은 검토 뒤에 더했다(4절). 잔기 12개의 탐색 예산도 검토 R15 뒤에 늘렸다(5절 (5)).

---

## C — minor / wording

- **C1.** md:117 "(두 개씩 곱하는 계산과 확률이 같은지 스크립트가 확인한다. 차 3.3×10⁻¹⁹)": the assert (py:319–320) is at N = 8, p = 3 only. Write "스크립트가 잔기 8개, p = 3에서 확인한다". Optional: add the same assert at N = 10 (0.1 s). I checked 10 / 14 / 18 / 22 qubits (below).
- **C2.** md:250 "(중앙값 12번)", out:83–84: the medians over 20 seeds are 11.5 (N = 8) and 12.5 (N = 10) in my replay; `{t_['seen']:.0f}` (py:456) rounds both halves to even and prints 12 twice. Print with `.1f`, then "중앙값 11.5번·12.5번".
- **C3.** md:225 "⟨C⟩로 고른 각도보다 γ가 훨씬 크다": by layer the factor is 2.2 / 1.8 / 1.2 (N = 8), 2.2 / 2.1 / 1.7 (N = 10), 2.9 / 2.8 / 1.5 (N = 12). Safer without a number: "⟨C⟩로 고른 각도보다 세 층 모두 γ가 크고 β의 크기도 크다."
- **C4.** md:30 heading "좋은 각도의 절반에도 못 미친다" holds for the medians (20 % / 43 % / 38 %, out:88); the best N = 12 seed reaches 20.1 / 35.5 = 57 %. → "**2,560샷의 훈련은 중앙값으로 좋은 각도의 절반에도 못 미친다.**"
- **C5.** md:273 "⟨C⟩를 줄이는 샷 훈련은 p ≤ 3에서 ⟨C⟩의 최적점에 닿아도 92.6배에 못 미친다": 92.6× is the threshold for 2,560 training shots (with ≤ 279 shots 78.8× suffices, md:383), and 78.8× is a found optimum. → "2,560샷을 쓰는 샷 훈련은 p ≤ 3에서 (가)가 찾은 ⟨C⟩의 최적점에 닿아도 92.6배에 못 미친다."
- **C6.** md:383 "0에서 천 개 남짓 사이": 16,495 − 15,321 = 1,174 (and 14,795 − 15,321 = −526) are not printed. Have out:136 print the two differences, or write "0에서 (16,495 − 15,321)개 사이".
- **C7.** md:349 quotes the first draft as "…담금질이 1.5배 짧다" and adds "첫 판의 숫자는 검토 R13의 기록에 남아 있다". R13 has no "1.5배"; it records "24,724 대 16,045" (R13:18, R13:73). → "훈련 비용을 0으로 쳐도 잔기 12개에서는 담금질이 짧다(24,724 대 16,045)".
- **C8.** md:252 "시드 10개도 4.4\~10.8배였다": the unit is p\*/random (R13:124: 6.7, 4.4, 9.9, 8.9; R15 B1: 10.81, 5.21, 5.13, 9.48, 7.41, 8.44). → "시드 10개도 p\*가 무작위의 4.4\~10.8배였다".
- **C9.** "⟨C⟩를 최소화한" / "⟨C⟩-최적" is still used for the N = 10 / 12 local-search values (md:16, md:32, md:207, md:254, md:270, md:359, md:398) although md:349 and CLAUDE.md:74 say such values are "찾은 값". One sentence in the §4 table would cover it: "이 장에서 '⟨C⟩-최적 각도'는 (가)에서 찾은 각도를 뜻한다."
- **C10.** py:248–253: the cache key `(id(T["C"]), T["sig"])` is correct in today's flow (see below) but an `id` can be recycled once the old array is freed, and in-place edits are invisible. Hold the array instead: `if T.get("_phC") is not T["C"] or T.get("_phs") != T["sig"]: …; T["_phC"] = T["C"]; T["_phs"] = T["sig"]`.
- **C11.** py:405: `N_SHORT` scans nine `MOVED` entries (the p = 1 transfer `("C", 10, 1)` set at py:344 is in the dict but not in the table) while out:80 says "여덟 가지". Harmless (its TTS is far above SA); filter on `k[2] >= 2`.
- **C12.** py:318: the comment "출발점이 '비용을 줄이는 쪽'인지 확인한다" now sits above the mixer check (py:319–321), not above py:322–323. py:309: the `else` branch of `nm2` cannot be reached (Nelder–Mead never returns worse than its start).
- **C13.** Regression block: nothing was removed and every tolerance was tightened or kept (see below). Not pinned among the new md numbers: `TRAIN[12]` ratio range 4.07–20.09 (md:30, md:248, md:273), `TRAIN3` ratio range 11.04–19.78, the count "3개" (md:252), the cum column (167,930), 20 % / 43 % / 38 %, the start ratios 3.42 / 4.61, the three `last_gain` values.
- **C14.** Plan md:115 (table of contents) still reads "HP 격자 모델을 QUBO로 바꿔 QAOA로 풀고"; the chapter states it does not build a QUBO (md:193). Older than this delta.
- **C15.** Figure panel (6): x ticks at 2.5, 7.5, … for an integer degree; `ax.set_xticks(range(2, 21, 2))`.

---

## R15 item → status

| R15 | Status | Where / remark |
| --- | --- | --- |
| A1 "검토 전에 정해 둔 방식" | resolved | md:15, md:127, md:129, md:268, md:277; out:103, out:114. Residual: B4 |
| A2 N = 12, p = 3 not converged | resolved | py:327 (400, two halves), out:58 (0.3863 / 0.4975 / 3.005e-04 / 78.78 / 15,321 = R13/R15 value), md:219–223, md:381–383; guard py:670. Residual: B2 |
| B1 five favourable seeds | resolved | py:448, out:85–86, md:15, md:248, md:252, md:351, md:401; reproduced with own code |
| B2 tuned-SA clause in §1 | resolved | md:32 (1.55배, out:107) |
| B3 verdict sentence in §1 | resolved | md:19 |
| B4 "거의 찾지 못한다" | resolved | md:30 (nit C4) |
| B5 p = 3 shot training | resolved | py:459, out:87, out:104, md:252, md:269, md:273. New wording issue B1; scope nit C5 |
| B6 SA schedule caveat | resolved | md:203, out:38 |
| B7 p\*-optimal angles | resolved | py:388, out:68, md:225; the printed N = 10 angles give 202.5× / 5,960 at N = 12 in my code (3-decimal angles; script 202.4× / 5,963) |
| C1 "G ≈ 100\~170" | resolved | md:139, out:113 (184배; ch14:198 is in §7 and says "CZ 수의 한계 ≈ 100 → 0.3 오류 예산") |
| C2 227.95 / 227.9 | resolved | md:379, out:131–133 |
| C3 "걸음" | resolved | py:203–204, py:616, md:83, md:291; figure label checked |
| C4 27·52·82 | resolved | md:221 |
| C5 "세 방식 모두" | resolved | md:273 |
| C6 readings stated as facts | resolved | md:254, md:30 |
| C7 unprinted numbers | mostly | printed now: 97.13 % (out:8), 2.87 (out:80), 1e-24 (out:113), 284,491. Left: "1.5배" (C7), new "천 개 남짓" (C6) |
| C8 rounding note | resolved | md:383, out:133 |
| C9 quotation marks | resolved | md:19 |
| C10 best-of-eight | resolved | md:240, out:80, chk py:672 (nit C11) |
| C11 sum 284,491 | resolved | md:24, out:80 (18.2배, 1/4.2) |
| C12 T ≈ 15,846 | resolved | md:139 |
| C13 "항 하나에 하나꼴" | resolved | md:131 (claim withdrawn) |
| C14 "이 서열에서는" | resolved | md:81, md:240, md:399 |
| C15 regression gaps | mostly | py:680–691 pin the eight transfer rows, ratios, ten sensitivity rows, warm start, p_L; py:1, md:421. Left: C13, B3 |

---

## Script review (task 3)

- **Grouped mixer** (py:270–274): `np.matmul(R, psi.reshape(2^q, 2^k, −1))` applies r^{⊗k} to the middle axis; groups are 4,4,4,2 (n = 14), 4,4,4,4,2 (18), 4×5,2 (22), 4,4,2 (10); an odd remainder would also work. The mixer is a symmetric tensor power, so it carries no bit order; the phase uses the table index unchanged, so big-endian is as before. Author's fast kernel vs my per-qubit implementation (own tables): max probability difference 3.5×10⁻¹⁸ (n = 10), 1.5×10⁻¹⁸ (14), 6.5×10⁻¹⁹ (18), 4.9×10⁻¹⁹ and 3.2×10⁻²¹ (22, p = 3 and p = 2). The `ref=True` path agrees equally.
- **Phase table**: `exp(−iγ·u)[inv]` with `u = unique(C)/σ` is the element-wise expression. No stale table in the λ sweep: each new `C` is allocated while the previous one is still referenced by `T10["C"]` (distinct ids; I got 5 distinct ids) and σ differs for every λ; after the sweep `keep[0]` is a copy with a new id. Replaying the sweep's assignment sequence with the copied functions: fast vs ref ≤ 1.2×10⁻¹⁹ for all five λ, and 1.6×10⁻¹⁹ vs my code after the restore. `TAB[6]` is a fresh dict (py:560) with `N_EVAL[6] = 0`; its result is also checked against the gate list (3.4×10⁻¹⁷).
- **N_EVAL**: `ref=True` does not count; py:319 resets, py:358 and the figure restore. N = 12: 3 (start) + 62 + 162 + 362 (가) + 381 (나) + 5 (other transfers) + 20×41 + 5×41 = 2,000 = out:143. N = 8 / 10 counts (11,452 / 6,252) are unchanged from before R15, i.e. the new kernel did not move any Nelder–Mead path there.
- **spsa_many**: seed lists unchanged (`[21, N, p, m, k, sd]`), so seeds 0–4 at N = 12 are the pre-R15 streams; the p = 3 run uses its own streams. Replay below.
- **New prints**: each recomputed from printed values — 10.9 (4.0–19.2), 7.1, 6.0 (4.1–7.2), 1.73, 0.60, 1.55, −279 / 850 / 0.33, 18 and 18.2, 1/4.2, 2.87, 0.0315, 20 % / 43 % / 38 %, 77.4 / 92.6, 184, e^(−55.16) = 1.1×10⁻²⁴, s.e. range 14,795–16,495 (p_L = 339/4,000). Labels say what is computed (nits C2, C11).
- **Regression block**: all 21 former asserts survive as `chk`/`near`; tolerances 1 % → 0.2–0.3 %, 0.5 % → 0.2 %, 0.2 % → 0.1 %; new checks for p_L, eight transfers, sensitivity rows, warm start, first-five median, p = 3 run, "(가) ≤ 0.3865", "only two transfers beat SA", "|15,321 − 15,601| < 850". Loosest pins: `atol=0.06` on 19.3 / 16.2 / 12.9 / 19.6 and `atol=2e-3` on the N = 8 means (the ratios pin those to 0.2 %). `near` has no default tolerance in use. See B3, C13.

## Independent recomputations (own code)

| Quantity | md / stdout | R16 | Match |
| --- | --- | --- | --- |
| Tables N = 6, 8, 10, 12: valid, E\_min, optima, mean, σ, invalid share | 284 / −2 / 8; 2,172 / 24; 16,268 / 376; 120,292 / −5 / 16; 3.562 / 3.073, 5.252 / 4.220, 6.348 / 4.990; 86.74 / 93.79 / 97.13 % | same | ✓ |
| Start `ramp(2)` N = 8 / 10 / 12; `ramp(3)` N = 12 | 2.008, 2.62×; 3.164, 2.58×; 3.927, 3.42×; 4.61× | 2.0076, 2.62; 3.1641, 2.58; 3.9266, 3.42; 4.61 | ✓ |
| Printed ⟨C⟩ angles, p = 3 (3 decimals) at N = 8 / 10 / 12 | 0.1031 / 14.48; 0.3578 / 20.53; 0.3863 / 0.4975 / 78.78 / 15,321 | 0.1031 / 14.44; 0.3578 / 20.54; 0.3863 / 0.4971 / 78.63 / 15,351 (fast and per-qubit mixer agree) | ✓ to angle rounding |
| Printed p\* angles, p = 3 at N = 8 / 10 / 12 | 22.9 / 135; 29.9 / 105; 262.8 / 4,592 | 22.86 / 135; 29.91 / 105; 262.76 / 4,592 | ✓ |
| Transfers to N = 12 (printed angles): C10, P10, C8, P8 | 70.5 / 17,114; 202.4 / 5,963; 62.1 / 19,453; 161.8 / 7,458 | 70.59 / 17,100; 202.48 / 5,960; 61.89 / 19,503; 162.07 / 7,446 | ✓ to angle rounding |
| Case N = 6 at γ = 0.8199, β = −0.2145 | 0.2904 / 0.0313 | 0.2904 / 0.03129 | ✓ |
| SPSA N = 8, 20 seeds | 4.21× (2.89–7.48); 3,305 (2,978–3,644) | 4.2070 (2.89–7.48); 3,305 (2,978–3,644) | ✓ |
| SPSA N = 10, 20 seeds | 5.10× (3.85–7.75); 3,187 (2,972–3,392) | 5.1025 (3.85–7.75); 3,187 (2,972–3,392) | ✓ |
| SPSA N = 12, p = 3, 5 seeds | 13.29× (11.04–19.78); 93,403 (63,591–111,906); seen 0 | 13.2887 (11.04–19.78); 93,403 (63,591–111,906); 0 | ✓ |
| SPSA N = 12, p = 2, 20 seeds | 7.24× (4.07–20.09); 169,302 (62,647–299,062); seen 0; first five 13.75× / 90,328; 3 shorter | 7.24 (4.07–20.09); 169,302 (62,647–299,062); 0; 13.75 / 90,328; 3 (seeds 1, 2, 15) | ✓ |

Per-seed TTS₉₉ (N = 12, p = 2, seeds 0–19): 152,825 · 84,121 · 79,854 · 216,686 · 90,328 · 202,621 · 247,939 · 156,615 · 129,027 · 105,112 · 168,232 · 170,373 · 212,570 · 193,554 · 97,555 · 62,647 · 175,621 · 225,577 · 262,078 · 299,062 (seeds 0–4 equal R15's replay of the pre-R15 run: same streams). p = 3 (seeds 0–4): 111,906 · 73,795 · 63,591 · 99,308 · 93,403.

## Verified OK

- **md delta vs stdout.** All 38 changed hunks read sentence by sentence; every number in md:15, 19, 24, 30, 32, 35, 117, 129, 139, 203, 217–225, 240, 248–254, 260–277, 351, 357, 361, 379–383, 399, 401, 421 equals the stdout (rounding included). Logic checked and sound: "스무 시드 모두 담금질보다 길다" (min 62,647 vs SA upper s.e. bound 16,495); "284,491 … 방식 1의 중앙값보다도 길다" (> 169,302; it is below the maximum 299,062, and the sentence says median); "좋은 각도의 절반에도 못 미친다" for the medians; "90,328보다 짧은 것은 셋뿐"; 0.33 / 850 / 279; "1.55배", "1.73배", "184배", "10⁻²⁴"; "근소하게 길다 → 짧다" (15,734 > 15,601 > 15,321).
- **Earlier-draft numbers.** 196,189, "12배", 16,045, 2,000 runs: R13:18, R13:128, R13:144, R15 A1. +0.3873 / 76.7배 / 15,734, 90,328 / 13.75배, "5.8배": `ch21_stdout.before_R15.txt`. 4.4–10.8배: R13:124 + R15 B1. 262.8배 / 4,592: R13:81. All attributed in the text ("첫 판 … R13의 기록", "이전 판", "둘째 판", "검토 R15").
- **Provenance of the stdout.** `scratchpad/c21/full15/run2.rc` = `rc=0`, `run2.err` empty, `run2.out` byte-identical to `ch21_stdout.txt`; the earlier run with old pins gave the same numbers (deterministic). mtimes: py < stdout / png / json < md < html. html contains the new text and none of eleven removed phrases. No stale chapter-21 number elsewhere in the repo (`grep` outside `ch21/` and `reviews/`).
- **Plan / CLAUDE.md / delivery note.** Plan md:235, md:253, md:293, md:315, CLAUDE.md:74 and `전송 목록.md` agree with the stdout (169,302; 62,647–299,062; 10.9배; 93,403; 15,601; 17,114; 5,963; 278,528; 76.7 → 78.8; 90,328; "약 20분"; R15 counts A 2 · B 7 · C 15; 5 decisions, 4 open questions; 14 new + 17 changed files).
- **Format.** Ten sections; no bare `~ _ * [`; no bare `|` outside tables; every table preceded by a blank line with constant column count; "실측" once (md:400, the ch14 rule); status lines md:111, md:144, md:310 follow CLAUDE.md; no `\uXXXX`.
- **Figure.** Regenerated with the run; panel (3) label "스텝 = 평가"; panel (4) shows 169,302 with the 62,647–299,062 bar and the two transfer markers.

## Could not verify

- The Claude Docs copies (chapter 21, plan) — no network; the md is the source.
- Execution of the two SDK snippets (no SDK).
- The first draft itself (numbers taken from R13) and R15's six independent seeds (taken from R15).
- Whether p, shots, iterations, gains and the SA schedule were fixed before any result was seen (the text now says this cannot be checked).
- Global optimality of any N ≥ 10 value (local searches; they agree with R13's independent searches to four digits).
