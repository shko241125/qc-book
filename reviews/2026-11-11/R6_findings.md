# R6 findings — ch18 "탐색 도구와 반례", second review (after the R5 edits)

Files reviewed (state of 2026-11-10 12:26): `/home/claude/qc-book/ch18/ch18_search_tools.md` (460 lines),
`ch18_search_tools.py` (785 lines), `ch18_stdout.txt`, `ch18_results.png`; the first review
`/home/claude/qc-book/reviews/2026-11-10/R5_findings.md`; the edit log
`/tmp/claude-0/-home-claude/fd63146d-511d-5ae9-8dab-3100fb85f0db/scratchpad/c18/r5_log.json` (54 edits);
cross-chapter lines `ch17/ch17_discovery.md:11, 52, 150, 192, 239` and `ch15/ch15_sources.md:43, 58, 60`.

Scratch code and outputs: `/tmp/claude-0/-home-claude/fd63146d-511d-5ae9-8dab-3100fb85f0db/scratchpad/mid/R6/`
(`c1_dla.py`, `c2_f4.py`, `c2b.py`, `c3_gram.py`, `c4_ccz.py`, `c5_f1.py`, `c6_hva.py`, `c7_cross.py`, `c8_bookf4.py`,
`c9_tail.py`, `c10_f2rep.py`, `run/rerun_stdout.txt`).

Re-run: a scratch copy of the script was run once (252 s on a loaded 2-core container; the book's run says 218 s).
Its stdout is identical to `ch18_stdout.txt` except the timing fields; the png is pixel-identical. Nothing under
`/home/claude/qc-book/` was modified.

Counts: **A 1 · B 7 · C 15**.

---

## Part 1 — status of every R5 finding

| R5 id | status | note |
| --- | --- | --- |
| A1 stable-rank table | addressed; author deviated (measured columns r = 1,000/4,000/16,000 instead of R5's upper-bound column). The deviation is better than R5's proposal. One cell is too optimistic → **B2** below | script header fixed (py:456, 471) |
| A2 exercise 3 | addressed correctly | 289,354 = 15 % ✓; "0.1초 남짓" ✓ (0.14 s) |
| A3 대칭 장벽 | addressed correctly in §3(나), trap table, §5 text, 문제 2 title, §9 | 3′ row wording still loose → C12 |
| A4 2×10⁻⁸ | addressed correctly | |
| B1 HVA baseline | **half addressed**: text now says the numbers are optimiser-limited; table column, figure (3) and script unchanged; the supporting numbers are attributed to an "독립 검토" that is not in the book → **B3, B7** | thresholds p = 4, 9 re-confirmed by me |
| B2 CCZ classes | addressed, but the new sentence about the 12-start re-run is wrong → **A1** | |
| B3 scope of the boundary table | addressed correctly ("0 이하" instead of R5's "음수" is the right deviation: the off-diagonals are −0.9045 or 0) | term "위상 문제" → C1 |
| B4 \|0…0⟩ vs \|+…+⟩ | addressed correctly (new paragraph, three Wick rows, ch17:150 changed) | the script function behind the new rows is fragile → **B1**; "다룬다" → **B6** |
| B5 central message | addressed correctly (§1, §7) | §9 last note not aligned → **B5** |
| B6 추천 질의 | addressed correctly | |
| B7 조건수 | addressed correctly (4.06×10⁵ printed by F2) | §2 analogy not aligned → C3 |
| B8 Cerezo caveat | addressed correctly (md:141, md:401) | |
| B9 커널 집중, track steps | addressed (definition + citation added; ch17:239 now lists four steps) | definition unconditional, wrong row → **B4** |
| C1, C2, C3, C5, C6, C7, C8, C12, C13, C14, C16, C17 | addressed correctly | |
| C4 references | (a) Valiant added but without "등급" → C7; (b), (c), (d) done | |
| C9 glosses | done for the nine listed terms; new unglossed terms → C8 | |
| C10 figure | panel (3) fixed (HVA clipped to the floor, axis label); panel (6) label still lies across the curve → C9 | |
| C11 unprinted numbers | done; code-cell comment still says "약 1초" → C6 | |
| C15 Majorana degree of Z_j | done in §3; §9 still says "4차 이상" → C2 | |

Left open from R5: B1 (table/figure), C10 (panel 6), and the KLM-exercise deviation noted under C14 (not touched; same as other chapters).

---

## A — factually wrong

### [A1] ch18 §5(A) — md:250 (new sentence from the R5-B2 edit); py:653–654
Current text:
> `시작점을 12개로 늘린 재실행(스크립트 F1)에서는 같은 4묶음이었지만, 독립 검토는 다른 난수로 두 묶음을 더 찾았다.`

What is wrong. F1 prints only the count ("70묶음 중 4묶음이 닿음 (시작점 4개일 때 4묶음)"). The four classes are **not** the same four.
I replicated F1's first call exactly (`c5_f1.py`: the book's own `Template`/`instantiate`/`canon`, `default_rng(1818)`, which F1 is the first
to use; result 70 classes, 4 solved, as in stdout) and listed them:

| run | solved classes |
| --- | --- |
| A3, 4 starts (stdout lines 29–32) | 01 01 02 12 02 12 · 01 02 01 12 02 12 · **01 02 01 12 12 02** · 01 02 12 02 12 01 |
| F1, 12 starts (seed 1818) | 01 01 02 12 02 12 · 01 02 01 12 02 12 · **01 02 12 01 12 02** · 01 02 12 02 12 01 |

The 12-start run found a class the 4-start run missed (01 02 12 01 12 02 — one of the "two more" the sentence credits to the independent
review) and lost a class the 4-start run had found (01 02 01 12 12 02; it ended at 0.1464). The book's own two runs already give five classes.
With my own template builder (`c4_ccz.py`, analytic gradient checked to 5×10⁻¹¹, 40 single starts per class) all six classes reach CCZ;
per-start success 0.72, 0.68, 0.12 (01 02 01 12 12 02), 0.70, 0.30 (01 02 12 01 12 02), 0.12 (01 02 12 12 01 02). So "시작점 하나가 성공할
확률은 높지 않아서" is true for three classes (12–30 %) and not for the other three (about 70 %).

Proposed replacement for the first five sentences of md:250 (up to "…두 묶음을 더 찾았다."):
```
주의: "닿지 못했다"는 수치 탐색의 결과이지 증명이 아니다. 닿을 수 있는 구조에서도 시작점 하나가 성공하는 비율은 구조에 따라 크게 달라서(이 장의 틀에서 10\~70%), 시작점 4개로는 닿는 구조를 놓칠 수 있다. 그래서 "4묶음"이라는 개수는 이 실행의 결과일 뿐이다. 시작점을 12개로 늘린 재실행(스크립트 F1)에서도 닿은 묶음은 4개였지만 같은 넷이 아니다. 처음 실행이 놓친 01·02·12·01·12·02를 찾았고, 처음 실행이 찾은 01·02·01·12·12·02는 놓쳤다. 두 실행을 합치면 닿는 묶음은 적어도 5개다.
```
Script: make F1 print the classes, e.g. after py:654
`for s in solved_classes: print("      닿은 구조:", " ".join(f"{a}{b}" for a, b in s))` (return the solved list from `count_solved`).
If the author wants the sixth class and the 10–70 % figure in the book, add to F1 a loop over the six named structures with 40 starts each
(about 20 s in total; my timing 3 s per class) and quote its output instead of "독립 검토" (see B3).

---

## B — misleading, incomplete, fragile, or contradicting another section/chapter

### [B1] script `gaussian_entropy` mishandles zero modes; the new GHZ rows are right only by rounding — py:527–531; md:341–342, 344; stdout F4
Current code: `nu = np.clip(np.linalg.eigvalsh(1j * G[:2 * l, :2 * l]), -1, 1); nu = nu[nu > 0]`

What is wrong. Eigenvalues of iΓ come in ±ν pairs and each pair is one mode. If a mode has ν = 0 exactly, the pair is two numbers at
rounding level, and `nu > 0` counts 0, 1 or 2 of them depending on their signs; only 1 is right (H₂(1/2) = 1 bit). The GHZ₋ quench has such
an exact zero mode in the half block at every time I looked at (t = 0, 1, 2, 3).
- Book's own code path (`c8_bookf4.py`): for GHZ₋ at t = 2 the two middle eigenvalues are −2.0×10⁻¹⁵ and +1.7×10⁻¹⁵ → counted once → 2.4979
  (the printed, correct value). The same function applied to GHZ± at t = 0 (eigenvalues 0.0, 0.0) returns 0.0000; the true entropy is 1.
- With Γ built by a different but equivalent route (`c2b.py`: my bit-operation Majoranas and a Taylor propagator) the same formula gives
  3.4979 for GHZ₋ at t = 2 (both rounding-level eigenvalues positive), and 2.6013 instead of 1.6013 at t = 1.
- Counting the upper half of the sorted spectrum gives S_가우스 = S_실제 for both sectors at all four times (|difference| below 10⁻¹²).

So the table rows `GHZ₊ 2.4977 | 2.4977 | 0` and `GHZ₋ 2.4979 | 2.4979 | 0` are correct, and the claim "|GHZ±⟩ 각각은 가우스 상태" stands
(it also follows analytically: GHZ± are stabilised by the n commuting bilinears i c₂ⱼ₊₁c₂ⱼ₊₂ and i c₀c₂ₙ₋₁). But on another BLAS the script
can print 3.4979 or 1.4979 for GHZ₋ and appear to refute the text. The n = 1000 run and D2/D4 are not affected (no zero mode there).

Proposed code (same output for every existing line of stdout):
```python
def gaussian_entropy(G, l):
    """공분산 Γ_ab = (i/2)⟨[c_a, c_b]⟩ 의 앞 2l×2l 블록에서 블록 엔트로피(비트). ±ν 쌍마다 모드 하나: 위쪽 l개만 센다(ν = 0 인 모드도 한 번)."""
    nu = np.clip(np.sort(np.linalg.eigvalsh(1j * G[:2 * l, :2 * l]))[l:], 0, 1)
    q = np.clip((1 + nu) / 2, 1e-15, 1 - 1e-15)
    return float(-(q * np.log2(q) + (1 - q) * np.log2(1 - q)).sum())
```
(checked in `c11_fix.py`: old and new agree on the |+…+⟩ row, the |0…0⟩ row and the n = 1000 values 2.001463 and 36.896295), and let F4
print the difference (py:713), so that the "0" in the table is a printed number: on the book's code path the fixed function gives
5.6×10⁻¹³ (GHZ₊) and 6.4×10⁻¹³ (GHZ₋). Table cells md:341–342: `6×10⁻¹³` instead of `0`.

### [B2] ch18 §5(C) stable-rank table, row 1 verdict — md:302 (and md:17 is fine)
Current text: `| 계수 5(추천 행렬형) | 3.0 | 0.136 | 0.049 | 0.031 | 0.016 | 1,000행 안에서 5% — 표본이 이득 |`

What is wrong.
1. 0.049 is a mean of three samples (py:669). My measurement on the same matrix type (`c3_gram.py`, own generator and seed): r = 1,000 →
   0.058 ± 0.015 (20 runs) and 0.055 ± 0.011 (20 more); r = 1,500 → 0.046; r = 2,000 → 0.041. R5 measured 0.059. The standard error of a
   3-run mean is about 0.008, so "5 % within 1,000 rows" is a favourable draw; 5 % needs roughly 1,200–1,500 rows — more than the matrix has.
   This also contradicts 문제 3 of the same chapter, which derives r ≈ 1,500 for 5 % at stable rank 3 (md:427).
2. The row is judged at 5 %, every other row at 10 %. At the common 10 % criterion the verdict is solid: 0.107 at r = 300, 0.089 at r = 400
   (mine), 0.136 → 0.049 between r = 200 and 1,000 in the book's own table.

Proposed replacement of the last cell (uses only numbers already in the table):
```
| 계수 5(추천 행렬형) | 3.0 | 0.136 | 0.049 | 0.031 | 0.016 | 10%에 200\~1,000행 — 행렬보다 적다, 표본이 이득 |
```
Recommended script change: 10 repetitions instead of 3 in F2 (py:669; about 20 s more). How noisy three are is visible in F2 itself: its r = 200
column (0.129, 0.894, 4.57, 8.13) disagrees with the C3 column printed in the book (0.136, 0.966, 4.41, 5.84); the tridiagonal 8.13 is the mean of
6.83, 8.85, 8.70 (`c10_f2rep.py`, exact replication of the rngF stream), three draws with a row sampled 4, 4 and 5 times, while the typical value is
5.4 ± 0.8 (300 runs, `c9_tail.py`). The other 14 cells of the book's table agree with my means within noise (see "verified").

### [B3] "독립 검토" is cited as evidence twice; the book and its script cannot back it — md:250, md:270, md:220
Current text:
> (md:250) `독립 검토는 다른 난수로 두 묶음을 더 찾았다.`
> (md:270) `독립 검토에서 시작점을 수백 개로 늘리자 n = 6의 p = 6\~8 값이 한 자릿수 이상 더 내려갔다.`
> (md:220) `모든 수치는 ch18_search_tools.py의 출력이다(seed = 18).`

A reader has no access to the review, and md:220 says every number comes from the script. Both claims are true as far as I can check — six
classes (A1 above); HVA n = 6 with my own implementation (`c6_hva.py`, adjoint gradient checked to 2×10⁻⁹, BFGS, random starts in [−π, π]):
p = 6: 1.34×10⁻⁴ in 250 starts (book 2.26×10⁻³, 17× lower); p = 7: 1.55×10⁻⁴ in 250 starts (book 7.35×10⁻⁴, 4.8× lower; R5 reached
1.8×10⁻⁵ with 600); p = 8: 4.3×10⁻⁷ in 200 starts (book 5.9×10⁻⁵, 136× lower); p = 9: below 10⁻¹³ (6 of 40 starts) — but they should
come from the book's own script or be worded without the reviewer.

Proposed. Add to section F of the script: (i) the six-structure loop of A1; (ii) `hva(Hd6, 6, 6, plus6, rngF, starts=250)` and the same for
p = 8 (26 s and 45 s here), printing the best error. Then md:270:
```
시작점을 250개로 늘리면(스크립트 F5) n = 6의 p = 6과 p = 8의 최선 오차는 한 자릿수 이상 더 내려간다.
```
(with the two printed values inserted; mine were 1.3×10⁻⁴ and 4×10⁻⁷). Until the script has it, a wording that needs no source:
```
이 값들은 시작점을 늘리면 더 내려간다(스크립트에는 넣지 않은 별도 확인에서 p = 6\~8이 대략 한 자릿수 이상).
```
md:220: `모든 수치는 ch18_search_tools.py의 출력이다(seed = 18, 보강 확인 F는 seed = 1818).`

### [B4] 커널 집중: stated unconditionally, and filed under the wrong step — md:141, md:146 (new paragraph), md:139
Current text:
> (md:146) `커널 집중(kernel concentration)은 큐비트 수가 늘수록 서로 다른 두 입력의 커널 값이 한 상수로 몰려, 그 차이를 읽는 데 필요한 샷 수가 지수로 느는 현상이다(Thanasilp et al. 2024). 불모 고원의 커널 판에 해당한다.`
> (md:141, 5′ 고전 대리 계산, 도구) `위 "고전 경계 검사" 표, 커널이면 커널 집중과 RBF 비교(22장)`

What is wrong.
1. The paper (arXiv:2208.11060v2 abstract, fetched today from arxiv.org/html/2208.11060v2) says kernel values "**can** be exponentially
   concentrated" and names four sources: expressivity of the data embedding, global measurements, entanglement, noise. The book's sentence
   reads as if every quantum kernel concentrates as n grows. (The rest matches Definition 1 of the paper: variance over pairs of inputs in
   O(1/bⁿ); with polynomially many shots the estimated Gram matrix carries no information about the inputs.)
2. The paragraph itself says "불모 고원의 커널 판" and "샷 수가 지수로" — that is trainability (3′) and measurement budget (4′). The table puts it
   under 5′, whose question is "훈련을 가능하게 한 구조가 고전 계산도 가능하게 하는가". A concentrated kernel is not classically surrogated; it is
   untrainable. The term also appears in the table (md:141) before it is defined (md:146).

Proposed replacements.
md:139:
```
| 3′ 훈련 가능성 | 틀이 답에 닿는가, 기울기(커널이면 커널 값의 차이)가 살아 있는가 | 궤도·DLA 차원과 기울기 분산(6·13장), 풀의 대칭 검사(이 장), 커널이면 커널 값의 분산(22장) | 불모 고원(기울기 분산이 n에 지수로 감소), 커널 집중(아래 설명), 풀이 보존하는 대칭에서 시작 상태의 구역이 답의 구역과 다름 |
```
md:141, third cell: `위 "고전 경계 검사" 표, 커널이면 고전 커널(RBF)과의 비교(22장)`
md:146:
```
커널 집중(kernel concentration)은 조건에 따라(표현력이 지나치게 큰 특성 맵, 전역 측정, 큰 얽힘, 잡음) 큐비트 수가 늘수록 서로 다른 두 입력의 커널 값이 한 상수로 지수적으로 몰리는 현상이다. 그러면 그 차이를 읽는 데 필요한 샷 수가 지수로 는다(Thanasilp et al. 2024). 불모 고원의 커널 판이어서 3′ 단계에서 본다. 검사는 n을 키우며 입력 쌍에 대한 커널 값의 분산을 재는 것이고, 22장에서 한다.
```
ch15:43 ("커널 집중 — 18장 휴리스틱 트랙") is then delivered as written.

### [B5] ch18 §9 last note still says what the revised §7 withdrew — md:436
Current text: `TFIM의 HVA가 불모 고원 없이 훈련되는 것과 그 동역학이 행렬 하나로 풀리는 것은 같은 사실의 두 얼굴이다.`

§7 (md:401) now says "불모 고원이 없다는 것과 잘 훈련된다는 것은 다르다. 5절의 HVA는 … 무작위 시작점 대부분이 국소 최소에 빠졌다" and that ADAPT's
success is not explained by the circuit being easy. The note keeps "훈련되는 것" and gives only 인수분해 as the place where the relation breaks.

Proposed replacement of the note:
```
- **도구와 반례는 한 쌍** — 탐색이 잘 되는 이유(구조)가 고전으로 쉬운 이유일 수 있다. TFIM의 HVA에 불모 고원이 없는 것과 그 회로가 2n×2n 행렬로 시뮬레이션되는 것은 같은 사실(Majorana 2차식, 다항 크기의 리 대수)의 두 얼굴이다. 비유: 시험이 쉬워서 점수가 잘 나온 것을 실력으로 착각하지 않기. 깨지는 곳: 이 관계는 경향이지 정리가 아니다. 불모 고원이 없어도 국소 최소 때문에 훈련이 안 될 수 있고(5절의 HVA), 이 장의 ADAPT처럼 회로는 자유 페르미온이 아닌데 표적 문제가 쉬워서 성공하는 경우도 있다. 구조가 있으면서 고전에 어려운 예(인수분해)가 책 전체의 출발점이다. 더 깊이: Cerezo et al. 2023.
```

### [B6] "|0…0⟩ 시작은 가우스 상태 둘의 중첩으로 다룬다" promises more than ch18 does — ch17:150, md:391; md:322
Current text:
> (ch17:150) `(18장 5절: 같은 모형을 |+…+⟩에서 급랭해 n = 1000, t = 40까지 약 10초. |0…0⟩ 시작은 가우스 상태 둘의 중첩으로 다룬다)`
> (md:391) `(|0…0⟩ 시작은 가우스 상태 둘의 중첩으로 다룬다)`
> (md:322) `가우스 계산 두 번과 두 구역 사이의 교차항(가우스 상태 사이의 행렬 원소)으로 풀린다. 교차항 계산은 이 장에서 구현하지 않았다.`

"다룬다" reads as "ch18 does it". ch18 does not: it shows that the two parity components are Gaussian and states that the cross term was
not implemented. The n = 1000, t = 40, "약 10초" part is delivered exactly (9.7 s; |+…+⟩ start).

The claim itself is correct, and I checked it by computation (`c2_f4.py`, `c7_cross.py`, n = 12). Since c₀ = Z₀ and Z₀|GHZ₊⟩ = |GHZ₋⟩, the state
at time t is (1 + γ(t))|GHZ₊(t)⟩/√2 with γ(t) = Σ_b (e^{At})_{b0} c_b, so every expectation value is a Wick expansion in the single Gaussian
state |GHZ₊(t)⟩: ⟨c_a⟩(t) = (e^{At})_{a0} (state vector vs 24×24 matrix: 7×10⁻¹⁵; this is the book's max|⟨c_a⟩| = 0.561), and
⟨Z_j⟩(t) for j = 0, 3, 6, 11 at t = 0.5, 1, 2, 3 from a Pfaffian of 24×24 data agrees with the state vector to 10 digits
(e.g. ⟨Z₆⟩(2) = 0.0964858613 both ways). This covers the observable of the CUDA-Q snippet (`spin.z(n // 2)`), which is exactly a cross term.
What does **not** follow from the Γ formula is the entropy that ch17 tabulated: for the |0…0⟩ quench S = 1.4980, while the sectors give
2.4977 and 2.4979 (not "minus one bit" exactly).

Proposed replacements.
ch17:150 parenthesis:
```
(18장 5절: 같은 모형을 |+…+⟩에서 급랭해 n = 1000, t = 40까지 약 10초. |0…0⟩ 시작은 가우스 상태 둘의 중첩이어서 같은 방법으로 다룰 수 있다 — 18장은 근거만 보이고 구현하지는 않았다)
```
md:391 parenthesis: `(|0…0⟩ 시작은 가우스 상태 둘의 중첩이어서 관측량 기댓값은 같은 방법으로 계산할 수 있다; 5절 (D))`
md:322, optional stronger version of the three sentences from "그래도 고전으로 쉽다." to "…구현하지 않았다.":
```
그래도 고전으로 쉽다. |GHZ±⟩ 각각은 가우스 상태이고(아래 Wick 검사) H가 패리티를 보존한다. 두 구역은 Majorana 하나로 이어진다(Z₀ = c₀이고 Z₀|GHZ₊⟩ = |GHZ₋⟩). 그래서 시각 t의 상태는 (1 + γ(t))|GHZ₊(t)⟩/√2, γ(t) = Σ\_b (e^{At})\_{b0} c\_b 이고, 관측량의 기댓값은 가우스 상태 하나 위의 Wick 전개가 된다(예: ⟨Z₀⟩(t) = (e^{At})\_{00}). 이 계산은 이 장의 스크립트에 넣지 않았다. 반쪽 엔트로피는 기댓값이 아니어서 3절의 식을 그대로 쓸 수 없다.
```
(The identities in this text are verified by my scratch code only; add a three-line check to F4 before printing them.)

### [B7] HVA table and figure (3) still show the optimiser-limited values (R5-B1, left open) — md:262–268, figure panel (3)
The text now warns, but the column `최선의 2배 이내로 간 시작점 비율` is still a ratio to a "best" that is 17× (p = 6) to 136× (p = 8) above
values that more starts find (B3), and the n = 6 triangles at 12–16 parameters in panel (3) sit one to two decades too high with no note in
the figure list. Minimum fix without re-running: rename the column and add one sentence to the figure description (md:357).
```
| n | HVA 층 수 p | 파라미터 | CNOT | 최선 오차(시작점 21개 중) | 그 최선값의 2배 이내로 간 시작점 비율 |
```
```
그림 (3)에서 n = 6의 HVA 삼각형 가운데 파라미터 12\~16개의 점은 시작점 21개 중 최선이라, 시작점을 늘리면 한두 자릿수 더 내려간다.
```

---

## C — minor

- **[C1] md:122** `몬테카를로에 위상 문제가 생기므로` — "위상 문제" has no gloss and reads as "topology problem". The usual name is the dynamical sign problem. → `몬테카를로 가중치가 복소 위상이 되어 부호 문제가 다시 생기므로(동역학적 부호 문제, dynamical sign problem)`. The scope statement itself is right: the ch17 model −ΣZZ − 0.9045ΣX − 0.8090ΣZ has off-diagonal elements −0.9045 or 0, and ch17:192 rates its quench "가능". ch15:58 lists stoquastic among the boundaries without this scope; consider adding "(바닥상태·열평형)" there.
- **[C2] md:435** `세로장·다음 이웃 결합·2차원처럼 4차 이상의 항이 하나라도 있으면` — the longitudinal field is odd-degree (Z_j: 1, 3, 5, 7 for n = 4, `c1_dla.py`), as §3 now says. → `세로장(홀수 차)·다음 이웃 결합·2차원(4차 이상)처럼 2차가 아닌 항이 하나라도 있으면 고차 상관이 독립적으로 움직여 쓸 수 없다. 초기 상태가 가우스 상태가 아닐 때도 Γ 하나로는 모자란다(5절의 |0…0⟩ 급랭).`
- **[C3] md:28** `계수가 높은 희소 행렬은 표본이 구조를 못 잡아서(5절) 가속이 남는다.` — after the B7 edit §5 says the surviving case also needs a good condition number. → `계수가 높은 희소 행렬은 표본이 구조를 못 잡는다(5절). 조건수까지 좋으면 가속이 남는다.`
- **[C4] md:250** `(구조 하나하나의 바닥은 0.375, 0.2714, 0.1464로 다르다)` — this is the CNOT-5 list from a 4-start run (stdout F1). → `(CNOT 5개에서 구조 하나하나가 멈춘 값은 0.375, 0.2714, 0.1464의 세 가지다)`.
- **[C5] md:270** `HVA는 n = 6의 p = 1\~9에 시작점 21개씩` — p = 1 has no warm start (py:302): 20 starts (stdout fraction 0.65 = 13/20). → `p = 1\~9에 시작점 20\~21개씩`.
- **[C6] md:151** code comment `(약 1초)` vs md:148 `(1\~2초)` — make the comment `(1~2초)`.
- **[C7] §10** md:449: the Thanasilp entry has no volume and no link, unlike its neighbours. Title, authors and arXiv id are right (arxiv.org/html/2208.11060v2); journal volume and DOI verified (Nature Communications 15 (2024), doi:10.1038/s41467-024-49287-w, via ideas.repec.org); article number 5200 is from memory, not re-fetched. → `Thanasilp, Wang, Cerezo, Holmes, "Exponential concentration in quantum kernel methods", *Nature Communications* 15, 5200 (2024), doi:10.1038/s41467-024-49287-w, [arXiv:2208.11060](https://arxiv.org/abs/2208.11060) — 커널 집중. 등급 A.` md:455: the Valiant entry lacks `등급 A.`
- **[C8] glosses on first use** `정류점` (md:270) → `정류점(stationary point)`; `기각 표본 추출` (md:296) → `기각 표본 추출(rejection sampling)`; `커널 집중` is used at md:141 before its gloss at md:146 (fixed by B4); `구역 가중치 보존` is the author's own label — say so once, e.g. `(이 책에서는 구역 가중치 보존이라 부른다)` at md:76.
- **[C9] figure** panel (6): the label "MPS라면 χ ≥ 2^S ≈ 1e+11" is still drawn across the curve near t = 35–39 (py:778; move it above the line, e.g. a positive y offset — not tried). Panel (3): "10⁻¹¹ 아래는 바닥에 표시" holds for the HVA markers only; the ADAPT lines leave the frame (harmless). Panel (5) shows the r = 200 column of the new table (0.136, 0.966, 4.41, 5.84) and its y-label says r = 200 — consistent with md:307.
- **[C10] md:220** `네 묶음으로 본다` and md:181 `스크립트는 여기에 네 가지를 더한다` do not mention sections (E) and (F), which the text cites seven times. Add to md:220: `(E)는 연습문제 검산, (F)는 본문 주장의 보강 확인이다.`
- **[C11] md:340–342** the difference column mixes precisions (`7×10⁻¹³`, `1.0`, `0`). After B1: `1.0` is 0.9999 (2.497924 − 1.498044); GHZ rows `6×10⁻¹³` (both sectors, `c11_fix.py`).
- **[C12] md:139** `풀이 보존하는 대칭이 답의 구역과 다름` — a symmetry is not a sector. → `풀이 보존하는 대칭에서 시작 상태의 구역이 답의 구역과 다름` (included in the B4 row).
- **[C13] md:135 vs ch17:52, ch17:239** ch18 says "17장의 4·5단계를 그대로 못 쓴다" and then replaces step 3 as well (3′–6′); ch17:52 says "4·5단계 대신". Add the reason in md:135: `1\~2단계(문제 정의, 구조)는 같다. 3단계(원시 요소 대응)도 증명된 부품이 없으므로 훈련 가능성 검사로 바꾸고, 그 뒤를 합해 네 단계로 한다.` ch17:239 (four steps, "19장이 이 트랙으로 평가한 뒤 21·22·23장이 실행한다") matches md:135–142.
- **[C14] md:401** `풀의 리 대수도 HVA보다 훨씬 크다(n = 4에서 56차원 대 16차원)` — correct but understated: the algebra is so(2ⁿ⁻¹) ⊕ so(2ⁿ⁻¹), dimension 2ⁿ⁻¹(2ⁿ⁻¹ − 1) = 12, 56, 240, 992 for n = 3–6 (Pauli-string closure, `c1_dla.py`), i.e. exponential, where HVA's is n². Optional: `(n = 3, 4에서 12, 56차원으로 n에 지수로 큰다. HVA는 9, 16차원)`. Only n = 3, 4 are in the script; the n = 5, 6 values are mine.
- **[C15] md:122** `하나라도 걸리면 5단계 탈락이다` followed by `부호 문제가 없다는 것도 다항 시간 수렴까지 보장하지는 않는다` — say which wins: `stoquastic 줄에 걸리면 탈락이 아니라 "양자 몬테카를로와 먼저 비교" 표시로 읽는다.` (unverified suggestion of wording; the two sentences are each correct.)

---

## Summary

1. A1 — the 12-start re-run did not find "the same four" CCZ classes: it found one the first run missed and lost one the first run had; the two runs of the book's own script already give five, and six exist.
2. B1 — `gaussian_entropy` counts a ν = 0 mode 0, 1 or 2 times depending on rounding signs; the new GHZ₋ row is correct only because the signs happened to be (−, +). One-line fix given.
3. B2 — "1,000행 안에서 5%" rests on a 3-sample mean (0.049); 20–40 runs give 0.055–0.058, and the chapter's own exercise needs 1,500 rows for 5 %. Judge the row at 10 % like the others.
4. B4/B5 — the kernel-concentration paragraph is unconditional and sits in the wrong step of the track; the last 학습 노트 still says "HVA가 불모 고원 없이 훈련되는 것".
5. B3/B6/B7 — "독립 검토" is used as a source inside the book; ch17 and §6 say ch18 "handles" the |0…0⟩ start, which it only argues (the argument is right — verified by a Pfaffian calculation); the HVA table/figure still carry optimiser-limited values.

Central message after the edits: §1 and §7 now say what the chapter shows (HVA: no barren plateau and classically simulable for the same reason;
ADAPT: non-free-fermion pool with an exponentially large algebra, success explained by the easy target). I found no remaining overclaim there
apart from the §9 note.

## Checked and found correct

- **Stdout vs text**, all new or changed numbers: F1 counts (70/4, 36/0, 0.1464; floors 0.1464/0.2714/0.375), F2 table cells (0.049, 0.031, 0.016; 0.417, 0.195, 0.096; 1.57, 0.84, 0.34; 1.87, 0.81, 0.39; condition number 4.06×10⁵), F3 (9, 16, 25; 12, 56; degrees 2, 4, 6), F4 (1.4980/2.4979/0.561; 2.4977; 2.4979), E3 (c = 0.65, r ≈ 1,507, 289,354 = 15 %, SVD 0.14 s), 2×10⁻⁸ with threshold 10⁻⁶ (`gtol=1e-6`), 60,622 / 2,091 / 960, fractions 5/21, 15/21, 4/21, 1/21, 8/21, the r = 200 column from C3 (5 runs) and the rest from F2 (3 runs) as md:307 says. Full re-run reproduces stdout and the figure.
- **DLA dimensions** (own closure with SVD rank check): {ΣZZ, ΣX} → 9, 16, 25 for n = 3, 4, 5 (the per-term generators give so(2n): 15, 28, 45); pool {Y_iZ_j} → 12, 56 for n = 3, 4. **Majorana degrees** of Y_iZ_j by explicit expansion in X-string Majoranas: 2, 4, 6 at distance 1, 2, 3; Z_j: 2j + 1; Z_jZ_{j+2}: 4.
- **F4** with an independent propagator and independent Majorana operators: S = 1.498044 / 2.497924 (|0…0⟩), 2.497720 (GHZ₊), 2.497890 (GHZ₋); max|⟨c_a⟩| = 0.5610; Γ(|0…0⟩) = average of the two sectors (10⁻¹⁵); |+…+⟩ row 2.001448 both ways; Γ(t) = e^{At}Γ(0)e^{−At} to 4×10⁻¹⁴; S(|0…0⟩) at t = 1, 3: 0.6013, 2.4364 (ch17's n = 16 values 0.601, 1.498 at t = 1, 2).
- **Stable-rank table, statistically** (20–30 runs at r = 200, 1,000; 12 at 4,000; 8 at 16,000): rank 40: 0.980 ± 0.065, 0.415 ± 0.031, 0.196 ± 0.010, 0.098 ± 0.004 (book 0.966, 0.417, 0.195, 0.096) — "10%에 약 16,000행" holds; rank 300: 4.56, 1.64, 0.77, 0.335 (book 4.41, 1.57, 0.84, 0.34); tridiagonal: 5.50, 2.04, 0.825, 0.394 (book 5.84, 1.87, 0.81, 0.39); rank 5: 0.130, 0.058, 0.029, 0.014 (book 0.136, 0.049, 0.031, 0.016). The upper-bound estimates 9×10², 1.6×10⁵, 9×10⁶, 1.4×10⁷ exceed the measured needs (about 330 and 16,000 rows) as the text says.
- **§3(나) and the trap table**: (i) g_P = 0 when Q|ψ⟩ = ±|ψ⟩, [H, Q] = 0, {P, Q} = 0 (⟨ψ|PH|ψ⟩ connects different sectors); (ii) sector weights conserved when every P commutes with Q; row 1 = (i) (Y_i anticommutes with ∏X, |+…+⟩ in the +1 sector), row 3 = (ii); Y_iX_j anticommutes with ∏X. Attribution of "symmetry roadblock" to case (i) agrees with the abstract of arXiv:2109.05340 as quoted by R5.
- **HVA**: thresholds p = 9 for n = 6 (below 10⁻¹³) and not p = 8 (4.3×10⁻⁷ in 200 starts); the warm start is a stationary point (both new derivatives equal ±⟨i[ΣX, ΣZZ]⟩, which is the last-β derivative).
- **Cross-chapter**: ch17:150 (n = 1000, t = 40, about 10 s, |+…+⟩ — delivered), ch17:239 (four steps, 19장 → 21·22·23장 — matches md:135–142), ch17:11, ch17:192 ("자유 페르미온 모형은 제외(18장)"), ch15:43 (Cerezo and 커널 집중 are both in ch18 §4), ch15:60 and ch15:111 (boundary tests exist in §4/§5).
- **Form**: no garbled Hangul (rare-syllable scan of the whole file; "놓치는" is fixed), every table has a consistent column count, figure references (1)–(6) and script references E1–E3, F1–F3 all resolve, ten sections present, six 학습 노트 with "깨지는 곳".

## Could not check

- PennyLane and CUDA-Q snippets (not runnable here); the `drain_pool` source-code statement was taken from R5.
- Whether a seventh 6-CNOT class or a 7-CNOT linear CCZ exists (my check covered the six named classes only; the linear case rests on the script's 12 starts and R5's 40).
- The article number 5200 of Thanasilp et al. (PMC page blocked by a captcha); the Shkolnikov paper beyond what R5 quoted.
- Whether the entropy of the |0…0⟩ quench (ch17's tabulated quantity) has a polynomial-time formula; I verified expectation values only.
- The rate "10–70 %" per start is from 40 starts per class (binomial error about ±0.07); the HVA p = 7 value needed more than my 250 starts to reach R5's 1.8×10⁻⁵.
- ch19 (a directory appeared during this review, script and figure only, no text yet) was not examined.
