# R15 — verification pass on the rewritten chapter 21 (`ch21/ch21_hp_folding.*`)

Reviewer: R15 · 2026-11-12 (UTC) · read-only on the repository. Own scratch code and outputs: `/tmp/claude-0/-home-claude/fd63146d-511d-5ae9-8dab-3100fb85f0db/scratchpad/r15/` (`lib.py`, `a_small.py`, `b_n12.py`, `c_depth.py`, `d_replica.py`, `*.out`). The chapter script was read, never imported or run.
Line numbers: `md:N` = `ch21_hp_folding.md` (416 lines, mtime 00:43:16), `py:N` = `ch21_hp_folding.py`, `out:N` = `ch21_stdout.txt`.

**Counts: A 2 · B 7 · C 15.**

## Verdict

The rewrite is sound in substance. Every number in the md that the script prints matches `ch21_stdout.txt` (no mismatch found in sections 1, 4, 5, 6, 7, 8, 9, 10), the R13 sign defect is fixed and guarded, and the three-protocol verdict structure is honest in direction: 방식 1 fails in every seed anyone has run (15 seeds: the chapter's 5, R13's 4, my 6), 방식 3 does meet the criterion literally (I reproduced 5,963 with my own p\*-maximisation at N = 10), and the four reasons for not calling it 통과 are fair.

Two things are still wrong or misleading:

1. "독립 검토 전에 정해 둔 방식" is not fully true of 방식 1 as run: the start point (sign), the number of seeds, the summary statistic, and the SA baseline's run count and seeds were all set after the review, and the md never says so (A1).
2. The N = 12, p = 3 "(가) ⟨C⟩ 최소화" value is not converged (80-evaluation cap). The converged value is 0.3863 / 78.8× / 15,321, not 0.3873 / 76.7× / 15,734, and that flips the sign of the exercise-3 answer "−134 … 훈련에 쓸 수 있는 샷이 없다" (A2).

Neither changes the registered verdict (불통과 under 방식 1). The 방식 1 headline "5.8배" rests on five seeds that turn out to be on the favourable side (B1).

---

## A — wrong number/claim, or misleading about the verdict

### A1. 방식 1 is labelled "독립 검토 전에 정해 둔 방식", but its start point, seed count and statistic (and the SA baseline's sample) were changed after the review; not disclosed

- Location: md:15 (table, "독립 검토 전에 정해 둔 방식"), md:127 (same phrase), md:129 ("p, 샷 수, 반복 수, 이득은 이 장을 쓰면서 정했고 결과를 보고 바꾸지 않았다"), md:266 ("검토 전에 정한 방식"), md:272, md:274 (lists "출발점" among the chapter's choices without saying it changed); out:95, out:105; py:374 (`ramp(p)` as SPSA start), py:404 (`SEEDS`), py:217 (`RUNS = 4000`).
- What is true: p = 2, 64 shots × 2 × 20 and the gains are the same as in the first draft (R13 B1 quotes the same a\_k, c\_k). What changed after R13, with results in hand: (i) the SPSA start `ramp(2)` (β sign; start p\* went from 0.24× to 3.42× random at N = 12), (ii) one seed → 5 seeds (N = 12) / 20 seeds, with the median as the verdict statistic, (iii) SA baseline 2,000 → 4,000 runs with new seeds, best L 1000 → 300, value 16,045 → 15,601. The first-draft 방식 1 number was 196,189 (one seed, 12× SA; R13 B1); it is now 90,328 (median of 5, 5.8×).
- The md carefully lists "p, 샷 수, 반복 수, 이득" (omitting 출발점) at md:129 and mentions the wrong start in §7, so nothing is literally false, but a reader of §1 and §5 (6) takes "검토 전에 정한 방식" to cover the whole procedure that produced 90,328. For a chapter whose main lesson is pre-registration hygiene this has to be explicit.
- The direction of the verdict is unaffected (first draft 12×; now 5.1–13.9×; all 15 known seeds fail).
- Proposed text, md:129, replacing the sentence "p, 샷 수, 반복 수, 이득은 … 확인할 수는 없다).":
  > p, 샷 수, 반복 수, 이득은 이 장을 쓰면서 정했고 결과를 보고 바꾸지 않았다(19장의 등록에는 없는 값이라 저장소의 기록으로 확인할 수는 없다). 검토 뒤에 바뀐 것은 셋이다. 출발점의 부호(첫 판은 β도 양수였다. 7절), 시드 수(첫 판은 시드 하나였고 잔기 12개의 값은 196,189였다. 지금은 5개·20개의 중앙값과 범위다), 담금질의 실행 수와 시드(2,000회 → 4,000회. 첫 판의 값은 16,045였다). 그러니 방식 1은 "검토 전에 정한 설정을, 검토 뒤에 고친 출발점에서 다시 돌린 것"이다. 불통과라는 방향은 첫 판(담금질의 12배)과 지금(5.1\~13.9배)이 같다.
- md:15, md:127, md:266: "독립 검토 전에 정해 둔 방식" → "설정(p, 샷, 반복, 이득)은 독립 검토 전에 정한 것. 출발점의 부호는 검토 뒤에 고쳤다(4절)". md:274: after "훈련 방법(SPSA, p = 2, 64샷 × 20회, 이득, 출발점)" add ", 시드 수와 중앙값으로 적는 것". (196,189 / 16,045 / 2,000 are R13's record of the first draft; I did not see the first draft.)

### A2. N = 12, p = 3 "(가) ⟨C⟩ 최소화" is not converged; the converged value flips the sign in exercise 3 and in the (6) comparison

- Location: md:219 (+0.3873, 0.4917, 2.926×10⁻⁴, 76.71, 15,734), md:30 ("⟨C⟩를 정확히 최소화하면 … p = 3이면 76.7배"), md:32, md:223, md:260, md:272 ("76.7배"), md:376 ("p = 3의 ⟨C⟩-최적 값(… 76.7배)"), md:378 ("12개는 음수(−134)다 … 훈련에 쓸 수 있는 샷이 없다"), md:394; out:56, out:98, out:126; py:294 (`BUDGET_12 = {…, 3: 80}`), assert py:608 (pins 76.71 to 1%).
- Evidence (`b_n12.py refine`): Nelder–Mead on ⟨C⟩ from the chapter's printed N = 12 angles, 185 further evaluations: ⟨C⟩ = 0.3863, P(valid) = 0.4975, p\* = 3.005×10⁻⁴ = 78.78× random, TTS₉₉ = 15,321, angles (+0.531 +1.032 +1.976 \| −0.158 −0.137 −0.116). This is exactly the value R13 reported as converged (R13 A1 table, A3) — the rewrite re-ran with an 80-evaluation cap (six parameters; seven evaluations go to the initial simplex) and landed short of it.
- Consequences:
  - Zero-training break-even vs SA is 77.4× (out:99). The chapter's 76.7× is just below, the converged 78.8× just above. 15,321 < 15,601 < 15,734.
  - Exercise 3 (md:378): SA − QAOA is −134 with the chapter's value and +280 with the converged one. "12개는 음수(−134)다 … 훈련에 쓸 수 있는 샷이 없다" is an artefact of the cap.
  - Either difference is far inside the SA sampling error: 134 is 0.16 s.e. and 280 is 0.34 s.e. of p\_L (s.e. range 14,795–16,495). So "표본 오차 안에서 같다" (md:32, md:223) is right; only the categorical sign statements are not.
  - The SA figure itself is the minimum over six noisy L estimates, so it is biased low: my own SA (4,000 runs, own seeds) gives 16,148 and 16,843 at L = 300 and 16,064 at L = 1000; R13 had 15,553–15,806; the chapter's own L = 1000 is 16,483.
- Fix (preferred): raise the p = 3 budget to about 270 evaluations (≈ 190 more state-vector evaluations, roughly 5–8 min), re-pin py:608, and rewrite from the new output. Minimum fix without a re-run of the search: say it. Proposed text for md:221, replacing "12개는 평가를 27·52·82회로 묶었다.":
  > 12개는 탐색을 25·50·80회에서 끊었고(표의 값은 찾은 값이지 최솟값이 아니다), p = 3은 아직 덜 내려갔다. 검토 R13과 R15가 같은 각도에서 더 돌린 값은 ⟨C⟩ = 0.3863, 무작위의 78.8배, TTS₉₉ 15,321이다. 담금질(15,601, 표준오차로 14,795\~16,495)과의 차는 어느 값으로도 표본 오차 안이고, 어느 쪽이 짧은지는 정해지지 않는다.
- Proposed text for md:378, replacing "담금질에 대해서는 8개 949, 10개 498이고 12개는 음수(−134)다. 즉 … 샷이 없다.":
  > 담금질에 대해서는 8개 949, 10개 498이다. 12개는 −134로 나오지만 이 차는 담금질의 표본 오차(약 ±850)보다 훨씬 작고, ⟨C⟩ 최소화를 끝까지 돌리면 +280이 된다(5절 (5)). 즉 잔기 12개에서는 이 각도를 훈련 없이 얻어도 담금질과 비기는 정도이고, 훈련에 쓸 수 있는 샷은 많아야 수백 개다.
  (the "+280", "±850" must be printed by the script before they go into the text.)
- md:30: "⟨C⟩를 정확히 최소화하면 35.4배, p = 3이면 76.7배다" → "정확한 ⟨C⟩를 줄여 찾은 값은 35.4배, p = 3이면 76.7배다" (or the converged number after a re-run). md:376: "⟨C⟩-최적 값" → "(가)에서 찾은 값".

---

## B — should fix

### B1. 방식 1's headline numbers (90,328, "5.8배", "13.8배") come from five seeds that are on the favourable side

- Location: md:15, md:30, md:248, md:266, md:270; out:80, out:95; py:404.
- Evidence: (i) `d_replica.py` replays the chapter's five seeds (`default_rng([21, 12, 2, 64, 20, sd])`, same call order, my own state-vector code): 8.03×, 14.80×, 15.62×, 5.64×, 13.75× → median 13.75× (5.64–15.62), TTS₉₉ 152,825 / 84,121 / 79,854 / 216,686 / 90,328 → median 90,328. The script output is reproduced exactly; optimum seen in training: 0 in all five. (ii) `b_n12.py spsa`: six independent seeds with the same gains, start and shots: 10.81×, 5.21×, 5.13×, 9.48×, 7.41×, 8.44× → median 7.93×, TTS₉₉ median 155,499 (114,237–237,983). R13's four seeds: 6.7×, 4.4×, 9.9×, 8.9× (125k–275k).
- Reading: three of the chapter's five seeds are above every one of the ten independent seeds. Pooled (15 seeds) the median is 8.4× → about 146,000 ≈ 9.3× SA. The code is right; the sample is small. "5.8배 길다" and "13.8배(중앙값)" understate the gap (in QAOA's favour) by roughly a factor 1.6. The verdict is unaffected: 15 of 15 seeds are longer than SA (minimum 79,854 vs the SA upper s.e. bound 16,495).
- Fix: run more seeds at N = 12 (15–20; 41 evaluations per seed, about 1 min each) or keep five and lead with the range. Proposed, md:15 column "담금질 대비": "5.1\~13.9배 길다(시드 5개. 중앙값 5.8배)". md:248 or md:250 add:
  > 시드 5개의 중앙값은 흔들린다. 검토 R13·R15가 같은 설정으로 따로 돌린 시드 10개는 4.4\~10.8배(중앙값 7.9배, TTS₉₉ 약 155,000)였다. 어느 시드도 담금질보다 짧지 않다.

### B2. §1 bullet "훈련 비용을 0으로 치면 QAOA는 담금질과 같거나 앞선다" drops the tuned-SA clause that §5 (4) and §7 supply

- Location: md:32 vs md:203, md:358.
- 15,734 (or 15,321) is 1.6× the tuned SA (9,894); 방식 2's 17,114 is 1.7×. Only 방식 3 and (나) are below 9,894. R13's proposed text for this bullet had the clause; the rewrite lost it. As it stands §1 says "같거나 앞선다" and §7 says "담금질의 일정을 맞추면 15,601이 9,894가 된다" without connecting them.
- Proposed (end of md:32): "… 12개에서는 15,734로 담금질(15,601)과 표본 오차 안에서 같다. 일정을 맞춘 담금질(9,894)보다는 1.6배 길다." (print 15,734/9,894 = 1.59 in the script.) Title: "**훈련 비용을 0으로 치면 QAOA는 일정을 고정한 담금질과 같거나 앞선다.**"

### B3. §1 never states the registered verdict in one sentence

- Location: md:11–19. The reader gets a table with one "**불통과**", one "불통과(근소한 차)" and one "통과", then "방식 2와 3은 … 사전 등록의 판정으로 쓰지 않는다". The sentence "사전 등록의 판정은 방식 1로 적는다: **불통과**" exists only at md:272.
- Proposed, first sentence of md:19: "사전 등록의 판정은 방식 1로 적는다: **불통과**. 방식 1은 무작위 추출보다 13.4배 짧지만 …".

### B4. "2,560샷의 훈련은 각도를 거의 찾지 못한다" overstates

- Location: md:30. Training raises p\* from 3.42× (start, out:80) to 13.75× (chapter median; 7.9× on independent seeds) at N = 12, and 2.6× → 4.2× / 5.1× at N = 8 / 10; the ⟨C⟩-optimal p = 2 values are 35.4× / 9.7× / 13.6×. So training gets a factor 1.6–4 and ends at 22–43% of the ⟨C⟩-optimal p\*. "거의 찾지 못한다" reads as "no better than the start".
- Proposed: "**2,560샷의 훈련은 좋은 각도의 절반에도 못 미친다.** 잔기 12개에서 한 샷 성공 확률은 출발점의 3.4배에서 훈련 뒤 13.8배(시드 5개의 중앙값, 5.6\~15.6배)로 오르지만, 같은 p = 2 회로에서 정확한 ⟨C⟩를 줄여 찾은 값은 35.4배, p = 3이면 76.7배다. 좋은 각도에서 출발해 샷으로 더 훈련하면 중앙값이 내려갔다(잔기 10개: 19.3배 → 16.2배)."

### B5. 방식 1 covers p = 2 only although the registered hypothesis says p ≤ 3; the stated reason is weak and a bound that closes the gap is available from printed numbers

- Location: md:129 ("잔기 12개에서 p = 3의 샷 훈련은 돌리지 않았다(평가 비용 때문이다 …)"), md:362.
- Five p = 3 seeds cost about 200 evaluations (the script already spends 456 at N = 12). More useful: with 2,560 training shots QAOA needs p\* ≥ 92.6× random to beat SA (out:99); the best ⟨C⟩ value found at p = 3 is 76.7× (78.8× converged). So no ⟨C⟩-minimising shot training at p ≤ 3 with this budget can pass, unless it lands far above the ⟨C⟩ optimum by luck. R13 ran three p = 3 seeds: 133k, 70k, 157k.
- Proposed (md:129): "잔기 12개에서 p = 3의 샷 훈련은 돌리지 않았다. 판정은 그래도 바뀌지 않는다. 2,560샷을 쓴 뒤 담금질을 이기려면 p\*가 무작위의 92.6배여야 하는데 p = 3에서 ⟨C⟩를 줄여 찾은 값이 76.7배다(5절 (6)). 검토 R13이 따로 돌린 p = 3의 세 시드는 70,000\~157,000이었다."

### B6. "온도 일정은 결과를 보기 전에 정했고 바꾸지 않았다" is still an unverifiable flat claim (R13 B4, half fixed)

- Location: md:203; out:38. The SPSA sentence got the caveat "(… 저장소의 기록으로 확인할 수는 없다)"; the SA sentence did not.
- Proposed: "온도 일정은 이 장을 쓰면서 정했고 결과를 보고 바꾸지 않았다(기록으로 확인할 수는 없다)."

### B7. The angles behind 방식 3 (the only "통과" cell) are printed nowhere

- Location: out:58 prints only the ⟨C⟩-optimal p = 3 angles; md:221 quotes N = 8 and N = 12 of those. The p\*-optimal N = 10 angles that give 5,963 are in neither file, so the chapter's most surprising number cannot be checked without the 21-minute run.
- My own search (Nelder–Mead on −p\* at N = 10 from the printed ⟨C⟩-optimal angles, 429 evaluations): (+1.059 +2.000 +2.982 \| −0.233 −0.175 −0.137); N = 10: 29.91×, TTS₉₉ 105.0; at N = 12 unchanged: ⟨C⟩ = 1.4992, P(valid) = 0.5409, p\* = 7.7205×10⁻⁴, 202.39×, TTS₉₉ 5,962.5 — identical to out:74.
- Fix: print `ORACLE[(10, 3)]["x"]` (and `ORACLE[(8, 3)]["x"]`) next to out:58 and quote them at md:240.

---

## C — minor / wording

- **C1.** md:139, out:104 "14장의 오류 예산(G ≈ 100\~170)": I find no 170 in `ch14/`. ch14:198 says "CZ 수의 한계 ≈ 100 → 0.3 오류 예산" (170 would be rG = 0.5). Either write "14장의 한계(CZ ≈ 100개, 오류 예산 0.3)의 184배" (18,388/100; have the script print it) or say where 170 comes from.
- **C2.** md:374 "227.95이고 TTS₉₉는 … 227.9다" — same number at two precisions (R13 C6, not fixed).
- **C3.** "걸음" for an SA step survives in the figure, panel (3) x-label "담금질 한 번의 길이 L (걸음 = 평가)" (py:561), and in py:6, py:201. md:83 uses "한 걸음" for a Trotter step — a third meaning; "한 단계" is safer.
- **C4.** md:221 "평가를 27·52·82회로 묶었다": the cap is 25·50·80 (py:294); 27·52·82 are the counted evaluations including the transferred-angle evaluation and the reporting evaluation.
- **C5.** md:270 "세 방식 모두 TTS₉₉가 예산 10⁷ 평가 안이다": out:101 prints 무작위 / 담금질 / QAOA 방식 1. Write "세 방법 모두(무작위 추출, 담금질, QAOA 방식 1)".
- **C6.** md:252 "이 샷 수의 훈련은 잡음이 기울기보다 커서 …": a reading, not shown by the script → "… 큰 것으로 보인다". md:30 "오히려 나빠졌다": it is the median (16.2×); the best of five seeds is 19.6× > 19.3×.
- **C7.** Numbers not printed by the script (all recomputed, none wrong): "97.13%" (md:356), "97%" (md:390), "세 배 가깝다" (md:240; 202.4/70.5 = 2.87), "1.5배" and the three first-draft claims (md:346; first-draft numbers, from R13), "수천 번" (md:24), "e^{−55}" (md:352), "262.8배 / 360회" (md:225; R13's number, labelled as such).
- **C8.** md:378: "(반올림 전 값으로 계산)" is attached to 2,927 only; 949 (1,163 − 215 = 948) and −134 (15,601 − 15,734 = −133) need it too. Put the note once at the start of the answer.
- **C9.** md:19 puts a paraphrase of the ch19 expectation in quotation marks. ch19:306 reads "불통과(담금질에 진다). 잔기 8·10개에서는 … 훈련에 수천 샷을 쓰는 순간 진다. 12개(무작위 1.2×10⁶)에서는 무작위보다 짧을 수 있다". Drop the quotation marks or quote verbatim. The claim "이 방식에서 맞았다" itself holds (N = 12: fails vs SA, 13.4× shorter than random; N = 8: 3,305 > 3,141; N = 10: 3,187 vs 3,208, a tie).
- **C10.** §5 (6) presents 방식 2 and 3 without saying they are the best of the eight transfer variants in the (전이) table. The p = 2 versions of 방식 3 fail (17,502 and 18,869). One sentence at md:272: "전이 여덟 가지 가운데 담금질보다 짧은 것은 p\*로 정한 p = 3의 둘(7,458, 5,963)뿐이다."
- **C11.** The pre-computation argument is fair as a lower bound ("표의 칸 수만으로도"); it would be tighter with the sum: 278,528 + 5,963 = 284,491 = 18.2× SA, 4.2× shorter than random, and longer than 방식 1 (90,328). Optional; needs a script print. Note for balance: SA's own L selection (six lengths × 4,000 runs ≈ 5.8×10⁷ evaluations, using the N = 12 answer) is uncharged too; the md says so at md:203, and SA at any L from 300 to 10,000 stays within 15.6k–20.6k, so the conclusion does not depend on it.
- **C12.** md:139 "(T ≈ 15,846)" is unexplained in this chapter (3·log₂(1/ε)·R with ε = 10⁻⁶; 12장 4절). Add "(12장의 어림)".
- **C13.** md:131 "항의 순서를 잘 잡으면 CX는 항 하나에 하나꼴까지 줄일 수 있고" has no source and no check. Add a reference or "…줄일 수 있다고 알려져 있고(이 장은 확인하지 않았다)".
- **C14.** md:81 and md:394 state angle transfer as a general property ("…거의 그대로 통한다"). Evidence is one sequence whose small instances are prefixes of the large one (§7 says so). Add "이 서열에서는".
- **C15.** Regression asserts (R13 B13, partly done): not pinned are the four p = 2 transfer rows, the 방식 1 median ratio 13.75, the ten sensitivity rows, the warm-start range, SA p\_L. py:608 pins the unconverged 76.71 (A2). py:1 says "실행 약 20분", md:415 "약 21분".

---

## R13 item → status

| R13 | Status | Where / remark |
| --- | --- | --- |
| A1 `ramp` sign | resolved | py:269–274; guard py:289–290; N = 8 multi-start py:304–306; p = 1 grid py:319–323; md:83, md:346. I confirm: correct-sign `ramp(2)` gives ⟨C⟩ 2.008, 2.62× at N = 8; wrong sign 5.742, 0.27× |
| A2 "p = 1 깎인다" | resolved | removed; md:225 gives 4.78× vs 5.65×, 6.36× vs 8.19×; listed as an error at md:346 |
| A3 "훈련 0이어도 담금질이 앞선다" | partly | md:32, md:223 rewritten; the tuned-SA clause is missing (B2); N = 12 value unconverged (A2) |
| A4 (나) at N = 12 | partly | re-run from the right basin: 255.6× with 80 evaluations (R13 ≥ 262.8× with 360); disclosed at md:225. Thresholds now 77.4× / 92.6× |
| A5 verdict per protocol | resolved | md:11–26, md:264–274; explicit decision (방식 1 = registered verdict, 2·3 post hoc). But see A1, B3 |
| A6 Perdomo-Ortiz 81 qubits | resolved | md:382, md:404 |
| B1 SPSA one seed, undisclosed start/gains | partly | seeds 20/20/5, start p\* and gains printed (md:129, md:244–248); p = 3 at N = 12 not run (B5); five seeds unrepresentative (B1); start change undisclosed (A1) |
| B2 "같은 회로 p = 3" | resolved | md:30 |
| B3 SA disclosures | resolved | md:101, md:195, md:201, md:203, md:286; seeds `default_rng([N, L])` py:224 |
| B4 unverifiable "실행 전에 정했다" | partly | SPSA caveat added (md:129); SA sentence not (B6); post-review changes not listed (A1) |
| B5 λ sweep start | resolved | py:471; values equal R13's (17.88 / 18.28 / 13.56 / 9.77 / 8.00); "때문으로 보인다" md:286 |
| B6 "주로 한 일은 벌점 피하기" | resolved | decomposition printed; md:223, md:225, md:356, md:390 |
| B7 symmetry sentence | resolved | md:193; py:170–179; reproduced (below) |
| B8 QUBO / gate wording / gate\_cost check | resolved | md:25, md:131, md:193, md:350; py:518–531 (2,050 / 2,277 / 265 / 3.4×10⁻¹⁷ reproduced) |
| B9 numbers not printed | mostly | remaining ones in C7, C1 |
| B10 noise remark as prediction | resolved | md:352 |
| B11 PERM / replica exchange | resolved | md:356, md:358, md:410 |
| B12 snippets | resolved | 255 terms (md:304), `m` rename (md:331), fallback text (md:342), 10-qubit-only note (md:294), lightning "확인하지 못했다" (md:301), ch19 reason (md:141) |
| B13 regression asserts | partly | re-pinned, guards added, `rng` now used; gaps in C15; pins the unconverged value (A2) |
| C1 exercise provenance | resolved | md:366; problem 5 added (N&C 6장) |
| C2 "걸음" | partly | md uses "스텝"; figure label and py docstrings not (C3) |
| C3 seeds | resolved | md:415 matches py:2, py:224, py:398 |
| C4 layer order | resolved | md:81 |
| C5 evaluation count | partly | C4 above |
| C6 227.95 / 227.9 | not | C2 above |
| C7 figure | resolved | panel (2) title and md:288; transfer points in panel (4) |
| C8 "예산 안에서 본다" | resolved | md:270 (wording nit C5) |
| C9 Boulebnane | resolved | md:348 |
| C10 Fingerhuth | resolved | md:356 |
| C11 P ≥ 0.99 rule | resolved | md:101 |
| C12 "거의 같아진다" | resolved | md:31 |
| C13 j − i ≥ 3 | resolved | md:58 |
| C14 "완벽했을 때" | resolved | md:124 ("…어디까지 가는가"의 참고값), out:44 ("찾은 최솟값") |
| C15 λ > 0 | resolved | md:183 |

No sentence from the first draft survives that contradicts the new results (searched for 이상적, 24,724, 16,045, 196,189, 절반쯤, 깎, 90배, 81.2, 14,866, 9,865, 상한, 완벽; the only hits are the deliberate account at md:346 and "증명된 상한도 아니다" at md:125).

---

## Independent recomputations (own code, `scratchpad/r15/`)

| Quantity | md / stdout | R15 | Match |
| --- | --- | --- | --- |
| N = 8 / 10 / 12 tables: valid, E\_min, optima, mean, σ | 2,172 / −2 / 24 / 3.562 / 3.073; 16,268 / −2 / 376 / 5.252 / 4.220; 120,292 / −5 / 16 / 6.348 / 4.990 | same (all-pairs table, no parity shortcut) | ✓ |
| N = 8, p = 3 at the printed angles: ⟨C⟩ / P(valid) / p\* | 0.1031 / 0.8016 / 2.121×10⁻² | 0.10313 / 0.80125 / 2.116×10⁻² at the 3-decimal angles; 0.10311 / 0.80164 / 2.121×10⁻² (14.478×, TTS 214.8) after refining from them; two mixer implementations agree | ✓ |
| N = 8, p = 1 global minimum | grid 1.2992, NM 1.2965 | 361 × 91 grid 1.2989 at (−0.912, +0.209) — the conjugate of (+, −), same probabilities | ✓ |
| Z terms, max degree, degree ≤ 2, CX (N = 8; 10) | 1,782 / 12 / 38 / 18,388; 28,671 / 16 / 70 / 409,602 | same; N = 8 histogram 38, 393, 804, 485, 60, 2 | ✓ |
| Symmetry: C invariant under all-bit flip and high-bit flip; both parities obeyed; allowed terms | 4^(N−2) − 1; alive 43.5% / 43.7% | invariances true (also low-bit flip); brute-force count of allowed masks 4,095 and 65,535; alive 43.52% / 43.75%; all-flip = right↔down, up↔left; high-flip = right↔left, up↔down | ✓ |
| N = 10, p = 3 printed ⟨C⟩ angles | 0.3578 / 0.6631 / 2.945×10⁻² / 20.53 / 154 | 0.35776 / 0.66327 / 2.946×10⁻² / 20.54 / 154.0 | ✓ |
| N = 10, p = 3 p\*-maximum | 4.290×10⁻², 29.9×, 105; 10.7 × 2.8 | 4.290×10⁻², 29.91×, 105.0; 10.7 × 2.8 | ✓ |
| 방식 3 (those angles at N = 12) | 1.4992 / 0.5409 / 7.720×10⁻⁴ / 202.4 / 5,963 | 1.4992 / 0.5409 / 7.7205×10⁻⁴ / 202.39 / 5,962.5 | ✓ |
| 방식 2 (printed N = 10 angles at N = 12) | 0.4178 / 0.4703 / 2.690×10⁻⁴ / 70.5 / 17,114 | 0.4181 / 0.4705 / 2.693×10⁻⁴ / 70.59 / 17,100 (3-decimal angles) | ✓ |
| N = 8 ⟨C⟩ angles at N = 12, p = 3 | 62.1× / 19,453 | 61.9× / 19,503 (3-decimal angles) | ✓ |
| N = 12, p = 3 at the printed angles | 0.3873 / 0.4917 / 2.926×10⁻⁴ / 76.71 / 15,734 | 0.3874 / 0.4912 / 2.917×10⁻⁴ / 76.47 / 15,784 (3-decimal angles; chunk and FWHT mixers agree) | ✓ as printed |
| N = 12, p = 3 ⟨C⟩ minimum (converged) | same row | 0.3863 / 0.4975 / 3.005×10⁻⁴ / 78.78 / 15,321 | ✗ (A2) |
| 방식 1, the chapter's five seeds replayed | 13.75× (5.64–15.62); 90,328 (79,854–216,686) | identical | ✓ |
| 방식 1, six independent seeds | — | 7.93× (5.13–10.81); 155,499 (114,237–237,983) | see B1 |
| SA N = 12 (4,000 runs, own seeds) | 15,601 (L = 300; 14,795–16,495) | L = 300: 16,148, 16,843; L = 1000: 16,064; T = 0.5 fixed, L = 300: 11,231 | ✓ within about 1.5 s.e.; chapter value on the low side |
| s.e. range of 15,601; thresholds | 14,795–16,495; 77.4× / 92.6× | same; ±2 s.e. 14,065–17,493; 17,114 sits at 1.63 s.e. | ✓ |
| Ratios | 5.8 (5.1–13.9), 13.4, 1.10, 0.38, 0.60, 18 | 5.79 (5.12–13.89), 13.36, 1.097, 0.382, 0.603, 17.85 | ✓ |
| Exercise 5 (Grover) | θ 0.001953, 402, 0.999998, 1/39, 1/3,003 | θ 0.0019531, π/4θ = 402.12 → 402 (round(π/4θ − ½) also 402), 0.999998, 38.8, 3,003.0 | ✓ |
| Exercise 1 (HPPH) | 64 / 36 / 8 / 28 / 34.5 | same | ✓ |
| Case file: table, σ, 8 optima, 256 terms | json | identical to my N = 6 table; γ = +0.8199, β = −0.2145 | ✓ |
| Case: ⟨C⟩ / p\* (diagonal) | 0.2904 / 0.0313 (random 0.0078) | 0.290386 / 0.031290 (0.0078) | ✓ |
| CUDA-Q kernel emulated exactly as written | CX 2,050, rz 255, rx 10; "0.29 근처", "±0.001" | CX 2,050, rz 255; same distribution to 3.0×10⁻¹⁷; 20,000-shot s.d. 0.0095 on ⟨C⟩, 0.0012 on p\* | ✓ |
| Gate list of the case (ASAP depth, T estimate) | 2,050 / 2,277 / 265 / T ≈ 15,846 | 2,050 / 2,277 / 265 / 15,846 | ✓ |

## Verified OK

- **Numbers.** Every table cell and in-text number in §1, §4, §5 (1)–(7), §6, §7, §8, §9, §10 that the script prints equals `ch21_stdout.txt` (manual pass plus an automatic token check of all numbers ≥ 3 digits or with decimals). The numbers the script does not print are listed in C7 and C1. No rounding inconsistency other than C2 and C8.
- **Stdout is current.** `ch21_stdout.txt` (00:37:46) is older than the last edit of the py (00:38:07), but the author's confirming re-run in `scratchpad/c21/rerun/` finished during this review with `rc=0`, `cmp=0`; I diffed its `out.txt` against `ch21_stdout.txt` myself (timing lines excluded): identical. All regression asserts therefore pass on the current py. html (00:43:17) is newer than the md (00:43:16) and contains the new numbers.
- **Verdict wording.** "표본 오차 안에서 같다" for 15,734 vs 15,601 (0.16 s.e.); "불통과(근소한 차)" for 17,114 (1.10×, 1.63 s.e. above, and above all ten SA estimates I know of — the chapter's two drafts, R13 and mine at L = 300 and 1000: 15,479–16,899); "5.8배" as arithmetic (but B1); "13.4배"; "예상은 이 방식에서 맞았다"; "방식 3이 일정을 맞춘 담금질(9,894)보다도 짧다"; the four reasons at md:21–26; the 18× pre-computation argument (a lower bound, stated as such; conclusion robust — C11). 방식 1 fails in 15 of 15 seeds.
- **Snippets (not executable here).** PennyLane: keys read from the case file exist; `RX(2β)`; expected 0.2904 / 0.0313 / 0.0078 correct. CUDA-Q: 255 terms, `m` rename, ladder / un-ladder order, `rz(2γw/σ)`, `rx(2β)`, big-endian mask, expected-output comments — consistent with the text and with the script's gate-list check. Status lines follow CLAUDE.md.
- **Format.** Ten sections in template order; `\~ \_ \*` escaped everywhere outside code; no bare `|` outside tables, `\|` used in running text (md:83, md:221, md:392); no paragraph starts with `|`; 13 tables with consistent column counts; "실측" appears once, referring to the ch14 rule (md:395), never for numpy values; "numpy 구현, 실행 완료" at md:111; one other SDK with a one-sentence reason (md:294).
- **Figure.** Six panels match md:288; panel (1) shows H at 0, 2, 5, 6, 8, 11 and five contacts; panel (2) title states the C ≤ 25 cut.
- **Environment line** md:415: Python 3.13.17, numpy 2.5.3, scipy 1.18.1, matplotlib 3.11.2 are what is installed; 1,286 s ≈ 21 min; 456 evaluations at N = 12.

## Could not verify

- Execution of either SDK snippet (no SDK).
- The first-draft values quoted in A1 (196,189; 16,045; 2,000 runs) — taken from R13, the first draft is no longer on disk.
- Whether the SA schedule and the SPSA p / shots / iterations / gains were fixed before any result was seen (B6); R13's quotation of the first draft shows they did not change between drafts.
- Global optimality of any N ≥ 10 QAOA value (local searches; my refinements only show the N = 12, p = 3 value was not converged).
