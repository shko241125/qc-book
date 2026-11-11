# R8 — ch00 "책 앞에: 기호와 비용 단위" (first review) and ch19 "후보 문제 선정" (second review)

Reviewed files
- Target 1: `/home/claude/qc-book/ch00/ch00_notation.md` (66 lines), checked cell by cell against `ch01`–`ch19` md files.
- Target 2: `/home/claude/qc-book/ch19/ch19_candidates.md`, `ch19_candidates.py`, `ch19_stdout.txt`, against the first review (`R7_findings.md`), the edit logs (`c19/r7_log.json`, 45 edits; `c19/q4_log.json`, 4 edits), ch17, ch18.

Scratch code and outputs: `/tmp/claude-0/-home-claude/fd63146d-511d-5ae9-8dab-3100fb85f0db/scratchpad/mid/R8/`
(`run/rerun.txt` fresh run of the ch19 script, `schw.py`, `kern.py`, `prereg23a.py` + `prereg23a_n20.txt`, `prereg20.py`, `pitch.py`, `levelstat.py`).

Not reported, as instructed: forward pointers to ch20–22; anything in section 6 or in a "게이트 비용 기록" paragraph of ch10–13.

Counts — Target 1: **A 3, B 8, C 9**. Target 2: **A 1, B 5, C 3**.

---

# Part 1 — Target 1: ch00_notation.md

## A — factually wrong

### T1-A1. ch00 §2, φ row — ch00_notation.md:25
Quote: "6·19장에서 φ(x)는 특성 맵이 만드는 상태"

Wrong chapter. ch06 does not contain the letter φ at all (`grep -c φ ch06_circuit_library.md` = 0). ch06 writes the feature state with a capital: "K(x, y) = |⟨Φ(x)|Φ(y)⟩|²" (ch06:45, ch06:235). Only ch19:66 writes "|⟨φ(x)|φ(y)⟩|²". So the two chapters the row cites disagree with each other, and the row hides it.
The clashes a reader actually meets are missing: Bloch-sphere longitude φ (ch01:40, ch01:48), the gate angle in U(θ, φ, λ) (ch02:60), pulse phase φ (ch07:39), the teleported state |φ⟩ (ch05:39), and the QSP phase angles φ₀…φ_d (ch16:64, 16 occurrences in the same chapter that also uses QPE).

Replacement row:

> | φ | QPE가 읽는 위상. 고유값이 e^{2πiφ}이고 φ는 0 이상 1 미만 | 10\~12장 | 1·2장에서는 Bloch 구의 경도와 U(θ, φ, λ)의 각, 7장에서는 펄스의 위상, 16장에서는 QSP 위상각 φ₀…φ\_d다. 5장에서는 순간이동으로 보낼 상태, 19장에서는 특성 맵이 만드는 상태 φ(x)의 이름이기도 하다(6장은 같은 상태를 Φ(x)로 적는다) |

(Alternatively change ch19:66 to Φ(x) so that ch06, ch19 and ch22 agree, and keep the row short.)

### T1-A2. ch00 §2 G row and §3 "rG" row — ch00_notation.md:37, :53
Quote (37): "| G | 회로의 2큐비트 게이트 수. rG는 오류의 기대 개수 | 14·19장 | — |" and (53) "…표본 비용이 지수로 든다(14·19장)"

ch14 uses neither G nor rG. It writes the count as N\_CZ and the product as "N\_CZ·p\_CZ" (ch14:38 formula, ch14:41, ch14:191). The only places with G and rG are ch19:94–100 and ch19:251, 309. In addition, in that ch14 formula p\_CZ is the spec-sheet error (ch14:41 "예산(사양서 r 그대로)"), i.e. the quantity ch00 calls r, while ch00:33 says p in ch14 is the depolarizing strength. A reader going from ch00 to ch14's budget formula is sent to the wrong symbol twice. The sentence "rG가 1을 크게 넘으면 오류 완화의 표본 비용이 지수로 든다" is in ch19:100 only; ch14:41 says "N\_CZ·p\_CZ ≳ 1이면 결과가 균등 분포에 묻힌다".

Replacement row (37):

> | G | 회로의 2큐비트 게이트 수. rG는 오류의 기대 개수 | 19장 | 14장은 같은 수를 N\_CZ로, 곱을 N\_CZ·p\_CZ로 적는다. 그 식의 p\_CZ는 사양서 오류율, 곧 이 표의 r이다 |

Replacement for the last cell of row 53:

> (1 − r)^G ≈ e^{−rG}. rG가 1을 크게 넘으면 오류 완화의 표본 비용이 지수로 든다(19장). 14장의 오류 예산 N\_CZ·p\_CZ가 같은 양이고, 거기서는 1을 넘으면 결과가 균등 분포에 묻힌다고 적었다

### T1-A3. ch00 §1 units bullet — ch00_notation.md:12
Quote: "장치를 다루는 7·14장만 ns, µs를 쓴다."

False. Counts of "<digit> µs / ns / ms" per chapter: ch17 has 11 + 3 (ch17:102 "논리 Toffoli 10 µs", ch17:163, ch17:175–180 "10 ns (1코어)"), ch19 has 5 + 1 + 2 (ch19:13 "약 10 µs", ch19:193 "t\_CZ = 68 ns … 2.2 ms", ch19:282 "170 µs"), ch16:306 ("170 µs"), ch09:215 ("170 μs"). These are wall-clock numbers that carry verdicts (break-even times in ch17, the FFT comparison in ch19).

Replacement for that sentence:

> 장치와 벽시계 시간을 다루는 곳에서만 초 단위(ns, µs, ms, s)를 쓴다. 7·14장의 게이트 시간과 T1·T2, 9·16·17장의 논리 게이트 시간과 손익분기, 19장의 FFT와 회로 실행 시간이 그렇다.

## B — misleading / incomplete / contradicts a chapter

### T1-B1. ch00 §2, t row — ch00_notation.md:24
Quote: "12장부터 t는 시간이다. 계수 큐비트 수는 12장에서 t\_b, 16장에서 m, 17장에서 b로 적는다"

(i) ch12 itself still uses t for the register size outside the cost paragraph: ch12:48 "해상도는 ΔE = 2π/(τ·2^t) — … t ≈ 20비트", and the table header ch12:171 "| t\_b | Q 호출 2^t−1 (BHMT의 M = 2^t) |", where t\_b and t appear in the same row for the same quantity. So "12장에서 t\_b" is only half true. (ii) "12장부터" is too late: ch07 already uses t as time throughout (ch07:33, 39, 62). (iii) ch15 uses t for the number of T gates (ch15:56 "T가 t개인 회로", ch15:85 "T가 t개면 정확 2^t 분기"), in a chapter that also has time steps t (ch15:11). m (ch16) and b (ch17) are correct (ch16:86 "m비트 QPE", ch17:57 "N\_W = 2^{b}").

Replacement row:

> | t | 계수(위상) 레지스터의 큐비트 수. 정밀도는 2⁻ᵗ | 10·11장 | 7장과 12장 이후에서 t는 시간이다. 계수 큐비트 수는 12장에서 t\_b, 16장에서 m, 17장에서 b로 적는다. 15장의 "T가 t개"에서 t는 T 게이트의 수다 |

And in ch12 (outside the paragraphs being edited): ch12:48 "ΔE = 2π/(τ·2^t)" → "ΔE = 2π/(τ·2^{t\_b})", "t ≈ 20비트" → "t\_b ≈ 20비트"; ch12:171 header "Q 호출 2^t−1 (BHMT의 M = 2^t)" → "Q 호출 2^{t\_b}−1 (BHMT의 M = 2^{t\_b})".

### T1-B2. ch00 §2, r row — ch00_notation.md:32
Quote: "| r | Trotter 단계 수 | 12·19장 | 14·19장에서는 게이트당 오류율(무작위 벤치마킹 값)이기도 하다. …"

(i) The most frequent r in the book is missing: the order in Shor's algorithm, ch11:29 "a의 위수 r = min{r ≥ 1 : a^r ≡ 1 (mod N)}" (used through ch11:153, 163 and in s/r). (ii) ch18 uses r for the number of samples (ch18:17 "표본 수 r", table ch18:288). (iii) ch19's text never uses r for Trotter steps: the table header is "Trotter 단계" (ch19:244) and r appears only as the error rate (ch19:94, 251, 309). r as a step count exists only as a variable in `ch19_candidates.py`.

Replacement row:

> | r | Trotter 단계 수 | 12장 | 11장에서는 위수(a^r ≡ 1 mod N인 가장 작은 r), 14·19장에서는 게이트당 오류율(무작위 벤치마킹 값), 18장에서는 표본 수다. 19장은 Trotter 단계 수를 글자 없이 "단계"로 적는다 |

### T1-B3. ch00 §2, p row (and §4) — ch00_notation.md:33, :62
Quote: "14장에서는 탈분극 채널의 강도. r에서 p = r(d+1)/d로 환산한다"

The formula is right (ch14:33 "p = r·(d+1)/d (d = 2^k) — 2큐비트 1.25배, 1큐비트 1.5배"), but d is not defined anywhere in ch00, so the conversion cannot be used from this page. Missing clash: ch17 uses p for the order of the speed-up (ch17:68 "제곱 가속(p = 2)", table ch17:175, ch17:250), in the chapter whose verdicts §4 summarises; ch11:29 uses p, q for the prime factors.

Replacement row:

> | p | QAOA·HVA의 층 수 | 13·18·21장 | 14장에서는 탈분극 채널의 강도(비항등 파울리 오류의 총확률)다. k큐비트 게이트면 d = 2ᵏ이고 r에서 p = r(d+1)/d로 환산한다(2큐비트 1.25배, 1큐비트 1.5배). 17장에서는 가속의 차수(p = 2가 제곱 가속), 11장에서는 N = pq의 소인수다 |

### T1-B4. ch00 §2, γ row — ch00_notation.md:27
Quote: "| γ | 시험 상태와 바닥상태의 겹침 확률 … | 12·17장 | 7장에서는 진폭 감쇠 채널의 감쇠 확률, 2장에서는 전역 위상의 각 |"

All three listed uses are correct (ch12:162, 214; ch17:102, 191, 236; ch07:62; ch02:142–153). Missing: the QAOA/HVA cost angle, which is the use a reader of ch13 sees ten times (ch13:34 "e^{−iγₗH\_ZZ} e^{−iβₗH\_X}", ch13:179 "γ ∈ \[0,π\]") and again in ch14:71 "(γ, β) = (0.2064, 0.3927)"; and ch18:58–64 where γ(U) is the Shende–Bullock–Markov invariant matrix. Also the overlap γ is introduced in ch10:261 ("겹침이 γ = |c\_0|² … 1/√γ번"), not ch12.

Replacement row:

> | γ | 시험 상태와 바닥상태의 겹침 확률. 반복 수는 1/γ, 진폭 증폭을 쓰면 1/√γ | 10·12·17장 | 13·14장에서는 QAOA·HVA의 비용 층 각도(γ, β), 18장에서는 2큐비트 합성의 판정 행렬 γ(U), 7장에서는 진폭 감쇠 채널의 감쇠 확률, 2장에서는 전역 위상의 각 |

### T1-B5. ch00 §2, S row — ch00_notation.md:34
Quote: "| S | 반쪽 사슬의 얽힘 엔트로피(비트) | 15·17\~19장 | — |"

The "—" (no clash) is contradicted by the page itself: line 10 defines the gate "S = P(π/2)". Other uses: ch16:94 "S = SWAP으로 W = S(2Π\_A − I)S(2Π\_A − I)"; the reflections S₀ and S\_f in Q = −AS₀A†S\_f (ch12:50, ch16:86, ten occurrences in ch09); and in ch19 the same chapter uses S for the entropy (ch19:237) and for a matrix (ch19:272 "I + 0.2·S", ch19:305 "I + γS"). Entanglement entropy S is also first used in ch03:220.

Replacement row:

> | S | 반쪽 사슬의 얽힘 엔트로피(비트) | 3·15·17\~19장 | 1절의 위상 게이트 S = P(π/2)와 글자가 같다. 9·12·16장의 S₀, S\_f는 반사이고, 16장 Szegedy 워크의 S는 SWAP이다. 19·20장의 I + γS에서 S는 무작위 희소 대칭 행렬이다 |

### T1-B6. ch00 §2, χ row — ch00_notation.md:35
Quote: "| χ | MPS의 결합 차원. 버린 가중치가 10⁻⁶ 이하가 되는 최소 값 | 15·17·19장 | 2^S는 χ의 하한일 뿐이다. 실제 χ는 몇 배 크다 |"

The 10⁻⁶ definition is ch17's and ch19's (ch17:79, ch17:132, ch19:138). ch15 does not use it: there χ is a truncation level that is set by hand (ch15:57 "슈미트 수를 χ로 제한", table ch15:162–172 with χ = 2, 8, 32, 128), and ch15:85 writes "필요 χ ≈ 2^{S}", which is exactly what the fourth column of this row says is not true. Other meanings of the letter: stabilizer rank "안정자 계수 χ(|T⟩^{⊗6}) = 7" (ch15:56) and the marking function in amplitude estimation "표식 χ", "S\_χ" (ch16:37, 86, 96).

Replacement row:

> | χ | MPS의 결합 차원(자름마다 남기는 Schmidt 계수의 수) | 15·17·19장 | 17·19장의 "필요한 χ"는 버린 가중치가 10⁻⁶ 이하가 되는 최소 값이다. 2^S는 그 하한일 뿐이고 실제 값은 몇 배 크다(15장의 "χ ≈ 2^S"는 이 하한을 가리킨다). 15장의 안정자 계수 χ와 16장 진폭 추정의 표식 χ는 다른 양이다 |

(ch15:85 "필요 χ ≈ 2^{S}" would be better as "필요 χ ≳ 2^{S}", matching ch15:57 and ch15:243.)

### T1-B7. ch00 §2, K row — ch00_notation.md:36
Quote: "19장 화자 예에서 K\_spk는 화자 수"

K\_spk is a variable name in `ch19_candidates.py:62`; it does not occur in any chapter text (`grep K_spk ch*/ch*.md` finds only this line of ch00). The clash that does exist in the text is the Kraus operator: ch05:192–202 "ρ → Σₖ KₖρKₖ†, … K₀ = |0⟩⟨0|", ch07:62 "ρ → Σ K\_i ρ K\_i†".

Replacement row:

> | K | 커널 값. 두 특성 상태의 내적 절댓값의 제곱 | 6·19·22장 | 5·7장의 Kₖ는 Kraus 연산자다 |

### T1-B8. ch00 §4, fidelity bullet — ch00_notation.md:61
Quote: "**충실도는 두 가지다.** 상태 충실도는 … (5·12·13·15·19장), 평균 게이트 충실도는 … (7·14장). 비충실도는 1에서 충실도를 뺀 값이다."

There is a third one, and it is the one the word "비충실도" is used for most often before ch19: ch18:37–41 defines "전역 위상을 무시한 공정 비충실도(process infidelity)" f = 1 − |Tr(U†V)|²/d² and uses it in all synthesis tables (ch18:222–248, ten occurrences). It is neither a state fidelity nor the average gate fidelity (F\_avg = (d·F\_pro + 1)/(d + 1)). ch04:55 uses yet another gate distance, 1 − |Tr(A†B)|/d. State fidelity is also used in ch10:173, 271.

Replacement bullet:

> - **충실도는 세 가지다.** 상태 충실도는 두 상태 내적의 절댓값 제곱이다(5·10·12·13·15·19장). 평균 게이트 충실도는 게이트 하나를 모든 입력 상태에 대해 평균한 값이다(7·14장). 18장의 회로 합성은 전역 위상을 무시한 공정 충실도 |Tr(U†V)|²/d²를 쓴다(d = 2ⁿ). 둘은 F\_avg = (d·F\_공정 + 1)/(d + 1)로 이어진다. 비충실도는 1에서 충실도를 뺀 값이다. 4장의 게이트 거리 1 − |Tr(A†B)|/d는 제곱하지 않은 값이어서 작은 오차에서는 공정 비충실도의 절반쯤이다.

## C — minor

### T1-C1. Global-phase pointer — ch00_notation.md:11
Quote: "그럴 때는 본문이 따로 밝힌다(4장)". The statement is made in ch02, not ch04: ch02:58 "제어 게이트의 목표로 쓰이면(3장) 그 "전역" 위상이 제어 큐비트와의 상대 위상이 되어 물리가 된다" and the §7 paragraph ch02:201 "전역 위상은 언제 물리가 되는가". ch04 has only the controlled-Rz exercise (ch04:219–221). (ch03, which ch02:58 points to, contains no such passage either.)
Replacement: "…그럴 때는 본문이 따로 밝힌다(2장 3절·7절. 제어 회전의 분해는 4장 연습문제 3)."

### T1-C2. M row — ch00_notation.md:26
Quote: "12·16장에서는 진폭 추정의 Q 연산자 호출 수(2의 거듭제곱)". In both chapters M = 2^(register bits) and the number of Q calls is M − 1 (ch12:171 "Q 호출 2^t−1 (BHMT의 M = 2^t)", ch16:86 "Q 호출 M − 1 = 2ᵐ − 1번"). ch19:221 also uses M for the number of spatial orbitals.
Replacement fourth cell: "12·16장에서는 진폭 추정의 격자 크기 M = 2^(계수 큐비트 수)이고 Q 호출 수는 M − 1이다. 14장에서는 판독 혼동 행렬, 19장 5절 (3)에서는 공간 궤도 수"

### T1-C3. N row — ch00_notation.md:23
The row lists ch19 under the base meaning 2ⁿ, but ch19 also uses N for the number of residues (ch19:253–262, "비트열 수 4^(N−1)") and for the number of electrons (ch19:221 "전자 N개 … C(M, N/2)²").
Replacement fourth cell: "11장에서는 인수분해할 수. 19장에서는 잔기 수(5절 (5))와 전자 수(5절 (3))로도 쓴다"

### T1-C4. α row — ch00_notation.md:28
The ch12 remark is correct (ch12:48 "H' = H + wI(w ≥ ‖H‖, 예를 들어 Σ|c\_k|)"). Missing: α, β as the amplitudes of a qubit state (ch01, ch03:59, ch05:39) and |α⟩ as the superposition of non-solutions (ch09:32, ch12:50).
Suggested addition to the fourth cell: "1\~5장의 α, β는 상태의 진폭이고, 9·12장의 α는 해가 아닌 상태들의 중첩에 붙인 이름이다"

### T1-C5. Which chapters use the query and block-encoding units — ch00_notation.md:48, :49
(48) "오라클 호출 수(질의 수) | … | 8·9·15장": the chapters say "쿼리" (ch08 39 times, ch09, ch11:228, ch15); "질의" in ch18 means sample-and-query access (ch18:78), a different thing. ch11 and ch16 (ch16:30 "호출 비용은 모두 입력 오라클(U 또는 A) 호출 수다") also count in this unit.
(49) "블록 인코딩 호출 수 | … | 12·16·17·20장": ch12 never counts block-encoding calls; it mentions them once as an outlook (ch12:210 "(16장)").
Replacements: first cell of row 48 "오라클 호출 수(쿼리 수)", third cell "8·9·11·15·16장. 호출 하나의 게이트 수는 따로 센다. 18장의 "표본·질의 접근"은 고전 자료 구조의 이름이다"; third cell of row 49 "16·17·20장(12장은 7절에서 예고만 한다). α에 비례한다. 호출 하나에는 SELECT와 PREP의 게이트 수가 든다".

### T1-C6. T-count row: which ε, and R used twice — ch00_notation.md:47, :52
The Ross–Selinger statement is correct (see "checked"). But ε here is the synthesis error of one rotation, while ch00:30 defines ε as "목표 정밀도" of the whole computation. With R rotations the errors add up to Rε (ch12's QPE row has R = 28,679, so ε = 10⁻⁶ per rotation is 0.03 in total). Also R means "회전 수" in row 47 and "한 번의 비용" in row 52.
Suggested addition to row 47: "여기서 ε은 회전 하나의 합성 오차다. 회전이 R개면 오차가 Rε까지 쌓이므로, 회로 전체의 오차를 ε\_전체로 묶으려면 ε = ε\_전체/R을 넣는다". Row 52: "한 번의 비용 R(T-count 줄의 회전 수 R과 다른 양이다)".

### T1-C7. Shots row — ch00_notation.md:50
Quote: "진폭 추정을 쓰면 1/ε". The row's unit is the shot; with amplitude estimation it is the number of calls to the preparation circuit inside one deep circuit that scales as 1/ε (ch16:37 "오차 ε에 O(1/ε)" is a call count), not the number of shots.
Replacement: "기댓값을 오차 ε로 알려면 샷 수가 1/ε²에 비례한다. 진폭 추정을 쓰면 샷 수가 아니라 상태 준비 회로의 호출 수가 1/ε에 비례하고, 회로가 그만큼 깊어진다"

### T1-C8. Source grades — ch00_notation.md:14
Quote: "B는 도구·해설 글로, 직관을 돕는 데만 쓴다." Two chapters use grade B for something else: ch11:268 "등급 B(2차 정리; 원논문 확인 필요)" and ch19:389, where a journal paper is graded B because its number (1.3×10¹²) was verified only through another paper, and that number is used as evidence in ch19:231.
Replacement: "B는 도구·해설 글과, 원문을 직접 확인하지 못해 2차 서술로 확인한 출처다. 직관을 돕거나 보조 수치를 대는 데 쓰고, 그 경우 본문이 확인 경로를 밝힌다."

### T1-C9. "2부와 4부의 회로에는 세 숫자를 함께 적는다" — ch00_notation.md:41
ch14 (2부) records only CZ and SWAP counts after routing (ch14:73), and ch08:149 gives no numbers for Simon. Suggested: "2부와 4부의 회로에는 가능한 곳마다 세 숫자를 함께 적는다("게이트 비용 기록"). 14장은 라우팅 뒤의 CZ 수만 적는다."

---

# Part 2 — Target 2: ch19, second review

Rerun: `run/rerun.txt` equals `ch19_stdout.txt` line for line except timings (FFT 14.1 → 10.6 µs, HP N = 12 0.28 → 0.22 s, N = 30 extrapolation 124 → 98 days, FFT ratio 157 → 210, total 18 → 17 s). The md's timing wording ("10 µs 안팎(실행마다 9\~14 µs)", "150배 이상", "3\~4달", "약 20초", "0.3초 안에") covers both runs.

Status of the R7 findings: A1, A2, B2, B9, B10, C1–C14 are resolved and the new numbers match stdout. B1, B3, B5–B8, B11 were edited, but the edits leave or introduce the problems below. B4 is resolved in rows 2 and the text but left one stale cell (T2-B5).

## A — factually wrong

### T2-A1. ch19 §5(3), §10 — ch19_candidates.md:231, :389
Quote (231): "N₂/cc-pVDZ의 1.4×10¹²은 지금까지 보고된 가장 큰 FCI 계산(1.3×10¹² 행렬식, Gao et al. 2024)과 같은 크기이고" and (389) "— 가장 큰 FCI 계산(1.3×10¹² 행렬식)."

The Gao et al. citation itself is right (authors, title, J. Chem. Theory Comput. 20, 1185–1192 (2024), C₃H₈/STO-3G, 1.3×10¹² determinants, called "by far, the largest FCI calculation reported to date" in arXiv:2402.13111, checked). But "지금까지" is out of date for a page dated 2026-11: Shayit, Liao, Upadhyay, Hu, Zhang, DePrince, Yang, Li, "Breaking the Quadrillion Determinant Barrier in Numerically Exact Configuration Interaction", arXiv:2505.20375 (May 2025), published as "Numerically exact configuration interaction at quadrillion-determinant scale", Nature Communications 16, 11016 (2025), reports numerically exact CI beyond 10¹⁵ determinants, three orders of magnitude above Gao et al., which it cites as the previous state of the art (https://www.emergentmind.com/papers/2505.20375 ; https://escholarship.org/content/qt3tw2m38x/qt3tw2m38x.pdf). Caveat on what I could read: the arXiv abstract speaks of a relativistic FCI of a magnesium atom with over 10¹⁵ determinants as "the largest FCI calculation reported to date"; the summary of the published version names an X2C-CASCI of HBrTe (1.05×10¹⁵ determinants) as "the largest CI calculation ever reported". Either way 1.3×10¹² is no longer the largest. The chapter's conclusion (N₂/cc-pVDZ is within reach of exact diagonalisation) gets stronger, not weaker.

Replacement for the sentence in md:231:

> N₂/cc-pVDZ의 1.4×10¹²은 2024년에 보고된 FCI 계산(C₃H₈/STO-3G, 1.3×10¹² 행렬식, Gao et al. 2024)과 같은 크기다. 2025년에는 10¹⁵ 행렬식 규모의 수치적으로 정확한 CI 계산까지 보고됐다(Shayit et al. 2025). 1s 궤도를 얼리면 4.3×10⁹으로 줄어 FCI 기준값이 오래전부터 있는 분자이기도 하다.

Replacement for md:389:

> - Gao, Imamura, Kasagi, Yoshida, "Distributed implementation of full configuration interaction for one trillion determinants", *J. Chem. Theory Comput.* 20, 1185 (2024) — 1.3×10¹² 행렬식의 FCI 계산. 수치는 [arXiv:2402.13111](https://arxiv.org/abs/2402.13111)의 서술로 확인했다. 등급 B. Shayit, Liao, Upadhyay, Hu, Zhang, DePrince, Yang, Li, "Numerically exact configuration interaction at quadrillion-determinant scale", *Nature Communications* 16, 11016 (2025), [arXiv:2505.20375](https://arxiv.org/abs/2505.20375) — 10¹⁵ 행렬식 규모. 초록만 확인했다. 등급 B.

## B — misleading / incomplete / inconsistent

### T2-B1. ch19 §5(7), verdict row 7 — ch19_candidates.md:264, :294; ch19_candidates.py:380
Quote (264): "12사이트 모형을 맨 진공에서 급랭해 보면(스크립트 R6) 전기장 항이 없을 때(J = 0) 반쪽 엔트로피가 t ≤ 3에서 4.58비트까지 자라지만, J = 2, m = 1에서는 0.91비트, J = 4에서는 0.66비트를 넘지 않는다. … 판정은 "가능(얽힘이 빠르게 자라는 매개변수 영역에 한함)"으로 적는다."

The code is a sound gauge-eliminated Schwinger model and matches its comment: H = Σ(XX + YY) + m Σ(−1)ⁿZₙ + J Σₙ Lₙ² with Lₙ = Σ_{k≤n}(Z_k + (−1)^k)/2 on links 0…10, and |1010…⟩ is the bare vacuum (it minimises the mass term, −12, and has Lₙ = 0 on every link). My independent operator-level build (`schw.py`) reproduces 4.58 / 0.91 / 0.66 bits and χ = 64 / 32 / 21. Three problems with how the text uses it:

1. The only fast-growing case shown, J = 0, is a free-fermion model. The Majorana degree-4 terms the same paragraph relies on come only from J Σ Lₙ²; with J = 0 the Hamiltonian is the XX chain (plus a staggered field), which ch18's test (D) catches with a 2n×2n matrix. So the chapter has no data point that is both interacting and fast-growing, and the verdict cell "급랭 엔트로피 4.58비트(J = 0)" supports "가능" with a classically trivial case.
2. The J = 0 run also has m = 0 (script:380 `(0.0, 0.0), (2.0, 1.0), (4.0, 1.0)`), so two parameters change between the compared runs; the text mentions only J.
3. H, J and m are not defined in the md (J otherwise means the spin coupling in this book, ch00:12).

An interacting, fast-growing regime does exist, so the verdict can stand with the right evidence (`schw.py`, same grid t = 0, 0.25, …, 3): J = 0.5, m = 0 → max 4.70 bits (χ max 64); J = 1, m = 0 → 4.33; J = 2, m = 0 → 3.14; J = 0, m = 1 → 3.35; J = 0.5, m = 1 → 2.50; J = 1, m = 1 → 1.49.

Fix. Script line 380: `for Jg, mg in ((0.0, 0.0), (0.5, 0.0), (2.0, 1.0), (4.0, 1.0)):` (the new line should print 4.70 bits, χ 64). Replacement for the two sentences in md:264:

> 12사이트 모형 H = Σ(XX + YY) + m Σ(−1)ⁿZₙ + J Σ Lₙ²(m은 질량, J는 전기장 에너지의 계수, Lₙ은 Gauss 법칙으로 정해지는 n번째 결합의 전기장)을 맨 진공 |1010…⟩에서 급랭해 보았다(스크립트 R6). 전기장 항을 끈 J = 0, m = 0은 자유 페르미온이라 견줄 대상이 아니고(4.58비트까지 자라지만 18장의 자유 페르미온 검사에 걸린다), 상호작용이 있는 J = 0.5, m = 0에서는 반쪽 엔트로피가 t ≤ 3에서 4.70비트까지 자란다. 반면 J = 2, m = 1에서는 0.91비트, J = 4, m = 1에서는 0.66비트를 넘지 않는다.

Verdict row 7 (md:294), fourth cell: "Majorana 4차, 급랭 엔트로피 4.70비트(J = 0.5, m = 0) 대 0.91비트(J = 2, m = 1)".

### T2-B2. ch19 §5(2) — ch19_candidates.md:206
Quote: "그러나 n과 함께 반드시 나빠지는 것은 커널 값의 크기와 정렬이지 1-NN이 아니다."

The new control (script R2) measures the 8- and 11-speaker data at n = 10 only, so it shows one thing: 1-NN is above chance there (0.37 vs 0.12, 0.70 vs 0.09; both match stdout). It does not show that 1-NN stops degrading with n. It does degrade. `kern.py` (the chapter's own `zz_feature_states`, `loo_1nn`, generator and noise 0.085; three seeds each, n = 2, 3, 4, 6, 8, 10):
- 11 speakers (all ten components carry speaker information): quantum 1-NN 0.72/0.64/0.76 → 0.81/0.75/0.85 → 0.84/0.77/0.79 → 0.80/0.67/0.75 → 0.67/0.66/0.62 → 0.65/0.63/0.63, while RBF rises to 1.00 from n = 4–6 on.
- 8 speakers: peak 0.90–0.92 at n = 3–4, then 0.73–0.79 (n = 6), 0.56–0.67 (n = 8), 0.39–0.50 (n = 10).
- 4 speakers: as in the chapter (0.91–0.99 at n = 2–3, 0.32–0.39 at n = 10).
So with informative features 1-NN still peaks near n = 3–4 and falls by 0.15–0.5 by n = 10. What the control separates is "falls to chance" (specific to four speakers) from "falls" (general here). The sentence over-corrects R7's B1.

Replacement for that sentence:

> 그러나 찍기 수준까지 무너지는 것은 화자가 4명뿐인 이 데이터의 사정이다. 특성 10개가 모두 화자 정보를 담는 11명 대조에서도 1-NN은 n = 10에서 0.70으로 RBF의 1.00에 못 미치지만, 찍기(0.09)와는 거리가 멀다. 어느 데이터에서나 n과 함께 나빠지는 것은 커널 값의 크기와 정렬이다.

(Optional: let R2 print the 11-speaker 1-NN at n = 4 and n = 10 so that "n ≥ 4에서 내려간다" has a source in the script.)

### T2-B3. Pre-registration, ch23a row, (b) and (c) — ch19_candidates.md:309
Quote: "(b) 포화 전 구간(S가 반쪽 사슬 큐비트 수의 1/3에 닿기 전)에서 log χ 대 t 직선 맞춤과 log χ 대 log t 직선 맞춤의 잔차 비교 … n = 12, 16, 20, t ≤ 3 … (b) 지수 맞춤의 잔차 제곱합이 거듭제곱 맞춤의 절반 이하이고, 기울기가 세 크기에서 10% 안에서 일치 (c) 세 크기에서 따로 구한 외삽값이 서로, 그리고 이 장의 값(…)과 n·t 모두 2배 안에서 일치하고, 그 크기에서 rG가 5보다 크다 | (a)(b)(c) 통과."

I ran the registered test (`prereg23a.py`, `prereg23a_n20.txt`; exact state vectors, both chains, n = 12, 16, 20; my n = 16 numbers reproduce the chapter's: S rate 2.91, log₂χ ≈ 1.76 + 3.58t, J₂ rate 3.01).

(b) is still not well-posed, and as written it fails:
- The time grid is not registered and decides the outcome. On the chapter's grid (step 0.25) the window "S < 반쪽의 1/3" holds 2 points at n = 12 (t ≤ 0.5) and 3 points at n = 16 (t ≤ 0.75): at n = 12 both two-parameter fits are exact and the ratio is 0/0; at n = 16 the ratio is 0.12 (Heisenberg) and 0.16 (J₂). On a 0.05 grid the same n = 16 data gives 0.49 and 0.45, right at the threshold 0.5.
- On the 0.05 grid the ratio SSE(지수)/SSE(거듭제곱) is: Heisenberg 1.02 (n = 12), 0.49 (n = 16), 0.22 (n = 20); J₁–J₂ 0.87, 0.45, 0.28. So n = 12 fails for both chains. The reason is the window, not the physics: at n = 12 it ends at t ≈ 0.6 with χ ≤ 20 (the maximum possible χ is 64), too short to tell the two laws apart.
- Slopes of log₂χ against t in the window: Heisenberg 4.55, 4.23, 4.02; J₁–J₂ 6.08, 5.41, 4.99. "10% 안에서 일치" is not defined (max/min is 1.13 and 1.22; deviation from the mean is ±7% and +11%/−9%). The J₁–J₂ chain fails under either reading.
- The row does not say whether all three sizes and both chains must pass.
So the registered expectation "(b) 통과" is contradicted by a two-minute calculation with the chapter's own definitions.

(c) depends on a recipe the row does not name. With the chapter's recipe (S slope on t = 0.25–1.25, log₂χ fit on t = 0.25–1.5, saturation ratio from t ≥ 2) the extrapolated (n, t) are: Heisenberg (36, 3.85), (30, 3.20), (32, 3.18) with rG = 12.4, 8.7, 9.1; J₁–J₂ (40, 4.98), (24, 2.84), (24, 2.56) with rG = 17.7, 6.1, 5.5. That passes, narrowly for J₁–J₂ (t ratio 1.95, rG 5.5). If the χ fit is taken on the window of (b) instead, the J₁–J₂ chain at n = 12 gives (16, 1.9) and rG = 2.6, which fails. At n = 12 the chapter's fit window lies inside the saturated region (χ max 64), so that size should not carry the verdict.

Replacement cells for the ch23a row.
지표 (b)(c):

> (b) 포화 전 구간(S가 반쪽 사슬 큐비트 수의 1/3에 닿기 전, 시간 간격 0.05)에서 log χ 대 t 직선 맞춤과 log χ 대 log t 직선 맞춤의 잔차 비교 (c) χ = 10⁴에 닿는 (n, t)의 외삽값과 그 크기 회로의 rG(r = 0.003). 외삽은 이 장 5절 (4)와 같은 절차로 한다(S의 기울기는 t = 0.25\~1.25, log₂χ 맞춤은 t = 0.25\~1.5, 포화 비율은 t ≥ 2의 평균). 베이스라인은 정확한 상태 벡터와 χ를 제한한 MPS

크기와 예산:

> n = 16, 20(판정용), n = 12(참고용. 구간이 t ≤ 0.65, χ ≤ 20에서 끝나 두 맞춤이 갈리지 않는다), t ≤ 3. Heisenberg 사슬과 J₂ = 0.5 사슬

통과 기준:

> (a) 지수 4 ± 0.4 (b) 두 사슬 모두 n = 20에서 지수 맞춤의 잔차 제곱합이 거듭제곱 맞춤의 절반 이하이고, n = 16과 20의 기울기가 둘의 평균에서 ±10% 안 (c) n = 16과 20에서 따로 구한 외삽값이 서로, 그리고 이 장의 값(n ≈ 30\~42, t ≈ 3.2\~4.6)과 n·t 모두 2배 안에서 일치하고, 그 크기에서 rG가 5보다 크다

(With these cells my scratch run gives (b) ratio 0.22 and 0.28, slopes within ±2.5% and ±4%; (c) passes with rG 8.7–9.1 and 5.5–6.1. These are reviewer numbers; the author should decide whether seeing them before ch23a is acceptable or whether the row should say that the n = 16 and 20 state-vector data were looked at when the thresholds were set.)

### T2-B4. Pre-registration, ch20 row — ch19_candidates.md:305
Quote: "성공 확률을 포함한 총 호출 수는 κ²에 가깝다(17장) | … QSVT 다항식 차수와 성공 확률 … | 차수의 κ 지수가 1 ± 0.2이고, 성공 확률을 포함한 총 호출 수의 κ 지수가 2 ± 0.3."

"총 호출 수" is not defined, and the three natural readings give three different verdicts. `prereg20.py` (the chapter's I + γS family with γ tuned to κ = 2, 4, 8, 16, 32; degree from the standard 1/x construction at ε = 10⁻³; success probability p = ‖A⁻¹b‖²/(2κ)² for ‖A‖ = 1):
- degree exponent 1.10 (N = 64 and 256): passes "1 ± 0.2";
- degree/√p (amplitude amplification, which is what ch17:195 "진폭 증폭 포함 총 κ²·log(κ/ε)" means), random b: exponent 1.61 (N = 64), 1.73 (N = 256). The first fails "2 ± 0.3", the second passes by 0.03;
- degree/p (repeat until success), random b: 2.12 (N = 64), 2.35 (N = 256). The second fails;
- b equal to the softest eigenvector: p = 0.25 for every κ, exponent 1.10 in both conventions.
κ² is the worst case over b; for a random b on this family most of the weight sits on well-conditioned eigenvectors, so the typical exponent is below 2. The row must fix the convention, the right-hand side and N, and the expectation "통과" should say for which of them.

Replacement cells.
가설:

> 잘 조건화된 희소 행렬에서 QSVT 역행렬 다항식의 차수 d(진폭 증폭을 뺀 블록 인코딩 호출 수)는 블록 인코딩의 유효 조건수 κ에 거의 선형이다. 진폭 증폭을 넣은 총 호출 수 d/√p\_s는 최악의 b에서 κ²이고(17장), 무작위 b에서는 κ¹과 κ² 사이다

지표와 고전 베이스라인:

> 해의 오차 10⁻³까지의 QSVT 다항식 차수 d, 성공 확률 p\_s = ‖A⁻¹b‖²/(2κ)²(‖A‖ = 1로 맞춘 뒤, b는 무작위 단위 벡터 10개의 평균), 총 호출 수 d/√p\_s. 베이스라인은 CG의 nnz × 반복 수

통과 기준:

> N = 256에서 차수의 κ 지수가 1 ± 0.2이고, 총 호출 수 d/√p\_s의 κ 지수가 1.3 이상 2.3 이하. 차수가 N과 무관한지는 검산 항목으로만 보고한다(다항식이 κ와 ε만으로 정해지므로 가설이 아니다)

### T2-B5. Verdict table row 6 left behind by the B4 fix — ch19_candidates.md:293 vs :289
Quote (293): "| 6 | 약물 활성 분류(QML) | 17장 1단계, 18장 5′ | 2번과 같음 | 탈락 | 싣지 않음 |"

Row 2 was changed from "18장 4′·5′" to "18장 3′·4′·6′(5′는 부분적)" and the text now says 5′ "는 방향만 보였다" (md:219). Row 6 still names 5′ as the stage that catches it and says "2번과 같음". The two rows now contradict each other.

Replacement third cell of row 6: "17장 1단계, 18장 3′·4′·6′"

## C — minor

### T2-C1. ch19 §5(1), frame of 4,096 samples — ch19_candidates.md:193; ch19_candidates.py:389–402
Quote: "N = 4,096(12큐비트)으로 4배 늘려도 필요한 샷 수는 그대로 32이고(스크립트 R1b: 16샷 0.987, 32샷 0.998)". Two sentences earlier the same paragraph says the threshold is 24 when sampled finely; "그대로 32" holds only on the powers-of-two grid that R1b uses (16, 32, 64). Finer runs (`pitch.py`): N = 1,024 crosses 99% at 23 shots (0.9906, 100,000 trials; 22 shots 0.9897), N = 4,096 at 20 shots (0.9906, 20,000 trials). Script details: the comment at py:399 says "±2 구간 = 같은 ±7.8 Hz 안", but the N = 1,024 run used ±1 bin = ±15.6 Hz (the result is the same with ±4 bins, checked); R1b is printed under the "(E) 연습문제 검산" banner although md:185 says the (E) block is exercise checks, and py:317–319 is a leftover duplicate "(E)" comment banner above the (R) block.
Replacement: "N = 4,096(12큐비트)으로 4배 늘려도 2의 거듭제곱으로 잰 필요 샷 수는 그대로 32이고(스크립트 R1b: 16샷 0.987, 32샷 0.998)". Script comment py:399: "구간 폭이 1/4 이 됐으므로 ±2 구간(±7.8 Hz)으로 잡았다. ±4 구간(±15.6 Hz)으로 잡아도 결과가 같다".

### T2-C2. What Buyens et al. is cited for — ch19_candidates.md:264
Quote: "Schwinger 모형은 가둠(confinement) 때문에 매개변수에 따라 얽힘이 느리게 자라고, 그런 영역은 MPS로 긴 시간까지 계산된다(Buyens et al. 2014)." The bibliographic data are right (Buyens, Haegeman, Van Acoleyen, Verschelde, Verstraete, PRL 113, 091601, arXiv:1312.6654). The paper simulates the real-time response to a background-field quench with gauge-invariant MPS up to about t = 25 with the bond dimension growing from 18 to 173, and says "It is this growth of the required bond dimensions, which can be traced back to the growth of entanglement, that makes the computations more costly at later times" (text of https://core.ac.uk/download/55887014.pdf). It does not attribute slow growth to confinement; that link is Kormos et al. 2017, which ch17 cites.
Replacement: "Schwinger 모형은 가둠(confinement)이 있는 이론이고, 가둠은 급랭 뒤 얽힘의 증가를 늦출 수 있다(17장의 Kormos et al. 2017). 실제로 배경 전기장을 켠 급랭은 MPS로 t ≈ 25까지 결합 차원 200 미만에서 계산됐다(Buyens et al. 2014)."

### T2-C3. Pre-registration, ch22 row: calibration target has no n — ch19_candidates.md:307
Quote: "임베딩 잡음은 RBF 커널의 NMI가 0.8\~0.9가 되도록 먼저 맞추고 고정한다". The RBF kernel is built on the first n principal components and its quality depends strongly on n: with 8 speakers and the chapter's noise, RBF 1-NN is 0.68–0.92 at n = 2 and 1.00 from n = 6 on (`kern.py`). One noise level cannot put it at 0.8–0.9 for every n in {2, 4, 6, 8, 10}.
Replacement: "임베딩 잡음은 n = 8에서 RBF 커널의 NMI가 0.8\~0.9가 되도록 먼저 맞추고 모든 n에 같은 값을 쓴다"

---

## Summary

ch00 has the gate definitions, the bit order, every conversion count, the Ross–Selinger rule, the depolarizing conversion and the verdict labels right. Its weak part is the symbol table's chapter claims. Three cells state something the chapters do not do (φ(x) in ch06, G and rG in ch14, "only ch07 and ch14 use ns/µs"), and for t, r, p, γ, S, χ, K the clash column leaves out the clash a reader meets first (order r in ch11, QAOA γ in ch13, the gate S, the Kraus Kₖ, speed-up order p in ch17) or lists one that is only a script variable (K\_spk). The fidelity bullet misses ch18's process fidelity.

In ch19 every changed number matches stdout and the new script block computes what the text says. One claim is out of date (largest FCI). Three of the fixes went too far or not far enough: the Schwinger "fast" case is the free-fermion point J = 0, the kernel paragraph now says 1-NN does not have to fall with n when the control (and a sweep) shows it does, and row 6 of the verdict table kept the old 5′ label. Two rewritten pre-registration rows are still not decidable as written: ch23a (b) depends on an unregistered time grid and fails at n = 12 for window reasons, and ch20's "총 호출 수" gives exponents 1.6, 2.1 or 1.1 depending on conventions the row does not fix.

## Checked and found correct

Target 1
- Bit order: big-endian for N&C, Cirq and the book's numpy, little-endian for Qiskit and Quirk, and the pointer to ch03 (ch03:22, 35–37, 200, 228).
- Rz(θ) = exp(−iθZ/2), P(λ) = diag(1, e^{iλ}), Z = P(π), S = P(π/2), T = P(π/4), "P와 Rz는 전역 위상만 다르다(2장)" (ch02:53, 58, 151).
- Units: entropy in bits (ch03:220, ch15:160, ch17:11, ch19:237); mHa and 1.6 mHa (ch12:48, ch13, ch17, ch19:231). Code labels "numpy 구현(실행 완료)" (19 hits), "작성됨, 미실행" in headings and on the first code line (every chapter), `import pennylane as qp` (no other alias found).
- Symbol rows found accurate as written: n; N for ch09/ch11; t for ch10/ch11 with m (ch16:86) and b (ch17:57); M for ch09 and ch14 (ch14:33); γ for ch12, ch17, ch07 (γ = 1 − e^{−t/T1}, ch07:62), ch02; α with w in ch12 (ch12:48); L for ch16/ch17 and ch13 (ch13:34); ε for ch16 (block-encoding error ch16:40, marked fraction ch16:94) and ch14 (ε₀₁, ε₁₀, ch14:33); κ for ch17/ch19 and the "steepness" use at ch16:201; p for ch13/ch18 (ch18:262) and p = r(d+1)/d (ch14:33, 161, 193); r as RB error rate in ch14.
- Cost conversions against ch04, ch08, ch09, ch18 and the cost paragraphs of ch10–12: CP(θ) = 2 CX (CP(π) = CZ counts 1), exp(−iθZZ) = 2, SWAP = 3, generic two-qubit gate = 3 (ch18:13, 58–64, 233), Toffoli = 6 CX and 7 T (ch04:41, 149, 158; ch18:248), 8-CX linear Toffoli and routed CZ counts in ch14 (ch14:43, 73).
- Depth convention: recomputed with an ASAP layer counter in which single-qubit gates take a layer and measurements do not. BV (n = 6, s = 101101) gives exactly 7 (ch08:149 "≈ 7"); one n = 3 Grover iteration gives 26 (13-gate CCZ) or 28 (H·CCX·H), against ch09:90 "≈ 30"; Toffoli alone 11.
- Ross & Selinger: abstract of arXiv:1403.2975 (read via https://alphaxiv.org/abs/1403.2975): "In the typical case, this yields circuit approximations of T-count 3log_2(1/epsilon) + O(log(log(1/epsilon)))"; 3·log₂(10⁶) = 59.8, so "약 60" is right; Quantum Inf. Comput. 16, 901–953 (2016) (rintonpress listing).
- 1/ε² shots for an expectation value and 1/ε calls with amplitude estimation; TTS₉₉ = R·ln(0.01)/ln(1 − P) with training added once (ch19:79–82, exercise 3); (1 − r)^G ≈ e^{−rG} (ch19:97); "19장 연습문제 3" pointer; MPS cost ∝ χ³ (ch15:57, ch19:339).
- §4: the three verdicts and "이 크기에서는 이득 없음" (ch17:104, 197; ch19:43), the three roles (ch19:33–37), state fidelity in ch05/12/13/15/19 and average gate fidelity in ch07/14 (ch07:164, ch14:33), r versus p (ch14:33), "예측"/"실측" (ch14:110, 195).

Target 2
- Fresh run equals `ch19_stdout.txt` except timings. Every number that is new or changed in the md against stdout: R1 (0.988/0.991/0.995, 24,528), R1b (0.987/0.998, 131,008, 49,152), R2 (0.984…0.001; 0.0020/0.37/0.12; 0.0022/0.70/0.09; RBF 1.00), R3 (4.33×10⁹), R4 (4.2–6.8; 1.76 + 3.58t; t ≈ 3.2; 9.4 bits; n ≳ 30; 2,916; 1.6×10⁻⁴; 17.2; 1.7; 0.18; J₂ 3.01), R5 (587 … 1,207,215), R6, E1 (37; 51,398). Arithmetic recomputed: 5,721, 3.4×10⁻⁸, 2^5.03 ≈ 33, 2.4–3.2×, 2.2 ms, 74.6 h, 4,928 for the rounded inputs.
- The new script block computes what the text says: R2 between-class variance is var(class means)/var(component) per PCA component; the 8- and 11-speaker controls reuse the generator, noise 0.085, PCA to 10, [0, π] scaling, the same kernel functions, and chance = (per − 1)/(K·per − 1); R4 fits log₂χ on t = 0.25–1.5 and S on 0.25–1.25, n = 2⌈S/0.646⌉, CNOT = 3((n/2)(r + 1) + (n/2 − 1)r); the J₁–J₂ Hamiltonian adds 0.5 Σ σᵢ·σᵢ₊₂ to the n = 16 chain; R5 uses n_opt/4^(N−1); R6 as described in T2-B1.
- Integrability: level-spacing ratio in sectors of fixed total S and reflection parity, n = 14 open chain (`levelstat.py`): J₂ = 0 gives ⟨r⟩ = 0.40 (Poisson 0.386), J₂ = 0.5 gives 0.545 (GOE 0.531). So "Bethe 가설풀이로 풀리는 적분 가능 모형" and "적분 가능성을 깬 J₁–J₂ 사슬" both hold. (J₂ = 0.5 is the Majumdar–Ghosh point, where only the ground state is known exactly; the spectrum is not integrable.)
- Citations: Wouters, De Nardis, Brockmann, Fioretto, Rigol, Caux, arXiv:1405.0172 (title and authors confirmed; the paper treats Δ ≥ 1 and gives the Δ = 1 details in a companion paper, so it covers the isotropic chain); Alba & Calabrese, arXiv:1608.00614, DOI 10.1073/pnas.1703516114 (XXZ chain, quasiparticle prediction checked against tDMRG/iTEBD); Buyens et al. authors/title/abstract; Gao, Imamura, Kasagi, Yoshida, JCTC 20, 1185–1192 (2024), C₃H₈/STO-3G, 1.3×10¹², as quoted in arXiv:2402.13111.
- Exercise 1: N&C Exercise 4.46 is the density-matrix parameter count (4ⁿ − 1); solution 2ⁿ + 2ⁿ(2ⁿ − 1) − 1 = 4ⁿ − 1; pure state 2·2ⁿ − 2; 2¹⁴/441 = 37.2; 2⁵⁶/1.402×10¹² = 51,398; no stale reference to the old exercise 1 remains (the 4,096-sample material is now in §5(1) and cites R1b).
- "18장 3′·4′·6′(5′는 부분적)" against ch18's track table (ch18:137–146): 3′ covers kernel concentration, 4′ the shot count, 5′ the classical surrogate, 6′ the same-budget comparison. §1, §3 table, §5(2) text, verdict row 2 and exercise 2(d) agree with each other.
- Consistency of the A1/B2/B10 edits across §1 item 4, §5(4), verdict row 4, pre-registration (c), §7 and §9: all use n ≈ 30–42, t ≈ 3.2–4.6, 2,900–5,700 CNOT, 10⁻⁴–10⁻⁸, "오류율 10배 개선 또는 오류 정정". ch17:192 does condition its "가능" on non-integrable models.
- Pre-registration rows ch21 (random-sampling TTS 3,141 / 3,208 / 1.2×10⁶ recomputed; pass size, TTS = ∞ rule and unit are stated), ch22 apart from T2-C3, ch23, and ch23a (a) (chapter data give exponent 4.02 on 10–40 steps) are checkable as written.
- "후보 10개에 드는 계산이 합쳐 10초였다" (§9): the blocks before (R) take 10.0 s here; the whole script 17–18 s.

## Could not check

- N&C Exercise 4.46 wording from the book itself (no copy here); confirmed from memory and from a public solution page whose heading and content are the n-qubit density-matrix count.
- Journal volume and page of Wouters et al. (PRL 113, 117202) and Alba & Calabrese (PNAS 114, 7947): arXiv listing pages did not show a journal reference line; the DOI prefix of the second confirms PNAS; numbers accepted from memory. arXiv abstract pages for 1403.2975 and 1312.6654 were not reachable directly and were read through mirrors.
- The exact system of the 10¹⁵-determinant calculation in T2-A1 (two summaries differ, see there); the replacement text avoids naming it.
- ch20–ch22 references in ch00 (forward pointers) and the cost paragraphs and section 6 of ch10–ch13 (being edited) were not reviewed. ch00's html build was not compared with the md.
- The exponents in T2-B4 use an ideal normalisation ‖A‖ = 1; a real LCU block encoding has a larger α, which changes p and the effective κ but not the point that the convention must be fixed.
- Whether UCCSD stays within 1.6 mHa for stretched LiH, and the sampling overhead of error mitigation at rG ≈ 17, as in R7.
