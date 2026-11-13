"""부록 A 5절 — BV '0000' 확정 실험의 설계·실행 전 예측·분석 (numpy 만 쓴다. QPU·네트워크 불필요)

  python appA/appA_bv_confirm.py                     설계 검산 + 실행 전 예측 + 판정 규칙 → appA/appA_bv_plan.json (제출 스크립트가 읽는다)
  python appA/appA_bv_confirm.py runs/<기록>.json     appA_bv_submit.py 가 남긴 기록을 읽어 같은 판정 규칙을 적용한다

무엇을 묻는가. 2026-11-11 ibm_fez 의 BV(s = 1011)에서 오답 '0000' 이 5.15% 나왔다. '첫 CZ 때 보조 큐비트(물리 3번)가 |1⟩ 이었다'로
설명되지만(부록 A 5절 ①②) 왜 |1⟩ 이었는지는 가르지 못했다. 이 스크립트는 그날 실제로 돈 ISA 게이트 목록(runs/ 의 고정 사본)에서
변형 회로를 만들고, 경우마다 '0000' 이 어떻게 나와야 하는지를 **실행 전에** 적어 둔다. 판정 규칙도 여기서 고정한다.

설계가 실행 순서에 기대지 않게 한다(독립 검토 R14 의 A1). IBM 문서("Repetition rate execution")는 한 작업의 회로들을 샷마다 한 번씩
돌아가며 실행한다고 적고(샷 0 에서 회로 0…k, 샷 1 에서 다시 0…k), PUB 의 순서와 한 묶음 실행은 보장하지 않는다고도 적는다.
'직전에 돈 회로'가 무엇인지 모르므로 첫 작업에는 회로 둘만 넣는다. 둘의 '0000' 만으로 네 경우가 갈리고 실행 순서도 함께 드러난다.

  첫 작업   P0  그날의 BV 그대로            P1  보조가 |0⟩ 으로 끝나게(마지막 H → X)
  경우      S  직전 회로의 들뜸이 남는다 + 같은 회로가 연달아 돈다      X  직전 회로의 들뜸이 남는다 + 두 회로가 번갈아 돈다
            I  직전 회로와 무관하다(시작부터 들떠 있거나 회로 안에서 들뜬다)   T  지금은 없다
  둘째 작업의 후보(첫 작업의 결과를 본 뒤에 정한다 — 여기서는 회로만 만들어 검산한다)
            P2 맨 앞에 보조 reset   P3 다른 허브 큐비트로 옮긴 BV   P4 끝에서 보조도 측정   PM 게이트 없는 큐비트(114)에 표지 X   Qm 첫 CZ 앞에서 보조를 측정

모형·보정값·맞춘 e 는 appA_qpu_runs.py 의 것을 그대로 쓴다(그 스크립트를 조용히 불러온다 — 약 6초).
"""
import contextlib, io, json, pathlib, sys
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
with contextlib.redirect_stdout(io.StringIO()):
    import appA_qpu_runs as A                       # Replay, cal(당일 보정값), gb(BV 의 ISA 게이트 목록), anc, e_fit, meas, isa, bell, real

SHOTS = 4000
ANC = A.anc[0]; GB = A.gb; CAL = A.cal
P_OBS = A.meas["bv1011"]["p"].get("0000", 0.0); N_OBS = A.meas["bv1011"]["n"]
K_RESET = [np.array([[1, 0], [0, 0]], complex), np.array([[0, 1], [0, 0]], complex)]      # 완벽한 reset: |0⟩⟨0| + |0⟩⟨1|


class R2(A.Replay):
    """A.Replay + reset(완벽한 것으로 친다. 걸리는 시간과 reset 자체의 오류는 넣지 않았다)."""
    def gate(self, name, phys, params):
        if name == "reset":
            self.kraus(K_RESET, phys); return
        super().gate(name, phys, params)


def run(gates, e=0.0, exc=None, **kw):
    """ISA 게이트 목록 → (끝의 측정 분포(Qiskit 문자열), 끝에서 보조가 |1⟩ 일 확률). e: 보조가 |1⟩ 에서 시작할 확률."""
    qs = {q for g in gates for q in g[1]}; meas = [(g[1][0], g[3][0]) for g in gates if g[0] == "measure"]
    init = exc if exc is not None else ({ANC: e} if e else None)
    d = R2(qs, CAL, init_exc=init, **kw).run(gates).dist(meas)
    end1 = R2(qs, CAL, init_exc=init, readout=False, **kw).run(gates).dist([(ANC, 0)]).get("1", 0.0)
    return d, end1


def p_data(dist, s):
    """데이터 비트 4개(문자열의 오른쪽 4자리 = c3c2c1c0)의 주변 확률."""
    return float(sum(v for k, v in dist.items() if k[-4:] == s))


def relabel(gates, mp):
    return [[g[0], [mp.get(q, q) for q in g[1]], g[2], g[3]] for g in gates]


# ---------------- 1) 회로 ----------------
ops = [g for g in GB if g[0] != "measure"]; ms = [g for g in GB if g[0] == "measure"]
first_cz = min(i for i, g in enumerate(ops) if g[0] == "cz" and ANC in g[1])
last_cz = max(i for i, g in enumerate(ops) if g[0] == "cz" and ANC in g[1])
tail = [i for i, g in enumerate(ops) if i > last_cz and g[1] == [ANC]]
assert [ops[i][0] for i in tail] == ["rz", "sx", "rz"], "보조의 마지막 블록이 rz·sx·rz(= H) 라고 전제한다"
assert not any(ANC in g[1] for g in ops[:first_cz]), "첫 CZ 앞에는 보조에 걸리는 게이트가 없다고 전제한다"
DATA = [g[1][0] for g in ops if g[0] == "cz"]; IDLE = list(A.idle)
P0 = [list(g) for g in GB]
P1 = [g for i, g in enumerate(ops) if i not in tail] + [["x", [ANC], [], []]] + ms            # 보조: |1⟩ → X → |0⟩ 으로 끝난다
JOB1 = {"P0": P0, "P1": P1}; ORDER = ["P0", "P1"]
P2 = [["reset", [ANC], [], []]] + P0
P4 = P0 + [["measure", [ANC], [], [4]]]
PM = ops + [["x", [IDLE[0]], [], []]] + ms                                                    # 게이트가 걸렸다는 표지: 114번의 비트가 1 이 된다
QM = ops[:first_cz] + [["measure", [ANC], [], [4]]] + ops[first_cz:] + ms                     # 첫 CZ 앞에서 보조를 읽는다(회로 중간 측정)
HUB_MAP = {ANC: 1003, **{q: 1000 + i for i, q in enumerate(DATA)}}                            # P3 검산용의 임의 번호(실제 큐비트는 장치를 보고 고른다)
CAND = {"P2": P2, "P4": P4, "PM": PM, "Qm": QM}

print("=" * 78 + "\n1) 회로 — 2026-11-11 의 ISA 게이트 목록(runs/20261111-063819Z_ibm_fez_ch14_isa_calib.json)에서 만든다\n" + "=" * 78)
print(f"    보조 큐비트 {ANC}, 보조와 CZ 하는 데이터 큐비트(순서대로) {DATA}, 게이트 없이 측정만 하는 큐비트 {IDLE}")
ideal = {}
for k, g in {**JOB1, "P2": P2, "P4": P4, "PM": PM}.items():
    d, end1 = run(g, noisy=False)
    # 동점(P4 는 보조 비트가 정확히 반반)이면 max 가 계산기마다 다른 쪽을 고르므로(2026-11-13: 기준 환경 11011, 맥·리눅스 01011) 1e-12 안의 출력을 사전순으로 모두 적고 첫째를 쓴다
    tops = sorted(s_ for s_ in d if d[s_] >= max(d.values()) - 1e-12); top = tops[0]; ideal[k] = top
    print(f"    {k:3s} cz {sum(1 for x in g if x[0] == 'cz')}, sx·x {sum(1 for x in g if x[0] in ('sx', 'x'))}, reset {sum(1 for x in g if x[0] == 'reset')}, "
          f"측정 {sum(1 for x in g if x[0] == 'measure')} | 잡음 없는 출력 {'·'.join(tops)} ({'각 ' if len(tops) > 1 else ''}확률 {d[top]:.4f}) | 끝에서 보조가 |1⟩ 일 확률 {end1:.3f}")
GOODS = ideal["P0"]
assert abs(run(P0, noisy=False)[0][GOODS] - 1) < 1e-9 and abs(run(P1, noisy=False)[0][GOODS] - 1) < 1e-9
assert abs(run(P1, noisy=False)[1]) < 1e-12 and abs(run(P0, noisy=False)[1] - 0.5) < 1e-12
assert abs(p_data(run(P0, e=1.0, noisy=False)[0], "0000") - 1) < 1e-9 and abs(p_data(run(P1, e=1.0, noisy=False)[0], "0000") - 1) < 1e-9, "보조가 |1⟩ 에서 시작하면 0000 (P0·P1 모두)"
assert abs(run(P1, e=1.0, noisy=False)[1] - 1) < 1e-9, "P1 에서 |1⟩ 로 시작한 보조는 |1⟩ 로 끝난다(들뜸이 다음 회로로 그대로 넘어간다)"
assert abs(p_data(run(P2, e=1.0, noisy=False)[0], GOODS) - 1) < 1e-9, "reset 이 완벽하면 P2 는 정답"
assert ideal["PM"] == "1111" and ideal["P4"][-4:] == GOODS
assert [g for g in QM if not (g[0] == "measure" and g[3] == [4])] == P0, "Qm 은 P0 에 보조의 중간 측정 하나만 더한 것"
d3, _ = (lambda g: (R2({q for x in g for q in x[1]}, CAL, noisy=False).run(g).dist([(x[1][0], x[3][0]) for x in g if x[0] == "measure"]), None))(relabel(P0, HUB_MAP))
assert abs(d3[GOODS] - 1) < 1e-9, "큐비트 번호만 바꾼 P3 도 같은 출력을 낸다"
print(f"    데이터 비트의 정답(Qiskit 문자열 c3c2c1c0) = {GOODS}. 보조가 |1⟩ 에서 시작하면 P0·P1 모두 0000 이고, P1 에서는 보조가 |1⟩ 로 끝난다.")
print("    둘째 작업의 후보 P2·P4·PM·Qm 과 번호만 바꾼 P3 도 잡음 없는 계산으로 검산했다(Qm 은 P0 에 중간 측정 하나를 더한 것인지만 확인).")

# ---------------- 2) 2026-11-11 의 자료를 실행 순서별로 다시 읽는다 ----------------
print("\n" + "=" * 78 + "\n2) 2026-11-11 의 자료 — 실행 순서에 따라 읽는 법이 달라진다\n" + "=" * 78)
pubs = A.isa[A.SJ]["pubs"]; uses = [i for i, g in enumerate(pubs) if any(ANC in x[1] for x in g)]
REP_DELAY = 250e-6                                  # 14장 5절의 가정(기록에는 없다). 샷 하나 ≈ 255 µs
T1_ANC = CAL["qubits"][str(ANC)]["T1_us"] * 1e-6
gap = len(pubs) * 255e-6
print(f"    그날의 Sampler 작업은 PUB {len(pubs)}개였고 보조 큐비트 {ANC} 을 쓰는 것은 PUB {uses}(BV) 뿐이다.")
print(f"    (가) 같은 회로의 샷이 연달아 돌았다면: BV 의 직전 샷은 BV 다. 보조는 절반이 |1⟩ 로 끝나므로 e = {A.e_fit:.4f} 는 '끝에서 |1⟩ 이던 것의 {A.e_fit / 0.5:.3f} 이 남는다'는 뜻이 된다.")
print(f"    (나) 회로를 번갈아 돌았다면(IBM 문서의 설명): BV 의 두 샷 사이는 약 {gap * 1e3:.2f} ms (반복 지연 250 µs 가정 × {len(pubs)}), T1 = {T1_ANC * 1e6:.0f} µs 의 {gap / T1_ANC:.1f}배다. "
      f"reset 이 전혀 없어도 남는 들뜸은 0.5 × exp(−{gap / T1_ANC:.1f}) = {0.5 * np.exp(-gap / T1_ANC):.4f} 이다 — e 의 {0.5 * np.exp(-gap / T1_ANC) / A.e_fit:.0%}. BV 자신의 잔류로는 5% 를 설명할 수 없다.")
base, _ = run(P0)
rows = []
for q in DATA:
    qs = {x for g in P0 for x in g[1]}; meas = [(g[1][0], g[3][0]) for g in P0 if g[0] == "measure"]
    dq = R2(qs, CAL, noisy=False, init_exc={q: 1.0}).run(P0).dist(meas)
    s = max(dq, key=dq.get); assert dq[s] > 0.999
    obs = A.meas["bv1011"]["p"].get(s, 0.0); rows.append((q, s, obs, base.get(s, 0.0)))
    print(f"    데이터 큐비트 {q} 가 |1⟩ 에서 시작하면 출력 {s}: 실측 빈도 {obs:.4f} (그 큐비트의 초기 들뜸의 상한), 당일 보정 재현 {base.get(s, 0.0):.4f}, 재현보다 많은 몫 {100 * (obs - base.get(s, 0.0)):+.2f}%p")
print(f"    세 데이터 큐비트는 매 샷 |1⟩ 로 읽히고 끝난다(s = 1011). (가)라면 이들의 잔류는 {100 * max(r[2] for r in rows):.1f}% 이하인데 보조는 {A.e_fit / 0.5:.1%} 다 → 3번 큐비트만의 문제이거나 측정 없이 끝난 큐비트만의 문제다.")
print("    (나)라면 데이터 큐비트도 약 1 ms 를 쉰 뒤라 잔류가 작은 것이 당연하고, 이 비교는 아무것도 가르지 못한다. 실제 순서는 기록에 없다.")

# ---------------- 3) 첫 작업 [P0, P1] 의 실행 전 예측 ----------------
print("\n" + "=" * 78 + "\n3) 첫 작업 [P0, P1] 의 실행 전 예측 — 보정값은 2026-11-11 의 것이다(실행하는 날의 값은 다를 것이다)\n" + "=" * 78)
_, end_e = run(P0, e=A.e_fit)
F = A.e_fit / end_e                                 # 끝에서 |1⟩ 이던 것 가운데 다음 회로의 시작까지 남는 비율(S·X 의 매개변수)
e_s1 = A.e_fit
for _ in range(200): e_s1 = F * run(P1, e=e_s1)[1]  # S: P1 의 직전 회로는 P1
e_x0, e_x1 = A.e_fit, A.e_fit
for _ in range(200): e_x0, e_x1 = F * run(P1, e=e_x1)[1], F * run(P0, e=e_x0)[1]    # X: P0 의 직전은 P1, P1 의 직전은 P0
E = {"S": (A.e_fit, e_s1), "X": (e_x0, e_x1), "I": (A.e_fit, A.e_fit), "T": (0.0, 0.0)}
NAME = {"S": "잔류 + 같은 회로가 연달아", "X": "잔류 + 두 회로가 번갈아", "I": "직전 회로와 무관", "T": "지금은 없다"}
PRED = {}
print(f"    잔류 비율 f = e / (P0 에서 보조가 끝에 |1⟩ 일 확률 {end_e:.3f}) = {F:.4f}. P1 은 보조를 |0⟩ 으로 끝내지만, |1⟩ 로 시작한 샷에서는 |1⟩ 로 끝난다(확률 {run(P1, e=1.0)[1]:.4f}).")
print("    경우                         보조의 초기 들뜸 e (P0, P1)    P(0000): P0      P1        성공 확률(게이트 중 열완화 넣은 모형 | 뺀 모형): P0              P1")
for h, (e0_, e1_) in E.items():
    r = {}
    for k, g, e in (("P0", P0, e0_), ("P1", P1, e1_)):
        d_in, _ = run(g, e=e); d_out, _ = run(g, e=e, gate_relax=False)
        r[k] = {"e": float(e), "p0000": p_data(d_in, "0000"), "success": p_data(d_in, GOODS), "success_no_gate_relax": p_data(d_out, GOODS)}
    PRED[h] = r
    print(f"    {h} {NAME[h]:22s}   ({e0_:.4f}, {e1_:.4f})              {100 * r['P0']['p0000']:5.2f}%    {100 * r['P1']['p0000']:5.2f}%     "
          f"{r['P0']['success']:.4f} | {r['P0']['success_no_gate_relax']:.4f}      {r['P1']['success']:.4f} | {r['P1']['success_no_gate_relax']:.4f}")
print(f"    (참고) 2026-11-11 실측 P(0000) = {P_OBS:.4f} ({N_OBS}샷, 표준오차 {np.sqrt(P_OBS * (1 - P_OBS) / N_OBS):.4f}), 성공 {A.meas['bv1011']['s']:.4f}")
print("    네 경우의 (P0, P1) 이 서로 다르다. S 면 실행은 연속이고 원인은 직전 회로의 들뜸이다. X 면 원인은 같고 실행이 교대다. I 면 실행 순서는 이 작업으로 알 수 없다.")
print("    성공 확률은 두 모형의 예측을 함께 적었다(부록 A 4.1절의 약속). 0000 의 판정에는 쓰지 않는다.")

# ---------------- 4) 판정 규칙 (실행 전에 고정) ----------------
RULES = {"shots": SHOTS, "gone_max": 0.003, "present_min": 0.01, "ratio_low": 1 / 3, "ratio_same": [0.6, 1 / 0.6], "ci": 0.95}


def ratio_ci(c_num, c_den, n):
    """p_num / p_den 과 그 95% 구간(로그 비의 델타 방법, 0 을 피하려 0.5 를 더한다)."""
    a, b = c_num + 0.5, c_den + 0.5
    r = a / b; se = np.sqrt(1 / a - 1 / (n + 0.5) + 1 / b - 1 / (n + 0.5))
    return float(r), float(r * np.exp(-1.96 * se)), float(r * np.exp(1.96 * se))


def verdict(c0, c1, n):
    """c0, c1: P0·P1 에서 데이터 비트가 0000 인 샷 수. 돌려주는 값: (경우 기호, 설명, 세부)."""
    p0, p1 = c0 / n, c1 / n; r, lo, hi = ratio_ci(c1, c0, n)
    det = {"p0": p0, "p1": p1, "ratio_p1_over_p0": r, "ci": [lo, hi]}
    if max(p0, p1) <= RULES["gone_max"]: return "T", "지금은 없다(그날의 일시적 상태였다). 가를 수 없다", det
    if max(p0, p1) < RULES["present_min"]: return "약함", "0000 이 있지만 약해졌다(1% 미만). 판정하지 않고 비율만 보고한다", det
    if p0 >= RULES["present_min"] and hi < RULES["ratio_low"]: return "S", "직전 회로의 들뜸이 남는다 + 같은 회로가 연달아 돈다", det
    if p1 >= RULES["present_min"] and lo > 1 / RULES["ratio_low"]: return "X", "직전 회로의 들뜸이 남는다 + 두 회로가 번갈아 돈다", det
    if min(p0, p1) >= RULES["present_min"] and lo > RULES["ratio_same"][0] and hi < RULES["ratio_same"][1]:
        return "I", "직전 회로와 무관하다(시작부터 들떠 있거나 회로 안에서 들뜬다)", det
    return "혼합", "어느 경우와도 맞지 않는다(원인이 섞였을 수 있다). 비율과 구간만 보고한다", det


print("\n" + "=" * 78 + "\n4) 판정 규칙 — 실행 전에 고정한다\n" + "=" * 78)
print("    p₀, p₁ = P0·P1 에서 데이터 비트가 0000 인 비율. r = p₁/p₀ 와 그 95% 구간 [lo, hi] 를 함께 적는다.")
print(f"    둘 다 {RULES['gone_max']:.1%} 이하 → T. 둘 다 {RULES['present_min']:.0%} 미만(그러나 T 는 아님) → '약함'(판정하지 않는다).")
print(f"    p₀ ≥ 1% 이고 hi < 1/3 → S.  p₁ ≥ 1% 이고 lo > 3 → X.  둘 다 1% 이상이고 구간이 ({RULES['ratio_same'][0]}, {RULES['ratio_same'][1]:.2f}) 안 → I.  나머지 → '혼합'(판정하지 않는다).")
rng = np.random.default_rng(1105)
SELF = {}
print("    자체 점검(합성 표본 500묶음씩 — 실측이 아니다. 샷 잡음이 판정을 바꾸는지만 본다):")
for h in ("S", "X", "I", "T"):
    ok = sum(verdict(rng.binomial(SHOTS, PRED[h]["P0"]["p0000"]), rng.binomial(SHOTS, PRED[h]["P1"]["p0000"]), SHOTS)[0] == h for _ in range(500))
    SELF[h] = ok / 500
    print(f"      예측값 그대로의 {h}: 규칙이 {h} 를 고른 비율 {SELF[h]:.3f}")
MIS = {"절반만 잔류(p₁ = p₀/2)": (P_OBS, P_OBS / 2, "혼합"), "약해진 경우(둘 다 0.7%)": (0.007, 0.007, "약함"), "P1 만 조금 낮음(p₁ = 0.7 p₀)": (P_OBS, 0.7 * P_OBS, None)}
for name, (a, b, want) in MIS.items():
    out = [verdict(rng.binomial(SHOTS, a), rng.binomial(SHOTS, b), SHOTS)[0] for _ in range(500)]
    tally = {k: out.count(k) / 500 for k in sorted(set(out))}
    SELF[name] = tally
    print(f"      예측에 없는 경우 — {name}: {tally}")
assert all(SELF[h] > 0.99 for h in ("S", "X", "I", "T")) and SELF["절반만 잔류(p₁ = p₀/2)"].get("혼합", 0) > 0.95 and SELF["약해진 경우(둘 다 0.7%)"].get("약함", 0) > 0.95
print("      이 점검은 예측 모형이 맞다는 확인이 아니다. 모형 밖의 결과는 '혼합'·'약함'으로 떨어지고 그때는 비율만 보고한다.")

# ---------------- 5) QPU 를 쓰지 않는 확인과 그 읽는 법, 둘째 작업의 후보 ----------------
print("\n" + "=" * 78 + "\n5) 먼저 할 일(QPU 불필요)과 둘째 작업의 후보\n" + "=" * 78)
print(f"    (a) init_error(직전 실험이 |1⟩ 이었을 때 초기화 뒤에 남는 |1⟩ 의 비율) — env/fetch_ibm_job.py 로 2026-11-11 작업 시각의 값을 받는다.")
print(f"        S 가 맞다면 e = init_error × (보조가 끝에 |1⟩ 일 확률 {end_e:.3f}) 이므로 3번 큐비트의 init_error ≈ {F:.3f} 이어야 한다. "
      f"데이터 큐비트 {DATA} 는 각각 {', '.join(f'{r[2]:.4f}' for r in rows)} 이하여야 한다.")
print(f"        3번 큐비트의 값이 0.01~0.02 수준이면 S 는 QPU 를 쓰기 전에 탈락한다. 0.10 근처면 S 와 맞는다.")
print(f"    (b) 2026-11-11 작업의 샷별 비트열(fetch_ibm_job.py --bitstrings, 미실행) — 0000 인 {round(P_OBS * N_OBS)}샷이 고르게 퍼져 있는지 한때에 몰려 있는지 본다. "
      "배열의 순서가 실행 순서라는 보장은 확인하지 못했다.")
print("    둘째 작업의 후보(첫 작업의 결과와 실행 순서를 안 뒤에 고른다. 읽는 법은 그때 실행 전에 고정한다):")
for k, why in (("P2", "맨 앞에 보조 reset — 0000 이 그대로면 초기화 reset 자체, 게이트가 걸리지 않은 샷, 회로 안의 들뜸, 누설 가운데 하나다"),
               ("P3", "같은 BV 를 다른 허브(차수 3) 큐비트로 — 판독 오류와 T1 이 3번과 비슷한 큐비트를 고른다. 큐비트 고유인지 회로 구조인지"),
               ("P4", "끝에서 보조도 측정 — 측정 없이 끝나는 것이 원인인지. 0000 인 샷에서 보조를 1 로 읽는 비율은 |1⟩ 시작이면 약 0.5(정상 샷과 같다), 게이트가 걸리지 않은 샷이면 0, 누설(|2⟩)이면 1 에 가깝다"),
               ("PM", "게이트 없는 114번에 표지 X — 게이트가 걸리지 않은 샷이면 표지 비트가 0 이다"),
               ("Qm", "첫 CZ 앞에서 보조를 측정 — 0000 인 샷에서 그 비트가 1 이면 '첫 CZ 때 |1⟩ 이었다'는 전제가 직접 확인된다. 회로 중간 측정의 지원 여부는 미확인")):
    print(f"      {k}: {why}")

# QPU 시간: 2026-11-11 의 Sampler 기록 둘(Bell, 14장 회로 넷)에 직선을 맞춘다. 기록은 초 단위 정수라 ±0.5초다.
n_b, t_b = A.bell["shots"], A.bell["usage_seconds"]; n_s, t_s = len(pubs) * A.real["shots"], A.real["qpu_seconds"][0]
slope = (t_s - t_b) / (n_s - n_b); icpt = t_b - slope * n_b
QPU1 = icpt + slope * len(ORDER) * SHOTS
print(f"    QPU 시간 어림: 2026-11-11 의 Sampler 작업 둘({n_b}샷 {t_b}초, {n_s}샷 {t_s}초)에 맞춘 직선 {icpt:.1f}초 + {slope * 1e6:.0f} µs × 샷 수 → 첫 작업({len(ORDER) * SHOTS}샷) 약 {QPU1:.0f}초.")

PLAN = {"note": "appA_bv_confirm.py 가 만든 실험 계획. 고치지 말고 스크립트를 다시 돌린다. appA_bv_submit.py 가 읽는다.",
        "source": "runs/20261111-063819Z_ibm_fez_ch14_isa_calib.json (작업 dbq0qmoqtr0c738f53u0 의 PUB 1 = BV)",
        "backend": A.isa[A.SJ]["backend"], "shots": SHOTS, "order": ORDER, "roles": {"ancilla": ANC, "data_cz": DATA, "idle_measured": IDLE},
        "good_data_bits": GOODS, "variants": JOB1, "job2_candidates": {**CAND, "P3": "P0 의 큐비트 번호를 다른 허브와 그 이웃으로 바꾼다"},
        "rules": RULES, "predicted": PRED, "f_carry": float(F), "qpu_seconds_estimate": float(QPU1),
        "observed_20261111": {"p0000": P_OBS, "shots": N_OBS, "e_fit": A.e_fit}}
if __name__ == "__main__":
    (HERE / "appA_bv_plan.json").write_text(json.dumps(PLAN, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    print(f"\n    계획 파일: appA/appA_bv_plan.json (첫 작업의 PUB 순서 {ORDER})")

assert abs(F - 0.1035) < 0.0005 and abs(PRED["S"]["P0"]["p0000"] - P_OBS) < 2e-4 and PRED["S"]["P1"]["p0000"] < 0.002
assert abs(PRED["X"]["P1"]["p0000"] - P_OBS) < 2e-3 and 0.004 < PRED["X"]["P0"]["p0000"] < 0.009
assert abs(PRED["I"]["P1"]["p0000"] - P_OBS) < 2e-3 and PRED["T"]["P0"]["p0000"] < 0.001 and 3.5 < QPU1 < 4.5

# ---------------- 6) 분석 — 기록 파일이 주어졌을 때만 ----------------
if __name__ == "__main__" and len(sys.argv) > 1:
    rec = json.loads(pathlib.Path(sys.argv[1]).read_text(encoding="utf-8"))
    print("\n" + "=" * 78 + f"\n6) 실측 분석 — {sys.argv[1]}\n" + "=" * 78)
    print(f"    {rec.get('backend')} · {rec.get('t_utc')} · 작업 {rec.get('job_id')} · 사용 {rec.get('qpu_seconds')}초")
    C = {}; n = None
    for k in rec["order"]:
        c = rec["pubs"][k].get("counts")
        if not c: print(f"    {k}: 결과 없음"); continue
        tot = sum(c.values()); n = n or tot
        c0000 = sum(v for s, v in c.items() if s[-4:] == "0000"); good = sum(v for s, v in c.items() if s[-4:] == rec["good_data_bits"]) / tot
        C[k] = c0000; p = c0000 / tot
        print(f"    {k:3s} 샷 {tot}: P(0000) = {p:.4f} ± {np.sqrt(p * (1 - p) / tot):.4f} | 성공 {good:.4f}")
        bs = rec["pubs"][k].get("bitstrings")
        if bs:                                       # 샷별 자료가 있으면: 구간별 빈도와 연달아 나오는 정도(배열 순서 = 실행 순서라는 보장은 확인하지 못했다)
            z = np.array([s[-4:] == "0000" for s in bs]); blocks = [int(b.sum()) for b in np.array_split(z, 8)]
            nxt = z[1:][z[:-1]]; print(f"        샷별: 8구간의 0000 횟수 {blocks}; 0000 다음 샷이 0000 인 비율 {nxt.mean() if len(nxt) else float('nan'):.3f} (전체 {z.mean():.4f})")
    if set(ORDER) <= set(C):
        h, text, det = verdict(C["P0"], C["P1"], n)
        print(f"    r = p₁/p₀ = {det['ratio_p1_over_p0']:.3f}, 95% 구간 [{det['ci'][0]:.3f}, {det['ci'][1]:.3f}]")
        print(f"    판정: {h} — {text}")
        if h in PRED:
            print(f"    그 경우의 실행 전 예측(2026-11-11 보정값): P(0000) = {PRED[h]['P0']['p0000']:.4f}, {PRED[h]['P1']['p0000']:.4f}; "
                  f"성공 P0 {PRED[h]['P0']['success']:.4f} | {PRED[h]['P0']['success_no_gate_relax']:.4f}, P1 {PRED[h]['P1']['success']:.4f} | {PRED[h]['P1']['success_no_gate_relax']:.4f}")
    else:
        print("    P0·P1 이 다 없어 판정하지 않는다")
