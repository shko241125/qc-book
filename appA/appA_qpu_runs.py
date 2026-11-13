"""부록 A — 첫 QPU 실행 기록의 해석 (numpy 구현, 네트워크 불필요)

읽는 파일(모두 저장소에 커밋됨):
  runs/20261111-062150Z_ibm_kingston_bell2.json        Bell 1개 (env/run_ibm.py)
  runs/20261111-063819Z_ibm_fez_ch14_record.json, runs/20261111-063700Z_fake_marrakesh_ch14_record.json
                                                       14장 회로 5개 (ch14_submit.py real / fake 기록의 고정 사본 — ch14/ 의 기록은 다음 실행이 덮어쓴다)
  runs/20261111-063819Z_ibm_fez_ch14_isa_calib.json    실제로 돈 ISA 회로와 그 시각의 보정값 (env/fetch_ibm_job.py)

하는 일:
  A1 Bell: TVD 와 샷 잡음 기준, 00/11 밖으로 샌 비율.
  A2 14장 회로: 실측 성공 확률(± 이항 표준오차)·TVD 를 14장 예측(ch14_hardware.py 출력)·fake 와 나란히.
  A3 '당일 보정값 재현': 실제 ISA 회로를 큐비트·연결별 보정값(게이트 오류, T1·T2, 비대칭 판독 오류)으로 다시 돌린다.
  A4 BV 의 '0000' 가설: 단일 파울리 결함을 모든 위치에 넣어 '0000' 을 만드는 결함을 찾고, 보조 큐비트 초기화 오류의 크기를 맞춘다.
  A5 HVA: 같은 재현으로 에너지와 판독 오류 몫, 회로가 만든 에너지 하강 중 남은 비율.
비트 순서: 기록의 counts 는 Qiskit 문자열(c_{m-1}…c0)이다. BV·Grover·QPE 는 측정 매핑으로 그것이 곧 책 순서지만,
GHZ 는 measure_all 이라 Qiskit 순서 그대로다 — GHZ 를 책 순서로 읽을 때만 뒤집는다.
"""
import itertools
import json
import pathlib

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
load = lambda p: json.loads((ROOT / p).read_text(encoding="utf-8"))
bell = load("runs/20261111-062150Z_ibm_kingston_bell2.json")
real, fake = load("runs/20261111-063819Z_ibm_fez_ch14_record.json"), load("runs/20261111-063700Z_fake_marrakesh_ch14_record.json")
isa = load("runs/20261111-063819Z_ibm_fez_ch14_isa_calib.json")
SJ, EJ = real["jobs"]["ghz5"]["job_id"], real["jobs"]["hva"]["job_id"]
NAMES = ["ghz5", "bv1011", "grover101", "qpe_t"]

# 14장 예측(ch14_hardware.py 출력 2) — 공개 중앙값 모형) — 비교 기준으로만 인용한다
PRED = {"ghz5": 0.918, "bv1011": 0.925, "grover101": 0.759, "qpe_t": 0.881}
PRED_HVA, PRED_HVA_NORO, E_HVA_IDEAL = -4.4329, -4.5551, -4.6462
th = np.arcsin(1 / np.sqrt(8)); P_GROVER = np.sin(5 * th) ** 2
IDEAL = {"ghz5": {"00000": .5, "11111": .5}, "bv1011": {"1011": 1.0}, "qpe_t": {"001": 1.0},
         "grover101": {**{format(i, "03b"): (1 - P_GROVER) / 7 for i in range(8)}, "101": P_GROVER}}
GOOD = {k: [b for b, v in IDEAL[k].items() if v > 0.2] for k in IDEAL}


def tvd(p, q):
    return 0.5 * sum(abs(p.get(k, 0) - q.get(k, 0)) for k in set(p) | set(q))


def norm(c):
    n = sum(c.values()); return {k: v / n for k, v in c.items()}, n


# ---------------- A1 Bell ----------------
cb = bell["counts_book_order"]; pb, nb = norm(cb)
leak = pb.get("01", 0) + pb.get("10", 0)
print(f"A1) Bell · {bell['backend']} · job {bell['job_id']} · {nb}샷 · CZ {bell['two_qubit_gates']} · 깊이 {bell['depth']} · 사용 {bell['usage_seconds']}초")
print(f"    책 순서 counts {cb} | TVD {tvd(pb, {'00': .5, '11': .5}):.4f}, 샷 잡음만일 때 평균 ≈ {np.sqrt(2 / np.pi) / (2 * np.sqrt(nb)):.4f} "
      f"| 00·11 밖 {leak:.4f} | 00:11 비 {pb['00'] / (pb['00'] + pb['11']):.3f} "
      f"(0.5 에서 {(cb['00'] / (cb['00'] + cb['11']) - 0.5) / np.sqrt(0.25 / (cb['00'] + cb['11'])):+.1f}σ)")

# ---------------- A2 실측 대 예측 ----------------
print(f"\nA2) 14장 회로 — 실측 {real['backend']} (sampler {SJ}, estimator {EJ}, {real['shots']}샷), fake {fake['backend']}")
meas = {}
for k in NAMES:
    p, n = norm(real["jobs"][k]["counts"]); pf, _ = norm(fake["jobs"][k]["counts"])
    s = sum(p.get(g, 0) for g in GOOD[k]); sf = sum(pf.get(g, 0) for g in GOOD[k])
    cz = sum(1 for g in isa[SJ]["pubs"][NAMES.index(k)] if g[0] == "cz")
    meas[k] = dict(p=p, n=n, s=s)
    print(f"    {k:10s} CZ {cz:3d} | 예측 {PRED[k]:.3f} | fake {sf:.4f} | 실측 {s:.4f} ± {np.sqrt(s * (1 - s) / n):.4f} "
          f"(실측−예측 {100 * (s - PRED[k]):+.1f}%p) | 실측 TVD {tvd(p, IDEAL[k]):.4f}")
for k in NAMES:
    top = sorted(real["jobs"][k]["counts"].items(), key=lambda kv: -kv[1])
    bad = [(b if k != "ghz5" else b[::-1], v) for b, v in top if b not in GOOD[k]][:3]
    print(f"    {k:10s} 오답 상위(책 순서): {bad}")
ev, sd = isa[EJ]["estimator"]["evs"], isa[EJ]["estimator"]["stds"]
print(f"    HVA ⟨H⟩ 실측 {ev:.4f} ± {sd:.4f} | 예측 {PRED_HVA} (판독 제외 {PRED_HVA_NORO}) | 차 {ev - PRED_HVA:+.3f} = {(ev - PRED_HVA) / sd:+.1f}σ "
      f"| fake {fake['jobs']['hva']['evs']:.4f} (Aer: 판독 오류 없는 정확값 + N(0, 0.02) 표본 잡음)")
print(f"    QPU 사용: Bell {bell['usage_seconds']}초 + sampler·estimator {real['qpu_seconds']}초 = {bell['usage_seconds'] + sum(real['qpu_seconds'])}초 "
      f"(Open 플랜 600초의 {100 * (bell['usage_seconds'] + sum(real['qpu_seconds'])) / 600:.1f}%)")

# ---------------- 재현 시뮬레이터 ----------------
T_1Q, T_CZ = 32e-9, 68e-9          # 14장 NOISE 의 Heron 게이트 시간(가정) — 유휴 열완화 시계용
I2 = np.eye(2, dtype=complex); X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1, -1]).astype(complex)
SX = 0.5 * np.array([[1 + 1j, 1 - 1j], [1 - 1j, 1 + 1j]])
RZ = lambda t: np.diag([np.exp(-0.5j * t), np.exp(0.5j * t)])
CZ = np.diag([1, 1, 1, -1]).astype(complex)
RX = lambda t: np.cos(t / 2) * I2 - 1j * np.sin(t / 2) * X
H1 = np.array([[1, 1], [1, -1]], complex) / np.sqrt(2)


class Replay:
    """ISA 게이트 목록을 밀도행렬로 다시 돌린다. 로컬 순서 = 물리 큐비트 번호 오름차순."""
    def __init__(self, qubits, cal, noisy=True, readout=True, init_exc=None):
        self.q = sorted(qubits); self.loc = {p: i for i, p in enumerate(self.q)}; n = self.n = len(self.q)
        self.cal, self.noisy, self.readout = cal, noisy, readout
        self.rho = np.zeros((2 ** n, 2 ** n), complex); self.rho[0, 0] = 1
        self.clock = dict.fromkeys(self.q, 0.0)
        for p, e in (init_exc or {}).items():     # 초기화 오류: 확률 e 로 |1⟩ 에서 시작
            self.kraus([np.sqrt(1 - e) * I2, np.sqrt(e) * X], [p])

    def kraus(self, ops, phys):
        n, k = self.n, len(phys); qs = [self.loc[p] for p in phys]
        r = self.rho.reshape([2] * 2 * n); out = 0
        for K in ops:
            Kt = K.reshape([2] * 2 * k)
            t = np.moveaxis(np.tensordot(Kt, r, axes=(list(range(k, 2 * k)), qs)), list(range(k)), qs)
            t = np.tensordot(t, Kt.conj(), axes=([n + q for q in qs], list(range(k, 2 * k))))
            out = out + np.moveaxis(t, list(range(2 * n - k, 2 * n)), [n + q for q in qs])
        self.rho = out.reshape(2 ** n, 2 ** n)

    def depol(self, r, phys):                      # RB 오류율 r → 탈분극 p = r(d+1)/d (14장과 같은 변환)
        k = len(phys); d = 2 ** k; p = r * (d + 1) / d
        ops = []
        for c in itertools.product([I2, X, Y, Z], repeat=k):
            P = c[0] if k == 1 else np.kron(c[0], c[1])
            w = 1 - p if all(m is I2 for m in c) else p / (d * d - 1)
            ops.append(np.sqrt(w) * P)
        self.kraus(ops, phys)

    def relax(self, t, p):
        c = self.cal["qubits"][str(p)]; T1, T2 = c["T1_us"] * 1e-6, c["T2_us"] * 1e-6
        g = 1 - np.exp(-t / T1); lam = 1 - np.exp(-t * max(1 / T2 - 1 / (2 * T1), 0))
        self.kraus([np.diag([1, np.sqrt(1 - g)]).astype(complex), np.array([[0, np.sqrt(g)], [0, 0]], complex)], [p])
        self.kraus([np.sqrt(1 - lam) * I2, np.sqrt(lam) * np.diag([1, 0]).astype(complex), np.sqrt(lam) * np.diag([0, 1]).astype(complex)], [p])

    def wait(self, phys, t0):
        for p in phys:
            if self.noisy and t0 > self.clock[p] + 1e-15: self.relax(t0 - self.clock[p], p)
            self.clock[p] = max(self.clock[p], t0)

    def gate(self, name, phys, params):
        if name == "rz": self.kraus([RZ(params[0])], phys); return    # 가상 게이트: 오류·시간 없음
        if name in ("h", "rx"):                   # A5 의 논리 회로 재구성용: 오류율·시간은 sx 하나로 친다
            U, dt, name = (H1 if name == "h" else RX(params[0])), T_1Q, "sx"
        else:
            U, dt = {"sx": (SX, T_1Q), "x": (X, T_1Q), "cz": (CZ, T_CZ)}[name]
        self.wait(phys, max(self.clock[p] for p in phys))
        self.kraus([U], phys)
        if self.noisy:
            self.depol(self.cal["gates"][f"{name}{phys}"], phys)
            for p in phys: self.relax(dt, p)
        for p in phys: self.clock[p] += dt

    def run(self, gates, upto=None):
        for g in gates[:upto]:
            if g[0] != "measure": self.gate(g[0], g[1], g[2])
        return self

    def dist(self, meas):
        """측정 (물리, 고전 비트) 목록 → Qiskit 문자열(c_{m-1}…c0) 분포."""
        if self.noisy: self.wait(self.q, max(self.clock.values()))
        p = np.real(np.diag(self.rho)).clip(0).reshape([2] * self.n)
        if self.readout and self.noisy:
            for ph, _ in meas:
                c = self.cal["qubits"][str(ph)]; a, b = c["p_meas1_prep0"], c["p_meas0_prep1"]
                M = np.array([[1 - a, b], [a, 1 - b]])
                p = np.moveaxis(np.tensordot(M, p, axes=([1], [self.loc[ph]])), 0, self.loc[ph])
        m = len(meas); out = {}
        for idx in itertools.product([0, 1], repeat=self.n):
            if p[idx] < 1e-15: continue
            bits = ["0"] * m
            for ph, cbit in meas: bits[m - 1 - cbit] = str(idx[self.loc[ph]])
            key = "".join(bits); out[key] = out.get(key, 0) + p[idx]
        return out


def replay(gates, cal, **kw):
    meas = [(g[1][0], g[3][0]) for g in gates if g[0] == "measure"]
    qs = {q for g in gates for q in g[1]}
    return Replay(qs, cal, **kw).run(gates).dist(meas)


cal = isa[SJ]["calibration"]
print(f"\nA3) 당일 보정값 재현 (보정 시각 {cal['last_update_date']}, 작업 {isa[SJ]['created']}) — 게이트 탈분극 + T1·T2(ASAP 유휴 포함) + 큐비트별 비대칭 판독")
czs = {k: [v for kk, v in cal["gates"].items() if kk.startswith("cz")] for k in ["all"]}["all"]
ro = [v["readout_error"] for v in cal["qubits"].values()]
print(f"    이 작업에 쓰인 연결의 CZ 오류 {min(czs):.2e} ~ {max(czs):.2e} (14장 중앙값 3.0e-03), 큐비트 {len(ro)}개의 판독 오류 {min(ro):.4f} ~ {max(ro):.4f} (14장 평균 0.013)")
rep = {}
for i, k in enumerate(NAMES):
    g = isa[SJ]["pubs"][i]
    p0 = replay(g, cal, noisy=False)
    s0 = sum(p0.get(b, 0) for b in GOOD[k])
    assert abs(s0 - sum(IDEAL[k][b] for b in GOOD[k])) < 1e-9, (k, s0)      # ISA 회로·비트 매핑 검산: 이상 성공 확률 재현
    pn, pr = replay(g, cal), replay(g, cal, readout=False)
    sn, sr = sum(pn.get(b, 0) for b in GOOD[k]), sum(pr.get(b, 0) for b in GOOD[k])
    rep[k] = pn; s, n = meas[k]["s"], meas[k]["n"]
    print(f"    {k:10s} 이상 {s0:.4f} | 재현 {sn:.4f} (판독 제외 {sr:.4f}) | 실측 {s:.4f} | 실측−재현 {100 * (s - sn):+.1f}%p = {(s - sn) / np.sqrt(s * (1 - s) / n):+.1f}σ")
g_ghz = rep["ghz5"]; q4flip = g_ghz.get("01111", 0) + g_ghz.get("10000", 0)
q4m = meas['ghz5']['p'].get('01111', 0) + meas['ghz5']['p'].get('10000', 0)
print(f"    GHZ q4 하나만 뒤집힘(책 순서 11110·00001): 재현 {q4flip:.4f} | 실측 {q4m:.4f} | 초과 {100 * (q4m - q4flip):.1f}%p")
import datetime as _dt
_gap = _dt.datetime.fromisoformat(isa[SJ]["created"]) - _dt.datetime.fromisoformat(cal["last_update_date"])
print(f"    보정 시각과 작업 생성 시각의 차: {_gap.seconds // 60}분 {_gap.seconds % 60}초")
gg = isa[SJ]["pubs"][0]; c2p = {g[3][0]: g[1][0] for g in gg if g[0] == "measure"}     # measure_all: 고전 비트 i = 논리 q_i
p4, p3 = c2p[4], c2p[3]; cq = cal["qubits"][str(p4)]
print(f"    GHZ 논리 q0…q4 = 물리 {[c2p[i] for i in range(5)]}; q4 = 물리 {p4}: 판독 오류 {cq['readout_error']:.4f}, T1 {cq['T1_us']:.0f} µs, T2 {cq['T2_us']:.0f} µs, "
      f"cz[{p3}, {p4}] 오류 {cal['gates'].get(f'cz[{p3}, {p4}]', cal['gates'].get(f'cz[{p4}, {p3}]')):.2e}; 물리 {p3} 의 T2 {cal['qubits'][str(p3)]['T2_us']:.0f} µs")

# ---------------- A4 BV '0000' ----------------
k = "bv1011"; gb = isa[SJ]["pubs"][1]; mb = [(g[1][0], g[3][0]) for g in gb if g[0] == "measure"]
anc = [q for q in {q for g in gb for q in g[1]} if q not in [m[0] for m in mb]]
data_cz = sorted({q for g in gb if g[0] == "cz" for q in g[1]} - set(anc))
idle = [m[0] for m in mb if not any(g[0] != "measure" and m[0] in g[1] for g in gb)]
print(f"\nA4) BV 보조 큐비트 가설 — ISA 회로: 물리 큐비트 {sorted({q for g in gb for q in g[1]})}, 보조(측정 안 함) {anc}, "
      f"보조와 CZ 하는 데이터 {data_cz}, 게이트 없이 측정만 하는 큐비트 {idle}(s 의 0 자리: H·H = I 로 지워짐)")
order = [(g[0], g[1]) for g in gb if g[0] != "measure"]
print("    ISA 게이트 순서(보조 관련):", [f"{n}{q}" for n, q in order if anc[0] in q])
qs = {q for g in gb for q in g[1]}
hits = []
for pos in range(-1, len(order)):              # pos 번째 게이트 '뒤'에 파울리 하나 (−1 = 시작 직후)
    for q in sorted(qs):
        for nm, P in (("X", X), ("Y", Y), ("Z", Z)):
            R = Replay(qs, cal, noisy=False)
            if pos == -1: R.kraus([P], [q])
            for i, g in enumerate([g for g in gb if g[0] != "measure"]):
                R.gate(g[0], g[1], g[2])
                if i == pos: R.kraus([P], [q])
            p0 = R.dist(mb).get("0000", 0)
            if p0 > 0.5: hits.append((pos, q, nm, p0))
first_cz = next(i for i, (n_, q_) in enumerate(order) if n_ == "cz" and anc[0] in q_)
assert hits and all(q == anc[0] and pos < first_cz and nm in "XY" and p0 > 0.999 for pos, q, nm, p0 in hits)
print(f"    단일 파울리 결함 {3 * len(qs) * (len(order) + 1)}가지 중 '0000' 을 만드는 것 {len(hits)}가지 — 모두 보조 {anc[0]} 의 X·Y, "
      f"첫 CZ {order[first_cz][0]}{order[first_cz][1]} 이전(= 보조가 |1⟩ 에서 시작한 것과 같다), 각각 P(0000) = 1")
p_obs = meas[k]["p"].get("0000", 0)
p_rep = rep[k].get("0000", 0)
c3 = cal["qubits"][str(anc[0])]
print(f"    P(0000): 실측 {p_obs:.4f} | 당일 보정 재현 {p_rep:.5f} (실측의 {p_rep / p_obs:.1%})")
es = np.linspace(0, 0.12, 241)
f = [replay(gb, cal, init_exc={anc[0]: e}).get("0000", 0) for e in es]
e_fit = float(np.interp(p_obs, f, es))
pr_fit = replay(gb, cal, init_exc={anc[0]: e_fit})
s_fit = sum(pr_fit.get(b, 0) for b in GOOD[k])
print(f"    보조 {anc[0]} 이 확률 e 로 |1⟩ 에서 시작한다고 두면 P(0000) 실측과 맞는 e = {e_fit:.4f} "
      f"→ 그때 성공 {s_fit:.4f} (실측 {meas[k]['s']:.4f}); 비교: 보정값 P(1|0 준비) {c3['p_meas1_prep0']:.4f}, 판독 오류 {c3['readout_error']:.4f} "
      f"(이 작업의 큐비트 {len(cal['qubits'])}개 중 판독 오류 최대: {c3['readout_error'] == max(v['readout_error'] for v in cal['qubits'].values())})")
for b in ("1000", "1010", "0011", "1001"):
    print(f"    {b}: 실측 {meas[k]['p'].get(b, 0):.4f} | 재현(e=0) {rep[k].get(b, 0):.4f} | 재현(e 맞춤) {pr_fit.get(b, 0):.4f} | 실측−재현 {100 * (meas[k]['p'].get(b, 0) - rep[k].get(b, 0)):+.1f}%p")
sb = sum(rep[k].get(b, 0) for b in GOOD[k])
print(f"    정리: 실측−재현 {100 * (meas[k]['s'] - sb):+.1f}%p 중 0000 초과분 {100 * (p_obs - p_rep):.1f}%p; e 를 넣은 뒤 남는 차 {100 * (meas[k]['s'] - s_fit):+.1f}%p; "
      f"e / P(1|0 준비) = {e_fit / c3['p_meas1_prep0']:.2f}")

# ---------------- A5 HVA ----------------
# Estimator 의 ISA 회로는 마지막 Rx(−2β) 층을 측정 기저 변경 블록(매개변수, 값은 작업 안의 samplex 에 있음)에 합쳐 버려서
# 숫자 부분만으로는 HVA 상태가 되지 않는다. 그래서 ISA 와 같은 물리 사슬·같은 CZ 6개로 13장 논리 회로를 다시 짠다
# (CX = H·CZ·H, 1큐비트 게이트 수는 ISA 와 다를 수 있다 — sx 오류 2.5e-4 수준이라 에너지에 주는 차는 작다).
gh = isa[EJ]["pubs"][0]; calh = isa[EJ]["calibration"]
mh = [(g[1][0], g[3][0]) for g in gh if g[0] == "measure"]
chain = [q for q, _ in sorted(mh, key=lambda m: m[1])]                              # 고전 비트 i ↔ 사슬 i 번째
assert [tuple(sorted(g[1])) for g in gh if g[0] == "cz"] == [tuple(sorted(e)) for e in zip(chain, chain[1:]) for _ in (0, 1)], "사슬·CZ 순서가 ISA 와 다르다"
GAMMA, BETA = 0.2064, 0.3927                                                        # ch14_submit.hva_tfim 과 같은 값


def hva_gates(basis):
    g = [["h", [q], []] for q in chain]
    for a, b in zip(chain, chain[1:]):
        cx = [["h", [b], []], ["cz", sorted([a, b]), []], ["h", [b], []]]
        g += cx + [["rz", [b], [-2 * GAMMA]]] + cx
    g += [["rx", [q], [-2 * BETA]] for q in chain]
    return g + ([["h", [q], []] for q in chain] if basis == "X" else [])


def energy(c, **kw):
    e = 0.0
    for basis in "ZX":
        d = Replay(set(chain), c, **kw).run(hva_gates(basis)).dist(mh)
        for s_, v in d.items():
            z = [1 - 2 * int(s_[len(s_) - 1 - i]) for i in range(4)]       # 고전 비트 i 의 고윳값
            e += v * (-sum(z[i] * z[i + 1] for i in range(3)) if basis == "Z" else -sum(z))
    return e


e0, en, enr = energy(calh, noisy=False), energy(calh), energy(calh, readout=False)
assert abs(e0 - E_HVA_IDEAL) < 1e-3, e0                                   # 재구성 회로가 13장 HVA 상태를 만든다는 검산
cal_no3 = json.loads(json.dumps(calh)); cal_no3["qubits"]["3"].update(p_meas1_prep0=0, p_meas0_prep1=0)
en_no3 = energy(cal_no3)
print(f"\nA5) HVA — 물리 사슬 {chain}, CZ 6 | 이상 {e0:.4f} | 당일 보정 재현 {en:.4f} (판독 제외 {enr:.4f}; 큐비트 3 판독만 제외 {en_no3:.4f}) "
      f"| 실측 {ev:.4f} ± {sd:.4f} → 실측−재현 {ev - en:+.3f} = {(ev - en) / sd:+.1f}σ")
print(f"    판독 오류 몫(재현 안에서): 전체 {enr - en:+.3f} 중 큐비트 3 {en_no3 - en:+.3f}")
print(f"    에너지 하강(|+⟩⁴ 의 −4 기준): 이상 {-4 - e0:.3f}, 재현 {-4 - en:.3f} ({(-4 - en) / (-4 - e0):.0%}), 실측 {-4 - ev:.3f} ({(-4 - ev) / (-4 - e0):.0%}), "
      f"fake {-4 - fake['jobs']['hva']['evs']:.3f} ({(-4 - fake['jobs']['hva']['evs']) / (-4 - e0):.0%})")
