"""14장 하드웨어 실행 1 — numpy 검증 스크립트
IBM Heron r2 급 공개 사양(중앙값)으로 노이즈 모델을 만들고, 8~13장 대표 회로 5개를
(1) heavy-hex 조각에 라우팅하고 (2) 밀도행렬 노이즈 에뮬레이션으로 QPU 결과를 예측한다.
규약: big-endian (q0 = 최상위 비트), seed=7.  실제 QPU 제출은 ch14_submit.py (미실행).
"""
import numpy as np
from collections import deque
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
matplotlib.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

rng = np.random.default_rng(7)

# ---------- 노이즈 파라미터: Heron r2 공개 사양 중앙값 (ibm_fez / ibm_marrakesh) ----------
# 사양서의 게이트 오류는 RB 오류율 r = 1 − F_avg. 탈분극 채널의 비항등 파울리 총확률 p 와는
# r = p·d/(d+1) 관계 (d = 2^k) → p = r·(d+1)/d: 2큐비트 1.25배, 1큐비트 1.5배 (DM.depolarize 직전에 변환).
# 판독 오류 ε 는 판독 중 T1 붕괴를 이미 포함한 측정값이므로 판독 전 열완화를 따로 넣지 않는다.
NOISE = dict(
    r_cz=3.0e-3,      # CZ RB 오류율 중앙값
    r_1q=2.5e-4,      # SX/X RB 오류율
    eps01=1.0e-2,     # 판독: 0 을 1 로 읽을 확률
    eps10=1.6e-2,     # 판독: 1 을 0 으로 읽을 확률 (평균 1.3e-2)
    T1=180e-6, T2=120e-6,
    t_1q=32e-9, t_cz=68e-9, t_meas=1.5e-6,
    rb_convert=True,  # False 면 p = r (이전 판의 단순화)
    idle=True,        # 유휴 큐비트도 기다린 시간만큼 열완화 (ASAP 일정)
    readout=True,     # False 면 판독 오류 없음 (TREX 로 완화된 Estimator 와 비교용)
)


def p_dep(r, k, noise):
    d = 2 ** k
    return r * (d + 1) / d if noise["rb_convert"] else r

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Z = np.diag([1.0, -1.0]).astype(complex)
H1 = np.array([[1, 1], [1, -1]]) / np.sqrt(2)


def rz(t): return np.diag([np.exp(-1j * t / 2), np.exp(1j * t / 2)])
def rx(t): return np.cos(t / 2) * I2 - 1j * np.sin(t / 2) * X


# ---------- 밀도행렬 시뮬레이터 ----------
class DM:
    """n 큐비트 밀도행렬. 게이트 뒤에 (탈분극 + T1/T2 열완화) 채널, 측정 때 판독 오류."""
    def __init__(self, n, noisy=True, noise=NOISE):
        self.n, self.noisy, self.p = n, noisy, noise
        self.rho = np.zeros((2 ** n, 2 ** n), complex); self.rho[0, 0] = 1.0
        self.count = dict(cz=0, sq=0)
        self.clock = np.zeros(n)                      # 큐비트별 현재 시각 (ASAP 일정)

    def _wait(self, qs, t0):
        """qs 를 시각 t0 까지 기다리게 한다 — 유휴 동안의 T1/T2 열완화."""
        for q in qs:
            if self.noisy and self.p["idle"] and t0 > self.clock[q] + 1e-15:
                self.relax(t0 - self.clock[q], [q])
            self.clock[q] = max(self.clock[q], t0)

    def _apply_ops(self, ops, qs):
        """크라우스 연산자 목록 ops (각각 2^k×2^k) 를 큐비트 qs 에 적용."""
        n = self.n; k = len(qs)
        rho = self.rho.reshape([2] * (2 * n))
        out = np.zeros_like(rho)
        for K in ops:
            Kt = K.reshape([2] * (2 * k))
            r = np.tensordot(Kt, rho, axes=(list(range(k, 2 * k)), qs))
            r = np.moveaxis(r, list(range(k)), qs)
            Kc = Kt.conj()
            r = np.tensordot(r, Kc, axes=([n + q for q in qs], list(range(k, 2 * k))))
            r = np.moveaxis(r, list(range(2 * n - k, 2 * n)), [n + q for q in qs])
            out += r
        self.rho = out.reshape(2 ** n, 2 ** n)

    def unitary(self, U, qs): self._apply_ops([U], qs)

    def depolarize(self, p, qs):
        k = len(qs); d = 2 ** k
        paulis = [I2, X, 1j * X @ Z, Z]
        ops = []
        import itertools
        for combo in itertools.product(range(4), repeat=k):
            P = np.array([[1.0 + 0j]])
            for c in combo: P = np.kron(P, paulis[c])
            w = (1 - p) if all(c == 0 for c in combo) else p / (d * d - 1)
            ops.append(np.sqrt(w) * P)
        self._apply_ops(ops, qs)

    def relax(self, t, qs):
        T1, T2 = self.p["T1"], self.p["T2"]
        g = 1 - np.exp(-t / T1)                       # 진폭 감쇠
        # 순수 위상 감쇠 (1/Tφ = 1/T2 − 1/2T1). 아래 Kraus 세트 {√(1−λ)I, √λ|0⟩⟨0|, √λ|1⟩⟨1|} 는
        # 비대각을 (1−λ) 배 하므로 λ = 1 − e^{−t/Tφ}. (7장의 {diag(1,√(1−λ)), diag(0,√λ)} 세트라면 지수에 2 가 붙는다.)
        lam = 1 - np.exp(-t * (1 / T2 - 1 / (2 * T1)))
        amp = [np.array([[1, 0], [0, np.sqrt(1 - g)]]), np.array([[0, np.sqrt(g)], [0, 0]])]
        ph = [np.sqrt(1 - lam) * I2, np.sqrt(lam) * np.diag([1.0, 0]), np.sqrt(lam) * np.diag([0, 1.0])]
        for q in qs:
            self._apply_ops(amp, [q]); self._apply_ops(ph, [q])

    # 네이티브 게이트: rz(가상, 무오류), sx/x/h (1큐비트, 오류 p_1q), cz (오류 p_cz)
    def g1(self, U, q, virtual=False):
        self.unitary(U, [q])
        if not virtual:
            self.count["sq"] += 1
            self.clock[q] += self.p["t_1q"]
            if self.noisy: self.depolarize(p_dep(self.p["r_1q"], 1, self.p), [q]); self.relax(self.p["t_1q"], [q])

    def h(self, q): self.g1(H1, q)
    def x(self, q): self.g1(X, q)
    def rz_(self, q, t): self.g1(rz(t), q, virtual=True)
    def rx_(self, q, t): self.g1(rx(t), q)

    def cz(self, a, b):
        U = np.diag([1, 1, 1, -1]).astype(complex)
        self._wait([a, b], max(self.clock[a], self.clock[b]))
        self.unitary(U, [a, b]); self.count["cz"] += 1
        self.clock[a] += self.p["t_cz"]; self.clock[b] += self.p["t_cz"]
        if self.noisy: self.depolarize(p_dep(self.p["r_cz"], 2, self.p), [a, b]); self.relax(self.p["t_cz"], [a, b])

    def cx(self, c, t): self.h(t); self.cz(c, t); self.h(t)
    def swap(self, a, b): self.cx(a, b); self.cx(b, a); self.cx(a, b)

    def probs(self, qs=None):
        if self.noisy: self._wait(list(range(self.n)), self.clock.max())   # 동시 판독: 먼저 끝난 큐비트는 기다린다
        p = np.real(np.diag(self.rho)).copy()
        if self.noisy and self.p["readout"]:
            # 판독 혼동 행렬을 큐비트마다 곱한다
            M = np.array([[1 - self.p["eps01"], self.p["eps10"]], [self.p["eps01"], 1 - self.p["eps10"]]])
            p = p.reshape([2] * self.n)
            for q in range(self.n):
                p = np.moveaxis(np.tensordot(M, p, axes=([1], [q])), 0, q)
            p = p.reshape(-1)
        if qs is None: return p
        # 부분 큐비트 주변 분포
        p = p.reshape([2] * self.n)
        keep = sorted(qs); drop = [q for q in range(self.n) if q not in keep]
        p = p.sum(axis=tuple(drop)) if drop else p
        return p.reshape(-1)


def bits(k, n): return format(k, f"0{n}b")
def tvd(p, q): return 0.5 * np.abs(p - q).sum()


# ---------- heavy-hex 조각과 라우팅 ----------
# 6개 물리 큐비트: 선 0-1-2-3-4 에 2 번이 가지 5 를 가진 "T" 조각 (차수 3 노드 = 2). heavy-hex 는 삼각형이 없다.
EDGES = [(0, 1), (1, 2), (2, 3), (3, 4), (2, 5)]
ADJ = {i: set() for i in range(6)}
for a, b in EDGES: ADJ[a].add(b); ADJ[b].add(a)


def shortest_path(a, b):
    prev = {a: None}; dq = deque([a])
    while dq:
        u = dq.popleft()
        if u == b: break
        for v in ADJ[u]:
            if v not in prev: prev[v] = u; dq.append(v)
    path = [b]
    while path[-1] != a: path.append(prev[path[-1]])
    return path[::-1]


class Routed:
    """논리 큐비트 → 물리 큐비트 사상을 유지하며, 인접하지 않은 2큐비트 게이트 앞에 SWAP 사슬을 넣는다."""
    def __init__(self, dm, layout):
        self.dm = dm; self.l2p = dict(layout); self.swaps = 0
    def _phys(self, q): return self.l2p[q]
    def _bring_adjacent(self, a, b):
        pa, pb = self._phys(a), self._phys(b)
        if pb in ADJ[pa]: return
        path = shortest_path(pa, pb)
        for u, v in zip(path[:-1], path[1:-1]):   # pa 를 pb 옆까지 밀어 간다 (마지막 변은 게이트 자체)
            self.dm.swap(u, v); self.swaps += 1
            p2l = {p: l for l, p in self.l2p.items()}
            lu, lv = p2l[u], p2l.get(v)
            self.l2p[lu] = v
            if lv is not None: self.l2p[lv] = u
    def h(self, q): self.dm.h(self._phys(q))
    def x(self, q): self.dm.x(self._phys(q))
    def rz(self, q, t): self.dm.rz_(self._phys(q), t)
    def rx(self, q, t): self.dm.rx_(self._phys(q), t)
    def cz(self, a, b): self._bring_adjacent(a, b); self.dm.cz(self._phys(a), self._phys(b))
    def cx(self, c, t): self._bring_adjacent(c, t); self.dm.cx(self._phys(c), self._phys(t))
    def cp(self, a, b, t):                       # 제어 위상 = Rz·CX·Rz·CX·Rz (CX 2개)
        self.rz(a, t / 2); self.rz(b, t / 2); self.cx(a, b); self.rz(b, -t / 2); self.cx(a, b)
    def ccz_linear(self, a, b, c):
        """선형 연결용 CCZ (Cirq CCZPowGate 분해: CX 8 + T/T† 7). 큐비트 b 가 a·c 양쪽과 인접해야 한다."""
        t = np.pi / 4
        sweep = lambda: (self.cx(a, b), self.cx(b, c))
        self.rz(a, t); self.rz(b, t); self.rz(c, t)
        sweep(); self.rz(b, -t); self.rz(c, t)
        sweep(); self.rz(c, -t)
        sweep(); self.rz(c, -t)
        sweep()
    def probs(self, logical_qs):
        p = self.dm.probs(None).reshape([2] * self.dm.n)
        phys = [self._phys(q) for q in logical_qs]
        drop = [q for q in range(self.dm.n) if q not in phys]
        p = p.sum(axis=tuple(drop)) if drop else p
        # 남은 축 순서는 물리 번호 순 → 논리 순서로 재배열
        remaining = sorted(phys)
        perm = [remaining.index(ph) for ph in phys]
        return np.transpose(p, perm).reshape(-1)


# ======================================================================
# 대표 회로 5개 (8~13장). 각 함수는 (Routed 회로 실행, 읽을 논리 큐비트, 이상적 '정답' 비트열 집합) 을 돌려준다.
# layout: 논리 → 물리. heavy-hex 조각의 차수 3 노드(물리 2)에 '허브' 큐비트를 둔다.
def circ_ghz(noisy):
    dm = DM(5, noisy); R = Routed(dm, {i: i for i in range(5)})
    R.h(0)
    for i in range(4): R.cx(i, i + 1)
    return R, [0, 1, 2, 3, 4], {"00000", "11111"}


def circ_bv(noisy, s="1011"):
    # 8장: 데이터 4 + 보조 1(|−>). 보조를 허브(물리 2)에 두면 s 의 1 이 3개여도 SWAP 없이 인접.
    # s 의 1 은 q0, q2, q3 — 이 셋을 허브의 세 이웃(물리 1, 3, 5)에, 0 인 q1 은 남는 자리(물리 0)에.
    dm = DM(6, noisy); R = Routed(dm, {0: 1, 1: 0, 2: 3, 3: 5, 4: 2})   # 논리 4 = 보조
    for q in range(4): R.h(q)
    R.x(4); R.h(4)
    for q, bit in enumerate(s):
        if bit == "1": R.cx(q, 4)
    for q in range(4): R.h(q)
    return R, [0, 1, 2, 3], {s}


def circ_grover(noisy, target="101", iters=2):
    # 9장: n=3, 표적 1개, k=2 회 반복. CCZ 는 선형 8-CX 판 → 가운데 큐비트(논리 1)가 양쪽과 인접.
    dm = DM(6, noisy); R = Routed(dm, {0: 1, 1: 2, 2: 3})
    for q in range(3): R.h(q)
    for _ in range(iters):
        for q, bit in enumerate(target):          # 오라클: 표적에서 위상 −1
            if bit == "0": R.x(q)
        R.ccz_linear(0, 1, 2)
        for q, bit in enumerate(target):
            if bit == "0": R.x(q)
        for q in range(3): R.h(q); R.x(q)          # 확산: H X CCZ X H
        R.ccz_linear(0, 1, 2)
        for q in range(3): R.x(q); R.h(q)
    return R, [0, 1, 2], {target}


def circ_qpe(noisy):
    # 10장: T 게이트(φ = 1/8)의 3비트 QPE. 표적(논리 3)을 허브에 두고 계수 0,1,2 를 이웃에.
    dm = DM(6, noisy); R = Routed(dm, {0: 1, 1: 3, 2: 5, 3: 2})
    R.x(3)                                        # 고유상태 |1>
    for q in range(3): R.h(q)
    # 계수 레지스터를 "리틀 엔디안"으로 쓴다: q0 가 U^1, q2 가 U^4. 그러면 뒤의 스왑 없는 QFT† 가
    # (스왑 포함 QFT)† 와 같아져 결과를 그대로 big-endian 으로 읽는다 (6장 Draper 덧셈기의 스왑 생략과 같은 요령; 10장은 스왑 포함 판).
    for q, power in zip(range(3), (1, 2, 4)):
        R.cp(q, 3, 2 * np.pi * power / 8)
    # 스왑 없는 QFT† (6장 QFT 순서의 역, 각도 부호 반전)
    R.h(2); R.cp(1, 2, -np.pi / 2); R.h(1); R.cp(0, 2, -np.pi / 4); R.cp(0, 1, -np.pi / 2); R.h(0)
    return R, [0, 1, 2], {"001"}


def circ_hva(noisy):
    # 13장: TFIM n=4 HVA p=1 (γ, β 최적 근사값) 상태의 에너지 측정 — Z 기저와 X 기저 두 회로를 합친다.
    # 에너지 ⟨H⟩ 를 돌려주므로 '정답'은 비트열이 아니라 에너지 값. 여기서는 Z 기저 회로만 돌리고 에너지는 별도 계산.
    gamma, beta = 0.2064, 0.3927  # 13장 HVA p=1 최적값 (mod π): E = -4.6462, E-E0 = 0.113
    def build(basis):
        dm = DM(6, noisy); R = Routed(dm, {i: i for i in range(4)})
        for q in range(4): R.h(q)
        for q in range(3):                        # e^{-iγ(-ZZ)} = CX Rz(-2γ) CX
            R.cx(q, q + 1); R.rz(q + 1, -2 * gamma); R.cx(q, q + 1)
        for q in range(4): R.rx(q, -2 * beta)     # e^{-iβ(-X)} = Rx(-2β)
        if basis == "X":
            for q in range(4): R.h(q)
        return R
    return build


def energy_from_probs(pz, px, n=4):
    zb = ((np.arange(2 ** n)[:, None] >> (n - 1 - np.arange(n))) & 1); zv = 1 - 2 * zb
    e = 0.0
    for i in range(n - 1): e += -1.0 * np.sum(pz * zv[:, i] * zv[:, i + 1])
    for i in range(n): e += -1.0 * np.sum(px * zv[:, i])
    return e


# ---------- 실행 ----------
from contextlib import contextmanager


@contextmanager
def noise_set(**kw):
    """전역 NOISE 를 잠시 바꾼다 (DM 의 기본 인자가 이 dict 객체를 가리킨다)."""
    old = dict(NOISE); NOISE.update(kw)
    try: yield
    finally: NOISE.clear(); NOISE.update(old)


def success(fn):
    R1, qs, good = fn(True); p1 = R1.probs(qs)
    return sum(p1[int(g, 2)] for g in good), R1, p1


print("1) 노이즈 모델 (Heron r2 공개 중앙값):", {k: v for k, v in NOISE.items()})
print(f"   탈분극 p = r·(d+1)/d: p_cz = {p_dep(NOISE['r_cz'], 2, NOISE):.2e}, p_1q = {p_dep(NOISE['r_1q'], 1, NOISE):.2e}")
tg = NOISE["t_cz"]; print(f"   CZ 한 번 동안의 T1 손실 {tg/NOISE['T1']:.1e}, 판독 1.5 µs 동안 {NOISE['t_meas']/NOISE['T1']:.1e} (판독 오류 ε10 에 이미 포함)")

CIRCS = [("GHZ-5 (3장)", circ_ghz), ("BV s=1011 (8장)", circ_bv), ("Grover n=3 k=2 (9장)", circ_grover), ("QPE 3비트 T (10장)", circ_qpe)]
results = {}
print("\n2) 회로별 이상 vs 노이즈 (CZ 수·SWAP, 성공 확률, TVD, 오류 예산)")
eps_avg = (NOISE["eps01"] + NOISE["eps10"]) / 2
for name, fn in CIRCS:
    R0, qs, good = fn(False); p0 = R0.probs(qs)
    succ1, R1, p1 = success(fn)
    n = len(qs)
    succ0 = sum(p0[int(g, 2)] for g in good)
    cz, sq = R1.dm.count["cz"], R1.dm.count["sq"]
    # 오류 예산: 사양서 숫자 그대로 (r_cz, r_1q, 평균 판독 오류) — 암산용
    budget = 1 - (1 - NOISE["r_cz"]) ** cz * (1 - NOISE["r_1q"]) ** sq * (1 - eps_avg) ** n
    T = R1.dm.clock.max()
    results[name] = dict(p0=p0, p1=p1, n=n, succ0=succ0, succ1=succ1, cz=cz, sq=sq, swaps=R1.swaps, budget=budget, T=T)
    print(f"   {name:22s}: CZ {cz:3d} (SWAP {R1.swaps}), 1q {sq:3d}, 길이 {T*1e6:.2f} µs | 성공 이상 {succ0:.3f} → 노이즈 {succ1:.3f} "
          f"| TVD {tvd(p0, p1):.3f} | 예산 {budget:.3f}, 실패−예산 {(1-succ1-budget)*100:+.1f}%p, (이상−노이즈)−이상·예산 {(succ0-succ1-succ0*budget)*100:+.1f}%p")

print("\n   모델 구성 요소별 기여 (성공 확률): 전체 | r→p 변환 끔 | 유휴 열완화 끔 | 판독 오류 끔")
for name, fn in CIRCS:
    row = [results[name]["succ1"]]
    for kw in (dict(rb_convert=False), dict(idle=False), dict(readout=False)):
        with noise_set(**kw): row.append(success(fn)[0])
    print(f"   {name:22s}: " + " | ".join(f"{v:.3f}" for v in row))

# HVA 에너지
build0 = circ_hva(False)
pz0 = build0("Z").probs([0, 1, 2, 3]); px0 = build0("X").probs([0, 1, 2, 3])
E0_exact = -4.758770
e_ideal = energy_from_probs(pz0, px0)
def hva_energy():
    b = circ_hva(True); return energy_from_probs(b("Z").probs([0, 1, 2, 3]), b("X").probs([0, 1, 2, 3]))
e_noisy = hva_energy()
with noise_set(readout=False): e_noro = hva_energy()
R = circ_hva(True)("Z")
print(f"   HVA p=1 TFIM n=4 (13장): CZ {R.dm.count['cz']}, 1q {R.dm.count['sq']} | ⟨H⟩ 이상 {e_ideal:.4f} → 노이즈 {e_noisy:.4f} ({e_noisy - e_ideal:+.3f}); "
      f"판독 오류 제외(TREX 완화 Estimator 와 비교할 값) {e_noro:.4f} ({e_noro - e_ideal:+.3f}); E0 = {E0_exact:.4f}")
results["HVA"] = dict(e_ideal=e_ideal, e_noisy=e_noisy, e_noro=e_noro)

# 샷 잡음을 얹은 '예상 히스토그램' (4000 샷)
print("\n3) 4000 샷 샘플링 (예상 히스토그램 상위 3개)")
shots = 4000
for name, r in results.items():
    if name == "HVA": continue
    counts = rng.multinomial(shots, r["p1"]); top = np.argsort(counts)[::-1][:3]
    print(f"   {name:22s}: " + ", ".join(f"{bits(k, r['n'])}:{counts[k]}" for k in top))

# 오류율 민감도: Grover 성공 확률 vs CZ 오류
print("\n4) CZ 오류율 민감도 (Grover n=3 k=2 성공 확률; 성공 − 1/8 의 감소 비율)")
sens = []
with noise_set(r_cz=0.0): s_zero = success(circ_grover)[0]
print(f"   r_cz = 0: 성공 {s_zero:.3f}")
prev = s_zero
for mult in (0.5, 1, 2, 4, 16):
    with noise_set(r_cz=3.0e-3 * mult): s1 = success(circ_grover)[0]
    sens.append((mult, s1))
    print(f"   r_cz × {mult:>4}: 성공 {s1:.3f}, 성공−1/8 = {s1-0.125:.3f}, (성공−1/8)/(r_cz=0 값) = {(s1-0.125)/(s_zero-0.125):.3f}")
k_per = np.log((sens[2][1] - 0.125) / (sens[1][1] - 0.125)) / (3e-3)
print(f"   r_cz 3e-3 증가당 (성공−1/8) 배율 ≈ {np.exp(k_per*3e-3):.3f}")

# 판독 오류 되돌리기 (연습문제 3): BV 예측 분포에 M^{-1} ⊗4
M = np.array([[1 - NOISE["eps01"], NOISE["eps10"]], [NOISE["eps01"], 1 - NOISE["eps10"]]]); Mi = np.linalg.inv(M)
pb = results["BV s=1011 (8장)"]["p1"].reshape([2] * 4)
for q in range(4): pb = np.moveaxis(np.tensordot(Mi, pb, axes=([1], [q])), 0, q)
pb = pb.reshape(-1)
print(f"\n   (연습문제 3) BV 판독 오류 되돌리기: 성공 {results['BV s=1011 (8장)']['succ1']:.4f} → {pb[int('1011', 2)]:.4f}, 최소 확률 {pb.min():.2e}")

# 연습문제 2: 잘못된 배정(q0 → U^4, q2 → U^1) + 스왑 없는 QFT† + 출력 뒤집기 (이상)
def qpe_wrong():
    dm = DM(4, False); R = Routed(dm, {i: i for i in range(4)})
    R.x(3)
    for q in range(3): R.h(q)
    for q, power in zip(range(3), (4, 2, 1)): R.cp(q, 3, 2 * np.pi * power / 8)
    R.h(2); R.cp(1, 2, -np.pi / 2); R.h(1); R.cp(0, 2, -np.pi / 4); R.cp(0, 1, -np.pi / 2); R.h(0)
    p = R.probs([0, 1, 2]).reshape(2, 2, 2).transpose(2, 1, 0).reshape(-1)   # 출력 뒤집어 읽기
    return p
pw = qpe_wrong()
print("   (연습문제 2) 잘못된 QPE 배정의 분포:", {bits(k, 3): round(float(v), 3) for k, v in enumerate(pw) if v > 1e-3})

# QPU 시간 예산 (Open Plan: 28일 10분)
print("\n5) QPU 시간 예산 추정")
rep_delay = 250e-6     # 기본 반복 지연(샷 사이 큐비트 리셋 대기) — IBM 기본값 수준
Tmax = max(r["T"] for k, r in results.items() if k != "HVA")
per_shot = rep_delay + NOISE["t_meas"] + Tmax
tot = per_shot * shots * 6
print(f"   가장 긴 회로 {Tmax*1e6:.1f} µs → 샷당 ≈ {per_shot*1e6:.0f} µs (반복 지연이 지배) → 4000 샷 × 6 회로 ≈ {tot:.1f} s = 10분의 {tot/600*100:.2f}%")

# ---------- 그림 ----------
fig, axs = plt.subplots(2, 2, figsize=(11, 7))
for ax, name in zip(axs.ravel(), [k for k in results if k != "HVA"]):
    r = results[name]; n = r["n"]; xs = np.arange(2 ** n)
    ax.bar(xs - 0.2, r["p0"], 0.4, label="이상", color="#4C6EF5")
    ax.bar(xs + 0.2, r["p1"], 0.4, label="노이즈 예측", color="#F03E3E")
    ax.set_xticks(xs); ax.set_xticklabels([bits(k, n) for k in xs], rotation=90 if n >= 4 else 0, fontsize=7)
    ax.set_title(f"{name} — CZ {r['cz']}, 성공 {r['succ0']:.2f}→{r['succ1']:.2f}"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig("ch14_results.png", dpi=130)
print("6) saved ch14_results.png")

# ---------- 7) QPU 실측 대조 (ch14_submit.py real 의 기록 — 파일이 없으면 건너뛴다) ----------
import json, pathlib
_here = pathlib.Path(__file__).resolve().parent
_rec = _here / "ch14_record_real.json"
if _rec.exists():
    rec = json.loads(_rec.read_text(encoding="utf-8"))
    # 2026-11-11 실행분은 GHZ 를 measure_all(Qiskit 순서)로 쟀다. bit_order = "book" 인 기록은 모두 책 순서.
    flip = {"ghz5"} if rec.get("bit_order") != "book" else set()
    sj = rec["jobs"]["ghz5"]["job_id"]
    isa_file = next((p for p in sorted((_here.parent / "runs").glob("*_ch14_isa_calib.json"))
                     if sj in json.loads(p.read_text(encoding="utf-8"))), None)
    isa_pubs = json.loads(isa_file.read_text(encoding="utf-8"))[sj]["pubs"] if isa_file else None
    print(f"\n7) QPU 실측 — {rec['backend']} {rec['t']}, {rec['shots']}샷, sampler {sj}, estimator {rec['jobs']['hva'].get('job_id')}, 사용 {rec.get('qpu_seconds')}초")
    keys = ["ghz5", "bv1011", "grover101", "qpe_t"]
    for i, (name, key) in enumerate(zip([c[0] for c in CIRCS], keys)):
        r = results[name]; c = rec["jobs"][key]["counts"]; N = sum(c.values())
        pm = np.zeros(2 ** r["n"])
        for b, v in c.items(): pm[int(b[::-1] if key in flip else b, 2)] += v / N
        good = [k for k in range(2 ** r["n"]) if r["p0"][k] > 0.2]
        s = pm[good].sum()
        cz = rec["jobs"][key].get("two_qubit_gates") or (sum(1 for g in isa_pubs[i] if g[0] == "cz") if isa_pubs else None)
        print(f"   {name:22s}: 실측 성공 {s:.3f} ± {np.sqrt(s * (1 - s) / N):.3f} (예측 {r['succ1']:.3f}, 차 {100 * (s - r['succ1']):+.1f}%p) "
              f"| TVD {tvd(r['p0'], pm):.3f} | 실기 트랜스파일 CZ {cz} (이 장 라우팅 {r['cz']})")
    h = rec["jobs"]["hva"]
    sd = h.get("stds")
    if sd is None and isa_file:
        sd = json.loads(isa_file.read_text(encoding="utf-8"))[h["job_id"]]["estimator"]["stds"]
    print(f"   HVA ⟨H⟩ 실측 {h['evs']:.3f} ± {sd:.3f} (예측 {results['HVA']['e_noisy']:.3f}, 차 {h['evs'] - results['HVA']['e_noisy']:+.3f} = {(h['evs'] - results['HVA']['e_noisy']) / sd:+.1f}σ)")
