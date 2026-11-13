"""21장 사례 2: 단백질 접힘 HP 격자 — numpy 검증 스크립트 (numpy·scipy 만 사용. 실행 시간은 출력 끝에)
난수: 여러 출발점은 default_rng(21), 담금질은 default_rng([N, L]), SPSA 는 default_rng([21, N, p, 샷, 반복, 시드 번호]).
2026-11-12 검토 R13 반영: QAOA 출발점의 부호(γ > 0, β < 0), 여러 출발점과 격자 점검, 각도 전이, SPSA 여러 시드.
2026-11-12 검토 R15 반영: N = 12 의 탐색을 더 돌리고(한 번 다시 출발) 마지막 구간의 개선폭을 찍는다, N = 12 의 SPSA 시드 5 → 20, p = 3 의 샷 훈련 추가, 믹서 계산을 네 큐비트씩 묶어 빠르게.
19장 5절의 사전 등록(21장 행)을 그대로 실행한다.
  가설: 같은 실행 예산에서 QAOA(p ≤ 3)가 담금질과 무작위 추출보다 짧은 TTS 로 최저 에너지 접힘을 찾는다.
  지표: TTS₉₉. 단위는 목적함수 평가 1회(회로 실행 1샷 = 담금질 1스텝 = 무작위 비트열 1개). 정답은 전수 열거.
  크기: 서열 HPHPPHHPHPPH 의 앞 8·10·12개 잔기(14·18·22큐비트). 예산 10⁷ 평가.
  통과 기준: 잔기 12개에서 QAOA 의 TTS(훈련 샷 포함)가 담금질과 무작위 추출 둘 다보다 짧다.
 1) 비용 표: 턴 인코딩의 모든 비트열에 대해 겹침 수와 H–H 접촉 수 (전수 열거와 대조)
 2) 벌점 λ: 최솟값이 유효한 접힘이 되는 가장 작은 λ
 3) 비용 해밀토니안의 Pauli-Z 전개(Walsh–Hadamard)와 QAOA 한 층의 게이트 비용
 4) 베이스라인: 균등 무작위 추출(해석식)과 담금질(벡터화한 4,000회 실행)
 5) QAOA: (가) 정확한 기댓값 ⟨C⟩ 최소화(훈련 비용 0), (나) 정답을 아는 각도 탐색, (전이) 작은 크기의 각도를 큰 크기에 그대로,
    (다) 샷으로 하는 실제 훈련(SPSA, 여러 시드)
 6) 사전 등록 기준에 대한 판정 — 훈련 방식별로
 7) 벌점 λ 에 따른 민감도
 (E) 연습문제 검산
규약: big-endian (비트열의 왼쪽 두 비트가 첫 걸음). 방향 0·1·2·3 = 오른쪽·위·왼쪽·아래.
"""
import json
import math
import os
import sys
import time
import numpy as np
from scipy.optimize import minimize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
matplotlib.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools"))
T_START = time.time()
rng = np.random.default_rng(21)
SEQ = "HPHPPHHPHPPH"
SIZES = (8, 10, 12)
BUDGET = 10 ** 7
DX = np.array([1, 0, -1, 0]); DY = np.array([0, 1, 0, -1])

# =====================================================================
# 1) 비용 표
# =====================================================================
print("=" * 78); print("1) 비용 표 — 턴 인코딩의 모든 비트열"); print("=" * 78)


def tables(seq):
    """모든 비트열 z (걸음마다 2비트, 4^(N−1) 개)에 대해 (겹치는 잔기 쌍의 수, 잇닿지 않은 H–H 격자 이웃 쌍의 수)."""
    N = len(seq); M = 4 ** (N - 1)
    z = np.arange(M, dtype=np.int64)
    X = np.zeros((N, M), np.int16); Y = np.zeros((N, M), np.int16)
    for k in range(N - 1):
        t = (z >> (2 * (N - 2 - k))) & 3                 # big-endian: 첫 걸음이 최상위 두 비트
        X[k + 1] = X[k] + DX[t]; Y[k + 1] = Y[k] + DY[t]
    ov = np.zeros(M, np.int16); con = np.zeros(M, np.int16)
    for i in range(N):
        for j in range(i + 1, N):
            d = np.abs(X[i] - X[j]) + np.abs(Y[i] - Y[j])
            if (j - i) % 2 == 0:
                ov += (d == 0)                           # 정사각 격자에서는 사슬 거리가 짝수인 쌍만 겹칠 수 있다
            elif j - i >= 3 and seq[i] == "H" and seq[j] == "H":
                con += (d == 1)                          # 사슬 거리가 홀수인 쌍만 격자 이웃이 될 수 있다
    return ov, con


def enumerate_hp(seq):
    """19장의 전수 열거(자기 회피 보행). 대조용."""
    L, pos, occ = len(seq), [(0, 0)], {(0, 0): 0}
    stat = {"count": 0, "best": 0, "nbest": 0}

    def energy():
        e = 0
        for i, (x, y) in enumerate(pos):
            if seq[i] == "H":
                for q in ((x + 1, y), (x, y + 1)):
                    j = occ.get(q)
                    if j is not None and seq[j] == "H" and abs(i - j) > 1:
                        e -= 1
        return e

    def grow(k):
        if k == L:
            stat["count"] += 1; e = energy()
            if e < stat["best"]: stat["best"], stat["nbest"] = e, 1
            elif e == stat["best"]: stat["nbest"] += 1
            return
        x, y = pos[-1]
        for q in ((x + 1, y), (x, y + 1), (x - 1, y), (x, y - 1)):
            if q not in occ:
                pos.append(q); occ[q] = k; grow(k + 1); pos.pop(); del occ[q]

    grow(1)
    return stat["count"], stat["best"], stat["nbest"]


TAB = {}
print("    잔기 N   큐비트 n   비트열 4^(N−1)   유효한 접힘   비율      최저 에너지   최저 접힘 수   전수 열거와 일치   (표 만드는 시간)")
for N in SIZES:
    t0 = time.time(); ov, con = tables(SEQ[:N]); dt = time.time() - t0
    valid = ov == 0; emin = -int(con[valid].max()); nopt = int((valid & (con == -emin)).sum())
    ref = enumerate_hp(SEQ[:N])
    TAB[N] = dict(ov=ov, con=con, valid=valid, emin=emin, nopt=nopt, n=2 * (N - 1), M=4 ** (N - 1), opt=valid & (con == -emin))
    ok = (int(valid.sum()), emin, nopt) == ref
    print(f"    {N:4d}      {2*(N-1):3d}     {4**(N-1):10d}     {int(valid.sum()):8d}    {valid.mean():6.2%}       {emin:3d}          {nopt:5d}           {'예' if ok else '아니오'}            {dt:.2f} s")
    assert ok
print("    겹침이 있는(무효) 비트열의 비율: " + ", ".join(f"N = {N}: {1 - TAB[N]['valid'].mean():.2%}" for N in SIZES))

# =====================================================================
# 2) 벌점 λ
# =====================================================================
print("\n" + "=" * 78); print("2) 벌점 λ — C(z) = −(H–H 접촉 수) + λ·(겹치는 쌍의 수)"); print("=" * 78)


def cost(N, lam):
    return (-TAB[N]["con"] + lam * TAB[N]["ov"]).astype(np.float64)


LAM_MIN = {}
for N in SIZES:
    T = TAB[N]; inv = ~T["valid"]
    # 겹침이 있는 비트열 가운데 (접촉 수 − 최저 접힘의 접촉 수)/겹침 수 의 최댓값보다 λ 가 커야 한다
    need = ((T["con"][inv] + T["emin"]) / T["ov"][inv]).max()
    LAM_MIN[N] = float(need)
    rows = []
    for lam in (0.5, 1.0, 2.0, 3.0):
        C = cost(N, lam); zmin = C == C.min()
        rows.append(f"λ = {lam}: 최솟값 {C.min():.1f}, 그 가운데 유효한 접힘 {int((zmin & T['valid']).sum())}/{int(zmin.sum())}")
    print(f"    N = {N}: 겹침 있는 비트열의 접촉 수 최대 {int(T['con'][inv].max())} (유효한 접힘의 최대 {-T['emin']}), 필요한 조건 λ > {need:.3f}.  " + ";  ".join(rows))
LAM = 2.0
print(f"    → 본 실험의 벌점: λ = {LAM} (세 크기 모두에서 최솟값이 정확히 최저 에너지 접힘들이다)")
for N in SIZES:
    C = cost(N, LAM); assert np.array_equal(C == C.min(), TAB[N]["opt"])
    TAB[N]["C"] = C; TAB[N]["sig"] = float(C.std()); TAB[N]["mean"] = float(C.mean())
print("    비용의 분포(λ = 2, 균등 분포 기준): " + "; ".join(f"N = {N}: 평균 {TAB[N]['mean']:.3f}, 표준편차 {TAB[N]['sig']:.3f}, 최댓값 {TAB[N]['C'].max():.0f}" for N in SIZES))

# =====================================================================
# 3) Pauli-Z 전개와 게이트 비용
# =====================================================================
print("\n" + "=" * 78); print("3) 비용 해밀토니안의 Pauli-Z 전개 (Walsh–Hadamard) 와 QAOA 한 층의 게이트 비용"); print("=" * 78)


def fwht(a):
    """빠른 Walsh–Hadamard 변환(정규화 1/2^n). 결과의 색인 s 는 Z 가 걸리는 큐비트의 집합(비트 마스크)."""
    a = a.astype(np.float64).copy(); n = int(math.log2(len(a))); h = 1
    while h < len(a):
        a = a.reshape(-1, 2, h); x = a[:, 0, :].copy(); a[:, 0, :] += a[:, 1, :]; a[:, 1, :] = x - a[:, 1, :]; a = a.reshape(-1); h *= 2
    return a / 2 ** n


def popcount(v):
    v = v.copy(); c = np.zeros(len(v), np.int16)
    while v.any():
        c += (v & 1).astype(np.int16); v >>= 1
    return c


GATE = {}
print("    N    큐비트   0 아닌 항 수 / 2^n        항의 최대 차수   차수 1·2 항   차수 3 이상   한 층의 2큐비트 게이트   한 층의 회전(Rz + Rx)")
for N in SIZES:
    T = TAB[N]; w = fwht(T["C"]); nz = np.abs(w) > 1e-9; nz[0] = False        # 항등 항(전역 위상)은 뺀다
    wt = popcount(np.nonzero(nz)[0].astype(np.int64))
    two_q = int((2 * (wt - 1)).sum()); rot = int(nz.sum()) + T["n"]             # 차수 k 항 = CX 2(k−1) + Rz 1, 믹서 = Rx n
    GATE[N] = dict(terms=int(nz.sum()), maxw=int(wt.max()), low=int((wt <= 2).sum()), high=int((wt >= 3).sum()), two_q=two_q, rot=rot,
                   hist=np.bincount(wt, minlength=T["n"] + 1))
    # 검산: 전개를 되돌리면 비용 표가 나온다
    back = fwht(w) * 2 ** T["n"]
    assert np.abs(back - T["C"]).max() < 1e-8
    print(f"    {N:2d}    {T['n']:3d}     {int(nz.sum()):9d} / {T['M']:8d} ({nz.mean():6.1%})     {int(wt.max()):3d}          {GATE[N]['low']:6d}      {GATE[N]['high']:9d}       {two_q:12d}            {rot:9d}")
print("    차수별 항 수(N = 8): " + ", ".join(f"{k}: {v}" for k, v in enumerate(GATE[8]["hist"]) if v))
# 대칭: 모든 비트를 뒤집으면 대각선(y = −x)에 대한 반사(00 ↔ 11, 01 ↔ 10), 걸음마다 상위 비트만 뒤집으면 180° 회전(00 ↔ 10, 01 ↔ 11).
# 두 대칭 때문에 상위 비트에 걸린 Z 와 하위 비트에 걸린 Z 가 각각 짝수 개인 항만 남는다 → 허용되는 항은 4^(N−2) − 1 개.
for N in SIZES:
    T = TAB[N]; z = np.arange(T["M"], dtype=np.int64); HI = sum(1 << (2 * k + 1) for k in range(N - 1)); LO = HI >> 1
    assert np.array_equal(T["C"], T["C"][z ^ (HI | LO)]) and np.array_equal(T["C"], T["C"][z ^ HI]), "비용이 두 뒤집기에 불변이어야 한다"
    idx = np.nonzero(np.abs(fwht(T["C"])) > 1e-9)[0].astype(np.int64); idx = idx[idx > 0]
    assert not (popcount(idx & HI) % 2).any() and not (popcount(idx & LO) % 2).any()
    GATE[N]["allowed"] = 4 ** (N - 2) - 1
print("    대칭(전체 뒤집기 = 반사, 상위 비트 뒤집기 = 180° 회전)이 허용하는 항 4^(N−2) − 1 개 가운데 살아 있는 것: "
      + ", ".join(f"N = {N}: {GATE[N]['terms']:,} / {GATE[N]['allowed']:,} ({GATE[N]['terms']/GATE[N]['allowed']:.1%})" for N in SIZES))
print("    (세는 법: 차수 k 인 Z 항 하나의 e^{−iγ w Z…Z} = CX 2(k−1)개 + Rz 1개, 전결합. 믹서는 Rx n개. 보조 큐비트로 차수를 낮추는 구성은 세지 않았다.)")

# =====================================================================
# 4) 베이스라인
# =====================================================================
print("\n" + "=" * 78); print("4) 베이스라인 — 균등 무작위 추출과 담금질 (단위: 목적함수 평가 1회)"); print("=" * 78)
LN = math.log(0.01)


def tts_from_p(p, per=1.0):
    """한 번의 시도(평가 per 회)의 성공 확률 p 에서 TTS₉₉. p = 0 이면 ∞, p ≥ 0.99 면 시도 한 번."""
    if p <= 0: return math.inf
    if p >= 0.99: return per
    return per * LN / math.log(1 - p)


RAND = {N: tts_from_p(TAB[N]["nopt"] / TAB[N]["M"]) for N in SIZES}
print("    균등 무작위 비트열: " + "; ".join(f"N = {N}: 한 개의 성공 확률 {TAB[N]['nopt']/TAB[N]['M']:.3e}, TTS₉₉ = {RAND[N]:,.0f}" for N in SIZES))


def anneal(N, steps, runs, t_hot=2.0, t_cold=0.2, seed=0):
    """담금질 runs 회를 한꺼번에. 한 스텝 = 걸음 하나의 방향을 다른 셋 중 하나로 바꾸는 제안 1회 = 평가 1회.
    돌려주는 것: 각 실행이 최저 접힘을 처음 본 스텝(못 봤으면 steps + 1)."""
    T = TAB[N]; C = T["C"]; opt = T["opt"]; g = np.random.default_rng(seed)
    z = g.integers(0, T["M"], runs); cz = C[z]; first = np.where(opt[z], 0, steps + 1)
    for s in range(1, steps + 1):
        temp = t_hot * (t_cold / t_hot) ** (s / steps)
        k = g.integers(0, N - 1, runs); sh = 2 * (N - 2 - k)
        new = (((z >> sh) & 3) + g.integers(1, 4, runs)) & 3
        z2 = (z & ~(np.int64(3) << sh)) | (new << sh); c2 = C[z2]
        acc = (c2 <= cz) | (g.random(runs) < np.exp(-(c2 - cz) / temp))
        z = np.where(acc, z2, z); cz = np.where(acc, c2, cz)
        hit = opt[z] & (first > steps)
        first[hit] = s
    return first


SA = {}; RUNS = 4000
LENGTHS = (30, 100, 300, 1000, 3000, 10000)
print(f"    담금질(온도 2.0 → 0.2 기하 감소, {RUNS}회 실행, 실행 길이별 성공률 p_L 과 TTS₉₉ = L·ln(0.01)/ln(1 − p_L); p_L ≥ 0.99 면 TTS₉₉ = L):")
print("      성공 = L 스텝 안에 최저 접힘을 한 번이라도 지난 것(무작위 추출·QAOA 와 같은 정의). 한 스텝 = 걸음 하나의 방향을 다른 셋 중 하나로 바꾸는 제안 = 평가 1회. 출발 비트열의 평가 1회는 세지 않았다.")
for N in SIZES:
    best = (math.inf, None, None); row = []
    for Lr in LENGTHS:
        first = anneal(N, Lr, RUNS, seed=[N, Lr]); pL = float((first <= Lr).mean()); tts = tts_from_p(pL, Lr)
        row.append((Lr, pL, tts))
        if tts < best[0]: best = (tts, Lr, pL)
    se = math.sqrt(best[2] * (1 - best[2]) / RUNS)                 # 성공률의 표준오차 → TTS 의 범위
    lo, hi = tts_from_p(best[2] + se, best[1]), tts_from_p(best[2] - se, best[1])
    SA[N] = dict(rows=row, tts=best[0], L=best[1], p=best[2], lo=lo, hi=hi)
    print(f"      N = {N:2d}: " + "  ".join(f"L={Lr}: {pL:.3f} → {tts:,.0f}" for Lr, pL, tts in row) + f"   ⇒ 가장 짧은 TTS₉₉ = {best[0]:,.0f} (L = {best[1]}; 성공률의 표준오차로 {lo:,.0f}~{hi:,.0f})")
SA_SENS = []
for th, tc in ((1.0, 0.2), (1.0, 0.1), (3.0, 0.3), (0.5, 0.5)):
    r_ = min(tts_from_p(float((anneal(12, Lr, RUNS, t_hot=th, t_cold=tc, seed=[5000, Lr, int(10 * th), int(100 * tc)]) <= Lr).mean()), Lr) for Lr in (300, 1000, 3000))
    SA_SENS.append((th, tc, r_))
print("    담금질의 온도 일정에 따른 민감도(N = 12, L = 300·1000·3000 중 가장 짧은 TTS₉₉): " + "; ".join(f"{th} → {tc}: {r_:,.0f}" for th, tc, r_ in SA_SENS)
      + f"  (본 실험의 일정 2.0 → 0.2 는 이 장을 쓰면서 정했고 결과를 보고 바꾸지 않았다 — 기록으로 확인할 수는 없다. 일정을 이 문제에 맞추면 {min(r_ for _, _, r_ in SA_SENS):,.0f} 까지 내려간다. 실행 길이 L 은 결과를 보고 가장 좋은 것을 골랐다 — 담금질에 유리한 선택이다)")

# =====================================================================
# 5) QAOA
# =====================================================================
print("\n" + "=" * 78); print("5) QAOA — |+⟩ⁿ 에서 시작, 층마다 e^{−iγ C/σ} 와 e^{−iβ ΣX} (l = 1 이 먼저 걸린다)"); print("=" * 78)
N_EVAL = {N: 0 for N in SIZES}
QUICK = bool(os.environ.get("CH21_QUICK"))                        # 개발용 빠른 점검(평가 횟수·시드를 줄인다). 이때는 회귀 검사를 건너뛴다


def _phase_table(T):
    """비용이 갖는 값은 몇십 가지뿐이다 — 값마다 위상을 한 번 계산해 색인으로 편다(원소별 exp 와 비트 단위로 같다)."""
    if T.get("_phC") is not T["C"] or T.get("_phs") != T["sig"]:    # 표가 바뀌면(7절의 λ 훑기) 다시 만든다. 배열 자체를 붙들어 둔다
        u, inv = np.unique(T["C"], return_inverse=True); T["_ph"] = (u / T["sig"], inv.reshape(-1)); T["_phC"] = T["C"]; T["_phs"] = T["sig"]
    return T["_ph"]


def qaoa_probs(N, params, ref=False):
    """params = (γ₁…γ_p, β₁…β_p). 비용은 표준편차 σ 로 나눈 값을 쓴다(각도를 크기 사이에서 옮길 수 있게).
    믹서 e^{−iβΣX} = Rx(2β)⊗ⁿ 는 큐비트 네 개씩 묶은 16×16 행렬을 곱한다(ref=True 는 두 개씩 einsum 으로 — 대조용)."""
    T = TAB[N]; n = T["n"]; p = len(params) // 2
    psi = np.full(T["M"], 2.0 ** (-n / 2), np.complex128)
    for l in range(p):
        c, s = math.cos(params[p + l]), -1j * math.sin(params[p + l]); r = np.array([[c, s], [s, c]])
        if ref:
            psi = psi * np.exp(-1j * params[l] * (T["C"] / T["sig"])); R = np.kron(r, r)
            for q in range(0, n, 2):                           # 큐비트 두 개씩 Rx⊗Rx
                psi = np.einsum("ab,xbz->xaz", R, psi.reshape(4 ** (q // 2), 4, -1), optimize=True)
            psi = psi.reshape(-1)
            continue
        u, inv = _phase_table(T); psi = psi * np.exp(-1j * params[l] * u)[inv]
        q = 0
        while q < n:                                           # 큐비트 k 개(최대 4)씩: (Rx)^{⊗k} 를 가운데 축에 곱한다
            k = min(4, n - q); R = r
            for _ in range(k - 1): R = np.kron(R, r)
            psi = np.matmul(R, psi.reshape(1 << q, 1 << k, -1)).reshape(-1); q += k
    if not ref: N_EVAL[N] += 1
    return (psi.real ** 2 + psi.imag ** 2)


def summary(N, params):
    T = TAB[N]; pr = qaoa_probs(N, params); p_rand = T["nopt"] / T["M"]
    s = dict(mean=float(pr @ T["C"]), ps=float(pr[T["opt"]].sum()), valid=float(pr[T["valid"]].sum()), x=np.array(params, float))
    s["ratio"] = s["ps"] / p_rand; s["tts"] = tts_from_p(s["ps"])
    s["g_valid"] = s["valid"] / float(T["valid"].mean())                          # 유효한 접힘의 확률이 균등 분포의 몇 배인가
    s["g_opt"] = (s["ps"] / s["valid"]) / (T["nopt"] / int(T["valid"].sum()))     # 유효한 접힘 안에서 최저 접힘의 비율이 몇 배인가
    return s


def ramp(p):
    """출발점: 선형 경사. e^{−iγC}·e^{−iβΣX} 에서 C 를 '줄이려면' γ 와 β 의 부호가 반대여야 한다.
    (|+⟩ⁿ 은 −ΣX 의 바닥상태다. H(s) = −(1−s)ΣX + sC 를 따라가는 한 단계는 e^{−i s C dt}·e^{+i(1−s)ΣX dt} — γ > 0, β < 0.
    13장의 MaxCut 은 비용을 키우는 문제라 둘 다 양수였다. 이 장의 첫 판은 그 부호를 그대로 써서 나쁜 국소 최솟값에 갇혔다 — 검토 R13.)"""
    f = (np.arange(p) + 0.5) / p
    return np.concatenate([0.6 * f, -0.6 * (1 - f)])


NM_LAST = {}


def nm(f, x0, maxfev, fatol=1e-7):
    r = minimize(f, np.array(x0, float), method="Nelder-Mead", options=dict(maxfev=maxfev, xatol=1e-4, fatol=fatol))
    NM_LAST.update(ok=bool(r.success), nfev=int(r.nfev))
    return r.x, float(r.fun)


def nm2(f, x0, budget, fatol=1e-7):
    """N = 12 용(평가가 느리다): 예산의 절반으로 한 번, 찾은 점에서 새 단체(simplex)로 한 번 더.
    돌려주는 것: (x, f, 둘째 탐색에서 f 가 내려간 폭, 둘째 탐색이 수렴 기준을 채웠는가). 폭이 작으면 '거의 다 내려왔다'는 뜻이지 증명은 아니다."""
    x1, f1 = nm(f, x0, budget // 2, fatol); x2, f2 = nm(f, x1, budget - budget // 2, fatol)
    return (x2, f2, f1 - f2, NM_LAST["ok"]) if f2 <= f1 else (x1, f1, 0.0, False)


def stretch(x):
    """p 층의 각도 → p + 1 층의 출발점(마지막 층을 한 번 더, 섞는 각은 절반)."""
    pp = len(x) // 2
    return np.concatenate([np.append(x[:pp], x[pp - 1]), np.append(x[pp:], x[-1] * 0.5)])


# 빠른 믹서 계산(네 큐비트씩)이 대조용 계산(두 큐비트씩 einsum)과 같은지 확인한다
_d = float(np.abs(qaoa_probs(8, ramp(3)) - qaoa_probs(8, ramp(3), ref=True)).max()); N_EVAL[8] = 0
assert _d < 1e-14, "빠른 믹서 계산이 대조용 계산과 다르다"
print(f"    점검: 네 큐비트씩 묶은 믹서 계산과 두 큐비트씩 einsum 계산의 확률 차(N = 8, p = 3): {_d:.1e}")
# 출발점이 '비용을 줄이는 쪽'인지 확인한다(부호가 틀리면 균등 분포의 평균보다 커진다)
START = {(N, p): summary(N, ramp(p)) for N in SIZES for p in (1, 2, 3)}
assert all(START[(N, p)]["mean"] < TAB[N]["mean"] for N in SIZES for p in (1, 2, 3)), "선형 경사 출발점이 비용을 줄이지 않는다 — 부호 확인"
print("    출발점(선형 경사, γ > 0·β < 0)에서 ⟨C⟩ / 균등 분포의 평균, p* 가 무작위의 몇 배: " + "; ".join(
    f"N = {N}, p = {p}: {START[(N, p)]['mean']:.3f} / {TAB[N]['mean']:.3f}, {START[(N, p)]['ratio']:.2f}배" for N in SIZES for p in (2,)))

BUDGET_12 = {1: 60, 2: 160, 3: 400} if not QUICK else {1: 4, 2: 4, 3: 4}    # N = 12 는 평가 하나가 0.3~1초라 횟수를 정해 둔다(절반씩 두 번: nm2)
# (가) 정확한 기댓값 ⟨C⟩ 을 최소화 (샷 잡음 없음, 훈련 비용은 세지 않는다)
IDEAL = {}; MOVED = {}
print("    (가) 정확한 기댓값 ⟨C⟩ 최소화 — 훈련 비용을 세지 않는 경우(찾은 최솟값. N = 10·12 는 지역 탐색)")
print("      N    p    ⟨C⟩ (균등: 평균)    유효한 접힘의 확률   한 샷 성공 확률 p*   무작위의 몇 배   TTS₉₉(훈련 0)     평가 횟수   분해: 유효 × 유효 안의 최저")
for N in SIZES:
    T = TAB[N]
    for p in (1, 2, 3):
        f = lambda x: float(qaoa_probs(N, x) @ T["C"])
        n0 = N_EVAL[N]
        if N == 8:                                            # 여러 출발점: 선형 경사, 한 층 적은 해를 늘린 것, 선형 경사를 흔든 것 6개
            starts = [ramp(p)] + ([stretch(IDEAL[(N, p - 1)]["x"])] if p > 1 else []) + [ramp(p) * (1 + 0.5 * rng.uniform(-1, 1, 2 * p)) for _ in range(6 if not QUICK else 1)]
            bx, bf = min((nm(f, x0, 600) for x0 in starts), key=lambda t: t[1])
        elif N == 10:
            starts = [IDEAL[(8, p)]["x"], ramp(p)]
            bx, bf = min((nm(f, x0, 300 if not QUICK else 60) for x0 in starts), key=lambda t: t[1])
        else:
            MOVED[("C", 10, p)] = summary(N, IDEAL[(10, p)]["x"])     # N = 10 의 각도를 그대로 옮겼을 때(N = 12 에서의 탐색 0회)
            bx, bf, gain_, ok_ = nm2(f, IDEAL[(10, p)]["x"], BUDGET_12[p])
            if MOVED[("C", 10, p)]["mean"] < bf: bx = IDEAL[(10, p)]["x"]
        s = summary(N, bx); s["evals"] = N_EVAL[N] - n0
        if N == 12: s["last_gain"] = gain_; s["nm_ok"] = ok_
        IDEAL[(N, p)] = s
        print(f"      {N:2d}   {p}    {s['mean']:+.4f} ({T['mean']:+.3f})         {s['valid']:.4f}             {s['ps']:.3e}          {s['ratio']:7.2f}        {s['tts']:12,.0f}       {s['evals']:5d}      {s['g_valid']:.1f} × {s['g_opt']:.1f}")
    print(f"           (균등 분포에서 유효한 접힘의 확률 {T['valid'].mean():.4f}, 한 샷 성공 확률 {T['nopt']/T['M']:.3e})")
print("      각도(γ₁…γ_p | β₁…β_p; 비용은 σ 로 나눈 값): " + ";  ".join(f"N = {N}, p = 3: " + " ".join(f"{v:+.3f}" for v in IDEAL[(N, 3)]["x"][:3]) + " | " + " ".join(f"{v:+.3f}" for v in IDEAL[(N, 3)]["x"][3:]) for N in SIZES))
print("      N = 12 의 탐색(예산 " + "·".join(str(BUDGET_12[p]) for p in (1, 2, 3)) + "회를 절반씩 두 번, 둘째는 첫째가 찾은 점에서 새로 출발): 둘째 탐색에서 ⟨C⟩ 가 더 내려간 폭 "
      + ", ".join(f"p = {p}: {IDEAL[(12, p)]['last_gain']:.1e} ({'수렴 기준 충족' if IDEAL[(12, p)]['nm_ok'] else '횟수 상한에서 끊김'})" for p in (1, 2, 3))
      + f".  N = 10 의 각도를 그대로 옮긴 값에서 탐색으로 더 내려간 폭(p = 3): {MOVED[('C', 10, 3)]['mean'] - IDEAL[(12, 3)]['mean']:.4f}")
# p = 1, N = 8 은 한 주기 전체를 격자로 훑어 Nelder–Mead 의 값이 전역 최솟값인지 확인한다(γ 의 주기는 2πσ: 비용이 정수라서. β 의 주기는 π/2)
_sig8 = TAB[8]["sig"]; _gs = np.linspace(-math.pi * _sig8, math.pi * _sig8, 241 if not QUICK else 41); _bs = np.linspace(-math.pi / 4, math.pi / 4, 61 if not QUICK else 21)
_n0 = N_EVAL[8]; GRID8 = min(float(qaoa_probs(8, [g_, b_]) @ TAB[8]["C"]) for g_ in _gs for b_ in _bs); N_EVAL[8] = _n0
print(f"      점검: N = 8, p = 1 의 격자({len(_gs)} × {len(_bs)}, 한 주기 전체) 최솟값 {GRID8:.4f}, Nelder–Mead {IDEAL[(8, 1)]['mean']:.4f}")
assert IDEAL[(8, 1)]["mean"] <= GRID8 + 1e-6, "p = 1 의 Nelder–Mead 가 격자 최솟값보다 나쁘다"

# (나) 정답을 아는 각도 탐색: 같은 회로에서 p* 자체를 키운다. 그 크기의 정답을 알아야 한다.
#      지역 탐색이므로 증명된 상한은 아니고, "⟨C⟩ 를 줄이는 대신 정답을 직접 겨냥했다면 어디까지 가는가"의 참고값이다.
ORACLE = {}
print("    (나) 정답을 아는 각도 탐색 — 같은 깊이의 회로에서 p* 를 직접 키운 값 (지역 탐색; 증명된 상한은 아니다)")
for N in SIZES:
    T = TAB[N]; row = []
    for p in ((1, 2, 3) if N < 12 else (3,)):
        g_ = lambda x: -float(qaoa_probs(N, x)[T["opt"]].sum())
        n0 = N_EVAL[N]
        if N == 8:
            starts = [IDEAL[(N, p)]["x"]] + ([stretch(ORACLE[(N, p - 1)]["x"])] if p > 1 else [])
            bx, bf = min((nm(g_, x0, 800, 1e-10) for x0 in starts), key=lambda t: t[1])
        elif N == 10:
            starts = [IDEAL[(N, p)]["x"], ORACLE[(8, p)]["x"]]
            bx, bf = min((nm(g_, x0, 400 if not QUICK else 60, 1e-10) for x0 in starts), key=lambda t: t[1])
        else:
            MOVED[("P", 10, p)] = summary(N, ORACLE[(10, p)]["x"])
            bx, bf, gain_, ok_ = nm2(g_, ORACLE[(10, p)]["x"], BUDGET_12[p], 1e-12)
            if MOVED[("P", 10, p)]["ps"] > -bf: bx = ORACLE[(10, p)]["x"]
        s = summary(N, bx); s["evals"] = N_EVAL[N] - n0; ORACLE[(N, p)] = s
        if N == 12: s["last_gain"] = gain_; s["nm_ok"] = ok_
        row.append(f"p = {p}: p* = {s['ps']:.3e} (무작위의 {s['ratio']:.1f}배, TTS₉₉ {s['tts']:,.0f}; 유효 {s['g_valid']:.1f} × 유효 안의 최저 {s['g_opt']:.1f}; 평가 {s['evals']}회)")
    print(f"      N = {N:2d}: " + ";  ".join(row))
print("      ⟨C⟩ 로 고른 각도의 p* 는 p* 로 고른 각도의 몇 %인가(p = 3): " + ", ".join(f"N = {N}: {IDEAL[(N, 3)]['ps']/ORACLE[(N, 3)]['ps']:.0%}" for N in SIZES)
      + "; p = 1: " + ", ".join(f"N = {N}: {IDEAL[(N, 1)]['ratio']:.2f}배 대 {ORACLE[(N, 1)]['ratio']:.2f}배" for N in (8, 10)))

print("      p* 로 고른 각도(γ₁…γ_p | β₁…β_p), p = 3: " + ";  ".join(f"N = {N}: " + " ".join(f"{v:+.3f}" for v in ORACLE[(N, 3)]["x"][:3]) + " | " + " ".join(f"{v:+.3f}" for v in ORACLE[(N, 3)]["x"][3:]) for N in SIZES)
      + f"   (N = 12 의 둘째 탐색에서 p* 가 더 오른 폭: 무작위의 {ORACLE[(12, 3)]['last_gain'] / (TAB[12]['nopt'] / TAB[12]['M']):.1f}배 — {'수렴 기준 충족' if ORACLE[(12, 3)]['nm_ok'] else '횟수 상한에서 끊김'})")
# (전이) 작은 크기에서 정한 각도를 N = 12 에 그대로 쓴다 — N = 12 에서의 훈련·탐색 0회
#        ⟨C⟩-최적은 정답을 몰라도 구할 수 있고, p*-최적은 그 작은 크기의 정답(전수 열거)을 써서 구한다. 어느 쪽도 N = 12 의 정답은 쓰지 않는다.
for p in (2, 3):
    MOVED[("C", 8, p)] = summary(12, IDEAL[(8, p)]["x"]); MOVED[("P", 8, p)] = summary(12, ORACLE[(8, p)]["x"])
MOVED[("P", 10, 2)] = summary(12, ORACLE[(10, 2)]["x"])
PRE_EVAL = TAB[8]["M"] + TAB[10]["M"]
print("    (전이) 작은 크기에서 정한 각도를 N = 12 에 그대로 — N = 12 에서의 훈련 샷 0, 탐색 0회")
print("      각도를 정한 방식           p    ⟨C⟩        유효한 접힘의 확률   p*          무작위의 몇 배   TTS₉₉(훈련 0)")
for kind, label in (("C", "⟨C⟩ 최소화"), ("P", "p* 최대화(작은 크기의 정답 사용)")):
    for src in (8, 10):
        for p in (2, 3):
            m_ = MOVED[(kind, src, p)]
            print(f"      N = {src:2d} 에서 {label:24s} {p}   {m_['mean']:+.4f}        {m_['valid']:.4f}          {m_['ps']:.3e}      {m_['ratio']:7.1f}       {m_['tts']:10,.0f}")
print(f"      이 각도를 미리 구하는 데 든 것(양자 장치 불필요): N = 8·10 의 비용 표 {TAB[8]['M']:,} + {TAB[10]['M']:,} = {PRE_EVAL:,}칸(작은 문제의 목적함수 평가)과 상태 벡터 계산. "
      f"N = 12 의 담금질 TTS₉₉ 의 {PRE_EVAL / SA[12]['tts']:.0f}배다 — 다른 서열에도 같은 각도가 통해야 나눠 낼 수 있는 비용이고, 이 장은 서열 하나만 봤다.")
_tp = MOVED[("P", 10, 3)]["tts"]; N_SHORT = [(k, MOVED[k]["tts"]) for k in sorted(MOVED) if k[2] >= 2 and MOVED[k]["tts"] < SA[12]["tts"]]
print(f"      여덟 가지 가운데 담금질({SA[12]['tts']:,.0f})보다 짧은 것: {len(N_SHORT)}가지 — " + ", ".join(f"N = {k[1]} 의 {'p*' if k[0] == 'P' else '⟨C⟩'} 각도 p = {k[2]} ({v:,.0f})" for k, v in N_SHORT)
      + f".  p* 전이 대 ⟨C⟩ 전이(N = 10, p = 3)의 p* 비 {MOVED[('P', 10, 3)]['ps'] / MOVED[('C', 10, 3)]['ps']:.2f}.  표 칸 수를 비용에 더하면 {PRE_EVAL:,} + {_tp:,.0f} = {PRE_EVAL + _tp:,.0f} "
      f"(담금질의 {(PRE_EVAL + _tp) / SA[12]['tts']:.1f}배, 무작위 추출의 1/{RAND[12] / (PRE_EVAL + _tp):.1f})")

# (다) 실제 훈련: 샷으로 추정한 ⟨C⟩ 을 SPSA 로 줄인다. 훈련에 쓴 샷도 평가로 센다.
P_TRAIN = 2; M_SHOTS = 64; K_ITERS = 20
SPSA_A, SPSA_C = 0.15, 0.15                                       # a_k = 0.15/(k+3)^0.602, c_k = 0.15/(k+1)^0.101 (Spall 의 표준 지수)


def train_spsa(N, p, m_shots, k_iters, seed, x0=None):
    """SPSA: 반복마다 회로 2개를 m_shots 번씩 실행해 표본 평균 비용의 기울기를 어림한다. 훈련에 쓴 샷을 모두 센다.
    '훈련 샷도 기회로 센 TTS₉₉' = 훈련 회로에서 나온 표본까지 쳐서, 최저 접힘을 한 번도 못 볼 확률이 1% 가 되는 누적 샷 수."""
    T = TAB[N]; g = np.random.default_rng(seed); x = np.array(ramp(p) if x0 is None else x0, float)
    st = dict(miss_log=0.0, shots=0, seen=0, cum=None)

    def run(xv):
        pr = qaoa_probs(N, xv); ps = float(pr[T["opt"]].sum())
        idx = np.minimum(np.searchsorted(np.cumsum(pr), g.random(m_shots)), T["M"] - 1)
        st["seen"] += int(T["opt"][idx].sum())
        if st["cum"] is None and ps > 0:
            need = (LN - st["miss_log"]) / math.log(1 - ps)       # 이 회로에서 몇 샷을 더 쓰면 누적 '못 볼 확률'이 1% 가 되는가
            if need <= m_shots: st["cum"] = st["shots"] + max(need, 0.0)
        st["miss_log"] += m_shots * math.log(1 - ps); st["shots"] += m_shots
        return float(T["C"][idx].mean())

    for k in range(k_iters):
        ak = SPSA_A / (k + 1 + 2) ** 0.602; ck = SPSA_C / (k + 1) ** 0.101
        delta = g.choice([-1.0, 1.0], 2 * p)
        x = x - ak * (run(x + ck * delta) - run(x - ck * delta)) / (2 * ck) * delta / T["sig"]
    s_ = summary(N, x)                                             # 보고용 평가(샷으로 세지 않는다)
    if st["cum"] is None:                                          # 훈련 중에는 1% 에 못 미쳤다 → 마지막 각도로 계속 추출
        st["cum"] = (st["shots"] + (LN - st["miss_log"]) / math.log(1 - s_["ps"])) if s_["ps"] > 0 else math.inf
    return dict(shots=st["shots"], ps=s_["ps"], ratio=s_["ratio"], tts=st["shots"] + s_["tts"], cum=st["cum"], seen=st["seen"], x=x, mean=s_["mean"], valid=s_["valid"])


def spsa_many(N, p, m, k, seeds, x0=None):
    runs = [train_spsa(N, p, m, k, [21, N, p, m, k, sd] + ([1] if x0 is not None else []), x0) for sd in range(seeds)]
    med = lambda key: float(np.median([r[key] for r in runs]))
    return dict(runs=runs, n=seeds, shots=runs[0]["shots"], ratio=med("ratio"), ratio_lo=min(r["ratio"] for r in runs), ratio_hi=max(r["ratio"] for r in runs),
                tts=med("tts"), tts_lo=min(r["tts"] for r in runs), tts_hi=max(r["tts"] for r in runs), cum=med("cum"), seen=med("seen"), mean=med("mean"))


SEEDS = {8: 20, 10: 20, 12: 20} if not QUICK else {8: 2, 10: 2, 12: 1}
TRAIN = {}
print(f"    (다) 샷으로 하는 훈련 — SPSA {K_ITERS}회 반복, 반복마다 회로 2개 × {M_SHOTS}샷 (p = {P_TRAIN}, 출발점은 선형 경사, "
      f"a_k = {SPSA_A}/(k+3)^0.602, c_k = {SPSA_C}/(k+1)^0.101). 시드마다 결과가 크게 달라 여러 시드의 중앙값과 범위를 적는다")
print("      N   시드 수  훈련 샷   출발점의 p*(무작위의 몇 배)   훈련 뒤 p*: 중앙값 (최소~최대)      TTS₉₉ = 훈련 샷 + 추출: 중앙값 (최소~최대)        훈련 샷도 기회로(중앙값)   훈련 중 최저 접힘을 본 샷(중앙값)")
for N in SIZES:
    TRAIN[N] = t_ = spsa_many(N, P_TRAIN, M_SHOTS, K_ITERS, SEEDS[N])
    print(f"      {N:2d}    {t_['n']:3d}     {t_['shots']:5d}          {START[(N, P_TRAIN)]['ratio']:.2f}배                   {t_['ratio']:6.2f}배 ({t_['ratio_lo']:.2f}~{t_['ratio_hi']:.2f})            "
          f"{t_['tts']:10,.0f} ({t_['tts_lo']:,.0f}~{t_['tts_hi']:,.0f})              {t_['cum']:10,.0f}                 {t_['seen']:.1f}")
_f5 = TRAIN[12]["runs"][:5]; F5 = (float(np.median([r["ratio"] for r in _f5])), float(np.median([r["tts"] for r in _f5])))
print(f"      N = 12 의 처음 5개 시드만 보면(검토 R15 전의 판이 쓴 표본): 중앙값 {F5[0]:.2f}배, TTS₉₉ {F5[1]:,.0f}; {TRAIN[12]['n']}개 가운데 TTS₉₉ 가 이보다 짧은 시드 {sum(r['tts'] < F5[1] for r in TRAIN[12]['runs'])}개")
TRAIN3 = spsa_many(12, 3, M_SHOTS, K_ITERS, 5 if not QUICK else 1)
print(f"      N = 12 에서 p = 3 으로 같은 훈련(시드 {TRAIN3['n']}개; 검토 R15 뒤에 더했다): 출발점 {START[(12, 3)]['ratio']:.2f}배 → 훈련 뒤 중앙값 {TRAIN3['ratio']:.2f}배 ({TRAIN3['ratio_lo']:.2f}~{TRAIN3['ratio_hi']:.2f}), "
      f"TTS₉₉ 중앙값 {TRAIN3['tts']:,.0f} ({TRAIN3['tts_lo']:,.0f}~{TRAIN3['tts_hi']:,.0f}), 훈련 중 최저 접힘을 본 샷(중앙값) {TRAIN3['seen']:.1f}. "
      f"p = 2 의 {TRAIN[12]['n']}개 시드 가운데 TTS₉₉ 가 p = 3 의 최대({TRAIN3['tts_hi']:,.0f})보다 짧은 것 {sum(r['tts'] < TRAIN3['tts_hi'] for r in TRAIN[12]['runs'])}개")
print(f"      훈련이 닿은 곳(N = 12, p = 2): 출발점 {START[(12, 2)]['ratio']:.1f}배 → 훈련 뒤 중앙값 {TRAIN[12]['ratio']:.1f}배 → 같은 회로의 (가) {IDEAL[(12, 2)]['ratio']:.1f}배 "
      f"(훈련 뒤 값은 (가)의 {TRAIN[12]['ratio'] / IDEAL[(12, 2)]['ratio']:.0%}; N = 8: {TRAIN[8]['ratio'] / IDEAL[(8, 2)]['ratio']:.0%}, N = 10: {TRAIN[10]['ratio'] / IDEAL[(10, 2)]['ratio']:.0%})")
TRAIN_SENS = []
for N, sd in ((8, 10), (10, 5)):
    for p_, m_, k_ in ((1, 64, 20), (2, 16, 20), (2, 256, 20), (2, 64, 60), (3, 64, 20)):
        t_ = spsa_many(N, p_, m_, k_, sd if not QUICK else 1); TRAIN_SENS.append((N, p_, m_, k_, t_["shots"], t_["ratio"], t_["tts"], t_["cum"]))
print("      훈련 설정에 따른 민감도(작은 크기에서만. N = 8 은 시드 10개, N = 10 은 5개의 중앙값; p, 샷, 반복 → 훈련 샷, 무작위의 몇 배, TTS₉₉, 훈련 샷도 기회로 센 TTS₉₉):")
for N in (8, 10):
    print(f"        N = {N:2d}: " + ";  ".join(f"p={p_}, {m_}샷, {k_}회 → {sh:,}, {ra:.1f}배, {tt:,.0f}, {cu:,.0f}" for n_, p_, m_, k_, sh, ra, tt, cu in TRAIN_SENS if n_ == N))
BEST_SENS10 = min((r for r in TRAIN_SENS if r[0] == 10), key=lambda r: r[6])
print(f"        N = 10 에서 가장 짧은 설정: p={BEST_SENS10[1]}, {BEST_SENS10[2]}샷, {BEST_SENS10[3]}회 → TTS₉₉ {BEST_SENS10[6]:,.0f} (담금질 {SA[10]['tts']:,.0f} 의 {BEST_SENS10[6]/SA[10]['tts']:.1f}배)")
WARM = spsa_many(10, 3, M_SHOTS, K_ITERS, 5 if not QUICK else 1, x0=IDEAL[(8, 3)]["x"]); WARM0 = summary(10, IDEAL[(8, 3)]["x"])
print(f"      좋은 각도에서 출발해 샷으로 더 훈련하면(N = 10, p = 3, N = 8 의 ⟨C⟩-최적 각도에서 출발, 시드 5개): 출발점 {WARM0['ratio']:.1f}배 → 훈련 뒤 중앙값 {WARM['ratio']:.1f}배 "
      f"({WARM['ratio_lo']:.1f}~{WARM['ratio_hi']:.1f}) — 2,560샷의 훈련은 좋은 각도를 {'더 좋게 만들지 못했다' if WARM['ratio'] <= WARM0['ratio'] else '조금 더 좋게 만들었다'}")

# =====================================================================
# 6) 판정
# =====================================================================
print("\n" + "=" * 78); print("6) 사전 등록 기준에 대한 판정 (19장 5절, 21장 행)"); print("=" * 78)
print("    N     무작위 추출     담금질        QAOA (다) 샷 훈련 p = 2: 중앙값      (가) ⟨C⟩ 최소화, 훈련 0 (가장 좋은 p)   (나) 정답을 아는 각도 탐색, 훈련 0")
BEST_IDEAL = {}; BEST_ORACLE = {}
for N in SIZES:
    bi = min(((IDEAL[(N, p)]["tts"], p) for p in (1, 2, 3))); bo = min(((ORACLE[(N, p)]["tts"], p) for p in (1, 2, 3) if (N, p) in ORACLE))
    BEST_IDEAL[N] = bi; BEST_ORACLE[N] = bo
    print(f"    {N:2d}    {RAND[N]:12,.0f}   {SA[N]['tts']:10,.0f}       {TRAIN[N]['tts']:14,.0f}                        {bi[0]:12,.0f} (p = {bi[1]})                    {bo[0]:12,.0f} (p = {bo[1]})")
N0 = 12; P_RAND0 = TAB[N0]["nopt"] / TAB[N0]["M"]; SA0 = SA[N0]; SA_TUNED = min(r_ for _, _, r_ in SA_SENS)
TC, TP = MOVED[("C", 10, 3)], MOVED[("P", 10, 3)]
PASS_SPSA = bool(TRAIN[N0]["tts"] < SA0["tts"] and TRAIN[N0]["tts"] < RAND[N0])
PASS_SPSA_ANY = any(r["tts"] < SA0["tts"] for r in TRAIN[N0]["runs"] + TRAIN3["runs"])
PASS_TC = bool(TC["tts"] < SA0["tts"] and TC["tts"] < RAND[N0]); PASS_TP = bool(TP["tts"] < SA0["tts"] and TP["tts"] < RAND[N0])
print(f"    기준: 잔기 12개에서 QAOA 의 TTS₉₉(훈련 샷 포함)가 담금질과 무작위 추출 둘 다보다 짧다. 담금질 {SA0['tts']:,.0f} (성공률의 표준오차로 {SA0['lo']:,.0f}~{SA0['hi']:,.0f}), 무작위 {RAND[N0]:,.0f}")
print(f"      방식 1 — N = 12 에서 샷으로 훈련(SPSA, p = 2, 2,560샷; p·샷·반복·이득은 검토 전에 정한 것, 출발점의 부호·시드 수·중앙값으로 적는 것은 검토 뒤): 중앙값 {TRAIN[N0]['tts']:,.0f} ({TRAIN[N0]['tts_lo']:,.0f}~{TRAIN[N0]['tts_hi']:,.0f}, 시드 {TRAIN[N0]['n']}개) "
      f"→ {'통과' if PASS_SPSA else '불통과'}. 담금질의 {TRAIN[N0]['tts']/SA0['tts']:.1f}배({TRAIN[N0]['tts_lo']/SA0['tts']:.1f}~{TRAIN[N0]['tts_hi']/SA0['tts']:.1f}배), 무작위 추출보다는 {RAND[N0]/TRAIN[N0]['tts']:.1f}배 짧다. 담금질보다 짧은 시드: {sum(r['tts'] < SA0['tts'] for r in TRAIN[N0]['runs'])}개")
print(f"        같은 훈련을 p = 3 으로(시드 {TRAIN3['n']}개; 검토 뒤에 더한 것): 중앙값 {TRAIN3['tts']:,.0f} ({TRAIN3['tts_lo']:,.0f}~{TRAIN3['tts_hi']:,.0f}) — 담금질의 {TRAIN3['tts']/SA0['tts']:.1f}배({TRAIN3['tts_lo']/SA0['tts']:.1f}~{TRAIN3['tts_hi']/SA0['tts']:.1f}배). "
      f"p = 2·3 의 {TRAIN[N0]['n'] + TRAIN3['n']}개 시드 가운데 담금질보다 짧은 것: {sum(r['tts'] < SA0['tts'] for r in TRAIN[N0]['runs'] + TRAIN3['runs'])}개")
print(f"      방식 2 — N = 10 의 ⟨C⟩-최적 각도(p = 3)를 그대로(훈련 0샷; 검토 뒤에 더한 사후 분석): {TC['tts']:,.0f} (무작위의 {TC['ratio']:.1f}배) → 글자대로면 {'통과' if PASS_TC else '불통과'}. 담금질의 {TC['tts']/SA0['tts']:.2f}배, 일정을 맞춘 담금질({SA_TUNED:,.0f})의 {TC['tts']/SA_TUNED:.2f}배")
print(f"      방식 3 — N = 10 의 p*-최적 각도(p = 3)를 그대로(훈련 0샷; 사후 분석. N = 10 의 정답을 쓴다): {TP['tts']:,.0f} (무작위의 {TP['ratio']:.1f}배) → 글자대로면 {'통과' if PASS_TP else '불통과'}. 담금질의 {TP['tts']/SA0['tts']:.2f}배, 일정을 맞춘 담금질({SA_TUNED:,.0f})의 {TP['tts']/SA_TUNED:.2f}배")
print(f"      참고 — N = 12 에서 더 다듬은 값(훈련 0 으로 침): (가) {IDEAL[(N0, 3)]['tts']:,.0f} ({IDEAL[(N0, 3)]['ratio']:.1f}배; 담금질과의 차 {IDEAL[(N0, 3)]['tts'] - SA0['tts']:+,.0f} — 담금질 값의 표준오차 범위 반폭 {(SA0['hi'] - SA0['lo']) / 2:,.0f} 의 {abs(IDEAL[(N0, 3)]['tts'] - SA0['tts']) / ((SA0['hi'] - SA0['lo']) / 2):.2f}배; 일정을 맞춘 담금질의 {IDEAL[(N0, 3)]['tts']/SA_TUNED:.2f}배), (나) {ORACLE[(N0, 3)]['tts']:,.0f} ({ORACLE[(N0, 3)]['ratio']:.1f}배; N = 12 의 정답을 쓰므로 실제로는 쓸 수 없다)")
NEED0 = -LN / SA0["tts"]; NEED = -LN / (SA0["tts"] - TRAIN[N0]["shots"])
print(f"      담금질을 이기는 데 필요한 p*: 훈련 0샷이면 무작위의 {NEED0/P_RAND0:.1f}배, 훈련 {TRAIN[N0]['shots']:,}샷을 쓴 뒤면 {NEED/P_RAND0:.1f}배. "
      f"찾은 값: 샷 훈련 p = 2 {TRAIN[N0]['ratio']:.1f}배(최대 {TRAIN[N0]['ratio_hi']:.1f}배), p = 3 {TRAIN3['ratio']:.1f}배(최대 {TRAIN3['ratio_hi']:.1f}배), (가) p = 3 {IDEAL[(N0, 3)]['ratio']:.1f}배, 전이(⟨C⟩) {TC['ratio']:.1f}배, 전이(p*) {TP['ratio']:.1f}배")
print("      작은 크기(훈련 0 으로 침): " + "; ".join(f"N = {N}: (가) {BEST_IDEAL[N][0]:,.0f} 대 담금질 {SA[N]['tts']:,.0f} ({SA[N]['tts']/BEST_IDEAL[N][0]:.1f}배 짧다)" for N in (8, 10))
      + "; 샷 훈련을 넣으면 " + "; ".join(f"N = {N}: {TRAIN[N]['tts']:,.0f} (무작위 {RAND[N]:,.0f}, 담금질 {SA[N]['tts']:,.0f})" for N in (8, 10)))
print(f"      예산 10⁷ 평가 안인가(TTS₉₉ 기준): 무작위 {'예' if RAND[N0] <= BUDGET else '아니오'}, 담금질 {'예' if SA0['tts'] <= BUDGET else '아니오'}, QAOA 방식 1 {'예' if TRAIN[N0]['tts_hi'] <= BUDGET else '아니오'}")
print("    보조 지표(훈련 비용을 뺀 한 샷 성공 확률이 무작위 추출의 몇 배인가, (가)): " + "; ".join(
    f"N = {N}: " + ", ".join(f"p={p} {IDEAL[(N, p)]['ratio']:.1f}배" for p in (1, 2, 3)) for N in SIZES))
print("    담금질이 무작위 추출을 줄이는 배율: " + ", ".join(f"N = {N}: {RAND[N]/SA[N]['tts']:.1f}배" for N in SIZES))
print(f"    단위: 양자 쪽 1샷과 고전 쪽 1스텝을 같은 값으로 쳤다. 가장 단순한 구성으로 N = 12 의 QAOA 한 층이 2큐비트 게이트 {GATE[12]['two_q']:,}개이고 p = 3 이면 {3*GATE[12]['two_q']:,}개다(3절). "
      f"14장의 오류율 r = 3×10⁻³ 을 곱한 rG 는 N = 8 의 한 층({GATE[8]['two_q']:,}개)에서 이미 {3e-3*GATE[8]['two_q']:.0f} 이다(14장 7절의 한계 'CZ ≈ 100개, 오류 예산 0.3' 의 {GATE[8]['two_q']/100:.0f}배; 오류 없이 끝날 확률의 어림 e^(−rG) = {math.exp(-3e-3*GATE[8]['two_q']):.0e}).")
print("    읽는 법: 판정은 훈련 방식에 달려 있다. 등록한 기준은 훈련 방식을 정하지 않았다. 설정을 검토 전에 정한 방식(1)으로는 불통과이고(출발점의 부호와 시드 수는 검토 뒤에 고쳤다), 사후에 더한 방식(2·3)은 담금질과 비기거나 앞선다.")

# =====================================================================
# 7) 벌점 λ 에 따른 민감도 (N = 10)
# =====================================================================
print("\n" + "=" * 78); print("7) 벌점 λ 에 따른 민감도 (N = 10, p = 2, 정확한 기댓값 최소화; 담금질은 같은 λ 의 비용으로 1,000회 실행)"); print("=" * 78)
SWEEP = []
print("     λ      비용의 σ    QAOA: 유효한 접힘의 확률   p*           무작위의 몇 배   분해: 유효 × 유효 안의 최저    담금질 TTS₉₉(가장 좋은 L)")
T10 = TAB[10]; keep = (T10["C"].copy(), T10["sig"], T10["mean"])
for lam in (0.5, 1.0, 2.0, 5.0, 20.0):
    C = cost(10, lam); T10["C"] = C; T10["sig"] = float(C.std()); T10["mean"] = float(C.mean())
    zmin = C == C.min(); okmin = bool(np.array_equal(zmin, T10["opt"]))
    f = lambda x: float(qaoa_probs(10, x) @ T10["C"])
    bx, bf = min((nm(f, x0, 200 if not QUICK else 40) for x0 in (IDEAL[(10, 2)]["x"], ramp(2))), key=lambda t: t[1])
    s = summary(10, bx)
    best = min(tts_from_p(float((anneal(10, Lr, 1000, seed=[77, int(10 * lam), Lr]) <= Lr).mean()), Lr) for Lr in (100, 300, 1000, 3000))
    SWEEP.append((lam, T10["sig"], s["valid"], s["ps"], s["ratio"], best, okmin, s["g_valid"], s["g_opt"]))
    print(f"    {lam:5.1f}    {T10['sig']:7.3f}        {s['valid']:.4f}              {s['ps']:.3e}      {s['ratio']:6.2f}        {s['g_valid']:5.1f} × {s['g_opt']:4.1f}              {best:10,.0f}        {'' if okmin else '(최솟값에 겹침 있는 비트열이 섞이거나 그쪽이 더 낮다)'}")
T10["C"], T10["sig"], T10["mean"] = keep

# =====================================================================
# (E) 연습문제 검산
# =====================================================================
print("\n" + "=" * 78); print("(E) 연습문제 검산"); print("=" * 78)
# 문제 1: N = 4, 서열 HPPH — 손으로 세는 예
ov4, con4 = tables("HPPH"); v4 = ov4 == 0
print(f"    문제 1: 서열 HPPH (N = 4, 6큐비트): 비트열 {len(ov4)}개, 유효한 접힘 {int(v4.sum())}개, 접촉 1개인 접힘 {int((v4 & (con4 == 1)).sum())}개, "
      f"겹침 있는 비트열 {int((~v4).sum())}개, 무작위 추출의 TTS₉₉ = {tts_from_p((v4 & (con4 == 1)).sum() / len(ov4)):.1f}")
E1 = (len(ov4), int(v4.sum()), int((v4 & (con4 == 1)).sum()))
# 문제 2: TTS 공식 — p 와 실행 길이
for p_, L_ in ((0.5, 100), (0.9, 1000), (0.02, 1)):
    print(f"    문제 2: 한 번의 시도가 평가 {L_}회, 성공 확률 {p_}: 시도 수 {LN/math.log(1-p_):.2f}, TTS₉₉ = {tts_from_p(p_, L_):,.2f}")
# 문제 3: 훈련 비용의 손익분기 — 훈련 뒤 p* 가 무작위의 r 배일 때, 훈련 샷이 몇 개 이하여야 무작위 추출을 이기는가
print("    문제 3 (차는 반올림 전 값으로 계산한다):")
for N in SIZES:
    pr = TAB[N]["nopt"] / TAB[N]["M"]; r = IDEAL[(N, 3)]["ratio"]
    print(f"    문제 3: N = {N}: 무작위 TTS₉₉ {RAND[N]:,.0f}, 담금질 {SA[N]['tts']:,.0f}; p = 3 의 (가)에서 찾은 p* 는 {r:.1f}배 → 추출만 {tts_from_p(r*pr):,.0f}; "
          f"훈련에 쓸 수 있는 샷: 무작위 추출을 이기려면 {RAND[N]-tts_from_p(r*pr):,.0f} 이하, 담금질을 이기려면 {SA[N]['tts']-tts_from_p(r*pr):,.0f} 이하"
          + (f" (담금질 값의 표준오차 범위 {SA[N]['lo']:,.0f}~{SA[N]['hi']:,.0f}, 반폭 {(SA[N]['hi'] - SA[N]['lo']) / 2:,.0f}; 범위의 두 끝에서 빼면 {SA[N]['lo']-tts_from_p(r*pr):,.0f} ~ {SA[N]['hi']-tts_from_p(r*pr):,.0f})" if N == 12 else ""))
# 문제 4: 벌점 항의 차수 — N = 8 에서 가장 큰 차수와 그 수
print(f"    문제 4: N = 8 (14큐비트) 의 Z 전개: 0 아닌 항 {GATE[8]['terms']}개, 최대 차수 {GATE[8]['maxw']} (그 차수의 항 {int(GATE[8]['hist'][GATE[8]['maxw']])}개), 차수 2 이하 {GATE[8]['low']}개")
# 문제 5 (N&C 6장 변형): 같은 탐색을 Grover 로 — 표적 16개 / 4^11 개
_M, _K = TAB[12]["M"], TAB[12]["nopt"]; _th = math.asin(math.sqrt(_K / _M)); _k = int(math.floor(math.pi / (4 * _th)))
E5 = (_k, math.sin((2 * _k + 1) * _th) ** 2)
print(f"    문제 5: N = 12 (표적 {_K}개 / {_M:,}개): θ = asin√(K/M) = {_th:.6f}, 반복 수 ⌊π/4θ⌋ = {_k}, 성공 확률 sin²((2k+1)θ) = {E5[1]:.6f} → 오라클 호출 {_k}회 "
      f"(담금질 TTS₉₉ {SA[12]['tts']:,.0f} 의 1/{SA[12]['tts']/_k:.0f}, 무작위 {RAND[12]:,.0f} 의 1/{RAND[12]/_k:,.0f}). 오라클 한 번은 비용 함수의 가역 계산이고 깊이는 오류 정정이 있어야 버틴다.")

# =====================================================================
# 스니펫용 사례 (N = 6, 10큐비트) — 6절의 SDK 스니펫이 읽는다
# =====================================================================
ov6, con6 = tables(SEQ[:6]); C6 = (-con6 + LAM * ov6).astype(np.float64); v6 = ov6 == 0; opt6 = v6 & (con6 == con6[v6].max())
TAB[6] = dict(ov=ov6, con=con6, valid=v6, emin=-int(con6[v6].max()), nopt=int(opt6.sum()), n=10, M=1024, opt=opt6, C=C6, sig=float(C6.std()), mean=float(C6.mean()))
N_EVAL[6] = 0
res6 = minimize(lambda x: float(qaoa_probs(6, x) @ C6), ramp(1), method="Nelder-Mead", options=dict(maxfev=300, xatol=1e-4, fatol=1e-7))
s6 = summary(6, res6.x); w6 = fwht(C6); nz6 = np.nonzero(np.abs(w6) > 1e-9)[0]
with open(os.path.join(HERE, "ch21_case.json"), "w", encoding="utf-8") as fh:
    json.dump({"note": "ch21_hp_folding.py 가 만든 스니펫용 사례. 서열 HPHPPH (N = 6, 10큐비트), 벌점 2, p = 1. 고치지 말고 스크립트를 다시 돌린다.",
               "seq": SEQ[:6], "n_qubits": 10, "lambda": LAM, "sigma": TAB[6]["sig"],
               "cost": [float(v) for v in C6], "optimal": [int(i) for i in np.nonzero(opt6)[0]],
               "gamma": float(res6.x[0]), "beta": float(res6.x[1]),
               "pauli_z_terms": [[int(i), float(f"{w6[i]:.17g}")] for i in nz6],
               "expect_mean_cost": s6["mean"], "expect_p_success": s6["ps"], "expect_p_valid": s6["valid"]}, fh, ensure_ascii=False)
import gate_cost as gc                                           # tools/gate_cost.py — 게이트 리스트가 표를 위상으로 건 회로와 같은지 수치로 확인
_g = [gc.G("h", q) for q in range(10)]
for i_ in nz6:
    if i_ == 0: continue
    qs_ = [q for q in range(10) if (int(i_) >> (9 - q)) & 1]                       # 빅 엔디안: 큐비트 0 이 최상위 비트
    lad = [gc.G("cx", a_, b_) for a_, b_ in zip(qs_, qs_[1:])]
    _g += lad + [gc.G("rz", qs_[-1], angle=2 * res6.x[0] * w6[i_] / TAB[6]["sig"])] + lad[::-1]
_g += [gc.G("rx", q, angle=2 * res6.x[1]) for q in range(10)]
_psi = np.zeros(1024, complex); _psi[0] = 1; _pr = np.abs(gc.apply(_g, _psi, 10)) ** 2
CASE_GATES = gc.cost(_g); CASE_DIFF = float(np.abs(_pr - qaoa_probs(6, res6.x)).max())
assert CASE_DIFF < 1e-12 and CASE_GATES["two_q"] == int(sum(2 * (bin(int(i_)).count("1") - 1) for i_ in nz6 if i_))
print(f"\n    스니펫용 사례(ch21_case.json): 서열 {SEQ[:6]}, 10큐비트, 유효한 접힘 {int(v6.sum())}개, 최저 에너지 {TAB[6]['emin']} ({TAB[6]['nopt']}개), Z 항 {len(nz6)}개(항등 포함); "
      f"p = 1 최적 각도 γ = {res6.x[0]:.4f}, β = {res6.x[1]:.4f} → ⟨C⟩ = {s6['mean']:.4f}, 유효 {s6['valid']:.4f}, p* = {s6['ps']:.4f} (무작위 {TAB[6]['nopt']/1024:.4f})")
print(f"    게이트로 푼 한 층(H 10 + Z 항 {len(nz6) - 1}개의 CX 사다리·Rz + Rx 10): {gc.fmt(_g)} — 표를 위상으로 건 계산과 확률의 최대 차 {CASE_DIFF:.1e}")

# =====================================================================
# 그림
# =====================================================================
fig, axs = plt.subplots(2, 3, figsize=(16, 9.6))
ax = axs[0, 0]                                                    # (1) 최저 접힘 하나 (N = 12)
z0 = int(np.nonzero(TAB[12]["opt"])[0][0]); xs, ys = [0], [0]
for k in range(11):
    t = (z0 >> (2 * (10 - k))) & 3; xs.append(xs[-1] + int(DX[t])); ys.append(ys[-1] + int(DY[t]))
ax.plot(xs, ys, "-", color="0.4", lw=2, zorder=1)
for i, (x_, y_) in enumerate(zip(xs, ys)):
    ax.scatter([x_], [y_], s=320, color="#E8590C" if SEQ[i] == "H" else "#74C0FC", edgecolor="k", zorder=2)
    ax.text(x_, y_, str(i), ha="center", va="center", fontsize=9, zorder=3)
for i in range(12):
    for j in range(i + 3, 12):
        if SEQ[i] == "H" and SEQ[j] == "H" and abs(xs[i] - xs[j]) + abs(ys[i] - ys[j]) == 1:
            ax.plot([xs[i], xs[j]], [ys[i], ys[j]], ":", color="#E8590C", lw=2.5, zorder=0)
ax.set_aspect("equal"); ax.margins(0.25); ax.set_xticks([]); ax.set_yticks([])
ax.set_title(f"(1) 최저 에너지 접힘 하나 (N = 12, E = {TAB[12]['emin']}; 주황 = H, 점선 = 접촉)")
ax = axs[0, 1]                                                    # (2) 비용 분포
vals, cnt = np.unique(TAB[12]["C"], return_counts=True)
ax.bar(vals, cnt / TAB[12]["M"], width=0.8, color="#868E96", label="균등 분포")
pr12 = qaoa_probs(12, IDEAL[(12, 3)]["x"]); N_EVAL[12] -= 1
ax.bar(vals, [pr12[TAB[12]["C"] == v].sum() for v in vals], width=0.4, color="#1971C2", label="QAOA p = 3 (⟨C⟩ 최소화)")
ax.set_yscale("log"); ax.set_xlim(-6, 25); ax.set_xlabel("비용 C = −접촉 + 2·겹침"); ax.set_ylabel("확률 (로그)"); ax.legend(fontsize=8)
ax.set_title(f"(2) N = 12 의 비용 분포: 균등 분포와 QAOA (C ≤ 25 만 그림; 최댓값 {int(TAB[12]['C'].max())})")
ax = axs[0, 2]                                                    # (3) 담금질 성공률
for N, mk in zip(SIZES, ("o-", "s-", "^-")):
    ax.semilogx([r[0] for r in SA[N]["rows"]], [r[1] for r in SA[N]["rows"]], mk, label=f"N = {N}")
ax.set_xlabel("담금질 한 번의 길이 L (스텝 = 평가)"); ax.set_ylabel("최저 접힘을 본 실행의 비율"); ax.legend(fontsize=8)
ax.set_title("(3) 담금질의 성공률")
ax = axs[1, 0]                                                    # (4) TTS 비교
xx = np.arange(len(SIZES)); wbar = 0.2
ax.bar(xx - 1.5 * wbar, [RAND[N] for N in SIZES], wbar, label="무작위 추출", color="#868E96")
ax.bar(xx - 0.5 * wbar, [SA[N]["tts"] for N in SIZES], wbar, label="담금질", color="#2F9E44")
ax.bar(xx + 0.5 * wbar, [TRAIN[N]["tts"] for N in SIZES], wbar, label="QAOA 샷 훈련 p = 2 (중앙값, 막대 = 최소~최대)", color="#1971C2",
       yerr=[[TRAIN[N]["tts"] - TRAIN[N]["tts_lo"] for N in SIZES], [TRAIN[N]["tts_hi"] - TRAIN[N]["tts"] for N in SIZES]], capsize=3)
ax.bar(xx + 1.5 * wbar, [BEST_IDEAL[N][0] for N in SIZES], wbar, label="QAOA ⟨C⟩ 최소화 (훈련 0)", color="#74C0FC")
ax.plot([2 + 1.5 * wbar], [TC["tts"]], "D", color="#E8590C", ms=7, label="N = 12: N = 10 의 ⟨C⟩-최적 각도 그대로")
ax.plot([2 + 1.5 * wbar], [TP["tts"]], "*", color="#E03131", ms=12, label="N = 12: N = 10 의 p*-최적 각도 그대로")
ax.axhline(BUDGET, color="k", ls=":", lw=1); ax.text(-0.45, BUDGET * 1.3, "예산 10⁷", fontsize=8)
ax.set_yscale("log"); ax.set_xticks(xx); ax.set_xticklabels([f"N = {N}" for N in SIZES]); ax.set_ylabel("TTS₉₉ (평가 횟수, 로그)"); ax.legend(fontsize=6.5)
ax.set_title("(4) TTS₉₉ 비교 — 낮을수록 좋다")
ax = axs[1, 1]                                                    # (5) p* 배율
for N, mk in zip(SIZES, ("o-", "s-", "^-")):
    ax.plot([1, 2, 3], [IDEAL[(N, p)]["ratio"] for p in (1, 2, 3)], mk, label=f"N = {N}: ⟨C⟩ 최소화")
    ps_ = [p for p in (1, 2, 3) if (N, p) in ORACLE]
    ax.plot(ps_, [ORACLE[(N, p)]["ratio"] for p in ps_], mk.replace("-", "--"), color="0.5", label=f"N = {N}: 정답을 아는 각도 탐색")
ax.set_yscale("log"); ax.set_xticks([1, 2, 3]); ax.set_xlabel("QAOA 층 수 p"); ax.set_ylabel("p* / 무작위 추출 (로그)"); ax.legend(fontsize=7)
ax.set_title("(5) 한 샷 성공 확률이 무작위의 몇 배인가")
ax = axs[1, 2]                                                    # (6) 항의 차수 분포
for N, mk in zip(SIZES, ("o-", "s-", "^-")):
    h = GATE[N]["hist"]; ks = np.nonzero(h)[0]
    ax.semilogy(ks, h[ks], mk, label=f"N = {N} (항 {GATE[N]['terms']:,}개)")
ax.set_xlabel("Z 항의 차수(걸리는 큐비트 수)"); ax.set_ylabel("항의 수 (로그)"); ax.legend(fontsize=8)
ax.set_xticks(range(2, 21, 2)); ax.set_title("(6) 비용 해밀토니안의 Z 전개: 차수별 항 수")
plt.tight_layout(); plt.savefig(os.path.join(HERE, "ch21_results.png"), dpi=110)

# =====================================================================
# 회귀 검사 — 본문이 인용하는 값 (어긋난 것을 모두 모아 한꺼번에 알린다)
# =====================================================================
_bad = []


def chk(cond, label):
    if not cond: _bad.append(label)


def near(actual, expected, label, rtol=0.0, atol=0.0):
    a = np.atleast_1d(np.array(actual, float)); e = np.atleast_1d(np.array(expected, float))
    if a.shape != e.shape or not np.allclose(a, e, rtol=rtol, atol=atol): _bad.append(f"{label}: 실제 {[float(f'{v:.6g}') for v in a]} / 적어 둔 값 {e.tolist()}")


chk([(int(TAB[N]["valid"].sum()), TAB[N]["emin"], TAB[N]["nopt"]) for N in SIZES] == [(2172, -2, 24), (16268, -2, 376), (120292, -5, 16)], "5절 (1)")
chk([LAM_MIN[N] for N in SIZES] == [0.0, 1.0, 1.0], "5절 (2): 필요한 λ")
chk([(GATE[N]["terms"], GATE[N]["maxw"], GATE[N]["two_q"]) for N in SIZES] == [(1782, 12, 18388), (28671, 16, 409602), (1048575, 20, 20971522)], "5절 (3): Z 전개")
chk([round(RAND[N]) for N in SIZES] == [3141, 3208, 1207215], "5절 (4): 무작위 추출")
if not QUICK:
    # 구조: 값이 바뀌어도 성립해야 하는 것
    near([IDEAL[(8, p)]["mean"] for p in (1, 2, 3)], [1.2965, 0.4740, 0.1031], "5절 (5가): N = 8 의 ⟨C⟩ 최솟값(검토 R13·R15 가 여러 출발점·격자로 확인한 값)", atol=2e-3)
    chk(all(IDEAL[(N, p)]["mean"] < IDEAL[(N, p - 1)]["mean"] for N in SIZES for p in (2, 3)), "층을 늘리면 ⟨C⟩ 가 내려가야 한다")
    chk(all(ORACLE[k]["ps"] >= IDEAL[k]["ps"] - 1e-12 for k in ORACLE), "(나)는 (가)보다 p* 가 작을 수 없다")
    chk(all(IDEAL[(N, 1)]["ratio"] > 1 for N in SIZES), "p = 1 의 ⟨C⟩-최적도 무작위보다 낫다")
    chk(IDEAL[(12, 3)]["mean"] <= 0.3865, "5절 (5가): N = 12, p = 3 의 ⟨C⟩ 가 검토 R13·R15 가 더 돌려 얻은 값(0.3863)까지 내려가지 않았다")
    chk(PASS_SPSA is False and not PASS_SPSA_ANY and max(TRAIN[12]["tts_hi"], TRAIN3["tts_hi"]) < RAND[12], "5절 (6): 방식 1 은 p = 2·3 의 모든 시드에서 담금질보다 길고 무작위보다 짧다")
    chk(PASS_TP is True and TP["tts"] < SA_TUNED and 0.8 < TC["tts"] / SA0["tts"] < 1.3, "5절 (6): 방식 3 은 담금질보다 짧고 방식 2 는 담금질과 비슷하다")
    chk(sorted(k for k, _ in N_SHORT) == [("P", 8, 3), ("P", 10, 3)], "5절 (5전이): 담금질보다 짧은 전이는 p* 로 정한 p = 3 의 둘뿐")
    chk(abs(IDEAL[(12, 3)]["tts"] - SA0["tts"]) < (SA0["hi"] - SA0["lo"]) / 2, "5절 (6): (가) p = 3 과 담금질의 차가 표본 오차 안")
    chk(WARM["ratio"] < WARM0["ratio"], "5절 (5다): 좋은 각도에서 출발한 샷 훈련은 중앙값을 올리지 못했다")
    chk(E5[0] == 402 and E5[1] > 0.99, "8절 문제 5")
    # 고정값: 스크립트를 바꾸지 않았는데 달라지면 무언가 깨진 것이다
    chk([round(SA[N]["tts"]) for N in SIZES] == [1163, 652, 15601] and [SA[N]["L"] for N in SIZES] == [300, 300, 300], "5절 (4): 담금질")
    near([SA[N]["p"] for N in SIZES], [0.695, 0.880, 0.085], "5절 (4): 담금질의 성공률 p_L", atol=6e-4)
    chk([round(r_) for _, _, r_ in SA_SENS] == [12989, 12908, 20281, 9894], "5절 (4): 담금질 일정 민감도")
    near([IDEAL[(N, p)]["ratio"] for N in (8, 10) for p in (1, 2, 3)], [4.78, 9.71, 14.48, 6.36, 13.56, 20.53], "5절 (5가): 무작위의 몇 배, N = 8·10", rtol=0.002)
    near([IDEAL[(12, p)]["ratio"] for p in (1, 2, 3)], [10.5514, 35.5271, 78.7821], "5절 (5가): 무작위의 몇 배, N = 12", rtol=0.002)
    near([ORACLE[k]["ratio"] for k in ((8, 3), (10, 3), (12, 3))], [22.8634, 29.9104, 262.766], "5절 (5나)", rtol=0.003)
    near([MOVED[(k, src, p)]["tts"] for k in ("C", "P") for src in (8, 10) for p in (2, 3)], [41130, 19453, 36154, 17114, 18869, 7458, 17502, 5963], "5절 (5전이)", rtol=0.002)
    near([TRAIN[N]["tts"] for N in SIZES] + [TRAIN[12]["tts_lo"], TRAIN[12]["tts_hi"]], [3305, 3187, 169302, 62647, 299062], "5절 (5다): TTS₉₉", rtol=0.001)
    near([TRAIN[N]["ratio"] for N in SIZES] + [F5[0], F5[1]], [4.207, 5.1025, 7.2402, 13.7543, 90328], "5절 (5다): 무작위의 몇 배, 처음 5개 시드", rtol=0.002)
    near([TRAIN3["ratio"], TRAIN3["tts"], TRAIN3["tts_lo"], TRAIN3["tts_hi"]], [13.2887, 93403, 63591, 111906], "5절 (5다): N = 12, p = 3", rtol=0.002)
    near([r[6] for r in TRAIN_SENS], [3337, 1399, 10888, 8187, 3010, 3333, 1292, 10776, 8074, 2920], "5절 (5다): 훈련 설정 민감도", rtol=0.001)
    near([WARM0["ratio"], WARM["ratio"], WARM["ratio_lo"], WARM["ratio_hi"]], [19.3, 16.2, 12.9, 19.6], "5절 (5다): 좋은 각도에서 출발", atol=0.06)
    near([r[4] for r in SWEEP], [17.88, 18.28, 13.56, 9.77, 8.00], "5절 (7): λ 민감도", rtol=0.002)
    near([TRAIN[12]["ratio_lo"], TRAIN[12]["ratio_hi"], TRAIN3["ratio_lo"], TRAIN3["ratio_hi"], START[(12, 2)]["ratio"], START[(12, 3)]["ratio"]], [4.07, 20.09, 11.04, 19.78, 3.42, 4.61], "5절 (5다): N = 12 의 범위와 출발점", atol=0.006)
    near([TRAIN[12]["cum"]], [167930], "5절 (5다): 훈련 샷도 기회로 센 TTS₉₉", rtol=0.001)
    near([TRAIN[N]["ratio"] / IDEAL[(N, 2)]["ratio"] for N in SIZES], [0.433, 0.376, 0.204], "5절 (5다): 훈련 뒤 p* 가 (가)의 몇 %인가", atol=0.002)
    chk(sum(r["tts"] < F5[1] for r in TRAIN[12]["runs"]) == 3 and sum(r["tts"] < TRAIN3["tts_hi"] for r in TRAIN[12]["runs"]) == 6, "5절 (5다): 처음 5개의 중앙값·p = 3 의 최대보다 짧은 p = 2 시드 수")
    chk(all(IDEAL[(12, p)]["last_gain"] < 1e-5 for p in (1, 2, 3)), "5절 (5가): N = 12 의 둘째 탐색에서 ⟨C⟩ 가 1e-5 넘게 더 내려갔다 — 예산이 모자라다")
chk(E1 == (64, 36, 8), "8절 문제 1")
chk(abs(s6["ps"] - 0.0313) < 5e-5 and abs(s6["mean"] - 0.2904) < 5e-5 and len(nz6) == 256 and CASE_GATES["two_q"] == 2050, "4·6절: 스니펫용 사례")

print(f"\n    QAOA 회로 평가 횟수(상태 벡터): " + ", ".join(f"N = {N}: {N_EVAL[N]}" for N in SIZES))
print(f"   (실행 시간: 총 {time.time()-T_START:.0f} s)")
print("    회귀 검사: " + ("모두 통과" if not _bad else "불일치 %d건" % len(_bad)) + (" (CH21_QUICK: 고정값 검사는 건너뜀)" if QUICK else ""))
assert not _bad, "회귀 검사 불일치 %d건:\n  " % len(_bad) + "\n  ".join(_bad)
