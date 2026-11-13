"""ch20_linear_systems.py — 20장 「사례 1: 선형계와 블록 인코딩」 검증 스크립트
QSVT 로 A x = b 를 푼다: 역행렬 다항식 → QSP 위상 → 블록 인코딩 위의 QSVT → 성공 확률·총 호출 수 → 고전 CG 와 비교.
19장의 사전 등록(20장 행)을 그대로 따른다: I + γS 족, κ = 2, 4, 8, 16, 32, N = 16~256, 해의 오차 10⁻³, 무작위 b 10개.
실행: python ch20_linear_systems.py   (numpy·scipy 만 사용, 약 25초, seed = 20)
큐비트 순서: 빅 엔디안. 블록 인코딩의 보조 큐비트가 앞(상위)이다.
"""
import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
from numpy.polynomial import chebyshev as C
from scipy.stats import binom

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools"))
import gate_cost as gc  # noqa: E402

t_start = time.time()
EPS = 1e-3                     # 사전 등록: 해의 오차 10⁻³
KAPPAS = (2, 4, 8, 16, 32)
SIZES = (16, 64, 256)
rng = np.random.default_rng(20)

# =====================================================================
# 1) 역행렬 다항식 — 1/x 를 [1/κ, 1] 에서 흉내 내는 홀수 다항식
# =====================================================================
print("=" * 78); print("1) 역행렬 다항식 P(x) ≈ 1/(2κx)  (|x| ∈ [1/κ, 1])"); print("=" * 78)


def inv_cheb(kappa, eps=EPS):
    """g_b(x) = (1 − (1 − x²)^b)/x 의 체비쇼프 계수(Childs–Kothari–Somma 2017, 보조정리 17~19).
    g_b = 4 Σ_{j=0}^{b−1} (−1)^j [2^{−2b} Σ_{i=j+1}^{b} C(2b, b+i)] T_{2j+1}.  b = ⌈κ² ln(2/ε)⌉ 이면 |x| ≥ 1/κ 에서 (1 − x²)^b ≤ ε/2.
    돌려주는 것: b, 그리고 P = g_b/(2κ) 의 T_{2j+1} 계수(끝없이 길어서 충분히 긴 앞부분만)."""
    b = int(np.ceil(kappa ** 2 * np.log(2 / eps)))
    j = np.arange(0, min(b, int(6 * np.sqrt(b * np.log(8 * b / eps))) + 10))
    tail = binom.sf(b + j, 2 * b, 0.5)                       # P(X ≥ b + j + 1), X ~ Bin(2b, 1/2)
    return b, 4 * (-1.0) ** j * tail / (2 * kappa)


def grid(kappa):
    xs = np.concatenate([np.linspace(1 / kappa, 1, 4001), (1 / kappa) * (1 + np.logspace(-6, 0, 200))])
    return xs[xs <= 1]


def cheb_full(c_odd):
    full = np.zeros(2 * len(c_odd)); full[1::2] = c_odd
    return full


def rel_err(kappa, c_odd):
    """sup_{x ∈ [1/κ,1]} |2κ·x·P(x) − 1| — 1/x 에 대한 상대 오차(격자 4,201점)."""
    xs = grid(kappa)
    return np.abs(2 * kappa * xs * C.chebval(xs, cheb_full(c_odd)) - 1).max()


def min_degree(kappa, eps=EPS):
    """상대 오차가 ε 이하인 가장 낮은 홀수 차수 d 와 그 계수. 꼬리 합 상한에서 출발해 격자 오차로 줄인다."""
    b, cj = inv_cheb(kappa, eps)
    a = np.abs(cj); tails = np.cumsum(a[::-1])[::-1]
    J = int(np.argmax(2 * kappa * tails <= eps / 2))           # 엄밀한 충분조건: Σ_{k>d}|c_k|·2κ ≤ ε/2
    d_rig = 2 * J - 1
    while J > 1 and rel_err(kappa, cj[:J - 1]) <= eps:
        J -= 1
    return b, 2 * J - 1, cj[:J].copy(), d_rig


POLY = {}
print("    κ        b    차수 d   (꼬리 합 상한 d)  (CKS 끊는 자리)   상대 오차     max|P|    ‖c‖₁    d/κ")
for kap in KAPPAS:
    b, d, c, d_rig = min_degree(kap)
    d_cks = 2 * int(np.ceil(np.sqrt(b * np.log(4 * b / EPS)))) + 1
    xx = np.linspace(-1, 1, 20001); pmax = np.abs(C.chebval(xx, cheb_full(c))).max()
    POLY[kap] = dict(b=b, d=d, c=c, d_rig=d_rig, d_cks=d_cks, err=rel_err(kap, c), pmax=pmax, l1=np.abs(c).sum())
    print(f"   {kap:3d}  {b:7d}    {d:5d}       {d_rig:5d}           {d_cks:5d}      {POLY[kap]['err']:.2e}    {pmax:.3f}    {POLY[kap]['l1']:.2f}   {d/kap:.2f}")
    assert pmax < 1 and POLY[kap]["err"] <= EPS
deg_slope = np.polyfit(np.log(KAPPAS), np.log([POLY[k]["d"] for k in KAPPAS]), 1)[0]
print(f"    차수의 κ 지수(κ = 2~32 로그-로그 맞춤): {deg_slope:.3f}   → 사전 등록 기준 1 ± 0.2")

print("    정밀도 의존(κ = 8):", end="")
EPS_ROWS = []
for e in (1e-2, 1e-3, 1e-4, 1e-6):
    _, d_e, c_e, _ = min_degree(8, e)
    EPS_ROWS.append((e, d_e, np.abs(C.chebval(np.linspace(-1, 1, 20001), cheb_full(c_e))).max()))
    print(f"  ε = {e:.0e} → d = {d_e} (max|P| {EPS_ROWS[-1][2]:.3f})", end="")
print("\n    읽는 법: 자릿수를 하나 늘릴 때 차수는 곱이 아니라 덧셈으로 는다(√(ln) 꼴). max|P| 가 1 을 넘으면 그대로는 QSVT 에 못 쓴다 — 전체를 줄여야 하고 그만큼 성공 확률이 준다.")


def law_degree(kappa):
    """κ 가 커서 직접 구하지 않는 곳에 쓰는 어림: d ≈ (d/κ at κ = 32)·κ."""
    return POLY[32]["d"] / 32 * kappa


# =====================================================================
# 2) QSP 위상 — 대칭 QSP, 표준(W_x) 규약, Newton
# =====================================================================
print("\n" + "=" * 78); print("2) QSP 위상 — 다항식 계수에서 위상 d + 1 개를 구한다"); print("=" * 78)


def qsp_wx(phis, x):
    """W_x 규약: U_Φ(x) = e^{iφ0 Z} Π_{j=1..d} [W(x) e^{iφj Z}],  W(x) = [[x, i√(1−x²)], [i√(1−x²), x]].  ⟨0|U_Φ|0⟩ 를 돌려준다."""
    x = np.atleast_1d(np.asarray(x, float)); s = np.sqrt(np.clip(1 - x * x, 0, None))
    a = np.exp(1j * phis[0]) * np.ones_like(x, complex); bb = np.zeros_like(x, complex)
    for p in phis[1:]:
        a, bb = a * x + bb * (1j * s), a * (1j * s) + bb * x
        a, bb = a * np.exp(1j * p), bb * np.exp(-1j * p)
    return a


def solve_qsp(c_odd, tol=1e-13, maxit=30):
    """Im⟨0|U_Φ(x)|0⟩ = Σ_j c_j T_{2j+1}(x) 인 대칭 위상을 Newton 으로 푼다(Dong–Lin–Ni–Wang 2024 의 축약 위상).
    미지수는 오른쪽 절반 위상 d̃ = (d+1)/2 개. 0 에서의 야코비안이 2I 이므로 출발점은 c/2 다."""
    dt = len(c_odd); th = (2 * np.arange(1, dt + 1) - 1) * np.pi / (4 * dt); nodes = np.cos(th)
    Ainv = np.linalg.inv(np.cos(np.outer(th, 2 * np.arange(dt) + 1)))       # 점 값 → 홀수 체비쇼프 계수
    full = lambda red: np.concatenate([red[::-1], red])
    F = lambda red: Ainv @ qsp_wx(full(red), nodes).imag
    red = c_odd / 2.0; hist = []
    for _ in range(maxit):
        r = F(red) - c_odd; hist.append(np.abs(r).max())
        if hist[-1] < tol: break
        Jm = np.empty((dt, dt)); h = 1e-7
        for k in range(dt):
            e = np.zeros(dt); e[k] = h; Jm[:, k] = (F(red + e) - F(red - e)) / (2 * h)
        red = red - np.linalg.solve(Jm, r)
    return full(red), hist


PHI = {}
xx = np.linspace(-1, 1, 4001)
print("    κ    차수 d   위상 수   Newton 갱신   계수 잔차    조밀 격자에서 |Im P_Φ − P| 최대    걸린 시간")
NEWTON = {}
for kap in KAPPAS:
    t0 = time.time(); phis, hist = solve_qsp(POLY[kap]["c"])
    dense = np.abs(qsp_wx(phis, xx).imag - C.chebval(xx, cheb_full(POLY[kap]["c"]))).max()
    PHI[kap] = phis; NEWTON[kap] = len(hist) - 1            # hist 는 잔차 평가 횟수(마지막은 수렴 확인만) → 갱신 횟수는 하나 적다
    print(f"   {kap:3d}   {POLY[kap]['d']:5d}    {len(phis):5d}       {len(hist) - 1:2d}        {hist[-1]:.1e}          {dense:.1e}                  {time.time()-t0:5.1f}초")
    assert dense < 1e-9 and np.allclose(phis, phis[::-1])
_dt = len(POLY[2]["c"]); _th = (2 * np.arange(1, _dt + 1) - 1) * np.pi / (4 * _dt)
_Ainv = np.linalg.inv(np.cos(np.outer(_th, 2 * np.arange(_dt) + 1)))
_F = lambda red: _Ainv @ qsp_wx(np.concatenate([red[::-1], red]), np.cos(_th)).imag
_J0 = np.column_stack([(_F(1e-6 * e) - _F(-1e-6 * e)) / 2e-6 for e in np.eye(_dt)])
JAC0 = np.abs(_J0 - 2 * np.eye(_dt)).max()
print(f"    출발점의 근거: 위상 0 에서의 야코비안과 2I 의 차(κ = 2, {_dt}×{_dt}) = {JAC0:.1e} → 첫 근사 Φ ≈ c/2.  ‖c‖₁ = {POLY[2]['l1']:.2f} 는 고정점 반복의 수렴 증명 범위(0.861) 밖이다")
assert JAC0 < 1e-8
_zz = np.exp(1j * np.linspace(0, 2 * np.pi, 20001))
CIRC_MAX = {k: float(np.abs(np.polynomial.polynomial.polyval(_zz, POLY[k]["c"])).max()) for k in KAPPAS}
print("    단위원 위에서 |Σ_j c_{2j+1} z^j| 의 최댓값(= ‖c‖₁; 계수의 부호가 번갈아서): " + ", ".join(f"κ={k}: {CIRC_MAX[k]:.4f}" for k in KAPPAS)
      + f" — 1 을 넘으면 뿌리 찾기 방식의 위상 해법은 그대로 쓸 수 없다(6절). 0.6 배 하면 {0.6 * CIRC_MAX[2]:.3f}")
print(f"    위상의 크기: κ = 8 에서 |φ| 최대 {np.abs(PHI[8]).max():.3f} rad, 가운데 두 개 {PHI[8][len(PHI[8])//2-1]:.4f}; 양 끝 {PHI[8][0]:.2e}")


# =====================================================================
# 3) 블록 인코딩과 QSVT — 행렬 족 I + γS
# =====================================================================
print("\n" + "=" * 78); print("3) 블록 인코딩 위의 QSVT — I + γS 족 (사전 등록 본 실험)"); print("=" * 78)


def sparse_sym(N, rng_):
    """19장과 같은 무작위 대칭 희소 행렬 S: 밀도 2/N 의 ±1 행렬 R 에서 (R + Rᵀ)/2."""
    R = sp.random(N, N, density=2.0 / N, random_state=rng_, format="csr", data_rvs=lambda k: rng_.choice([-1.0, 1.0], size=k))
    return ((R + R.T) * 0.5).tocsr()


def family(S, kappa):
    """A = (I + γS)/‖I + γS‖ 에서 γ 를 이분법으로 맞춰 κ(A) = kappa 로 만든다. A 는 양의 정부호, 스펙트럼 ⊂ [1/κ, 1].
    사전 등록대로 S 는 크기마다 하나로 고정하고 γ 만 바꾼다."""
    N = S.shape[0]; ev = np.linalg.eigvalsh(S.toarray()); lo, hi = ev[0], ev[-1]
    assert lo < 0 < hi
    f = lambda g: (1 + g * hi) / (1 + g * lo) - kappa
    a, b = 0.0, -1 / lo * (1 - 1e-12)
    for _ in range(200):
        m = 0.5 * (a + b)
        if f(m) > 0: b = m
        else: a = m
    g = 0.5 * (a + b); A = sp.identity(N, format="csr") + g * S
    return (A / (1 + g * hi)).tocsr(), g, S


def dilation(Ad):
    """1큐비트 보조의 블록 인코딩(유니터리 팽창): U = [[A, √(I−A²)], [√(I−A²), −A]].  에르미트이면서 유니터리(반사). α = 1.
    회로가 아니라 행렬이다 — 호출 수를 세는 실험에는 이것으로 충분하고, 호출 하나의 게이트 수는 5)·6) 에서 따로 본다."""
    w, V = np.linalg.eigh(Ad); s = np.sqrt(np.clip(1 - w * w, 0, None)); B = (V * s) @ V.T
    return np.block([[Ad, B], [B, -Ad]])


def qsvt_apply(U, n_sys, phis_w, vec, proj_diag):
    """QSVT 열을 벡터에 적용한다. U 는 에르미트 블록 인코딩(각 고유값 λ 의 2차원 부분공간에서 [[λ, s], [s, −λ]]).
    W = i·e^{−iπ/4 Z_Π} U e^{−iπ/4 Z_Π} (Z_Π = 2Π − I) 이므로
      U_Φ = i^d · e^{i(φ0−π/4)Z_Π} U e^{i(φ1−π/2)Z_Π} U ⋯ U e^{i(φd−π/4)Z_Π}.
    proj_diag = Z_Π 의 대각(±1). 돌려주는 것: U_Φ·vec 전체."""
    d = len(phis_w) - 1
    sh = np.full(d + 1, np.pi / 2); sh[0] = sh[-1] = np.pi / 4
    ang = phis_w - sh
    v = np.exp(1j * ang[-1] * proj_diag) * vec
    for k in range(d - 1, -1, -1):
        v = U @ v
        v = np.exp(1j * ang[k] * proj_diag) * v
    return (1j) ** d * v


def qsvt_solve(U, N, phis_w, bvec, proj_diag):
    """Im 부분을 뽑는다: (U_Φ − U_{−Φ})/(2i) 의 Π 블록 = P(A).  부호 보조 큐비트 하나를 더 쓰는 LCU 에 해당한다."""
    pad = np.zeros(U.shape[0], complex); pad[:N] = bvec
    vp = qsvt_apply(U, N, phis_w, pad, proj_diag); vm = qsvt_apply(U, N, -phis_w, pad, proj_diag)
    return ((vp - vm) / 2j)[:N]


def cg_error_iters(A, b, xtrue, tol=EPS, maxit=100000):
    """켤레 기울기법을 직접 돌려 해의 상대 오차 ‖x_k − x‖/‖x‖ 가 tol 이하가 되는 첫 반복 수."""
    x = np.zeros_like(b); r = b.copy(); p = r.copy(); rs = r @ r; nx = np.linalg.norm(xtrue)
    for k in range(1, maxit + 1):
        Ap = A @ p; al = rs / (p @ Ap); x += al * p; r -= al * Ap
        if np.linalg.norm(x - xtrue) / nx <= tol: return k
        rs2 = r @ r; p = r + (rs2 / rs) * p; rs = rs2
    return maxit


ROWS = {}
print("     N    κ    차수 d   nnz    QSVT−다항식   해의 오차(최대)  p_s 표본 평균   공식 평균    p_s 최악 b    총 호출 d/√p_s   최악 b    CG 반복   CG nnz×반복")
SMAT = {}
for N in SIZES:
    SMAT[N] = sparse_sym(N, rng)
    for kap in KAPPAS:
        A, gam, S = family(SMAT[N], kap); Ad = A.toarray(); w = np.linalg.eigvalsh(Ad)
        assert abs(w[-1] - 1) < 1e-9 and abs(w[-1] / w[0] - kap) < 1e-6 * kap
        U = dilation(Ad); pd = np.concatenate([np.ones(N), -np.ones(N)])
        assert np.abs(U @ U - np.eye(2 * N)).max() < 1e-9 and np.abs(U - U.T).max() < 1e-12
        cfull = cheb_full(POLY[kap]["c"]); d = POLY[kap]["d"]
        evals, evecs = np.linalg.eigh(Ad); PA = (evecs * C.chebval(evals, cfull)) @ evecs.T        # P(A) 를 고유분해로
        q_err, s_err, ps, ps_f, cg_it = [], [], [], [], []
        for _ in range(10):
            b = rng.normal(size=N); b /= np.linalg.norm(b)
            xq = qsvt_solve(U, N, PHI[kap], b, pd)
            xt = np.linalg.solve(Ad, b)
            q_err.append(np.abs(xq - PA @ b).max())
            s_err.append(np.linalg.norm(xq.real / np.linalg.norm(xq) - xt / np.linalg.norm(xt)))
            ps.append(np.linalg.norm(xq) ** 2); ps_f.append(np.linalg.norm(xt) ** 2 / (2 * kap) ** 2)
            cg_it.append(cg_error_iters(A, b, xt))
        vmax = evecs[:, -1]; ps_worst = np.linalg.norm(qsvt_solve(U, N, PHI[kap], vmax, pd)) ** 2    # λ = 1 인 고유벡터: ‖A⁻¹b‖ = 1
        assert A.nnz == np.count_nonzero(Ad), "저장된 0 이 nnz 에 섞이면 CG 비용이 부풀려진다"
        ROWS[(N, kap)] = dict(d=d, nnz=A.nnz, ps_list=np.array(ps), psf_list=np.array(ps_f), q_err=max(q_err), s_err=max(s_err), ps=np.mean(ps), ps_f=np.mean(ps_f), ps_worst=ps_worst,
                              tot=d / np.sqrt(np.mean(ps)), tot_worst=d / np.sqrt(ps_worst), cg=np.mean(cg_it), cgcost=A.nnz * np.mean(cg_it), gamma=gam, amax=np.abs(Ad).max(), srow=int((Ad != 0).sum(axis=1).max()))
        r = ROWS[(N, kap)]
        print(f"   {N:4d}  {kap:3d}   {d:5d}  {r['nnz']:5d}    {r['q_err']:.1e}      {r['s_err']:.2e}       {r['ps']:.5f}      {r['ps_f']:.5f}     {r['ps_worst']:.2e}     {r['tot']:9.0f}    {r['tot_worst']:8.0f}    {r['cg']:5.1f}   {r['cgcost']:9.0f}")
        assert r["q_err"] < 1e-8 and r["s_err"] <= EPS
        assert abs(ps_worst * (2 * kap) ** 2 - 1) < 2.1 * EPS



def snippet_check(Ad, bvec, phis):
    """4절 Qiskit 스니펫과 같은 게이트 열을 numpy 로 돌린다. 큐비트(상위부터): 부호 보조, 블록 인코딩 보조, 시스템.
    Qiskit 의 정의 RZZ(θ) = exp(−iθ/2·Z⊗Z), RZ(θ) = exp(−iθ/2·Z) 를 그대로 쓴다. 돌려주는 것: (부호 = 1, 보조 = 0) 블록."""
    N = len(bvec); d = len(phis) - 1
    zs = np.repeat([1.0, -1.0], 2 * N); za = np.tile(np.repeat([1.0, -1.0], N), 2)
    rzz = lambda th: np.exp(-1j * th / 2 * zs * za); rz = lambda th: np.exp(-1j * th / 2 * za)
    Ub = np.kron(np.eye(2), dilation(Ad)); Hs = np.kron(np.array([[1.0, 1.0], [1.0, -1.0]]) / np.sqrt(2), np.eye(2 * N))
    sh = np.full(d + 1, np.pi / 2); sh[0] = sh[-1] = np.pi / 4
    v = np.zeros(4 * N, complex); v[:N] = bvec; v = Hs @ v
    v = rz(2 * sh[d]) * (rzz(-2 * phis[d]) * v)
    for k in range(d - 1, -1, -1):
        v = Ub @ v; v = rz(2 * sh[k]) * (rzz(-2 * phis[k]) * v)
    return (Hs @ v)[2 * N:3 * N]


# 스니펫용 사례 파일: N = 16, κ = 2 (d = 25). 4절의 Qiskit 스니펫과 6절의 PennyLane 스니펫이 이 파일을 읽는다.
A_c, _, _ = family(SMAT[16], 2); a_c = A_c.toarray(); b_c = np.random.default_rng(2020).normal(size=16); b_c /= np.linalg.norm(b_c)   # 본 실험의 난수 열을 건드리지 않게 따로 뽑는다
x_c = np.linalg.solve(a_c, b_c); blk_c = snippet_check(a_c, b_c, PHI[2])
ev_c, V_c = np.linalg.eigh(a_c); PA_c = (V_c * C.chebval(ev_c, cheb_full(POLY[2]["c"]))) @ V_c.T
CASE_DIFF = np.abs(blk_c - (1j) ** (1 - POLY[2]["d"]) * (PA_c @ b_c)).max()
CASE_PS = float(np.linalg.norm(blk_c) ** 2)
CASE_ERR = float(np.linalg.norm((blk_c / (1j) ** (1 - POLY[2]["d"])).real / np.linalg.norm(blk_c) - x_c / np.linalg.norm(x_c)))
mono = C.cheb2poly(cheb_full(POLY[2]["c"])); mono[::2] = 0.0
MONO_MAX = float(np.abs(mono).max())
MONO_ERR = float(np.abs(np.polynomial.polynomial.polyval(grid(2), mono) - C.chebval(grid(2), cheb_full(POLY[2]["c"]))).max())
with open(os.path.join(HERE, "ch20_case.json"), "w", encoding="utf-8") as fh:
    json.dump({"note": "ch20_linear_systems.py 가 만든 스니펫용 사례. N = 16, kappa = 2, 차수 25. 고치지 말고 스크립트를 다시 돌린다.",
               "N": 16, "kappa": 2, "degree": POLY[2]["d"], "A": [[float(f"{v:.17g}") for v in row] for row in a_c],
               "b": [float(f"{v:.17g}") for v in b_c], "phis_wx": [float(f"{v:.17g}") for v in PHI[2]],
               "cheb_odd": [float(f"{v:.17g}") for v in POLY[2]["c"]], "x_normalized": [float(f"{v:.17g}") for v in x_c / np.linalg.norm(x_c)],
               "p_success": CASE_PS, "solution_error": CASE_ERR}, fh, ensure_ascii=False, indent=1)
print(f"    스니펫용 사례(N = 16, κ = 2, d = 25 → ch20_case.json): 스니펫과 같은 게이트 열의 (부호 1, 보조 0) 블록 vs i^(1−d)·P(A)b: {CASE_DIFF:.1e}; p_s = {CASE_PS:.4f}, 해의 오차 {CASE_ERR:.2e}")
print(f"    같은 다항식의 단항식 계수(PennyLane 입력): 최대 크기 {MONO_MAX:.3g}, [1/κ, 1] 에서 체비쇼프 꼴과의 차 {MONO_ERR:.1e}")
CASE_INF = float(np.abs(a_c @ a_c.T).sum(axis=1).max()); CASE_2 = float(np.linalg.norm(a_c, 2))
print(f"    같은 행렬의 노름: ‖A‖₂ = {CASE_2:.6f}, ‖A·Aᵀ‖_∞(행 합) = {CASE_INF:.4f} — 블록 인코딩 도구가 2-노름이 아닌 노름으로 정규화하면 A 가 줄어든다(6절)")
assert CASE_DIFF < 1e-10 and CASE_ERR <= EPS

# =====================================================================
# 4) 사전 등록 기준에 대한 판정
# =====================================================================
print("\n" + "=" * 78); print("4) 사전 등록 기준(19장 5절)에 대한 결과"); print("=" * 78)
lk = np.log(KAPPAS)
SLOPE = {}
for N in SIZES:
    s_deg = np.polyfit(lk, np.log([ROWS[(N, k)]["d"] for k in KAPPAS]), 1)[0]
    s_tot = np.polyfit(lk, np.log([ROWS[(N, k)]["tot"] for k in KAPPAS]), 1)[0]
    s_wst = np.polyfit(lk, np.log([ROWS[(N, k)]["tot_worst"] for k in KAPPAS]), 1)[0]
    s_cg = np.polyfit(lk, np.log([ROWS[(N, k)]["cg"] for k in KAPPAS]), 1)[0]
    SLOPE[N] = (s_deg, s_tot, s_wst, s_cg)
    print(f"    N = {N:3d}: 차수 ∝ κ^{s_deg:.2f},  총 호출(무작위 b) ∝ κ^{s_tot:.2f},  총 호출(최악 b) ∝ κ^{s_wst:.2f},  CG 반복 ∝ κ^{s_cg:.2f}")
ok_deg = abs(SLOPE[256][0] - 1) <= 0.2; ok_tot = 1.3 <= SLOPE[256][1] <= 2.3
print(f"    기준 1 (N = 256 에서 차수의 κ 지수가 1 ± 0.2): {SLOPE[256][0]:.2f} → {'통과' if ok_deg else '불통과'}")
print(f"    기준 2 (N = 256 에서 총 호출 d/√p_s 의 κ 지수가 1.3 이상 2.3 이하): {SLOPE[256][1]:.2f} → {'통과' if ok_tot else '불통과'}")
_fit = lambda ys, ks=KAPPAS: float(np.polyfit(np.log(ks), np.log(ys), 1)[0])
ROB = {"차수: 꼬리 합 상한으로": _fit([POLY[k]["d_rig"] for k in KAPPAS]), "차수: CKS 끊는 자리로": _fit([POLY[k]["d_cks"] for k in KAPPAS]),
       "총 호출: b 마다 d/√p_s 를 구해 평균": _fit([np.mean(ROWS[(256, k)]["d"] / np.sqrt(ROWS[(256, k)]["ps_list"])) for k in KAPPAS]),
       "총 호출: 공식 p_s 로": _fit([ROWS[(256, k)]["d"] / np.sqrt(ROWS[(256, k)]["ps_f"]) for k in KAPPAS]),
       "총 호출: κ = 2~16 만": _fit([ROWS[(256, k)]["tot"] for k in KAPPAS[:-1]], KAPPAS[:-1]), "총 호출: κ = 4~32 만": _fit([ROWS[(256, k)]["tot"] for k in KAPPAS[1:]], KAPPAS[1:])}
print("    판정의 견고성(N = 256; 등록 뒤에 이 장이 정한 세부를 바꿨을 때의 지수): " + "; ".join(f"{k} {v:.2f}" for k, v in ROB.items()))
PS_AGREE = max(float(np.abs(r_["ps"] / r_["ps_f"] - 1)) for r_ in ROWS.values())
print(f"    p_s 표본 평균과 공식 평균의 상대 차(15개 설정의 최대): {PS_AGREE:.2%}")
print(f"    검산 (차수가 N 과 무관): κ = 8 의 차수는 N = 16, 64, 256 에서 {[ROWS[(N, 8)]['d'] for N in SIZES]}")
r = ROWS[(256, 32)]
print(f"    N = 256, κ = 32: 양자 총 호출 {r['tot']:.0f}번(진폭 증폭의 상수 π/2 를 넣으면 {np.pi/2*r['tot']:.0f}) 대 CG {r['cg']:.1f}회 반복 × nnz {r['nnz']} = 곱셈 {r['cgcost']:.0f}번")
SDD = {}
for kap in KAPPAS:
    Ak = family(SMAT[256], kap)[0].toarray(); off = np.abs(Ak).sum(axis=1) - np.abs(np.diag(Ak))
    SDD[kap] = (int((np.diag(Ak) >= off - 1e-12).sum()), float((off / np.diag(Ak)).max()))
print("    대각 지배 여부(N = 256; 행마다 a_ii ≥ Σ_{j≠i}|a_ij| 인 행의 수 / 256, 그리고 비율 Σ|a_ij|/a_ii 의 최댓값): " + ", ".join(f"κ={k}: {SDD[k][0]}행, {SDD[k][1]:.2f}" for k in KAPPAS))
print(f"    무작위 b 의 성공 확률이 최악의 b 보다 큰 정도(N = 256): " + ", ".join(f"κ={k}: {ROWS[(256,k)]['ps']/ROWS[(256,k)]['ps_worst']:.1f}배" for k in KAPPAS))

# 보강(사후 확인, 사전 등록 기준 밖): 무작위 단위 벡터 b 에서 E[b_i²] = 1/N 이므로 E[p_s] = (1/N) Σ_i 1/(2κλ_i)² 를 표본 없이 계산한다.
print("    보강(사후): 표본 대신 기댓값 E[p_s] = (1/N)Σ 1/(2κλ_i)² 로 다시 잰 총 호출의 κ 지수, 그리고 바닥 E[p_s] ≥ 1/(4N)")
EXACT = {}
rng_x = np.random.default_rng(2021)                              # 본 실험의 난수 열과 분리
for N in (16, 64, 256, 1024, 4096):
    mu = np.linalg.eigvalsh((SMAT[N] if N in SMAT else sparse_sym(N, rng_x)).toarray())
    eps_ = []
    for kap in KAPPAS:
        g_ = (kap - 1) / (mu[-1] - kap * mu[0]); lam = (1 + g_ * mu) / (1 + g_ * mu[-1])
        assert abs(lam[-1] / lam[0] - kap) < 1e-8 * kap
        eps_.append(np.mean(1 / (2 * kap * lam) ** 2))
    sl = np.polyfit(lk, np.log([POLY[k]["d"] / np.sqrt(e_) for k, e_ in zip(KAPPAS, eps_)]), 1)[0]
    EXACT[N] = (sl, eps_[-1], 1 / (4 * N))
    print(f"      N = {N:4d}{'' if N in SMAT else ' (새로 뽑은 S)'}: 총 호출 ∝ κ^{sl:.2f};  κ = 32 에서 E[p_s] = {eps_[-1]:.5f}, 바닥 1/(4N) = {1/(4*N):.5f} (E[p_s] 의 {1/(4*N)/eps_[-1]*100:.0f}%)")


# =====================================================================
# 5) 유효 조건수 — 실제 블록 인코딩은 α 만큼 손해를 본다
# =====================================================================
print("\n" + "=" * 78); print("5) 유효 조건수 κ_eff = α·κ — 블록 인코딩이 A/α 를 담을 때 다항식이 보는 조건수"); print("=" * 78)
print("    희소 접근 블록 인코딩(Gilyén et al. 2019 보조정리 48)은 α = s·max|a_ij| (s = 행당 비영 원소 수의 최댓값). ‖A‖ = 1 로 맞췄으므로 κ_eff = α·κ.")
print("     N    κ     s    max|a_ij|    α      κ_eff    차수 d(κ) → d(κ_eff)    배율")
EFF = {}
for kap in KAPPAS:
    r = ROWS[(256, kap)]; alpha = r["srow"] * r["amax"]; keff = alpha * kap
    d_eff = min_degree(keff)[1] if keff <= 70 else int(round(law_degree(keff)))
    EFF[kap] = (r["srow"], r["amax"], alpha, keff, d_eff)
    print(f"    256  {kap:3d}   {r['srow']:3d}     {r['amax']:.4f}    {alpha:5.2f}   {keff:7.1f}     {r['d']:5d} → {d_eff:6d}{'*' if keff > 70 else ' '}        {d_eff/r['d']:.1f}배")
print("    (* 는 κ_eff 가 커서 d ≈ 13.2·κ_eff 로 어림한 값)")


def pauli_l1(Ad):
    """16×16 실대칭 행렬의 Pauli 전개 계수 절댓값의 합 Σ|Tr(P·A)|/16 — LCU 블록 인코딩의 α."""
    P1 = {"I": np.eye(2), "X": np.array([[0, 1], [1, 0]]), "Y": np.array([[0, -1j], [1j, 0]]), "Z": np.diag([1.0, -1.0])}
    tot = 0.0; nz = 0
    for a_ in "IXYZ":
        for b_ in "IXYZ":
            Pab = np.kron(P1[a_], P1[b_])
            for c_ in "IXYZ":
                Pabc = np.kron(Pab, P1[c_])
                for d_ in "IXYZ":
                    cf = np.trace(np.kron(Pabc, P1[d_]) @ Ad) / 16
                    tot += abs(cf); nz += abs(cf) > 1e-12
    return tot, nz


A16, _, _ = family(SMAT[16], 8); a16 = A16.toarray()
al_p, n_p = pauli_l1(a16); al_s = int((a16 != 0).sum(axis=1).max()) * np.abs(a16).max()
print(f"    N = 16, κ = 8 의 같은 행렬을 세 가지로 담으면: 이상적 팽창 α = 1 (κ_eff = 8.0) | 희소 접근 α = {al_s:.2f} (κ_eff = {8*al_s:.1f}) | Pauli LCU α = {al_p:.2f}, 항 {n_p}개 (κ_eff = {8*al_p:.1f})")


# =====================================================================
# 6) 게이트 수준의 끝에서 끝까지 — n = 4, A = I + γ(ΣZZ + ΣX), LCU 블록 인코딩
# =====================================================================
print("\n" + "=" * 78); print("6) 게이트 리스트로 푼 작은 사례 — A = I + γ(Σ Z_i Z_{i+1} + Σ X_i), n = 4, κ = 4"); print("=" * 78)
I2 = np.eye(2); X2 = np.array([[0.0, 1.0], [1.0, 0.0]]); Z2 = np.diag([1.0, -1.0])


def kron_all(ops):
    out = np.array([[1.0]])
    for o in ops: out = np.kron(out, o)
    return out


n_s = 4
TERMS = [("I", ())] + [("ZZ", (i, i + 1)) for i in range(n_s - 1)] + [("X", (i,)) for i in range(n_s)]      # 8항 → 색인 3큐비트
Smat = sum(kron_all([Z2 if q in qs else I2 for q in range(n_s)]) for name, qs in TERMS if name == "ZZ") \
    + sum(kron_all([X2 if q in qs else I2 for q in range(n_s)]) for name, qs in TERMS if name == "X")
evS = np.linalg.eigvalsh(Smat); kap6 = 4.0
gam6 = (kap6 - 1) / (evS[-1] - kap6 * evS[0])                    # (1 + γλmax)/(1 + γλmin) = κ
A6 = np.eye(2 ** n_s) + gam6 * Smat; w6 = np.linalg.eigvalsh(A6)
alpha6 = 1 + 7 * gam6; keff6 = alpha6 / w6[0]
print(f"    S 의 고유값 범위 [{evS[0]:.4f}, {evS[-1]:.4f}], γ = {gam6:.4f} → A 의 고유값 [{w6[0]:.4f}, {w6[-1]:.4f}], κ = {w6[-1]/w6[0]:.3f}")
print(f"    LCU: 항 8개(I 1개 + ZZ 3개 + X 4개), 계수 (1, γ×7), α = 1 + 7γ = {alpha6:.4f}, ‖A‖ = {w6[-1]:.4f} → κ_eff = α/λ_min = {keff6:.3f} (κ 의 {keff6/kap6:.2f}배)")

# 큐비트 배치(빅 엔디안): 0 = 부호 보조 c | 1,2,3 = 색인 a0 a1 a2 | 4,5 = 작업 w1 w2 | 6..9 = 시스템 q0..q3
QC, QA, QW, QS = 0, (1, 2, 3), (4, 5), (6, 7, 8, 9); NQ = 10


def ccry(c1, c2, t, th):
    """이중 제어 Ry(θ) = CRy(θ/2)(c2,t)·CX(c1,c2)·CRy(−θ/2)(c2,t)·CX(c1,c2)·CRy(θ/2)(c1,t)."""
    if abs(th) < 1e-12: return []
    return [gc.G("cry", c2, t, angle=th / 2), gc.G("cx", c1, c2), gc.G("cry", c2, t, angle=-th / 2), gc.G("cx", c1, c2), gc.G("cry", c1, t, angle=th / 2)]


def prep_gates(amp):
    """실수·양수 진폭 8개를 색인 3큐비트에 싣는 Ry 나무. PREP|000⟩ = Σ amp_i |i⟩."""
    p = amp ** 2; a0, a1, a2 = QA
    ang = lambda p1, p0: 2 * np.arcsin(np.sqrt(p1 / (p0 + p1))) if p0 + p1 > 0 else 0.0
    g = [gc.G("ry", a0, angle=ang(p[4:].sum(), p[:4].sum()))]
    t1 = [ang(p[2:4].sum(), p[0:2].sum()), ang(p[6:8].sum(), p[4:6].sum())]                     # a0 = 0, 1 일 때 a1 의 각
    g += [gc.G("ry", a1, angle=t1[0]), gc.G("cry", a0, a1, angle=t1[1] - t1[0])]
    t2 = [[ang(p[1], p[0]), ang(p[3], p[2])], [ang(p[5], p[4]), ang(p[7], p[6])]]              # (a0, a1) 일 때 a2 의 각
    g += [gc.G("ry", a2, angle=t2[0][0]), gc.G("cry", a0, a2, angle=t2[1][0] - t2[0][0]), gc.G("cry", a1, a2, angle=t2[0][1] - t2[0][0])]
    g += ccry(a0, a1, a2, t2[1][1] - t2[1][0] - t2[0][1] + t2[0][0])
    return g


def and3(flag_bits):
    """색인이 주어진 비트열일 때만 w2 = 1 이 되게 한다(Toffoli 2개, 작업 큐비트 w1 경유). 되돌리려면 뒤집어 쓴다."""
    a0, a1, a2 = QA; w1, w2 = QW
    xs = [gc.G("x", q) for q, bit in zip(QA, flag_bits) if bit == 0]
    return xs + [gc.G("ccx", a0, a1, w1), gc.G("ccx", w1, a2, w2)], [gc.G("ccx", w1, a2, w2), gc.G("ccx", a0, a1, w1)] + xs


def select_gates():
    """SELECT = Σ_i |i⟩⟨i| ⊗ P_i. 항마다 색인을 w2 에 모으고 w2 로 제어한 Pauli 를 건 뒤 되돌린다(단순 구성)."""
    g = []
    for i, (name, qs) in enumerate(TERMS):
        if name == "I": continue
        bits = [(i >> 2) & 1, (i >> 1) & 1, i & 1]
        comp, uncomp = and3(bits)
        ctl = [gc.G("cz", QW[1], QS[q]) for q in qs] if name == "ZZ" else [gc.G("cx", QW[1], QS[qs[0]])]
        g += comp + ctl + uncomp
    return g


amp6 = np.sqrt(np.array([1.0] + [gam6] * 7) / alpha6)
PREP = prep_gates(amp6); U6 = PREP + select_gates() + gc.inverse(PREP)


def relabel(gates, shift):
    return [(nm, tuple(q - shift for q in qq), an) for nm, qq, an in gates]


U6mat = gc.unitary(relabel(U6, 1), NQ - 1)                         # 부호 보조를 뺀 9큐비트 유니터리(512×512)
prep_col = gc.unitary(relabel(PREP, 1), NQ - 1)[:, 0].reshape(8, -1)[:, 0]
blk = U6mat.reshape(8, 4, 16, 8, 4, 16)[0, 0, :, 0, 0, :]          # 색인 000, 작업 00 블록
print(f"    PREP 검산: ‖PREP|000⟩ − √(c/α)‖ = {np.abs(prep_col - amp6).max():.1e};  블록 검산: ‖⟨0|U|0⟩ − A/α‖ = {np.abs(blk - A6/alpha6).max():.1e};  U 는 유니터리: {np.abs(U6mat.conj().T @ U6mat - np.eye(512)).max():.1e}")
assert np.abs(prep_col - amp6).max() < 1e-12 and np.abs(blk - A6 / alpha6).max() < 1e-12
leak = np.abs(U6mat.reshape(8, 4, 16, 8, 4, 16)[:, 1:, :, :, 0, :]).max()
U_sec = U6mat.reshape(8, 4, 16, 8, 4, 16)[:, 0, :, :, 0, :].reshape(128, 128)       # 작업 = 00 인 부분공간
HERM6 = float(np.abs(U_sec - U_sec.conj().T).max()); REFL6 = float(np.abs(U_sec @ U_sec - np.eye(128)).max())
print(f"    작업 = 00 부분공간에서 U 는 에르미트이고 반사다: ‖U − U†‖ = {HERM6:.1e}, ‖U² − I‖ = {REFL6:.1e} (3절의 W 변환이 쓰는 성질)")
assert HERM6 < 1e-12 and REFL6 < 1e-12
print(f"    작업 큐비트는 매 호출 뒤 00 으로 돌아온다(작업 ≠ 00 블록의 최대 성분 {leak:.1e})")
print(f"    U 한 번의 게이트 비용: {gc.fmt(U6)}")


def proj_phase(phi, shift):
    """e^{i(sφ − shift)Z_Π}: Π = 색인 000. s = +1(c = 0), −1(c = 1). Z_Π = −Z_{w2} 이므로 e^{−iφ Z_c Z_w2}·e^{+i·shift·Z_w2}."""
    comp, uncomp = and3([0, 0, 0])
    return comp + [gc.G("rzz", QC, QW[1], angle=phi), gc.G("rz", QW[1], angle=-2 * shift)] + uncomp


b6, d6, c6, _ = min_degree(keff6); phi6, hist6 = solve_qsp(c6)
sh6 = np.full(d6 + 1, np.pi / 2); sh6[0] = sh6[-1] = np.pi / 4
circ = [gc.G("h", QC)] + proj_phase(phi6[-1], sh6[-1])
for k in range(d6 - 1, -1, -1):
    circ += U6 + proj_phase(phi6[k], sh6[k])
circ += [gc.G("h", QC)]
ONE_PHASE = proj_phase(phi6[1], sh6[1])
print(f"    κ_eff = {keff6:.3f} 의 다항식: b = {b6}, 차수 d = {d6}, 위상 {len(phi6)}개(Newton 갱신 {len(hist6) - 1}회, 잔차 {hist6[-1]:.1e})")
print(f"    사영 위상 한 번의 게이트 비용: {gc.fmt(ONE_PHASE)}")
t6 = time.time()
b_vec = rng.normal(size=16); b_vec /= np.linalg.norm(b_vec)
psi0 = np.zeros(2 ** NQ, complex); psi0[:16] = b_vec                # c = 0, 색인 000, 작업 00, 시스템 = b
psi = gc.apply(circ, psi0, NQ)
out = psi.reshape(2, 8, 4, 16)[1, 0, 0, :]                          # 부호 보조 1, 색인 000, 작업 00
ev6, V6 = np.linalg.eigh(A6 / alpha6); PA6 = (V6 * C.chebval(ev6, cheb_full(c6))) @ V6.T
want = PA6 @ b_vec; x_true = np.linalg.solve(A6, b_vec)
ph = np.vdot(want, out) / abs(np.vdot(want, out))                   # 전역 위상(i^{d+1} 을 회로에서 뺐다)
ps6 = np.linalg.norm(out) ** 2
err6 = np.linalg.norm(out / ph / np.linalg.norm(out) - x_true / np.linalg.norm(x_true))
print(f"    회로 전체: 게이트 {len(gc.expand(circ))}개(전개 후), 10큐비트 상태 벡터에 적용 {time.time()-t6:.1f}초")
assert abs(ph - (1j) ** (1 - d6)) < 1e-9                             # 회로에서 뺀 전역 위상 i^d 와 Im 추출의 i: 블록 = i^{1−d}·P(A/α)b
print(f"    끝 상태의 (c = 1, 색인 000, 작업 00) 블록 vs i^(1−d)·P(A/α)b: ‖차‖ = {np.abs(out/ph - want).max():.1e}")
print(f"    성공 확률 p_s = {ps6:.5f} (공식 ‖A⁻¹b‖²·λ_min²/4 = {np.linalg.norm(x_true)**2 * w6[0]**2 / 4:.5f}), 해의 오차 {err6:.2e}")
assert np.abs(out / ph - want).max() < 1e-8 and err6 <= EPS
cost6 = gc.cost(circ)
print(f"    QSVT 한 번(진폭 증폭 없음)의 게이트 비용: {gc.fmt(circ)}")
G6 = gc.two_qubit_count(circ); GU = gc.two_qubit_count(U6); GP = gc.two_qubit_count(ONE_PHASE)
print(f"    검산: d·G_U + (d+1)·G_Π = {d6}·{GU} + {d6+1}·{GP} = {d6*GU + (d6+1)*GP} (리스트에서 센 값 {G6})")
assert G6 == d6 * GU + (d6 + 1) * GP
for r_cz, tag in ((3.0e-3, "14장 공개 중앙값"), (2.0e-3, "부록 A: ibm_fez 2026-11-11 실행에 쓰인 연결의 범위 1.58~2.46e-3 의 가운데")):
    print(f"    CZ 오류율 r = {r_cz:.1e} ({tag}): rG = {r_cz*G6:.1f}, 오류 없이 끝날 확률 e^(−rG) = {np.exp(-r_cz*G6):.1e}")
print(f"    14장의 오류 예산 rG ≈ 0.3~0.5 에 견주면: {2.0e-3*G6/0.5:.0f}배 ~ {3.0e-3*G6/0.3:.0f}배")
AA6 = np.pi / 2 / np.sqrt(ps6); T_CZ = 68e-9                          # 14장: Heron r2 의 CZ 길이 68 ns
M6 = int(np.floor(np.pi / (4 * np.arcsin(np.sqrt(ps6))))); G6_AA = (2 * M6 + 1) * G6
print(f"    진폭 증폭까지 넣으면 반복 m = {M6}회, QSVT 회로 2m + 1 = {2*M6+1}번(어림 (π/2)/√p_s = {AA6:.1f}배): 2큐비트 게이트 {G6_AA}개 (b 준비와 반사의 게이트는 뺀 값)")
cg6 = cg_error_iters(sp.csr_matrix(A6), b_vec, x_true); nnz6 = int((A6 != 0).sum())
print(f"    같은 문제의 CG: 반복 {cg6}회 × nnz {nnz6} = 곱셈 {cg6*nnz6}번. 양자 쪽은 오류가 없다고 치고 CZ 를 차례로만 세어도 {G6_AA} × 68 ns = {G6_AA*T_CZ*1e3:.1f} ms")

# =====================================================================
# 7) 조건이 맞는 족과 맞지 않는 족 — Poisson 대 I + γS
# =====================================================================
print("\n" + "=" * 78); print("7) 조건이 맞지 않는 족 — Poisson 행렬은 κ 가 N 과 함께 자란다"); print("=" * 78)


def poisson1d(N): return sp.diags([-1.0, 2.0, -1.0], [-1, 0, 1], shape=(N, N), format="csr")
def poisson2d(L): I_ = sp.identity(L, format="csr"); T = poisson1d(L); return (sp.kron(I_, T) + sp.kron(T, I_)).tocsr()


FAM = []
print("    행렬 족                N       κ        차수 d      E[p_s]   (표본 10개 평균)   양자 총 호출 d/√E[p_s]    CG 반복(오차 1e-3)   nnz    CG nnz×반복")
for N in SIZES:
    L = int(round(np.sqrt(N)))
    for name, A in (("1D Poisson", poisson1d(N)), ("2D Poisson", poisson2d(L)), ("I + γS (γ = 0.2)", (sp.identity(N, format="csr") + 0.2 * SMAT[N]).tocsr())):
        Ad = A.toarray(); w = np.linalg.eigvalsh(Ad); kap = w[-1] / w[0]; An = Ad / w[-1]
        nnz = int(np.count_nonzero(Ad))                   # A.nnz 는 저장된 0 까지 센다(sp.kron 이 작은 행렬에서 BSR 을 돌려준다) — 검토 R12
        direct = kap <= 70
        d = min_degree(kap)[1] if direct else int(round(law_degree(kap)))
        ps_exact = float(np.mean(1 / (2 * kap * (w / w[-1])) ** 2))      # 무작위 단위 벡터의 기댓값. 표본 10개의 평균은 20~40% 흔들린다
        ps, its = [], []
        for _ in range(10):
            b = rng.normal(size=N); b /= np.linalg.norm(b); xt = np.linalg.solve(An, b)
            ps.append(np.linalg.norm(xt) ** 2 / (2 * kap) ** 2); its.append(cg_error_iters(sp.csr_matrix(An), b, xt))
        FAM.append((name, N, kap, d, ps_exact, d / np.sqrt(ps_exact), np.mean(its), nnz * np.mean(its), direct, np.mean(ps), nnz))
        print(f"    {name:18s} {N:4d}   {kap:9.1f}   {d:8d}{' ' if direct else '*'}    {ps_exact:.5f}     ({np.mean(ps):.4f})        {d/np.sqrt(ps_exact):12.0f}            {np.mean(its):6.1f}         {nnz:5d}    {nnz*np.mean(its):9.0f}")
print("    (* 는 d ≈ 13.2·κ 로 어림. E[p_s] = (1/N)Σ 1/(2κλ_i)² — 3)·4) 에서 확인한 공식의 기댓값. 표본 평균은 참고로만 적는다.)")
FAM_SLOPE = {}
for name in ("1D Poisson", "2D Poisson", "I + γS (γ = 0.2)"):
    rows = [f for f in FAM if f[0] == name]; ln = np.log([f[1] for f in rows])
    FAM_SLOPE[name] = tuple(float(np.polyfit(ln, np.log([f[i] for f in rows]), 1)[0]) for i in (2, 5, 7))
    print(f"    {name:18s}: κ ∝ N^{FAM_SLOPE[name][0]:.2f},  양자 총 호출 ∝ N^{FAM_SLOPE[name][1]:.2f},  CG 곱셈 ∝ N^{FAM_SLOPE[name][2]:.2f}")

# =====================================================================
# 8) 읽어 내기 — 해 벡터가 아니라 기댓값 하나
# =====================================================================
print("\n" + "=" * 78); print("8) 읽어 내기 비용 — 기댓값 ⟨x|M|x⟩ 하나를 오차 δ 로"); print("=" * 78)
A8, _, _ = family(SMAT[256], 8); a8 = A8.toarray(); b8 = rng.normal(size=256); b8 /= np.linalg.norm(b8)
x8 = np.linalg.solve(a8, b8); x8 /= np.linalg.norm(x8)
Mdiag = np.where(np.arange(256) < 128, 1.0, -1.0)                 # M = 첫 큐비트의 Z
mval = float(x8 @ (Mdiag * x8)); var = 1 - mval ** 2
r8 = ROWS[(256, 8)]
for delta in (0.1, 0.01):
    shots = var / delta ** 2
    print(f"    N = 256, κ = 8, M = Z_0: ⟨M⟩ = {mval:+.4f}, 분산 {var:.3f} → δ = {delta}: 샷 {shots:,.0f}번 × 한 번의 호출 {np.pi/2*r8['tot']:.0f} = 블록 인코딩 호출 {shots*np.pi/2*r8['tot']:,.0f}번")
print(f"    해 벡터 전체(성분 256개)를 같은 오차로 읽으려면 진폭마다 이만큼이 든다. CG 는 곱셈 {r8['cgcost']:.0f}번으로 벡터 전체를 준다.")

# =====================================================================
# (E) 연습문제 검산
# =====================================================================
print("\n" + "=" * 78); print("(E) 연습문제 검산"); print("=" * 78)
# 문제 1: A = diag(1, 1/2), κ = 2, 이상적인 f(x) = 1/(2κx)
fA = np.array([1 / (2 * 2 * 1.0), 1 / (2 * 2 * 0.5)])
E1 = {}
for tag, bv in (("b = (1, 1)/√2", np.array([1.0, 1.0]) / np.sqrt(2)), ("b = (1, 0)", np.array([1.0, 0.0])), ("b = (0, 1)", np.array([0.0, 1.0]))):
    out_ = fA * bv; E1[tag] = float(out_ @ out_)
    print(f"    문제 1: {tag}: f(A)b = {out_.round(4)}, p_s = {E1[tag]:.5f}, (π/2)/√p_s = {np.pi/2/np.sqrt(E1[tag]):.2f}, x ∝ {(out_/np.linalg.norm(out_)).round(4)}")
assert abs(E1["b = (1, 1)/√2"] - 5 / 32) < 1e-15 and abs(E1["b = (1, 0)"] - 1 / 16) < 1e-15 and abs(E1["b = (0, 1)"] - 1 / 4) < 1e-15
# 문제 2: b ≥ κ² ln(2/ε) 이면 |x| ≥ 1/κ 에서 (1 − x²)^b ≤ ε/2
for kap in (8, 32):
    b_ = POLY[kap]["b"]
    print(f"    문제 2: κ = {kap}: κ²·ln(2/ε) = {kap**2*np.log(2/EPS):.2f} → b = {b_}; (1 − 1/κ²)^b = {(1 - 1/kap**2)**b_:.3e} ≤ e^(−b/κ²) = {np.exp(-b_/kap**2):.3e} ≤ ε/2 = {EPS/2:.1e}")
    assert (1 - 1 / kap ** 2) ** b_ <= np.exp(-b_ / kap ** 2) <= EPS / 2
# 문제 3: 진폭 증폭의 정확한 반복 수 (N = 256, κ = 32 의 평균 p_s)
E3 = {}
for tag, p_, d_ in (("N = 256, κ = 32, 무작위 b", ROWS[(256, 32)]["ps"], 423), ("N = 256, κ = 32, 최악의 b", ROWS[(256, 32)]["ps_worst"], 423), ("6) 의 n = 4 사례", ps6, d6)):
    th_ = np.arcsin(np.sqrt(p_)); m_ = int(np.floor(np.pi / (4 * th_))); pf = np.sin((2 * m_ + 1) * th_) ** 2
    E3[tag] = (m_, pf, (2 * m_ + 1) * d_)
    print(f"    문제 3: {tag}: p_s = {p_:.5f}, θ = {th_:.5f} rad, 반복 m = ⌊π/(4θ)⌋ = {m_}, 증폭 뒤 성공 확률 sin²((2m+1)θ) = {pf:.4f}, 블록 인코딩 호출 (2m+1)·d = {(2*m_+1)*d_}  (어림 (π/2)·d/√p_s = {np.pi/2*d_/np.sqrt(p_):.0f})")
# 문제 4: α 는 성공 확률을 바꾸지 않고 차수만 바꾼다 — 6) 의 행렬을 α = 1.8826 과 α' = ‖A‖ = 1.6 으로 담았을 때
for al_ in (alpha6, w6[-1]):
    ke_ = al_ / w6[0]; _, d_al, c_al, _ = min_degree(ke_)
    ev_, V_ = np.linalg.eigh(A6 / al_); ps_al = np.linalg.norm(((V_ * C.chebval(ev_, cheb_full(c_al))) @ V_.T) @ b_vec) ** 2
    print(f"    문제 4: α = {al_:.4f}: κ_eff = {ke_:.3f}, 차수 d = {d_al}, p_s = ‖P(A/α)b‖² = {ps_al:.5f}  (공식 ‖A⁻¹b‖²λ_min²/4 = {np.linalg.norm(x_true)**2*w6[0]**2/4:.5f})")

# =====================================================================
# 9) 그림
# =====================================================================
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
matplotlib.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
from matplotlib.ticker import NullFormatter, NullLocator


def plain_ticks(ax, vals):
    """로그 x축의 눈금을 실제로 쓴 값(κ 또는 N)에만 두고 보통 숫자로 적는다."""
    ax.set_xticks(list(vals)); ax.set_xticklabels([str(v) for v in vals])
    ax.xaxis.set_minor_locator(NullLocator()); ax.xaxis.set_minor_formatter(NullFormatter())


fig, axs = plt.subplots(3, 2, figsize=(12, 12.6))
ax = axs[0, 0]; xg = np.linspace(-1, 1, 2001); kk = 8
ax.plot(xg, C.chebval(xg, cheb_full(POLY[kk]["c"])), label=f"P(x), 차수 {POLY[kk]['d']}")
xr = np.linspace(1 / kk, 1, 200); ax.plot(xr, 1 / (2 * kk * xr), "k--", lw=1, label="1/(2κx)"); ax.plot(-xr, -1 / (2 * kk * xr), "k--", lw=1)
ax.axvspan(-1 / kk, 1 / kk, color="0.9"); ax.set_ylim(-1.05, 1.05); ax.set_xlabel("x (A/α의 고유값)"); ax.set_ylabel("P(x)"); ax.legend(fontsize=8)
ax.set_title("(1) 역행렬 다항식(κ = 8): 회색 띠 안에서는 1/x를 따라가지 않는다")
ax = axs[0, 1]
ax.loglog(KAPPAS, [POLY[k]["d"] for k in KAPPAS], "o-", label="차수 d (오차 10⁻³)")
ax.loglog(KAPPAS, [POLY[k]["d_cks"] for k in KAPPAS], "s--", label="CKS 끊는 자리")
ax.loglog(KAPPAS, [ROWS[(256, k)]["cg"] for k in KAPPAS], "^-", label="CG 반복 수 (N = 256)")
ax.set_xlabel("조건수 κ"); ax.set_ylabel("횟수 (로그)"); ax.legend(fontsize=8)
plain_ticks(ax, KAPPAS); ax.set_title(f"(2) 차수 ∝ κ^{deg_slope:.2f}, CG 반복 ∝ κ^{SLOPE[256][3]:.2f}")
ax = axs[1, 0]
for N in SIZES: ax.loglog(KAPPAS, [ROWS[(N, k)]["ps"] for k in KAPPAS], "o-", label=f"무작위 b, N = {N}")
ax.loglog(KAPPAS, [ROWS[(256, k)]["ps_worst"] for k in KAPPAS], "k--", label="최악의 b: 1/(2κ)²")
ax.set_xlabel("조건수 κ"); ax.set_ylabel("성공 확률 p_s (로그)"); ax.legend(fontsize=8)
plain_ticks(ax, KAPPAS); ax.set_title("(3) 성공 확률: 무작위 b는 최악보다 낫지만 κ와 함께 준다")
ax = axs[1, 1]
ax.loglog(KAPPAS, [ROWS[(256, k)]["tot"] for k in KAPPAS], "o-", label="양자 총 호출 d/√p_s (무작위 b)")
ax.loglog(KAPPAS, [ROWS[(256, k)]["tot_worst"] for k in KAPPAS], "s--", label="양자 총 호출 (최악의 b)")
ax.loglog(KAPPAS, [ROWS[(256, k)]["cgcost"] for k in KAPPAS], "^-", label="CG: nnz × 반복 수")
ax.set_xlabel("조건수 κ"); ax.set_ylabel("횟수 (로그, 단위가 서로 다르다)"); ax.legend(fontsize=8)
plain_ticks(ax, KAPPAS); ax.set_title(f"(4) N = 256: 총 호출 ∝ κ^{SLOPE[256][1]:.2f} (최악 κ^{SLOPE[256][2]:.2f})")
ax = axs[2, 0]
for name, mk in (("1D Poisson", "o-"), ("2D Poisson", "s-"), ("I + γS (γ = 0.2)", "^-")):
    rows = [f for f in FAM if f[0] == name]
    ax.loglog([f[1] for f in rows], [f[5] for f in rows], mk, label=f"{name}: 양자 총 호출")
    ax.loglog([f[1] for f in rows], [f[7] for f in rows], mk.replace("-", ":"), color="0.5", label="CG 곱셈 수(같은 표식)" if name == "1D Poisson" else None)
ax.set_xlabel("행렬 크기 N"); ax.set_ylabel("횟수 (로그, 단위가 서로 다르다)"); ax.legend(fontsize=8)
plain_ticks(ax, SIZES); ax.set_title("(5) 크기 의존: 양자 총 호출(실선)과 CG 곱셈 수(회색 점선)")
ax = axs[2, 1]
lab = ["이상적 팽창", "희소 접근", "Pauli LCU"]; vals = [8.0, 8 * al_s, 8 * al_p]
ax.bar(lab, vals, color=["C0", "C1", "C2"]); ax.set_ylabel("유효 조건수 κ_eff")
for i, v in enumerate(vals): ax.text(i, v, f"{v:.1f}", ha="center", va="bottom")
ax.set_title("(6) 같은 행렬(N = 16, κ = 8)을 담는 방법에 따른 유효 조건수")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "ch20_results.png"), dpi=120)

# ── 회귀 검사: 본문이 인용하는 값 ──
assert [POLY[k]["d"] for k in KAPPAS] == [25, 51, 105, 211, 423], "5절 (1): 차수"
assert [POLY[k]["b"] for k in KAPPAS] == [31, 122, 487, 1946, 7784], "5절 (1): b"
assert [POLY[k]["d_cks"] for k in KAPPAS] == [41, 81, 169, 353, 735], "5절 (1): CKS 상한"
assert [r_[1] for r_ in EPS_ROWS] == [67, 105, 139, 211], "5절 (1): 정밀도 의존"
assert abs(deg_slope - 1.021) < 5e-4, "5절 (1): 차수의 κ 지수"
assert ok_deg and ok_tot, "5절 (4): 사전 등록 기준"
assert abs(SLOPE[256][1] - 1.64) < 0.005 and abs(SLOPE[256][2] - 2.02) < 0.005 and abs(SLOPE[256][3] - 0.54) < 0.005, "5절 (4): N = 256 의 지수"
assert abs(SLOPE[16][1] - 1.38) < 0.005 and abs(SLOPE[64][1] - 1.53) < 0.005, "5절 (4): N = 16, 64 의 지수"
assert round(ROWS[(256, 32)]["tot"]) == 6659 and round(ROWS[(256, 32)]["cgcost"]) == 24560, "5절 (3): N = 256, κ = 32"
assert [EFF[k][4] for k in KAPPAS] == [251, 463, 889, 1736, 3432], "5절 (5): 유효 조건수의 차수"
assert abs(al_s - 4.82) < 0.005 and abs(al_p - 2.84) < 0.005 and n_p == 88, "5절 (5): N = 16 의 α"
assert d6 == 61 and GU == 206 and GP == 26 and G6 == 14178, "5절 (6): 게이트 수"
assert gc.depth(circ) == 22675 and gc.t_count(circ) == (13692, 1892), "5절 (6): 깊이와 T"
assert abs(ps6 - 0.04679) < 5e-6 and abs(keff6 - 4.706) < 5e-4, "5절 (6): 성공 확률과 κ_eff"
assert [f[3] for f in FAM] == [1539, 125, 25, 22626, 425, 33, 353839, 1539, 33], "5절 (7): 족별 차수"
assert [E3[k][2] for k in E3] == [10575, 42723, 427] and [E3[k][0] for k in E3] == [12, 50, 3], "8절 문제 3: 진폭 증폭 반복 수"
assert abs(CASE_PS - 0.1249) < 5e-5 and abs(MONO_MAX - 4.94e5) < 5e2, "4·6절: 스니펫용 사례"
assert [round(v[2], 3) for v in EPS_ROWS] == [0.740, 0.881, 1.005, 1.216], "5절 (1): 정밀도에 따른 max|P|"
assert NEWTON == {2: 5, 4: 4, 8: 4, 16: 4, 32: 4} and len(hist6) - 1 == 4, "5절 (2)·(6): Newton 갱신 횟수"
assert abs(CIRC_MAX[2] - 1.5639) < 5e-5, "6절: 단위원 위 최댓값"
assert [round(ROWS[(N, k)]["tot"]) for N in SIZES for k in KAPPAS] == [71, 214, 641, 1257, 3535, 72, 212, 703, 1812, 5038, 73, 224, 706, 2380, 6659], "5절 (3): 총 호출"
assert [round(ROWS[(N, k)]["cgcost"]) for N in SIZES for k in KAPPAS] == [303, 444, 592, 651, 725, 1311, 2056, 2682, 3308, 3785, 5317, 8862, 12660, 17724, 24560], "5절 (3): CG nnz × 반복"
assert np.allclose([ROWS[(N, k)]["ps"] for N in SIZES for k in KAPPAS],
                   [0.12482, 0.05686, 0.02680, 0.02820, 0.01432, 0.11960, 0.05788, 0.02230, 0.01356, 0.00705, 0.11712, 0.05182, 0.02209, 0.00786, 0.00403], atol=6e-6), "5절 (3): p_s 평균"
# QSVT − 다항식의 차는 차수 d 번의 곱에서 쌓인 반올림이라 끝자리가 계산기마다 다르다(2026-11-13: 기준 환경 3.8e-14, 맥 4.2e-14).
# 상한은 배정밀도 ε × 최대 차수(2.2e-16 × 423 ≈ 9.4e-14)에서 잡은 1e-13 이다 — 본문은 "10⁻¹³ 미만(반올림 수준)"으로 적는다.
Q_ERR_MAX = max(r_["q_err"] for r_ in ROWS.values())
print(f"    오차 상한 점검: 해의 오차 최대 {max(r_['s_err'] for r_ in ROWS.values()):.2e}, QSVT − 다항식 최대 {Q_ERR_MAX:.1e} (ε·d_max = {np.finfo(float).eps * max(r_['d'] for r_ in ROWS.values()):.1e}), p_s 표본·공식 차 최대 {PS_AGREE:.2%}")
assert max(r_["s_err"] for r_ in ROWS.values()) < 7.5e-4 and Q_ERR_MAX < 1e-13 and PS_AGREE < 2e-3, "5절 (3): 오차 상한"
assert [round(EXACT[N][0], 2) for N in (16, 64, 256, 1024, 4096)] == [1.38, 1.58, 1.64, 1.69, 1.71], "5절 (4) 사후 확인: 기댓값으로 잰 지수"
assert [round(v, 2) for v in ROB.values()] == [1.05, 1.05, 1.65, 1.64, 1.67, 1.64], "5절 (4): 판정의 견고성"
assert [SDD[k][0] for k in KAPPAS] == [256, 238, 204, 204, 167], "7절: 대각 지배인 행의 수"
assert [round(ROWS[(256, k)]["ps"] / ROWS[(256, k)]["ps_worst"], 1) for k in KAPPAS] == [1.9, 3.3, 5.6, 8.1, 16.5], "5절 (4): 무작위 b 대 최악의 b"
assert G6_AA == 99246 and M6 == 3, "5절 (6): 진폭 증폭 포함 2큐비트 게이트 수"
assert [f[10] for f in FAM] == [46, 64, 74, 190, 288, 298, 766, 1216, 1266] and round(FAM[1][7]) == 448, "5절 (7): nnz(저장된 0 제외)와 2D Poisson N = 16 의 CG 비용"
assert [tuple(round(v, 2) for v in FAM_SLOPE[n_]) for n_ in ("1D Poisson", "2D Poisson", "I + γS (γ = 0.2)")] == [(1.96, 2.46, 1.99), (0.91, 1.41, 1.58), (0.09, 0.16, 1.06)], "5절 (7): N 지수"
assert round(np.pi / 2 * r8["tot"]) == 1110 and abs(mval - (-0.0361)) < 5e-5, "5절 (8): 읽어 내기"
print(f"\n   (실행 시간: 총 {time.time()-t_start:.0f} s)")
