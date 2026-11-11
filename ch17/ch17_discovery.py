"""17장 발견 절차 — numpy 검증 스크립트
문제 → 구조 → 원시 요소 → 자원 추정 → 고전 베이스라인 의 다섯 단계 체크리스트를
(A) 고전 베이스라인 검사(얽힘으로 MPS 비용 재기), (B) 큐비트화 QPE 자원 계산기, (C) 제곱·고차 가속의 손익분기 계산기로
실제 숫자로 돌려 본다. 규약: big-endian. 난수를 쓰지 않아 결과는 결정론적이다.
"""
import os, time
import numpy as np
from scipy.sparse import kron as skron, identity as sid, csr_matrix
from scipy.sparse.linalg import eigsh, expm_multiply
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
matplotlib.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

T0 = time.time()
X = np.array([[0, 1], [1, 0]], float); Z = np.diag([1.0, -1.0])


def op_at(o, i, n):
    mats = [sid(2, format="csr") for _ in range(n)]; mats[i] = csr_matrix(o)
    out = mats[0]
    for m_ in mats[1:]: out = skron(out, m_, format="csr")
    return out


def tfim(n, J=1.0, h=1.0, hz=0.0):
    """열린 사슬 TFIM H = −J Σ Z_i Z_{i+1} − h Σ X_i (12~16장과 같은 모형), 희소 행렬.
    hz ≠ 0 이면 세로장 −hz Σ Z_i 를 더한 비적분(mixed-field) Ising."""
    H = sum(-J * op_at(Z, i, n) @ op_at(Z, i + 1, n) for i in range(n - 1))
    H = H + sum(-h * op_at(X, i, n) for i in range(n))
    return H + sum(-hz * op_at(Z, i, n) for i in range(n)) if hz else H


def schmidt(psi, n, cut):
    s = np.linalg.svd(psi.reshape(2 ** cut, 2 ** (n - cut)), compute_uv=False)
    return s ** 2


def entropy_bits(p):
    p = p[p > 1e-15]; return float(-(p * np.log2(p)).sum())


def chi_needed(p, infid=1e-6):
    """슈미트 가중치 p(내림차순)에서 버린 가중치 ≤ infid 가 되는 최소 결합 차원(가운데 자름 하나)."""
    c = np.cumsum(np.sort(p)[::-1]); return int(np.searchsorted(c, 1 - infid) + 1)


# ======================================================================
# (A) 5단계 '고전 베이스라인': 1D 바닥상태(면적 법칙) vs 급랭(quench) 동역학(부피 법칙 성장)
print("A) 고전 베이스라인 검사 — 가운데 자름의 얽힘 엔트로피와 필요한 MPS 결합 차원 χ (버린 가중치 ≤ 1e-6)")
gs_rows = []
for n in (6, 8, 10, 12, 14, 16, 18):
    e0, v = eigsh(tfim(n), k=1, which="SA")
    p = schmidt(v[:, 0], n, n // 2)
    gs_rows.append((n, float(e0[0]), entropy_bits(p), chi_needed(p)))
    print(f"   바닥상태 n = {n:2d}: E0 = {e0[0]:11.6f}, S = {entropy_bits(p):.4f} 비트, χ(1e-6) = {chi_needed(p):3d}")
ns_g = np.array([r[0] for r in gs_rows]); S_g = np.array([r[2] for r in gs_rows])
slope_g = np.polyfit(np.log2(ns_g), S_g, 1)[0]
print(f"   S vs log2 n 기울기 = {slope_g:.3f} (열린 사슬 임계점 이론값 c/6 = 1/12 = 0.083, c = 1/2)")

n = 16; Hs = tfim(n)
psi = np.zeros(2 ** n, complex); psi[0] = 1.0                       # |0…0⟩ 에서 급랭
q_rows = []
ts = np.arange(0, 6.01, 0.5)
states = expm_multiply(-1j * Hs, psi, start=0, stop=ts[-1], num=len(ts), endpoint=True)
for t, ph in zip(ts, states):
    p = schmidt(ph, n, n // 2)
    q_rows.append((t, entropy_bits(p), chi_needed(p)))
print(f"   급랭 동역학 n = {n} (|0…0⟩, J = h = 1):")
for t, S, c in q_rows:
    print(f"     t = {t:3.1f}: S = {S:.3f} 비트, χ(1e-6) = {c:4d}")
tq = np.array([r[0] for r in q_rows]); Sq = np.array([r[1] for r in q_rows])
lin = (tq >= 1) & (tq <= 3)
vS = np.polyfit(tq[lin], Sq[lin], 1)[0]
print(f"   t = 1~3 구간 S 증가 속도 ≈ {vS:.3f} 비트/단위시간 (선형 성장 → χ ~ 2^S 는 t 에 지수)  |  최대 S = n/2 = {n//2}")

# 주의: 1D TFIM 은 자유 페르미온 모형이라 위 급랭은 2n×2n 행렬로도 정확히 풀린다(18장). 위 표는 "MPS 의 비용"이다.
# 고전으로 어렵다는 근거는 비적분 모형에서 재야 한다: Kim–Huse(2013) 매개변수의 mixed-field Ising, 같은 n = 16.
Hni = tfim(n, J=1.0, h=0.9045, hz=0.8090)
y1 = np.array([1, 1j]) / np.sqrt(2); psi_y = np.array([1.0 + 0j])
for _ in range(n): psi_y = np.kron(psi_y, y1)                      # |+y…+y⟩
ni_rows, conf_rows = [], []
ts_ni = np.arange(0, 6.01, 1.0)
for init, rows in ((psi_y, ni_rows), (psi, conf_rows)):
    for t, ph in zip(ts_ni, expm_multiply(-1j * Hni, init, start=0, stop=ts_ni[-1], num=len(ts_ni), endpoint=True)):
        p = schmidt(ph, n, n // 2); rows.append((t, entropy_bits(p), chi_needed(p)))
print(f"   비적분 급랭 n = {n} (J = 1, h = 0.9045, 세로장 0.8090):")
print("     |+y…+y⟩ 에서: " + ", ".join(f"t={t:.0f}: S={S:.3f} χ={c}" for t, S, c in ni_rows))
print("     |0…0⟩ 에서 : " + ", ".join(f"t={t:.0f}: S={S:.3f} χ={c}" for t, S, c in conf_rows))
vni = np.polyfit([r[0] for r in ni_rows[1:5]], [r[1] for r in ni_rows[1:5]], 1)[0]
print(f"     |+y…+y⟩: t = 1~4 구간 S 증가 속도 ≈ {vni:.3f} 비트/단위시간;  |0…0⟩: 최대 S = {max(r[1] for r in conf_rows):.3f} 비트(가둠, confinement)")

# ======================================================================
# (B) 4단계 '자원 추정': 큐비트화 + QPE 로 바닥상태 에너지 (16장 부품의 조립)
#   LCU 항 L 개, α = Σ|c|. W 한 번의 Toffoli ≈ SELECT(L − 1, unary iteration) + PREP 2회(일반 계수: QROM·별칭 표본 ≈ L 씩)
#   + 반사(m − 1). 정밀도 ε 의 QPE: 위상을 가장 가까운 격자점으로 반올림하는 모형에서 오차 ≤ π/2^b → 에너지 오차 ≤ α·π/2^b
#   → W 호출 ≈ 2^b ≥ πα/ε. (표준 QPE 는 가장 가까운 두 격자점 중 하나를 확률 ≥ 8/π² 로 준다 — 보장하려면 1비트 + 반복)
print("\nB) 자원 추정 계산기 — 큐비트화 QPE (16장 부품) 의 Toffoli 수와 논리 큐비트 (근사 모델, 상수 2배 안팎)")


def qpe_qubitization_cost(L, alpha, n_sys, eps, mu=10):
    m = int(np.ceil(np.log2(L)))
    b = int(np.ceil(np.log2(np.pi * alpha / eps)))
    calls = 2 ** b
    tof_W = (L - 1) + 2 * L + (m - 1)
    qubits = n_sys + m + mu + m + 1                          # 시스템 + 색인 + keep 비트 μ + alt 색인 + 위상 1(반복 QPE); unary iteration 보조 m−1·균등 비교 레지스터는 빼서 하한 쪽
    return dict(m=m, b=b, calls=calls, tof_W=tof_W, tof=calls * tof_W, qubits=qubits)


# 공식 검산(반올림 모형): n = 4..8 TFIM 의 정확한 E0 로 위상 θ0 = arccos(E0/α) 를 b 비트로 반올림했을 때의 에너지 오차 ≤ ε 인가
print("   공식 검산 (TFIM, ε = 0.01): b = ⌈log2(πα/ε)⌉ 비트 반올림 오차")
for n in (4, 6, 8):
    e0 = eigsh(tfim(n), k=1, which="SA")[0][0]
    alpha = (n - 1) + n; eps = 0.01
    c = qpe_qubitization_cost(2 * n - 1, alpha, n, eps)
    th = np.arccos(e0 / alpha); M = 2 ** c["b"]
    th_hat = np.round(th / (2 * np.pi) * M) / M * 2 * np.pi
    print(f"     n = {n}: α = {alpha}, b = {c['b']} (W 호출 {c['calls']}), |ΔE| = {abs(alpha*np.cos(th_hat) - e0):.2e} ≤ ε = {eps}: {abs(alpha*np.cos(th_hat) - e0) <= eps}")

models = []
for name, n_sys, L, alpha in [
    ("TFIM 1D n=100", 100, 199, 199.0),
    ("TFIM 1D n=1000", 1000, 1999, 1999.0),
    ("TFIM 2D 10×10", 100, 2 * 10 * 9 + 100, 2 * 10 * 9 + 100.0),
    ("Heisenberg 2D 10×10", 100, 3 * 2 * 10 * 9, 3 * 2 * 10 * 9.0),
    ("Heisenberg 2D 20×20", 400, 3 * 2 * 20 * 19, 3 * 2 * 20 * 19.0),
]:
    eps = 1e-3 * n_sys                                        # 에너지 밀도 정밀도 10⁻³ (사이트당)
    c = qpe_qubitization_cost(L, alpha, n_sys, eps)
    t_run = c["tof"] * 10e-6                                  # 가정: 논리 Toffoli 10 µs (직렬)
    models.append((name, n_sys, L, alpha, eps, c, t_run))
    print(f"   {name:20s}: L = {L:5d}, α = {alpha:7.0f}, ε = {eps:5.2f} → b = {c['b']:2d}, W 호출 {c['calls']:6d}, "
          f"Toffoli/W {c['tof_W']:5d}, 총 Toffoli {c['tof']:.2e}, 논리 큐비트 ≈ {c['qubits']}, 실행 ≈ {t_run/3600:.1f} 시간 (Toffoli 10 µs)")
print("   → 에너지 밀도 정밀도가 고정이면 α/ε 가 n 에 무관해 W 호출 수는 일정하고, 비용은 W 하나의 크기(∝ L)로 늘어난다")

# ======================================================================
# (C) 5단계 '고전 베이스라인'의 시간 비교: N^{1/p} 양자 호출 vs N 고전 평가 의 손익분기
#   고전: N·t_c,  양자: N^{1/p}·C·τ  (C = 오라클 한 번의 Toffoli 수, τ = 논리 Toffoli 시간)
#   N* = (Cτ/t_c)^{p/(p−1)},  손익분기 실행 시간 T* = N*·t_c
print("\nC) 손익분기 계산기 (Babbush et al. 2021 의 논리를 그대로 숫자로)")
be_rows = []
for C_, tau in ((1e3, 1e-6), (1e4, 1e-5)):
    for tc, label in ((1e-8, "1코어 10 ns"), (1e-11, "1000코어 → 10 ps")):
        for p_ in (2, 3, 4):
            Nstar = (C_ * tau / tc) ** (p_ / (p_ - 1)); Tstar = Nstar * tc
            # 수치 검산: N t_c = N^{1/p} C τ 를 이분법으로
            lo, hi = 0.0, 300.0
            for _ in range(200):
                mid = (lo + hi) / 2
                if 10 ** mid * tc < 10 ** (mid / p_) * C_ * tau: lo = mid
                else: hi = mid
            be_rows.append((C_, tau, tc, label, p_, Nstar, Tstar, 10 ** lo))
            print(f"   C = {C_:.0e} Toffoli, τ = {tau*1e6:4.0f} µs, 고전 {label:14s}, p = {p_}: N* = {Nstar:.2e} (이분법 {10**lo:.2e}), T* = {Tstar:.2e} s = {Tstar/86400:.2e} 일")

# ---------- 그림 ----------
fig, axs = plt.subplots(2, 2, figsize=(12, 8.4))
ax = axs[0, 0]
ax.plot(ns_g, S_g, "o-", color="#4C6EF5", label="TFIM 바닥상태 S (가운데 자름)")
ax.plot(ns_g, S_g[0] + slope_g * (np.log2(ns_g) - np.log2(ns_g[0])), "k:", lw=1, label=f"적합 {slope_g:.3f}·log₂n (이론 1/12 = 0.083)")
ax.set_xscale("log", base=2); ax.set_xticks(ns_g); ax.set_xticklabels([str(v) for v in ns_g]); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xlabel("사슬 길이 n"); ax.set_ylabel("얽힘 엔트로피 S (비트)"); ax.set_title("(1) 1D 임계 바닥상태: S ∝ log n (χ = 4~6)"); ax.legend(fontsize=8)

ax = axs[0, 1]
ax.semilogy(tq, [r[2] for r in q_rows], "s-", color="#F03E3E", label="TFIM 급랭 χ(1e-6), n = 16 (자유 페르미온 — 18장)")
ax.semilogy([r[0] for r in ni_rows], [r[2] for r in ni_rows], "o-", color="#343A40", label="비적분 Ising 급랭 χ, |+y…+y⟩")
ax.semilogy([r[0] for r in conf_rows], [r[2] for r in conf_rows], "o:", color="#868E96", label="비적분 Ising 급랭 χ, |0…0⟩ (가둠)")
ax.axhline(gs_rows[5][3], color="#4C6EF5", ls="--", lw=1, label=f"바닥상태 χ(1e-6), n = 16: {gs_rows[5][3]}")
ax.set_xlabel("시간 t"); ax.set_ylabel("필요 결합 차원 χ (로그)"); ax.set_title("(2) 동역학은 MPS의 χ 가 t 에 지수 — 모형·초기 상태에 따라 다르다"); ax.legend(fontsize=7)

ax = axs[1, 0]
names = [r[0] for r in models]; tofs = [r[5]["tof"] for r in models]
ax.barh(range(len(models)), tofs, color=["#74C0FC", "#74C0FC", "#868E96", "#4C6EF5", "#4C6EF5"])
ax.set_yticks(range(len(models))); ax.set_yticklabels(names, fontsize=9); ax.invert_yaxis()
ax.set_xscale("log"); ax.set_xlim(1e6, 1e9); ax.set_xlabel("총 Toffoli 수 (큐비트화 QPE, 에너지 밀도 10⁻³)")
ax.set_title("(3) 자원 추정: 비용은 W 크기(∝ 항 수 L)로 는다")

ax = axs[1, 1]
for (C_, tau, tc, lab), col in zip([(1e3, 1e-6, 1e-8, "C=10³, τ=1µs vs 1코어"), (1e3, 1e-6, 1e-11, "C=10³, τ=1µs vs 1000코어"),
                                    (1e4, 1e-5, 1e-8, "C=10⁴, τ=10µs vs 1코어"), (1e4, 1e-5, 1e-11, "C=10⁴, τ=10µs vs 1000코어")],
                                   ("#74C0FC", "#4C6EF5", "#F59F00", "#F03E3E")):
    rows = [r for r in be_rows if r[0] == C_ and r[1] == tau and r[2] == tc]
    ax.semilogy([r[4] for r in rows], [r[6] for r in rows], "o-", color=col, label=lab)
for yv, lab in ((86400, "1일"), (86400 * 365, "1년")):
    ax.axhline(yv, color="#868E96", ls=":", lw=1); ax.text(3.9, yv * 1.4, lab, fontsize=8, color="#868E96", ha="right")
ax.set_xticks([2, 3, 4]); ax.set_xticklabels(["제곱 (p=2)", "3차 (p=3)", "4차 (p=4)"])
ax.set_ylabel("손익분기 실행 시간 T* (s, 로그)"); ax.set_title("(4) 손익분기 시간: 가정에 따라 100초~32년, 차수가 높을수록 짧다"); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "ch17_results.png"), dpi=120)
print("\n6) saved ch17_results.png")
print(f"   (실행 시간 {time.time()-T0:.1f} s)")
