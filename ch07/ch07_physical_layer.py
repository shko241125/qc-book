"""7장 물리층 한눈에 — 게이트 = 해밀토니안의 시간 진화, T1/T2 감쇠, Lindblad 방정식을 numpy로 검증한다.
실행: python ch07_physical_layer.py
환경: Python 3.13, numpy 2.5 (SDK 없이 동작). 단위: ħ=1, 시간 ns, 주파수 GHz·2π.
"""
import os
import numpy as np
from scipy.linalg import expm

np.set_printoptions(precision=4, suppress=True, linewidth=120)
rng = np.random.default_rng(seed=7)

I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
ket0 = np.array([1, 0], dtype=complex); ket1 = np.array([0, 1], dtype=complex)
sm = np.array([[0, 1], [0, 0]], dtype=complex)        # σ- = |0><1| : 감쇠 연산자 (|1> → |0>)
def Rx(t): return expm(-1j * t/2 * X)
def fidelity_op(U, V):                                  # 과정 충실도(2큐비트 아님, 단일): |Tr(U†V)|²/d²
    return abs(np.trace(U.conj().T @ V))**2 / 4
def bloch(rho): return np.real([np.trace(rho @ X), np.trace(rho @ Y), np.trace(rho @ Z)])

# ── 1. 게이트 = 해밀토니안의 시간 진화: Rabi 진동 ─────────────────────
# 회전 틀(rotating frame)에서 공진 구동: H = (Ω/2) X.  U(t) = exp(-i H t) = Rx(Ω t)
Omega = 2*np.pi * 0.05        # Rabi 주파수 50 MHz → 2π·0.05 rad/ns
ts = np.linspace(0, 40, 401)  # ns
p1 = []
for t in ts:
    U = expm(-1j * (Omega/2) * X * t)
    psi = U @ ket0
    p1.append(abs(psi[1])**2)
p1 = np.array(p1)
t_pi = np.pi / Omega          # π 펄스 길이
U_pi = expm(-1j * (Omega/2) * X * t_pi)
print(f"1) Rabi 진동: Ω/2π = 50 MHz → π 펄스 = {t_pi:.1f} ns.  P(1) 최대 {p1.max():.4f}, 처음 1 에 닿는 t = {ts[np.argmax(p1 > 1-1e-9)]:.1f} ns")
print("   U(t_π) =\n", U_pi, "\n   = -i·X ?", np.allclose(U_pi, -1j*X), "  (전역 위상 -i 를 빼면 X 게이트)")
U_2pi = expm(-1j * (Omega/2) * X * 2*t_pi)
print("   U(2t_π) = -I ?", np.allclose(U_2pi, -I), " (스핀-1/2 의 2π 회전 = -1, 연습문제 1)")
U_half = expm(-1j * (Omega/2) * X * t_pi/2)
print("   U(t_π/2) = Rx(π/2) ?", np.allclose(U_half, Rx(np.pi/2)))
assert abs(t_pi - 10.0) < 1e-9, "5절 Rabi: Ω/2π = 50 MHz 에서 π 펄스 = 10.0 ns"
assert abs(p1.max() - 1) < 1e-9 and abs(ts[np.argmax(p1 > 1-1e-9)] - 10.0) < 1e-9, "5절 Rabi: P(1) 이 10 ns 에서 정확히 1 에 닿는다"
assert np.allclose(U_pi, -1j*X), "5절 펄스 표: U(t_π) = −i·X"
assert np.allclose(U_2pi, -I), "5절 펄스 표·8절 문제 1: U(2t_π) = −I"
assert np.allclose(U_half, Rx(np.pi/2)), "5절 펄스 표: U(t_π/2) = Rx(π/2)"

# 구동 축을 바꾸면 Ry: H = (Ω/2)(cos φ X + sin φ Y)  — 펄스의 위상 φ 가 회전축을 정한다
phi = np.pi/2
U_y = expm(-1j * (Omega/2) * (np.cos(phi)*X + np.sin(phi)*Y) * t_pi)
print("   구동 위상 φ=π/2 인 π 펄스 = -i·Y ?", np.allclose(U_y, -1j*Y), " → 펄스 위상 = 회전축 (가상 Rz 의 근거, 4장)")
assert np.allclose(U_y, -1j*Y), "5절 펄스 표: 위상 φ = π/2 인 π 펄스 = −i·Y"

# 이조(detuning) Δ 가 있으면: H = (Δ/2) Z + (Ω/2) X → 회전축이 기울고 P(1) 최대가 1 에 못 미친다
for Delta_MHz in (0, 25, 50):
    Delta = 2*np.pi * Delta_MHz * 1e-3
    Hd = (Delta/2)*Z + (Omega/2)*X
    pmax = max(abs((expm(-1j*Hd*t) @ ket0)[1])**2 for t in ts)
    print(f"   이조 Δ/2π = {Delta_MHz:>2d} MHz: P(1) 최대 = {pmax:.4f}  (이론 Ω²/(Ω²+Δ²) = {Omega**2/(Omega**2+Delta**2):.4f})")
    assert abs(pmax - {0: 1.000, 25: 0.800, 50: 0.500}[Delta_MHz]) < 5e-4, "5절 이조: P(1) 최대 = 1.000, 0.800, 0.500"
    assert abs(pmax - Omega**2/(Omega**2+Delta**2)) < 5e-5, "5절 이조: 이론값 Ω²/(Ω²+Δ²) 과 네 자리 일치"

# ── 2. Lindblad 방정식: T1 (진폭 감쇠), T2 (위상 감쇠) ────────────────
def lindblad_rhs(rho, Hs, cops):
    d = -1j * (Hs @ rho - rho @ Hs)
    for L in cops:
        LdL = L.conj().T @ L
        d += L @ rho @ L.conj().T - 0.5 * (LdL @ rho + rho @ LdL)
    return d
def mesolve(rho0, Hs, cops, ts):
    """4차 Runge–Kutta 로 Lindblad 방정식 적분 (QuTiP mesolve 의 최소 대체)"""
    rho = rho0.copy(); out = [rho.copy()]
    for i in range(1, len(ts)):
        h = ts[i] - ts[i-1]
        k1 = lindblad_rhs(rho, Hs, cops); k2 = lindblad_rhs(rho + h/2*k1, Hs, cops)
        k3 = lindblad_rhs(rho + h/2*k2, Hs, cops); k4 = lindblad_rhs(rho + h*k3, Hs, cops)
        rho = rho + h/6*(k1 + 2*k2 + 2*k3 + k4); out.append(rho.copy())
    return np.array(out)

T1, Tphi = 100.0, 150.0                    # μs 단위로 생각해도 되고 ns 로 생각해도 된다 — 비율만 중요
gamma1 = 1/T1; gamma_phi = 1/Tphi
T2 = 1/(1/(2*T1) + 1/Tphi)                 # 1/T2 = 1/(2T1) + 1/Tφ
tt = np.linspace(0, 300, 601)
# (a) |1> 에서 시작, T1 감쇠만
rho_1 = np.outer(ket1, ket1.conj())
r_T1 = mesolve(rho_1, 0*Z, [np.sqrt(gamma1)*sm], tt)
p1_t = np.real(r_T1[:, 1, 1])
# (b) |+> 에서 시작, T1 + 순수 위상감쇠
plus = H @ ket0; rho_plus = np.outer(plus, plus.conj())
r_T2 = mesolve(rho_plus, 0*Z, [np.sqrt(gamma1)*sm, np.sqrt(gamma_phi/2)*Z], tt)
coh_t = np.abs(r_T2[:, 0, 1])              # 비대각 원소 |ρ01| = Bloch 벡터 x 성분 / 2
idx = np.searchsorted(tt, T1)
print(f"\n2) Lindblad: T1={T1:.0f}, Tφ={Tphi:.0f} → T2 = {T2:.1f}  (T2 ≤ 2T1 = {2*T1:.0f} ✓, 연습문제 2)")
print(f"   (a) |1> 의 P(1) at t=T1: {p1_t[idx]:.4f}  (이론 e^-1 = {np.exp(-1):.4f})")
idx2 = np.searchsorted(tt, T2)
print(f"   (b) |+> 의 결맞음 |ρ01| at t=T2: {coh_t[idx2]:.4f}  (이론 0.5·e^-1 = {0.5*np.exp(-1):.4f})")
print("   (b) t=300 의 Bloch 벡터:", bloch(r_T2[-1]), " → x 가 사라지고 z 가 +1 로 (|0> 로 떨어짐)")
assert abs(T2 - 85.7) < 5e-2 and T2 <= 2*T1 == 200, "5절 T1·T2 감쇠: T2 = 85.7 (≤ 2T1 = 200)"
assert abs(p1_t[idx] - 0.3679) < 5e-5 and abs(np.exp(-1) - 0.3679) < 5e-5, "5절 T1·T2 표: |1> 의 P(1) = 0.3679 = e⁻¹ (t = T1)"
assert tt[idx2] == 86.0 and abs(coh_t[idx2] - 0.1833) < 5e-5 and abs(0.5*np.exp(-1) - 0.1839) < 5e-5, "5절 T1·T2 표: |+> 의 |ρ01| = 0.1833 (격자점 t = 86), 이론 0.1839"
assert np.allclose(bloch(r_T2[-1]), (0.030, 0, 0.950), atol=5e-4), "5절 T1·T2 표: t = 300 의 Bloch 벡터 (0.030, 0, 0.950)"

# ── 3. Kraus 표현: 진폭 감쇠·위상 감쇠 채널 = Lindblad 의 유한 시간 해 ──
def amp_damp(rho, g):
    K0 = np.array([[1, 0], [0, np.sqrt(1-g)]]); K1 = np.array([[0, np.sqrt(g)], [0, 0]])
    return K0 @ rho @ K0.conj().T + K1 @ rho @ K1.conj().T
def phase_damp(rho, lam):
    K0 = np.array([[1, 0], [0, np.sqrt(1-lam)]]); K1 = np.array([[0, 0], [0, np.sqrt(lam)]])
    return K0 @ rho @ K0.conj().T + K1 @ rho @ K1.conj().T
t = 70.0; i_t = np.searchsorted(tt, t)
g = 1 - np.exp(-gamma1*t)
print(f"\n3) Kraus 검산 at t={t:.0f}: 진폭감쇠 γ = 1-e^(-t/T1) = {g:.4f}")
print("   amp_damp(|1><1|) == Lindblad 해?", np.allclose(amp_damp(rho_1, g), r_T1[i_t], atol=1e-6))
# 위상감쇠 + 진폭감쇠 합성: 비대각 원소가 e^{-t/T2} 로 줄어드는지
lam = 1 - np.exp(-2*gamma_phi*t)          # 위상감쇠 Kraus 의 λ: 비대각이 √(1-λ) 배 → √(1-λ) = e^(-t/Tφ)
rho_k = phase_damp(amp_damp(rho_plus, g), lam)
print("   phase∘amp (|+><+|) == Lindblad 해?", np.allclose(rho_k, r_T2[i_t], atol=1e-6),
      f"  |ρ01| = {abs(rho_k[0,1]):.4f} vs 0.5·e^(-t/T2) = {0.5*np.exp(-t/T2):.4f}")
# 연습문제 3: 위상감쇠 λ=1/2 후 |+> 의 순도
rho_e3 = phase_damp(rho_plus, 0.5)
print(f"   연습문제 3: 위상감쇠 λ=1/2 후 |+> 의 ρ =\n{rho_e3}\n   순도 Tr ρ² = {np.real(np.trace(rho_e3 @ rho_e3)):.4f}, Bloch 길이 = {np.linalg.norm(bloch(rho_e3)):.4f}")
assert abs(g - 0.5034) < 5e-5, "5절 Kraus 검산: t = 70 의 진폭 감쇠 γ = 0.5034"
assert np.allclose(amp_damp(rho_1, g), r_T1[i_t], atol=1e-6), "5절 Kraus 검산: 진폭 감쇠 = Lindblad 해 (10⁻⁶ 안)"
assert np.allclose(rho_k, r_T2[i_t], atol=1e-6), "5절 Kraus 검산: 위상∘진폭 감쇠 = Lindblad 해 (10⁻⁶ 안)"
assert abs(abs(rho_k[0,1]) - 0.2210) < 5e-5 and abs(0.5*np.exp(-t/T2) - 0.2210) < 5e-5, "5절 Kraus 검산: |ρ01| = 0.2210 = 0.5e^{−t/T2}"
assert abs(abs(phase_damp(amp_damp(rho_plus, g), 1 - np.exp(-gamma_phi*t))[0, 1]) - 0.279) < 5e-4, "10절 기록: λ 를 1 − e^{−t/Tφ} 로 두면 비대각 0.279 (vs 0.221)"
assert np.allclose(rho_e3, [[0.5, 0.3536], [0.3536, 0.5]], atol=5e-5), "8절 문제 3: ρ = [[0.5, 0.3536], [0.3536, 0.5]]"
assert abs(np.real(np.trace(rho_e3 @ rho_e3)) - 0.75) < 1e-9 and abs(np.linalg.norm(bloch(rho_e3)) - 0.707) < 5e-4, "8절 문제 3: 순도 0.75, Bloch 길이 0.707"
rho_e3b = phase_damp(rho_e3, 0.5)                    # 같은 채널을 한 번 더 (출력 없음, 검사만)
assert abs(rho_e3b[0, 1] - 0.25) < 1e-9 and abs(np.linalg.norm(bloch(rho_e3b)) - 0.5) < 1e-9 \
    and abs(np.real(np.trace(rho_e3b @ rho_e3b)) - 0.625) < 1e-9, "8절 문제 3: 한 번 더 걸면 ρ01 = 0.25, Bloch 길이 0.5, 순도 0.625"

# ── 4. 게이트 충실도 vs 펄스 길이: 느린 게이트는 T1/T2 에 먹힌다 ─────
# π 펄스를 길이 t_g 로 늘리며(Ω = π/t_g) Lindblad 로 적분, 결과를 이상적 X 와 비교
def avg_gate_fidelity_X(t_g, T1, Tphi, nsteps=400):
    Om = np.pi / t_g
    Hs = (Om/2) * X
    cops = [np.sqrt(1/T1)*sm, np.sqrt(1/Tphi/2)*Z]
    ts_g = np.linspace(0, t_g, nsteps)
    # 평균 게이트 충실도 ≈ 6개 축 상태(±x,±y,±z)에 대한 상태 충실도 평균 (단일 큐비트 2-design)
    states = [ket0, ket1, H@ket0, H@ket1, (ket0+1j*ket1)/np.sqrt(2), (ket0-1j*ket1)/np.sqrt(2)]
    fs = []
    for s in states:
        rho0 = np.outer(s, s.conj()); rf = mesolve(rho0, Hs, cops, ts_g)[-1]
        target = X @ s
        fs.append(np.real(target.conj() @ rf @ target))
    return float(np.mean(fs))
print("\n4) X 게이트 평균 충실도 vs 펄스 길이 (T1=100, Tφ=150, 같은 단위)")
tgs = [1, 5, 10, 20, 40]
Fg = []
for tg in tgs:
    F = avg_gate_fidelity_X(tg, T1, Tphi); Fg.append(F)
    print(f"   t_g = {tg:>2d}: F = {F:.4f}   1-F = {1-F:.4f}")
for tg, F, F_md, inf_md, lin_md in zip(tgs, Fg, (0.9945, 0.9730, 0.9474, 0.9004, 0.8206), (0.0055, 0.0270, 0.0526, 0.0996, 0.1794),
                                       (0.0056, 0.0278, 0.0556, 0.1111, 0.2222)):
    assert abs(F - F_md) < 5e-5 and abs((1 - F) - inf_md) < 5e-5, "5절 게이트 충실도 표: F, 1 − F (t_g = 1, 5, 10, 20, 40)"
    assert abs((tg/3)*(1/T1 + 1/Tphi) - lin_md) < 5e-5, "5절 게이트 충실도 표: 1차 공식 (t_g/3)(1/T1 + 1/Tφ) 열"
assert abs((40/3)*(1/100e3 + 1/100e3) - 2.7e-4) < 5e-6, "5절: t_g = 40 ns, T1 = Tφ = 100 μs 에서 1 − F ≈ 2.7×10⁻⁴"
# 검산: "이상적 X 게이트 뒤에 Kraus 채널 한 번" 근사 (Aer·PennyLane 방식) vs Lindblad 적분, t_g = 40
def avg_fid_kraus_after_gate(t_g, T1, Tphi):
    g = 1 - np.exp(-t_g/T1); lam = 1 - np.exp(-2*t_g/Tphi)
    states = [ket0, ket1, H@ket0, H@ket1, (ket0+1j*ket1)/np.sqrt(2), (ket0-1j*ket1)/np.sqrt(2)]
    fs = []
    for s_ in states:
        target = X @ s_
        rho_after = phase_damp(amp_damp(np.outer(target, target.conj()), g), lam)
        fs.append(np.real(target.conj() @ rho_after @ target))
    return float(np.mean(fs))
Fk = avg_fid_kraus_after_gate(40, T1, Tphi)
print(f"   게이트 뒤 Kraus 근사 (t_g=40): 1-F = {1-Fk:.5f}  vs Lindblad 적분 1-F = {1-Fg[-1]:.5f}  (차이 {abs(Fk-Fg[-1]):.1e})")
assert abs((1 - Fk) - 0.17925) < 5e-6 and abs((1 - Fg[-1]) - 0.17937) < 5e-6, "6절: t_g = 40 에서 1−F 가 0.17925(Kraus) vs 0.17937(Lindblad)"
assert abs(abs(Fk - Fg[-1]) - 1.2e-4) < 5e-6, "6절: 두 방식의 차이 1.2×10⁻⁴"

# ── 5. 그림 ───────────────────────────────────────────────────────
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Noto Sans CJK JP"; plt.rcParams["axes.unicode_minus"] = False
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 3.4))
ax1.plot(ts, p1, color="#4C6EF5"); ax1.axvline(t_pi, ls="--", color="#868E96"); ax1.text(t_pi+0.5, 0.5, "π 펄스", color="#868E96")
ax1.set_xlabel("시간 (ns)"); ax1.set_ylabel("P(1)"); ax1.set_title("Rabi 진동: Ω/2π = 50 MHz")
ax2.plot(tt, p1_t, color="#4C6EF5", label="P(1), |1> 시작 (T1)"); ax2.plot(tt, 2*coh_t, color="#F03E3E", label="2|ρ01|, |+> 시작 (T2)")
ax2.axvline(T1, ls=":", color="#4C6EF5"); ax2.axvline(T2, ls=":", color="#F03E3E")
ax2.set_xlabel("시간"); ax2.set_title("T1 감쇠와 T2 위상감쇠"); ax2.legend(frameon=False, fontsize=8)
ax3.semilogy(tgs, 1-np.array(Fg), 'o-', color="#4C6EF5"); ax3.set_xlabel("π 펄스 길이 t_g"); ax3.set_ylabel("1 − F"); ax3.set_title("느린 게이트일수록 오류 ↑")
plt.tight_layout(); plt.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "ch07_results.png"), dpi=150); print("5) saved ch07_results.png")
