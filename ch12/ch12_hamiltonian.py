"""12장 해밀토니안 시뮬레이션 — Trotter 분해로 스핀 사슬을 진화시키고, QPE 로 바닥상태 에너지를, 진폭 추정으로 확률을 읽는다.
실행: python ch12_hamiltonian.py
환경: Python 3.13, numpy 2.5, scipy (expm). 큐비트 순서: 빅 엔디안(3장 규약). ħ = 1.
"""
import os
import numpy as np
from scipy.linalg import expm

np.set_printoptions(precision=4, suppress=True, linewidth=120)

I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
def kron(*ops):
    out = np.array([[1]], dtype=complex)
    for o in ops: out = np.kron(out, o)
    return out
def op_on(g, q, n): ops = [I]*n; ops[q] = g; return kron(*ops)
def op2_on(g1, q1, g2, q2, n): ops = [I]*n; ops[q1] = g1; ops[q2] = g2; return kron(*ops)
def ket0n(n): v = np.zeros(2**n, dtype=complex); v[0] = 1; return v

# ── 1. 모델: 횡자기장 이징 사슬  H = -J Σ Z_i Z_{i+1} - h Σ X_i  (열린 경계) ─────
def tfim_terms(n, J, h):
    """(계수, 파울리 문자열) 리스트 — Qiskit SparsePauliOp 와 같은 표현"""
    terms = []
    for i in range(n-1): terms.append((-J, ''.join('Z' if k in (i, i+1) else 'I' for k in range(n))))
    for i in range(n): terms.append((-h, ''.join('X' if k == i else 'I' for k in range(n))))
    return terms
PAULI = {'I': I, 'X': X, 'Y': Y, 'Z': Z}
def pauli_matrix(s): return kron(*[PAULI[c] for c in s])
def hamiltonian(terms): return sum(c * pauli_matrix(s) for c, s in terms)
n, J, h = 4, 1.0, 1.0
terms = tfim_terms(n, J, h); Hm = hamiltonian(terms)
evals, evecs = np.linalg.eigh(Hm)
print(f"1) TFIM n={n}, J={J}, h={h}: 항 {len(terms)}개, 바닥 에너지 E0 = {evals[0]:.6f}, 첫 들뜸 E1 = {evals[1]:.6f}, 스펙트럼 폭 = {evals[-1]-evals[0]:.4f}")
print("   ZZ 항끼리는 교환, X 항끼리는 교환, ZZ 와 X 는 비교환 → [H_ZZ, H_X] ≠ 0:",
      not np.allclose(hamiltonian(terms[:n-1]) @ hamiltonian(terms[n-1:]) - hamiltonian(terms[n-1:]) @ hamiltonian(terms[:n-1]), 0))

# ── 2. Trotter 분해: 1차·2차, 오차 vs 스텝 수 ─────────────────────────────
def exp_term(c, s, dt):
    """e^{-i c P dt} 를 게이트로: P 는 파울리 문자열. (= Qiskit PauliEvolutionGate 의 한 항)"""
    P = pauli_matrix(s); return np.cos(c*dt) * np.eye(len(P)) - 1j * np.sin(c*dt) * P   # P² = I 이므로
def trotter1(terms, t, steps):
    dt = t/steps; U = np.eye(2**n, dtype=complex)
    step = np.eye(2**n, dtype=complex)
    for c, s in terms: step = exp_term(c, s, dt) @ step
    for _ in range(steps): U = step @ U
    return U
def trotter2(terms, t, steps):
    dt = t/steps; step = np.eye(2**n, dtype=complex)
    for c, s in terms: step = exp_term(c, s, dt/2) @ step
    for c, s in reversed(terms): step = exp_term(c, s, dt/2) @ step
    U = np.eye(2**n, dtype=complex)
    for _ in range(steps): U = step @ U
    return U
t = 2.0
U_exact = expm(-1j * Hm * t)
psi0 = ket0n(n)                                                   # |0000> = 모든 스핀 up (Z 고유상태), X 장이 뒤섞는다
print(f"\n2) Trotter 오차 (t = {t}, 초기 |0000>): ‖U_trot − U_exact‖₂ 와 상태 충실도")
rows = []
for steps in (1, 2, 4, 8, 16, 32):
    U1 = trotter1(terms, t, steps); U2 = trotter2(terms, t, steps)
    e1 = np.linalg.norm(U1 - U_exact, 2); e2 = np.linalg.norm(U2 - U_exact, 2)
    f1 = abs(np.vdot(U_exact @ psi0, U1 @ psi0))**2; f2 = abs(np.vdot(U_exact @ psi0, U2 @ psi0))**2
    rows.append((steps, e1, e2, f1, f2))
    print(f"   steps={steps:>2d} (dt={t/steps:.3f}): 1차 오차 {e1:.4f}, 2차 오차 {e2:.5f};  충실도 1차 {f1:.4f}, 2차 {f2:.6f}")
r = np.array(rows)
print(f"   스텝 2배 → 1차 오차 비 {r[-2,1]/r[-1,1]:.2f} (이론 2), 2차 오차 비 {r[-2,2]/r[-1,2]:.2f} (이론 4)")
# 게이트 수: 항당 ZZ = CNOT·Rz·CNOT (2 CNOT), X 항 = Rx 1개
n_zz, n_x = n-1, n
print(f"   스텝당 게이트: 1차 = ZZ {n_zz}개×(CNOT 2 + Rz 1) + Rx {n_x} = CNOT {2*n_zz}, 단일 {n_zz+n_x};  2차 = CNOT {4*n_zz} (양 끝 항 합치면 {4*n_zz - 2*(n_zz>0)})")

# ── 3. 관측량 동역학: 자화 <Z_0>(t) — 정확 vs Trotter(steps 고정) ────────────
ts = np.linspace(0, 6, 61)
Z0 = op_on(Z, 0, n)
mz_exact = [np.real(np.vdot(expm(-1j*Hm*tt) @ psi0, Z0 @ (expm(-1j*Hm*tt) @ psi0))) for tt in ts]
mz_tr = {}                                                        # 고정 dt 로 순차 진화 (실제 회로가 하는 일)
for dt in (1.0, 0.25):
    step = np.eye(2**n, dtype=complex)
    for c, s_ in terms: step = exp_term(c, s_, dt/2) @ step
    for c, s_ in reversed(terms): step = exp_term(c, s_, dt/2) @ step
    psi = psi0.copy(); times = [0.0]; vals = [np.real(np.vdot(psi, Z0 @ psi))]
    for k in range(int(round(6/dt))):
        psi = step @ psi; times.append((k+1)*dt); vals.append(np.real(np.vdot(psi, Z0 @ psi)))
    mz_tr[dt] = (np.array(times), np.array(vals))
print(f"\n3) <Z_0>(t): t=6 에서 정확 {mz_exact[-1]:.4f}, 2차 Trotter dt=1 {mz_tr[1.0][1][-1]:.4f}, dt=0.25 {mz_tr[0.25][1][-1]:.4f}")

# ── 4. 바닥 에너지를 QPE 로: U = e^{+i(H+wI)τ}, φ = (E+w)τ/2π, 시험 상태의 겹침이 확률 ──
# 에너지를 [0, 1) 위상으로: φ = (E - E_min)/(E_max - E_min) 이 되도록 τ 와 오프셋을 잡는다 (스펙트럼 폭은 Gershgorin 으로 상한)
w = sum(abs(c) for c, _ in terms)                                   # ‖H‖ ≤ Σ|c_i| = 7
tau = 2*np.pi / (2*w) * 0.99                                       # 위상 φ = (E + w)·τ/2π ∈ (0, 1)
U = expm(+1j * (Hm + w*np.eye(2**n)) * tau)                         # 오프셋 w 를 더해 위상 φ = (E+w)τ/2π ∈ (0,1); 부호 +i 는 φ 를 양으로 읽기 위한 편의 (e^{-iH'τ} 의 역 — QPE 결과는 k ↔ 2^t-k 로 대응)
def dft_matrix(T):
    wq = np.exp(2j*np.pi/T); j = np.arange(T); return wq**np.outer(j, j)/np.sqrt(T)
def qpe_energy_probs(U, trial, tb):
    T = 2**tb; M = len(trial)
    psi = np.zeros((T, M), dtype=complex); psi[:] = trial/np.sqrt(T)
    Upows = {}; cur = U.copy()
    for j in reversed(range(tb)):                                   # q[j] 가 U^{2^{tb-1-j}} 제어 — 거듭제곱을 반복 제곱으로
        Upows[j] = cur; cur = cur @ cur
    idx = np.arange(T)
    for j in range(tb):
        mask = ((idx >> (tb-1-j)) & 1) == 1
        psi[mask] = psi[mask] @ Upows[j].T
    psi = dft_matrix(T).conj().T @ psi
    return (abs(psi)**2).sum(axis=1)
tb = 8
gs = evecs[:, 0]
# 시험 상태 1: 모든 스핀 |+> (X 장의 바닥상태, J→0 극한), 시험 상태 2: |0000>
trial_plus = kron(*([H]*n)) @ ket0n(n); trial_0 = ket0n(n)
for name, trial in (("|++++>", trial_plus), ("|0000>", trial_0)):
    gamma = abs(np.vdot(gs, trial))**2
    p = qpe_energy_probs(U, trial, tb)
    ks = np.argsort(p)[::-1][:3]
    E_from_k = [k/2**tb * 2*np.pi/tau - w for k in ks]
    print(f"\n4) QPE 에너지 (t_b={tb}, 시험 {name}, 바닥상태 겹침 γ={gamma:.3f}): 상위 피크 E ≈ {np.round(E_from_k, 3)} (확률 {np.round(p[ks], 3)});  정확 E0 = {evals[0]:.4f}, E1 = {evals[1]:.4f}")
    # 바닥 피크의 확률 질량 (±1 빈)
    k0 = int(round((evals[0] + w)*tau/(2*np.pi) * 2**tb)); mass = p[[(k0-1) % 2**tb, k0, (k0+1) % 2**tb]].sum()
    print(f"   E0 빈 k={k0} ±1 의 확률 {mass:.3f} ≈ γ (고유값 샘플러: 10장)  → 1/γ ≈ {1/gamma:.1f} 회 반복 또는 진폭 증폭 1/√γ ≈ {1/np.sqrt(gamma):.1f} 회")
print(f"   에너지 분해능 = 2π/(τ·2^{tb}) = {2*np.pi/(tau*2**tb):.4f}  (t_b 한 비트당 2배, 제어 U 는 2^{tb}-1 = {2**tb-1} 번)")

# ── 5. 진폭 추정: QPE(Grover 연산자) 로 a = sin²θ 를 읽는다 — 샷 노이즈 1/√N vs Heisenberg 1/N ──
# 상태 준비 A: n_a 큐비트, '좋은' 상태 = |1...1>, a = sin²θ 를 Ry 로 직접 설정
na = 2; theta_true = 0.3; a_true = np.sin(theta_true)**2
A = kron(np.array([[np.cos(theta_true), -np.sin(theta_true)], [np.sin(theta_true), np.cos(theta_true)]]), np.eye(2))   # q0 에 Ry(2θ): P(q0=1) = sin²θ
good = np.zeros(2**na, dtype=complex); good[2] = 1; good[3] = 1      # q0 = 1 인 상태들
S_f = np.eye(2**na) - 2*np.diag(abs(good))                          # 좋은 상태에 -1
S_0 = np.eye(2**na); S_0[0, 0] = -1
Q = -A @ S_0 @ A.conj().T @ S_f                                      # 9장 진폭 증폭 연산자, 고유값 e^{±2iθ}
def amp_est_probs(Q, A, tb):
    T = 2**tb; M = Q.shape[0]
    psi = np.zeros((T, M), dtype=complex); psi[:] = (A @ ket0n(na))/np.sqrt(T)
    Upows = {}; cur = Q.copy()
    for j in reversed(range(tb)): Upows[j] = cur; cur = cur @ cur
    idx = np.arange(T)
    for j in range(tb):
        mask = ((idx >> (tb-1-j)) & 1) == 1; psi[mask] = psi[mask] @ Upows[j].T
    psi = dft_matrix(T).conj().T @ psi
    return (abs(psi)**2).sum(axis=1)
print(f"\n5) 진폭 추정: a = sin²θ = {a_true:.4f} (θ = {theta_true})")
for tb in (3, 5, 7):
    p = amp_est_probs(Q, A, tb); k = int(np.argmax(p))
    theta_hat = np.pi * k / 2**tb                                   # 고유값 e^{±2iθ} → φ = ±θ/π → k/2^t = θ/π 또는 1 - θ/π
    a_hat = np.sin(theta_hat)**2
    oracle_calls = 2**tb - 1
    # 같은 수의 '오라클 호출'(A 실행) 을 고전 샘플링에 쓰면: 표준오차 √(a(1-a)/N)
    se_classical = np.sqrt(a_true*(1-a_true)/oracle_calls)
    bound = 2*np.pi*np.sqrt(a_true*(1-a_true))/2**tb + np.pi**2/4**tb       # BHMT 2002 정리 12: |â-a| ≤ 2π√(a(1-a))/M + π²/M², M = 2^t
    print(f"   t_b={tb}: 최빈 k={k:>3d} → â = {a_hat:.4f} (오차 {abs(a_hat-a_true):.4f}, P={p[k]:.3f});  Q 호출 {oracle_calls:>3d}회 → Heisenberg 오차 상한 {bound:.4f} vs 고전 {oracle_calls}샷 표준오차 {se_classical:.4f}")

# ── 6. 그림 ───────────────────────────────────────────────────────
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Noto Sans CJK JP"; plt.rcParams["axes.unicode_minus"] = False
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 3.4))
ax1.loglog(r[:, 0], r[:, 1], 'o-', color="#4C6EF5", label="1차 Trotter"); ax1.loglog(r[:, 0], r[:, 2], 's-', color="#F03E3E", label="2차 Trotter")
ax1.loglog(r[:, 0], r[-1, 1]*r[-1, 0]/r[:, 0], '--', color="#868E96", lw=0.8); ax1.loglog(r[:, 0], r[-1, 2]*(r[-1, 0]/r[:, 0])**2, ':', color="#868E96", lw=0.8)
ax1.set_xlabel("스텝 수"); ax1.set_ylabel("‖U_trot − U_exact‖"); ax1.set_title("Trotter 오차: 1/steps, 1/steps²"); ax1.legend(frameon=False, fontsize=8)
ax2.plot(ts, mz_exact, color="#1f2328", label="정확"); ax2.plot(*mz_tr[1.0], 's--', ms=4, color="#F03E3E", label="2차, dt=1"); ax2.plot(*mz_tr[0.25], '.:', ms=4, color="#4C6EF5", label="2차, dt=0.25")
ax2.set_xlabel("t"); ax2.set_ylabel("<Z_0>"); ax2.set_title("TFIM 자화 동역학 (n=4)"); ax2.legend(frameon=False, fontsize=8)
p = qpe_energy_probs(U, trial_plus, 8); Es = np.arange(2**8)/2**8 * 2*np.pi/tau - w
ax3.plot(Es, p, color="#4C6EF5");
for e in evals: ax3.axvline(e, ls=":", color="#868E96", lw=0.6)
ax3.set_xlabel("E (k 로부터 환산)"); ax3.set_ylabel("확률"); ax3.set_title("QPE 에너지 분포, 시험 |++++>"); ax3.set_xlim(evals[0]-0.5, evals[-1]+0.5)
plt.tight_layout(); plt.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "ch12_results.png"), dpi=150); print("6) saved ch12_results.png")

# ── 게이트 비용 기록 (2큐비트 게이트 / 깊이 / T-count) ──
# 위에서 행렬로 계산한 Trotter·QPE·진폭 추정 회로를 게이트 리스트로 다시 짜서 세고, 위의 행렬·분포와 같은지(전역 위상 제외) 확인한다.
# 규약(tools/gate_cost.py): 기저 CX + 1큐비트 게이트, 전결합. exp(−iθ Z⊗Z) = 2 CX + Rz(2θ), 제어 Rz/Ry/Rx = 2 CX + 회전 2,
# CP(θ) = 2 CX + 위상 회전 3 (CP(π) = CZ), SWAP = 3 CX, CCZ = 6 CX + 7 T. 깊이 = 게이트를 가능한 한 앞에 놓은(ASAP) 층 수.
# π/4 의 배수가 아닌 회전은 개수 R 로 세고 T ≈ 3·log2(1/ε)·R (ε = 1e-6).
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools"))
import gate_cost as gc
_g = gc.G
CHAIN = [0, 1, 2]; BRICK = [0, 2, 1]                               # ZZ 결합 (i, i+1) 을 거는 순서: 본문(01, 12, 23) / 짝·홀(01, 23, 12)
def zz_layer(dt_, order): return [_g("rzz", i, i+1, angle=-J*dt_) for i in order]      # e^{-i(-J) Z_i Z_{i+1} dt}
def x_layer(dt_): return [_g("rx", q, angle=-2*h*dt_) for q in range(n)]                # e^{-i(-h) X dt} = Rx(-2h dt)
def trotter1_list(t_, steps, order=CHAIN): return (zz_layer(t_/steps, order) + x_layer(t_/steps)) * steps
def trotter2_list(t_, steps, order=CHAIN):                         # 본문 trotter2 와 같은 순서: 항 순서대로 dt/2, 역순으로 dt/2
    d = t_/steps; return (zz_layer(d/2, order) + x_layer(d/2) + x_layer(d/2)[::-1] + zz_layer(d/2, order[::-1])) * steps
def trotter2_grouped_list(t_, steps, order=BRICK):                 # ZZ 끼리 교환 → ZZ 묶음 전체를 이웃 스텝과 합친 꼴: ZZ(dt/2) [X(dt) ZZ(dt)]^{r-1} X(dt) ZZ(dt/2)
    d = t_/steps; g = zz_layer(d/2, order)
    for k in range(steps): g += x_layer(d) + zz_layer(d if k < steps - 1 else d/2, order)
    return g
print("\n7) 게이트 비용 기록 — TFIM n = 4. 기저 CX + 1큐비트, 전결합; e^{−iθZZ} = 2 CX + Rz, 제어 회전 = 2 CX + 회전 2, CP = 2 CX + 위상 회전 3, CCZ = 6 CX + 7 T; 깊이 = ASAP 층 수; T ≈ 3·log2(1/ε)·R, ε = 1e-6")
for r_ in (1, 2, 4, 8, 16, 32):                                    # 2) 의 모든 스텝 수에서 게이트 리스트 = 본문 행렬
    for order in (CHAIN, BRICK):
        assert gc.same_up_to_phase(gc.unitary(trotter1_list(t, r_, order), n), trotter1(terms, t, r_))
        assert gc.same_up_to_phase(gc.unitary(trotter2_list(t, r_, order), n), trotter2(terms, t, r_))
        assert gc.same_up_to_phase(gc.unitary(gc.merge_adjacent(trotter2_list(t, r_, order)), n), trotter2(terms, t, r_))
    assert gc.same_up_to_phase(gc.unitary(trotter2_grouped_list(t, r_), n), trotter2(terms, t, r_))
print("   [검증] 위 2) 의 스텝 수 1…32 전부에서 아래 게이트 리스트의 유니터리 = trotter1/trotter2 행렬: True")
print(f"   1차 1스텝, 본문 항 순서(ZZ 01·12·23 → X): {gc.fmt(trotter1_list(t, 32, CHAIN)[:n_zz+n_x])}")
print(f"   1차 1스텝, ZZ 를 01·23·12 순으로(서로 교환 → 같은 유니터리): {gc.fmt(trotter1_list(t, 32, BRICK)[:n_zz+n_x])}")
print(f"   2차 1스텝, 본문 순서 그대로: {gc.fmt(trotter2_list(t, 32)[:2*(n_zz+n_x)])};  이어 붙은 회전을 합치면(가운데 Rx 8 → 4): {gc.fmt(gc.merge_adjacent(trotter2_list(t, 32)[:2*(n_zz+n_x)]))}")
for t_, r_, note in ((2.0, 32, "위 2) 의 마지막 행"), (6.0, 24, "위 3) 의 dt = 0.25")):
    raw = trotter2_list(t_, r_); mg = gc.merge_adjacent(raw); gr = trotter2_grouped_list(t_, r_)
    if t_ != t: assert all(gc.same_up_to_phase(gc.unitary(x, n), trotter2(terms, t_, r_)) for x in (raw, mg, gr))
    assert gc.two_qubit_count(raw) == 12*r_ and gc.two_qubit_count(mg) == 10*r_ + 2 and gc.two_qubit_count(gr) == 6*r_ + 6
    print(f"   2차 {r_}스텝 — t = {t_:g}, {note}:")
    print(f"     본문 순서 그대로 [12r]: {gc.fmt(raw)}")
    print(f"     이어 붙은 회전 병합(가운데 Rx, 스텝 경계의 ZZ 01) [10r + 2]: {gc.fmt(mg)}")
    print(f"     ZZ 를 01·23·12 순으로, 병합 없이 [12r]: {gc.fmt(trotter2_list(t_, r_, BRICK))}")
    print(f"     ZZ 묶음 전체를 스텝 경계에서 병합(ZZ 끼리 교환, 01·23·12 순) [6r + 6]: {gc.fmt(gr)}")
print(f"   1차 32스텝 (t = 2), 본문 순서: {gc.fmt(trotter1_list(t, 32))};  01·23·12 순: {gc.fmt(trotter1_list(t, 32, BRICK))}")

# 에너지 QPE: 본문은 U = expm(+i(H+wI)τ) 를 행렬로 두었다(게이트 없음). 게이트로 짜려면 U 를 Trotter 로 바꿔야 한다:
# U_s = e^{+iwτ}·(2차 Trotter, 시간 −τ, s 스텝). 제어 U_s^{2^j} = 시간 2^j·τ 를 2^j·s 스텝으로 (ZZ 묶음 병합 꼴), 회전만 제어.
def ctrl_evolution_list(c, sys_q, m, d):                           # 제어 [ZZ(d/2) {X(d) ZZ(d)}^{m-1} X(d) ZZ(d/2)], 시간 부호는 d 에 포함
    def czz(dd): return [x for i in BRICK for x in (_g("cx", sys_q[i], sys_q[i+1]), _g("crz", c, sys_q[i+1], angle=-2*J*dd), _g("cx", sys_q[i], sys_q[i+1]))]
    def cx_(dd): return [_g("crx", c, q, angle=-2*h*dd) for q in sys_q]
    g = czz(d/2)
    for k in range(m): g += cx_(d) + czz(d if k < m - 1 else d/2)
    return g
def iqft_list(t_):
    g = []
    for j in range(t_):
        g.append(_g("h", j))
        for k in range(j+1, t_): g.append(_g("cp", k, j, angle=2*np.pi/2**(k-j+1)))
    g += [_g("swap", j, t_-1-j) for j in range(t_//2)]
    return gc.inverse(g)
def energy_qpe_list(tb_, s):                                        # 계수 q[0:tb_] (q[j] 가 U^{2^{tb_-1-j}}), 계 q[tb_:tb_+n]; 시험 상태 준비는 뺀다
    sys_q = list(range(tb_, tb_ + n)); g = [_g("h", j) for j in range(tb_)]
    for j in range(tb_):
        pw = 2**(tb_-1-j); g += [_g("p", j, angle=w*tau*pw)] + ctrl_evolution_list(j, sys_q, pw*s, -tau/s)
    return g + iqft_list(tb_)
U_s1 = np.exp(1j*w*tau) * trotter2(terms, -tau, 1)
ok_c = gc.same_up_to_phase(gc.unitary([_g("p", 0, angle=w*tau)] + ctrl_evolution_list(0, [1, 2, 3, 4], 1, -tau), 5), np.kron(np.diag([1, 0]), np.eye(16)) + np.kron(np.diag([0, 1]), U_s1)); assert ok_c
print(f"   제어 U_s (s = 1: 2차 Trotter 1스텝, 회전만 제어 — 제어 ZZ = 4 CX, 제어 Rx = 2 CX): {gc.fmt(ctrl_evolution_list(0, [1, 2, 3, 4], 1, -tau))}  [= 제어(e^{{iwτ}}·trotter2(−τ, 1)): {ok_c}]")
print(f"   에너지 QPE t_b = 8 (12큐비트): 본문의 U = expm 은 게이트 리스트가 없다 → U 를 U_s 로 바꾼 회로만 센다. 2큐비트 = 20·255·s + 12·8 + 68 = 5100s + 164")
for s in (1, 8):
    eg = energy_qpe_list(8, s); U_s = np.exp(1j*w*tau) * trotter2(terms, -tau, s)
    psi_g = gc.apply([_g("h", q) for q in range(8, 8 + n)] + eg, np.eye(2**(8 + n), dtype=complex)[0], 8 + n).reshape(2**8, 2**n)   # 시험 |++++>
    p_g = (abs(psi_g)**2).sum(axis=1); p_ref = qpe_energy_probs(U_s, trial_plus, 8); ok_e = np.allclose(p_g, p_ref); assert ok_e
    assert gc.two_qubit_count(eg) == 5100*s + 164
    ks = np.argsort(p_ref)[::-1][:2]
    print(f"     s = {s}: {gc.fmt(eg)}  [분포 = qpe_energy_probs(U_s, |++++>): {ok_e}]  ‖U_s − U‖ = {np.linalg.norm(U_s - U, 2):.4f}, 상위 피크 E = {np.round([k/2**8*2*np.pi/tau - w for k in ks], 3)} (확률 {np.round(p_ref[ks], 3)}) — 본문(expm) 은 −4.735 (0.43), −4.79 (0.24)")

# 진폭 추정: Q = −A S_0 A† S_f, A = Ry(2θ) ⊗ I.  게이트: S_f = Z(q0), A† = Ry(−2θ), S_0 = X X·CZ·X X, A = Ry(2θ), 전역 −1
def Q_list(a0, a1): return [_g("z", a0), _g("ry", a0, angle=-2*theta_true), _g("x", a0), _g("x", a1), _g("cz", a0, a1), _g("x", a0), _g("x", a1), _g("ry", a0, angle=2*theta_true)]
def cQ_list(c, a0, a1):                                             # 제어 Q: 제어 Z, 제어 Ry, X X·CCZ·X X, 제어 Ry, 그리고 −1 을 제어 큐비트의 Z 로
    return [_g("cz", c, a0), _g("cry", c, a0, angle=-2*theta_true), _g("x", a0), _g("x", a1), _g("ccz", c, a0, a1), _g("x", a0), _g("x", a1),
            _g("cry", c, a0, angle=2*theta_true), _g("z", c)]
def amp_est_list(tb_):
    a0, a1 = tb_, tb_ + 1; g = [_g("ry", a0, angle=2*theta_true)] + [_g("h", j) for j in range(tb_)]
    for j in range(tb_): g += cQ_list(j, a0, a1) * 2**(tb_-1-j)
    return g + iqft_list(tb_)
ok_q = gc.same_up_to_phase(gc.unitary(Q_list(0, 1), 2), Q); ok_cq = gc.same_up_to_phase(gc.unitary(cQ_list(0, 1, 2), 3), np.kron(np.diag([1, 0]), np.eye(4)) + np.kron(np.diag([0, 1]), Q)); assert ok_q and ok_cq
print(f"   진폭 추정: Q 한 번(제어 없음) {gc.fmt(Q_list(0, 1))}  [= Q (전역 위상 제외): {ok_q}];  제어 Q 한 번 {gc.fmt(cQ_list(0, 1, 2))}  [= 제어 Q (부호 포함): {ok_cq}]")
for tb_ in (3, 5, 7):
    ag = amp_est_list(tb_); p_g = (abs(gc.apply(ag, np.eye(2**(tb_ + 2), dtype=complex)[0], tb_ + 2).reshape(2**tb_, 4))**2).sum(axis=1)
    ok_a = np.allclose(p_g, amp_est_probs(Q, A, tb_)); assert ok_a
    print(f"     t_b = {tb_} ({tb_+2}큐비트, 제어 Q {2**tb_-1}번 + 역 QFT): {gc.fmt(ag)}  [분포 = amp_est_probs: {ok_a}]")
