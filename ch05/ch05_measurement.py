"""5장 측정·고전 제어·동적 회로 — 사영 측정, 붕괴, 조건부 게이트, 텔레포테이션을 numpy로 검증한다.
실행: python ch05_measurement.py
환경: Python 3.13, numpy 2.5 (SDK 없이 동작). 큐비트 순서: 빅 엔디안(3장 규약).
"""
import os
import numpy as np

np.set_printoptions(precision=4, suppress=True)
rng = np.random.default_rng(seed=7)

I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
P0 = np.diag([1, 0]).astype(complex); P1 = np.diag([0, 1]).astype(complex)
ket0 = np.array([1, 0], dtype=complex); ket1 = np.array([0, 1], dtype=complex)

def kron(*ops):
    out = np.array([[1]], dtype=complex) if ops[0].ndim == 2 else np.array([1], dtype=complex)
    for o in ops: out = np.kron(out, o)
    return out
def single(g, q, n): ops = [I]*n; ops[q] = g; return kron(*ops)
def cnot(c, t, n):
    a = [I]*n; a[c] = P0; b = [I]*n; b[c] = P1; b[t] = X
    return kron(*a) + kron(*b)
def proj(q, bit, n): return single(P1 if bit else P0, q, n)

def measure(psi, q, n):
    """큐비트 q를 계산 기저로 사영 측정. (결과, 붕괴된 정규화 상태) 반환"""
    p1 = np.vdot(psi, proj(q, 1, n) @ psi).real
    bit = int(rng.random() < p1)
    post = proj(q, bit, n) @ psi
    return bit, post / np.linalg.norm(post)

def fidelity(a, b): return abs(np.vdot(a, b))**2

# ── 1. 사영 측정과 붕괴 ───────────────────────────────────────────
psi = kron(H @ ket0, ket0)                      # |+>|0>
bit, post = measure(psi, 0, 2)
bit_again = measure(post, 0, 2)[0]
print("1) |+>|0> 의 q0 측정 →", bit, " 붕괴 후:", post, " (다시 재면 같은 결과)", bit_again)
assert bit == 0 and bit_again == 0, "5절 붕괴와 반복 측정: q0 측정 0, 다시 재도 0"
assert np.allclose(post, [1, 0, 0, 0]), "5절 붕괴와 반복 측정: 붕괴 후 상태 (1, 0, 0, 0) = |00>"

# ── 2. 부분 측정: Bell 상태에서 한 쪽을 재면 다른 쪽이 정해진다 ────
bell = cnot(0, 1, 2) @ single(H, 0, 2) @ kron(ket0, ket0)
stats = {}
for _ in range(1000):
    b0, post = measure(bell, 0, 2)
    b1, _ = measure(post, 1, 2)
    stats[(b0, b1)] = stats.get((b0, b1), 0) + 1
print("2) Bell 상태 순차 측정 1000회:", dict(sorted(stats.items())), "→ 항상 같다")
assert stats == {(0, 0): 489, (1, 1): 511}, "5절: Bell 순차 측정 1000회 = (0,0) 489회, (1,1) 511회, 다른 조합 0회"

# ── 3. 텔레포테이션 ───────────────────────────────────────────────
def random_state():
    v = rng.normal(size=2) + 1j*rng.normal(size=2); return v / np.linalg.norm(v)

def teleport(phi, record=None):
    """q0: 보낼 상태 φ, q1·q2: Bell 쌍. 반환: q2의 최종 상태(2차원), 측정 비트"""
    n = 3
    psi = kron(phi, ket0, ket0)
    psi = single(H, 1, n) @ psi; psi = cnot(1, 2, n) @ psi          # Bell 쌍 (q1, q2)
    psi = cnot(0, 1, n) @ psi; psi = single(H, 0, n) @ psi          # Alice: Bell 측정 준비
    m0, psi = measure(psi, 0, n); m1, psi = measure(psi, 1, n)      # 중간 측정
    if m1: psi = single(X, 2, n) @ psi                              # 고전 제어: m1 → X
    if m0: psi = single(Z, 2, n) @ psi                              #           m0 → Z
    # q2 만 꺼낸다: q0,q1 은 |m0 m1> 로 확정되어 있으므로 해당 블록을 읽는다
    idx = (m0 << 2) | (m1 << 1)
    out = psi[idx:idx+2]
    return out / np.linalg.norm(out), (m0, m1)

print("3) 텔레포테이션 (난수 상태 200개)")
fids, bits = [], {}
for _ in range(200):
    phi = random_state()
    out, mb = teleport(phi)
    fids.append(fidelity(phi, out)); bits[mb] = bits.get(mb, 0) + 1
print(f"   충실도 최소 {min(fids):.6f}, 평균 {np.mean(fids):.6f}  (1.0 = 완벽)")
print("   측정 비트 분포:", dict(sorted(bits.items())), "→ 4가지가 균등, 즉 비트만 보면 φ 정보가 없다")
assert abs(min(fids) - 1.000000) < 5e-7 and abs(np.mean(fids) - 1.000000) < 5e-7, "5절 텔레포테이션: 충실도 최소·평균 1.000000 (200개)"
assert bits == {(0, 0): 67, (0, 1): 46, (1, 0): 40, (1, 1): 47}, "5절 텔레포테이션: 측정 비트 분포 67, 46, 40, 47"

# 보정을 생략하면?
def teleport_nocorrection(phi):
    n = 3
    psi = kron(phi, ket0, ket0)
    psi = single(H, 1, n) @ psi; psi = cnot(1, 2, n) @ psi
    psi = cnot(0, 1, n) @ psi; psi = single(H, 0, n) @ psi
    m0, psi = measure(psi, 0, n); m1, psi = measure(psi, 1, n)
    idx = (m0 << 2) | (m1 << 1); out = psi[idx:idx+2]
    return out / np.linalg.norm(out), (m0, m1)
f_nc = [fidelity(phi, teleport_nocorrection(phi)[0]) for phi in [random_state() for _ in range(200)]]
print(f"   보정(X,Z) 생략 시 충실도 평균 {np.mean(f_nc):.4f} → 고전 비트 2개 없이는 상태가 안 넘어간다")
assert abs(np.mean(f_nc) - 0.4431) < 5e-5, "5절 텔레포테이션: 보정 생략 시 평균 충실도 0.4431 (1절 요약의 0.44)"

# ── 4. 지연 측정 원리: 측정+고전제어 = 제어 게이트 + 끝에서 측정 ────
def teleport_unitary(phi):
    """측정을 끝으로 미루고 X, Z 를 CNOT, CZ 로 바꾼 버전. q2 의 환산 밀도행렬을 돌려준다."""
    n = 3
    psi = kron(phi, ket0, ket0)
    psi = single(H, 1, n) @ psi; psi = cnot(1, 2, n) @ psi
    psi = cnot(0, 1, n) @ psi; psi = single(H, 0, n) @ psi
    psi = cnot(1, 2, n) @ psi                                       # m1 → X  ⇒  CNOT(q1→q2)
    CZ02 = kron(P0, I, I) + kron(P1, I, Z); psi = CZ02 @ psi        # m0 → Z  ⇒  CZ(q0, q2)
    rho = np.outer(psi, psi.conj()).reshape(2,2,2,2,2,2)
    return np.einsum('abcabf->cf', rho)                             # q0,q1 추적 → q2
phi = random_state()
rho2 = teleport_unitary(phi)
print("4) 지연 측정: 유니터리 버전의 q2 밀도행렬이 |φ><φ| 와 같은가:", np.allclose(rho2, np.outer(phi, phi.conj())))
assert np.allclose(rho2, np.outer(phi, phi.conj())), "5절 지연 측정: 유니터리 버전의 q2 밀도행렬 = |φ><φ|"

# ── 5. 측정 기저 변경: X 기저 측정 = H 후 Z 측정 ────────────────────
plus = H @ ket0
cnt = sum(measure(H @ plus, 0, 1)[0] for _ in range(1000))
print("5) |+> 를 X 기저로 측정(H 후 Z 측정) 1000회 중 결과 1:", cnt, "회 → 확정적으로 0")
assert cnt == 0, "5절 기저 변경: |+> 를 X 기저로 1000회 재면 1이 0회"

# ── 6. 반복 측정으로 상태 준비: 측정 후 조건부 X 로 |0> 만들기 (reset) ──
def reset(psi, q, n):
    b, psi = measure(psi, q, n)
    if b: psi = single(X, q, n) @ psi
    return psi
psi = kron(random_state(), random_state())
psi = reset(psi, 0, 2)
print("6) reset 후 q0 가 |0> 인 확률:", np.vdot(psi, proj(0, 0, 2) @ psi).real)
assert abs(np.vdot(psi, proj(0, 0, 2) @ psi).real - 1.0) < 1e-9, "5절 reset: q0 가 |0> 일 확률 1.0"

# ── 7. 측정 오류 모델: 읽기 오류 p 가 텔레포테이션 충실도에 미치는 영향 ──
def teleport_noisy(phi, p):
    n = 3
    psi = kron(phi, ket0, ket0)
    psi = single(H, 1, n) @ psi; psi = cnot(1, 2, n) @ psi
    psi = cnot(0, 1, n) @ psi; psi = single(H, 0, n) @ psi
    m0, psi = measure(psi, 0, n); m1, psi = measure(psi, 1, n)
    r0 = m0 ^ int(rng.random() < p); r1 = m1 ^ int(rng.random() < p)   # 읽기 오류: 비트가 뒤집혀 보고됨
    if r1: psi = single(X, 2, n) @ psi
    if r0: psi = single(Z, 2, n) @ psi
    idx = (m0 << 2) | (m1 << 1); out = psi[idx:idx+2]
    return out / np.linalg.norm(out)
print("7) 읽기 오류율 p 에 따른 평균 충실도 (난수 상태 300개)")
ps = [0.0, 0.01, 0.03, 0.05, 0.1]
fid_p = []
for p in ps:
    f = np.mean([fidelity(phi, teleport_noisy(phi, p)) for phi in [random_state() for _ in range(300)]])
    fid_p.append(f); print(f"   p={p:.2f}: F={f:.4f}")
assert ps == [0.0, 0.01, 0.03, 0.05, 0.1], "5절 판독 오류 표: p 열"
for f, f_md in zip(fid_p, (1.0000, 0.9860, 0.9612, 0.9386, 0.8583)):
    assert abs(f - f_md) < 5e-5, "5절 판독 오류 표: 평균 충실도(난수 상태 300개)"
assert abs(fid_p[3] - 0.94) < 5e-3, "1절 요약: p = 0.05 에서 충실도 0.94"

# ── 8. 그림 ───────────────────────────────────────────────────────
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]; plt.rcParams["axes.unicode_minus"] = False
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.4))
bins = np.linspace(0, 1, 21)
ax1.hist(np.clip(f_nc, 0, 1), bins=bins, color="#868E96", alpha=0.7, label="보정 생략"); ax1.hist(np.clip(fids, 0, 1), bins=bins, color="#4C6EF5", label="보정 포함")
ax1.set_xlabel("충실도 |<φ|ψ_out>|²"); ax1.set_ylabel("횟수"); ax1.set_title("텔레포테이션 200회 충실도"); ax1.legend(frameon=False)
ax2.plot(ps, fid_p, 'o-', color="#4C6EF5"); ax2.set_xlabel("읽기 오류율 p"); ax2.set_ylabel("평균 충실도"); ax2.set_title("측정 오류 → 충실도")
plt.tight_layout(); plt.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "ch05_results.png"), dpi=150); print("8) saved ch05_results.png")

# ── 연습문제 검산 (8절) ──────────────────────────────────────────
# 8절 해답이 인용하는 수치·행렬을 전부 다시 계산한다. 위쪽 코드·그림은 건드리지 않고,
# 난수는 위의 rng 대신 이 블록 전용 ex_rng 를 쓴다(앞 출력이 바뀌지 않게).
print("── 연습문제 검산 (8절) ──")
ex_rng = np.random.default_rng(seed=508)
def ex_state():
    v = ex_rng.normal(size=2) + 1j*ex_rng.normal(size=2); return v / np.linalg.norm(v)

# 문제 1: Bell 쌍이 |Ψ-> 일 때의 보정표
ex_psi_m = (kron(ket0, ket1) - kron(ket1, ket0)) / np.sqrt(2)          # |Ψ-> = (|01> - |10>)/√2
XZ = X @ Z
print(f"[연습 1] <Ψ-|(I⊗XZ)|Φ+> = {np.vdot(ex_psi_m, single(XZ, 1, 2) @ bell).real:+.4f} → |Ψ-> = (I⊗XZ)|Φ+>")
ex_A = ex_rng.normal(size=(2, 2)) + 1j*ex_rng.normal(size=(2, 2))
print("[연습 1] (A⊗I)|Φ+> = (I⊗Aᵀ)|Φ+> (난수 A):", np.allclose(single(ex_A, 0, 2) @ bell, single(ex_A.T, 1, 2) @ bell),
      "  (XZ)ᵀ = -XZ :", np.allclose(XZ.T, -XZ))
print(f"[연습 1] <Ψ-|(XZ⊗I)|Φ+> = {np.vdot(ex_psi_m, single(XZ, 0, 2) @ bell).real:+.4f} → 어느 쪽에 걸든 전역 위상(-1)만 다른 같은 상태")
def ex_bob(phi, m0, m1, pair):
    """teleport() 와 같은 회로. pair='Ψ-' 이면 Bell 쌍 생성 뒤 q1 에 Z, X 를 차례로 건다.
    측정은 난수 대신 결과 (m0, m1) 로 사영한다. 반환: 보정 전 q2 상태, 그 결과의 확률"""
    n = 3
    psi = kron(phi, ket0, ket0)
    psi = single(H, 1, n) @ psi; psi = cnot(1, 2, n) @ psi
    if pair == 'Ψ-': psi = single(Z, 1, n) @ psi; psi = single(X, 1, n) @ psi
    psi = cnot(0, 1, n) @ psi; psi = single(H, 0, n) @ psi
    post = proj(1, m1, n) @ proj(0, m0, n) @ psi
    idx = (m0 << 2) | (m1 << 1); out = post[idx:idx+2]
    return out / np.linalg.norm(out), np.vdot(post, post).real
ex_paulis = {'I': I, 'X': X, 'Z': Z, 'XZ': XZ}
ex_phis = [ex_state() for _ in range(20)]
ex_table = {}
for ex_pair in ('Φ+', 'Ψ-'):
    for ex_m in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        ex_minf = {name: min(fidelity(phi, G @ ex_bob(phi, *ex_m, ex_pair)[0]) for phi in ex_phis) for name, G in ex_paulis.items()}
        ex_ps = [ex_bob(phi, *ex_m, ex_pair)[1] for phi in ex_phis]
        ex_table[(ex_pair, ex_m)] = [name for name, f in ex_minf.items() if np.isclose(f, 1.0)]
        if ex_pair == 'Ψ-':
            print(f"[연습 1] |Ψ-> 쌍, (m0,m1)={ex_m}: 충실도 1 이 되는 보정 = {'/'.join(ex_table[(ex_pair, ex_m)])}"
                  f"   (난수 상태 20개 최소 충실도 " + " ".join(f"{k}:{v:.4f}" for k, v in ex_minf.items())
                  + f";  결과 확률 {min(ex_ps):.4f}~{max(ex_ps):.4f})")
print("[연습 1] 비교: |Φ+> 쌍의 보정 X^m1 Z^m0 =",
      ", ".join(f"{m}→{'/'.join(ex_table[('Φ+', m)])}" for m in [(0, 0), (0, 1), (1, 0), (1, 1)]))
ex_flip = {'I': 'XZ', 'XZ': 'I', 'X': 'Z', 'Z': 'X'}                    # XZ 를 한 번 더 곱한 것(전역 위상 제외)
print("[연습 1] |Ψ-> 표 = |Φ+> 표에 XZ 를 합성한 것(I↔XZ, X↔Z 로 뒤집힘):",
      all(ex_table[('Ψ-', m)] == [ex_flip[ex_table[('Φ+', m)][0]]] for m in [(0, 0), (0, 1), (1, 0), (1, 1)]))
assert abs(abs(np.vdot(ex_psi_m, single(XZ, 1, 2) @ bell)) - 1) < 1e-9, "8절 문제 1: |Ψ-> = (I⊗XZ)|Φ+> (전역 위상 제외)"
assert np.allclose(single(ex_A, 0, 2) @ bell, single(ex_A.T, 1, 2) @ bell) and np.allclose(XZ.T, -XZ), "8절 문제 1: (A⊗I)|Φ+> = (I⊗Aᵀ)|Φ+>, (XZ)ᵀ = −XZ"
assert [ex_table[('Ψ-', m)] for m in [(0, 0), (0, 1), (1, 0), (1, 1)]] == [['XZ'], ['Z'], ['X'], ['I']], "8절 문제 1: |Ψ-> 쌍의 보정 (0,0)→XZ, (0,1)→Z, (1,0)→X, (1,1)→I"
assert [ex_table[('Φ+', m)] for m in [(0, 0), (0, 1), (1, 0), (1, 1)]] == [['I'], ['X'], ['Z'], ['XZ']], "8절 문제 1: |Φ+> 쌍의 보정 X^m1 Z^m0"

# 문제 2: 복제 불가 — |0> 과 |+> 의 내적
ex_00 = kron(ket0, ket0)
print(f"[연습 2] <0|+> = {np.vdot(ket0, plus).real:.4f},  입력 쪽 <00|+0> = {np.vdot(ex_00, kron(plus, ket0)).real:.4f},"
      f"  복제됐다면 출력 쪽 <00|++> = {np.vdot(ex_00, kron(plus, plus)).real:.4f}  → 내적이 보존되지 않는다")
ex_cn = cnot(0, 1, 2) @ kron(plus, ket0)
print("[연습 2] CNOT|+>|0> = |Φ+> :", np.allclose(ex_cn, bell), f"  |<++|CNOT|+0>|² = {fidelity(kron(plus, plus), ex_cn):.4f}  (복사였다면 1)")
assert abs(np.vdot(ex_00, kron(plus, ket0)).real - 0.7071) < 5e-5, "8절 문제 2: 입력 쪽 내적 1/√2 ≈ 0.7071"
assert abs(np.vdot(ex_00, kron(plus, plus)).real - 0.5) < 1e-9, "8절 문제 2: 복제됐다면 출력 쪽 내적 1/2"
assert np.allclose(ex_cn, bell) and not np.allclose(ex_cn, kron(plus, plus)), "8절 문제 2: CNOT|+>|0> = |Φ+> ≠ |+>|+>"

# 문제 3: reset 의 Kraus 표현 K0 = |0><0|, K1 = |0><1|
K0 = np.outer(ket0, ket0.conj()); K1 = np.outer(ket0, ket1.conj())
print("[연습 3] K0 = P0, K1 = X·P1 (reset() 의 두 갈래: 결과 0 이면 사영만, 1 이면 사영 뒤 X):",
      np.allclose(K0, proj(0, 0, 1)) and np.allclose(K1, single(X, 0, 1) @ proj(0, 1, 1)))
print("[연습 3] K0†K0 + K1†K1 = I (추적 보존):", np.allclose(K0.conj().T @ K0 + K1.conj().T @ K1, I))
def ex_reset_channel(rho): return K0 @ rho @ K0.conj().T + K1 @ rho @ K1.conj().T
ex_dev = max(np.abs(ex_reset_channel(np.outer(s, s.conj())) - np.diag([1, 0])).max() for s in [ex_state() for _ in range(20)])
print(f"[연습 3] 난수 순수 상태 20개: 출력과 diag(1, 0) 의 최대 차이 = {ex_dev:.1e}")
ex_mix = 0.3 * np.outer(ex_phis[0], ex_phis[0].conj()) + 0.7 * np.outer(ex_phis[1], ex_phis[1].conj())   # 혼합 상태
ex_out = ex_reset_channel(ex_mix)
print("[연습 3] 혼합 상태 ρ: 출력 = (ρ00 + ρ11)|0><0| :", np.allclose(ex_out, (ex_mix[0, 0] + ex_mix[1, 1]) * K0),
      f"  출력 대각 = ({ex_out[0, 0].real:.4f}, {ex_out[1, 1].real:.4f})")
ex_o0 = ex_reset_channel(np.outer(ket0, ket0.conj())); ex_o1 = ex_reset_channel(np.outer(ket1, ket1.conj()))
print(f"[연습 3] 유니터리가 아니다: 직교 입력 |0>, |1> (겹침 {abs(np.vdot(ket0, ket1))**2:.4f}) 의 출력 겹침 Tr(ρ0'ρ1') = {np.trace(ex_o0 @ ex_o1).real:.4f},"
      f"  K0, K1 이 유니터리? {np.allclose(K0.conj().T @ K0, I)}, {np.allclose(K1.conj().T @ K1, I)}")
assert np.allclose(K0.conj().T @ K0 + K1.conj().T @ K1, I), "8절 문제 3: K0†K0 + K1†K1 = I (추적 보존)"
assert ex_dev < 1e-12, "8절 문제 3: 난수 상태를 넣으면 출력이 항상 diag(1, 0)"
assert np.allclose(ex_out, (ex_mix[0, 0] + ex_mix[1, 1]) * K0), "8절 문제 3: 채널 출력 = (ρ00 + ρ11)|0><0|"
assert abs(np.trace(ex_o0 @ ex_o1).real - 1) < 1e-9 and not np.allclose(K0.conj().T @ K0, I) and not np.allclose(K1.conj().T @ K1, I), "8절 문제 3: reset 은 유니터리가 아니다(모든 입력을 한 점으로)"
