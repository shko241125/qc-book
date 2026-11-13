"""6장 회로 라이브러리 해부 — n_local, zz_feature_map, QFT, Draper 덧셈기를 게이트부터 numpy로 조립해 검증한다.
실행: python ch06_circuit_library.py
환경: Python 3.13, numpy 2.5 (SDK 없이 동작). 큐비트 순서: 빅 엔디안(3장 규약).
"""
import os
import numpy as np
from itertools import product

np.set_printoptions(precision=4, suppress=True, linewidth=120)
rng = np.random.default_rng(seed=7)

I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
P0 = np.diag([1, 0]).astype(complex); P1 = np.diag([0, 1]).astype(complex)
def Ry(t): return np.array([[np.cos(t/2), -np.sin(t/2)], [np.sin(t/2), np.cos(t/2)]], dtype=complex)
def Pgate(lam): return np.diag([1, np.exp(1j*lam)])

def kron(*ops):
    out = np.array([[1]], dtype=complex) if ops[0].ndim == 2 else np.array([1], dtype=complex)
    for o in ops: out = np.kron(out, o)
    return out
def single(g, q, n): ops = [I]*n; ops[q] = g; return kron(*ops)
def cnot(c, t, n):
    a = [I]*n; a[c] = P0; b = [I]*n; b[c] = P1; b[t] = X; return kron(*a) + kron(*b)
def cphase(c, t, lam, n):
    a = [I]*n; a[c] = P0; b = [I]*n; b[c] = P1; b[t] = Pgate(lam); return kron(*a) + kron(*b)
def swap(a_, b_, n): return cnot(a_, b_, n) @ cnot(b_, a_, n) @ cnot(a_, b_, n)
def ket(bits):
    v = np.array([1], dtype=complex)
    for b in bits: v = np.kron(v, np.array([0, 1]) if b == '1' else np.array([1, 0]))
    return v
def seq(gates):              # 시간 순서 리스트 → 전체 유니터리
    n = gates[0].shape[0]; M = np.eye(n, dtype=complex)
    for g in gates: M = g @ M
    return M

# ── 1. QFT: 게이트 조립 vs DFT 행렬 ──────────────────────────────
def qft_gates(n):
    """N&C 그림 5.1: H 와 제어-위상, 끝에 SWAP 으로 순서 뒤집기"""
    g = []
    for j in range(n):
        g.append(single(H, j, n))
        for k in range(j+1, n):
            g.append(cphase(k, j, 2*np.pi / 2**(k-j+1), n))   # 제어 k → 표적 j, 각 2π/2^(k-j+1)
    for j in range(n//2):
        g.append(swap(j, n-1-j, n))
    return g
def dft_matrix(N):
    w = np.exp(2j*np.pi/N)
    return np.array([[w**(j*k) for k in range(N)] for j in range(N)]) / np.sqrt(N)
for n in (2, 3, 4):
    Uq = seq(qft_gates(n))
    ok = np.allclose(Uq, dft_matrix(2**n))
    n_cp = n*(n-1)//2
    print(f"1) QFT n={n}: 게이트 조립 = DFT 행렬? {ok}   H {n}개, 제어위상 {n_cp}개, SWAP {n//2}개")
    assert ok, "5절 QFT = DFT: n = 2, 3, 4 에서 게이트 조립 행렬 = DFT 행렬"
    assert len(qft_gates(n)) == n + n_cp + n//2, "5절 QFT = DFT: 게이트 수 H n + 제어 위상 n(n−1)/2 + SWAP ⌊n/2⌋"
assert len(qft_gates(4)) == 12, "5절 QFT = DFT: n = 4 에서 4 + 6 + 2 = 12개"

# QFT 는 |x> 를 '위상이 x 에 비례해 도는' 곱상태로 보낸다
n = 3; Uq = seq(qft_gates(n))
for x in (0, 1, 5):
    psi = Uq @ ket(format(x, '03b'))
    phases = np.angle(psi) / (2*np.pi)        # 각 진폭의 위상 / 2π
    print(f"   QFT|{x}>: 진폭 크기 모두 {abs(psi[0]):.4f}, 위상/2π = {np.round(phases % 1, 3)}  (k·x/8 의 소수부)")
    assert np.allclose(abs(psi), 0.3536, atol=5e-5), "5절 QFT 위상 표: 진폭 크기 모두 1/√8 = 0.3536"
    ph_md = {0: [0, 0, 0, 0, 0, 0, 0, 0], 1: [0, .125, .25, .375, .5, .625, .75, .875], 5: [0, .625, .25, .875, .5, .125, .75, .375]}[x]
    assert np.allclose(np.exp(2j*np.pi*phases), np.exp(2j*np.pi*np.array(ph_md))), "5절 QFT 위상 표: 진폭 위상/2π (x = 0, 1, 5)"

# ── 2. Draper QFT 덧셈기: |a>|b> → |a>|a+b mod 2^n> ───────────────
def draper_adder_gates(n):
    """레지스터 a = q[0:n], b = q[n:2n] (각각 빅 엔디안). QFT(b) → 제어 위상(a 가 제어) → QFT†(b)"""
    N = 2*n
    g = []
    # QFT on b (SWAP 없이: 비트 순서를 뒤집힌 채로 두고 역 QFT 도 똑같이 뒤집어 쓰면 SWAP 이 필요 없다)
    for j in range(n):
        g.append(single(H, n+j, N))
        for k in range(j+1, n):
            g.append(cphase(n+k, n+j, 2*np.pi/2**(k-j+1), N))
    # 위상 더하기: b 의 j 번째(상위부터) 큐비트에 a 의 k 번째 큐비트가 제어하는 각 2π/2^(k-j+1)  (k>=j)
    for j in range(n):
        for k in range(j, n):
            g.append(cphase(k, n+j, 2*np.pi/2**(k-j+1), N))
    # 역 QFT on b
    inv = []
    for j in range(n):
        inv.append(single(H, n+j, N))
        for k in range(j+1, n):
            inv.append(cphase(n+k, n+j, 2*np.pi/2**(k-j+1), N))
    for gate in reversed(inv): g.append(gate.conj().T)
    return g
n = 3
Uadd = seq(draper_adder_gates(n))
ok_all = True
for a, b in product(range(2**n), repeat=2):
    inp = ket(format(a, f'0{n}b') + format(b, f'0{n}b'))
    out = Uadd @ inp
    expect = ket(format(a, f'0{n}b') + format((a+b) % 2**n, f'0{n}b'))
    if not np.isclose(abs(np.vdot(expect, out)), 1.0): ok_all = False; print("   실패:", a, b)
print(f"2) Draper 덧셈기 n={n}: 64가지 (a,b) 전수 검사 모두 통과? {ok_all}   게이트 {len(draper_adder_gates(n))}개, 보조 큐비트 0")
print("   예: a=5, b=6 → 5+6 mod 8 = 3;  |101>|110> → 측정 b 레지스터 =",
      format(int(np.argmax(abs(Uadd @ ket('101110')))) & 0b111, '03b'))
assert ok_all, "5절 Draper 덧셈기: (a, b) 64가지 전부 |a>|a+b mod 8>"
assert len(draper_adder_gates(n)) == 18 and Uadd.shape[0] == 2**(2*n), "5절 Draper 덧셈기: 게이트 18개(6 + 6 + 6), 보조 큐비트 0"
assert format(int(np.argmax(abs(Uadd @ ket('101110')))) & 0b111, '03b') == '011', "5절 Draper 덧셈기: |101>|110> → b 레지스터 011"

# ── 3. zz_feature_map: 데이터 → 회로 → 커널 ─────────────────────
def zz_feature_map_gates(x, reps=1):
    """Qiskit 정의: H^n → P(2x_i) on i → (CX, P(2(π-x_i)(π-x_j)), CX) on pairs, 반복"""
    n = len(x); g = []
    for _ in range(reps):
        for i in range(n): g.append(single(H, i, n))
        for i in range(n): g.append(single(Pgate(2*x[i]), i, n))
        for i in range(n):
            for j in range(i+1, n):
                g.append(cnot(i, j, n)); g.append(single(Pgate(2*(np.pi-x[i])*(np.pi-x[j])), j, n)); g.append(cnot(i, j, n))
    return g
def feature_state(x, reps=1):
    n = len(x); return seq(zz_feature_map_gates(x, reps)) @ ket('0'*n)
def kernel(x, y, reps=2): return abs(np.vdot(feature_state(x, reps), feature_state(y, reps)))**2
x = np.array([0.3, 1.2]); y = np.array([0.35, 1.1]); z = np.array([2.5, -0.8])
print(f"3) ZZ feature map (2큐비트, reps=2): K(x,x)={kernel(x,x):.4f}  K(x,y)={kernel(x,y):.4f}  K(x,z)={kernel(x,z):.4f}")
print("   → 가까운 점은 커널이 크고 먼 점은 작다. 커널 행렬은 대칭·양의 준정부호:",
      end=" ")
pts = [x, y, z, np.array([0.0, 0.0])]
K = np.array([[kernel(p, q) for q in pts] for p in pts])
print(np.allclose(K, K.T), (np.linalg.eigvalsh(K) > -1e-10).all())
assert abs(kernel(x,x) - 1.0000) < 5e-5 and abs(kernel(x,y) - 0.9620) < 5e-5 and abs(kernel(x,z) - 0.7247) < 5e-5, "5절 ZZ feature map 커널 표: 1.0000, 0.9620, 0.7247"
assert abs(np.linalg.norm(x - y) - 0.11) < 5e-3 and abs(np.linalg.norm(x - z) - 2.97) < 5e-3, "5절 ZZ feature map 커널 표: 거리 0.11, 2.97"
assert np.allclose(K, K.T) and (np.linalg.eigvalsh(K) > -1e-10).all(), "5절: 네 점의 커널 행렬은 대칭·양의 준정부호"

# ── 4. n_local (real_amplitudes): Ry 층 + CX 사다리 ────────────────
def real_amplitudes_gates(n, reps, theta):
    """Qiskit real_amplitudes(entanglement='reverse_linear'): Ry 층, CX 사다리, ... , Ry 층. 파라미터 n*(reps+1)개"""
    g = []; p = 0
    for r in range(reps + 1):
        for q in range(n): g.append(single(Ry(theta[p]), q, n)); p += 1
        if r < reps:
            for q in range(n-2, -1, -1): g.append(cnot(q, q+1, n))   # reverse_linear
    return g
n, reps = 4, 2
theta = rng.uniform(0, 2*np.pi, n*(reps+1))
Ua = seq(real_amplitudes_gates(n, reps, theta))
psi = Ua @ ket('0'*n)
print(f"4) real_amplitudes n={n}, reps={reps}: 파라미터 {len(theta)}개, CX {reps*(n-1)}개, 출력 진폭이 전부 실수? {np.allclose(psi.imag, 0)}")
print("   진폭:", psi.real)
assert len(theta) == 12 and reps*(n-1) == 6, "5절 real_amplitudes: 파라미터 12개, CX 6개"
assert len(psi) == 16 and np.allclose(psi.imag, 0), "5절 real_amplitudes: 출력 진폭 16개가 전부 실수"

# 바렌 고원 맛보기: 큐비트 수를 늘리며 ∂<Z0>/∂θ0 의 분산
def grad_var(n, reps, samples=60):
    gs = []
    for _ in range(samples):
        th = rng.uniform(0, 2*np.pi, n*(reps+1))
        def expZ(th):
            s = seq(real_amplitudes_gates(n, reps, th)) @ ket('0'*n)
            return np.vdot(s, single(Z, 0, n) @ s).real
        d = th.copy(); d[0] += np.pi/2; e1 = expZ(d); d[0] -= np.pi; e2 = expZ(d)
        gs.append((e1 - e2)/2)                   # parameter-shift 규칙
    return np.var(gs)
print("5) 바렌 고원: Var[∂<Z0>/∂θ0], reps=n (깊이가 큐비트 수에 비례)")
vars_ = []
for nq in (2, 4, 6, 8):
    v = grad_var(nq, nq, samples=80); vars_.append(v); print(f"   n={nq}: {v:.5f}")
for v, v_md in zip(vars_, (0.2414, 0.0560, 0.0179, 0.0047)):
    assert abs(v - v_md) < 5e-5, "5절 불모 고원 표: Var[∂<Z0>/∂θ0] (n = 2, 4, 6, 8)"
assert abs(vars_[0] - 0.24) < 5e-3 and abs(vars_[3] - 0.005) < 5e-4, "1절 요약·9절: n = 2 → 8 에서 0.24 → 0.005"
assert all(1/4.5 < vars_[i+1] / vars_[i] < 1/3 + 0.01 for i in range(3)), "5절: 2큐비트 늘어날 때마다 약 1/3~1/4"

# ── 6. 그림 ───────────────────────────────────────────────────────
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]; plt.rcParams["axes.unicode_minus"] = False
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.4))
Uq = seq(qft_gates(3))
ax1.imshow(np.angle(Uq), cmap="twilight", vmin=-np.pi, vmax=np.pi); ax1.set_title("QFT(3큐비트) 행렬의 위상"); ax1.set_xlabel("열 k"); ax1.set_ylabel("행 j")
ax2.semilogy([2, 4, 6, 8], vars_, 'o-', color="#4C6EF5"); ax2.set_xlabel("큐비트 수 n"); ax2.set_ylabel("Var[d<Z0>/dθ0]"); ax2.set_title("기울기 분산이 지수적으로 줄어든다")
plt.tight_layout(); plt.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "ch06_results.png"), dpi=150); print("6) saved ch06_results.png")

# ── 연습문제 검산 (8절) ──────────────────────────────────────────
# 8절 해답(과 7절이 인용한 reps = 1 커널값)의 수치를 전부 다시 계산한다. 위쪽 코드·난수(rng)·그림은 건드리지 않는다.
print("── 연습문제 검산 (8절) ──")

# 문제 1: Draper 덧셈기를 거꾸로 돌리면 뺄셈기
ex_n = 3
ex_g = draper_adder_gates(ex_n)
ex_Uadd = seq(ex_g)
ex_Usub = seq([gate.conj().T for gate in reversed(ex_g)])              # 모든 게이트의 역을 역순으로
print("[연습 1] 거꾸로 돌린 회로 = U† :", np.allclose(ex_Usub, ex_Uadd.conj().T))
ex_ok = all(np.isclose(abs(np.vdot(ket(format(a, '03b') + format((b - a) % 8, '03b')), ex_Usub @ ket(format(a, '03b') + format(b, '03b')))), 1.0)
            for a, b in product(range(8), repeat=2))
print("[연습 1] |a>|b> → |a>|b-a mod 8> : 64가지 (a,b) 전수 검사 모두 통과?", ex_ok)
for ex_a, ex_b in [(2, 5), (5, 2)]:
    ex_out = ex_Usub @ ket(format(ex_a, '03b') + format(ex_b, '03b'))
    ex_idx = int(np.argmax(abs(ex_out)))
    print(f"[연습 1] U†|{ex_a:03b}>|{ex_b:03b}> → a = {ex_idx >> 3:03b}, b = {ex_idx & 0b111:03b}"
          f"  ({ex_b} - {ex_a} = {ex_b - ex_a} ≡ {(ex_b - ex_a) % 8} mod 8, 확률 {abs(ex_out[ex_idx])**2:.4f})")
ex_nq = ex_n + ex_n*(ex_n-1)//2; ex_np = ex_n*(ex_n+1)//2             # QFT 게이트 6개, 위상 더하기 게이트 6개
ex_struct = seq(ex_g[:ex_nq] + [gate.conj().T for gate in ex_g[ex_nq:ex_nq+ex_np]] + ex_g[ex_nq+ex_np:])
print("[연습 1] 시간 순서 'QFT → 위상 빼기(각의 부호만 반대) → QFT⁻¹' = U† :", np.allclose(ex_struct, ex_Usub))
assert np.allclose(ex_Usub, ex_Uadd.conj().T) and ex_ok, "8절 문제 1: 거꾸로 돌린 회로는 |a>|b> → |a>|b−a mod 8> (64가지)"
assert np.allclose(ex_struct, ex_Usub), "8절 문제 1: QFT → 위상 빼기 → QFT⁻¹ = U†"

# 문제 2: 주기 r 인 빗 상태의 QFT 는 N/r 의 배수에서만 피크
ex_Uq = seq(qft_gates(3))
def ex_comb(r, offset):
    v = sum(ket(format(k, '03b')) for k in range(offset, 8, r)); return v / np.linalg.norm(v)
ex_p4 = np.abs(ex_Uq @ ex_comb(4, 0))**2                               # (|0> + |4>)/√2
ex_p2 = np.abs(ex_Uq @ ex_comb(2, 1))**2                               # (|1> + |3> + |5> + |7>)/2
print("[연습 2] QFT (|0>+|4>)/√2 의 P(k), k=0..7 =", np.round(ex_p4, 4) + 0.0, " 피크 k =", [int(k) for k in np.flatnonzero(ex_p4 > 1e-9)])
print("[연습 2] QFT (|1>+|3>+|5>+|7>)/2 의 P(k), k=0..7 =", np.round(ex_p2, 4) + 0.0, " 피크 k =", [int(k) for k in np.flatnonzero(ex_p2 > 1e-9)])
print("[연습 2] 오프셋을 바꿔도 확률 분포는 같다: 주기 4 (오프셋 0~3)",
      all(np.allclose(np.abs(ex_Uq @ ex_comb(4, o))**2, ex_p4) for o in range(4)),
      " 주기 2 (오프셋 0, 1)", all(np.allclose(np.abs(ex_Uq @ ex_comb(2, o))**2, ex_p2) for o in range(2)))
ex_amp = ex_Uq @ ex_comb(2, 1)
print("[연습 2] 오프셋 1 은 위상에만 남는다: 주기 2 상태의 피크 진폭 위상/2π =",
      {int(k): float(np.round((np.angle(ex_amp[k]) / (2*np.pi)) % 1, 3)) for k in (0, 4)}, "(오프셋 0 이면 둘 다 0)")
assert np.allclose(ex_p4, [0.25, 0, 0.25, 0, 0.25, 0, 0.25, 0]), "8절 문제 2: 주기 4 → k = 0, 2, 4, 6 에서 각 0.25"
assert np.allclose(ex_p2, [0.5, 0, 0, 0, 0.5, 0, 0, 0]), "8절 문제 2: 주기 2 → k = 0, 4 에서 각 0.5"
assert all(np.allclose(np.abs(ex_Uq @ ex_comb(4, o))**2, ex_p4) for o in range(4)) \
    and all(np.allclose(np.abs(ex_Uq @ ex_comb(2, o))**2, ex_p2) for o in range(2)), "8절 문제 2: 오프셋은 확률 분포를 바꾸지 않는다"

# 문제 3 (과 7절): ZZ feature map 의 reps = 1, 2, 3 커널
ex_x = np.array([0.3, 1.2]); ex_y = np.array([0.35, 1.1]); ex_z = np.array([2.5, -0.8])   # 3절(본문 5절 표)과 같은 세 점
ex_Kxy = [kernel(ex_x, ex_y, r) for r in (1, 2, 3)]
ex_Kxz = [kernel(ex_x, ex_z, r) for r in (1, 2, 3)]
print("[연습 3] K(x,y), reps = 1, 2, 3 :", " → ".join(f"{v:.4f}" for v in ex_Kxy))
print("[연습 3] K(x,z), reps = 1, 2, 3 :", " → ".join(f"{v:.4f}" for v in ex_Kxz))
print(f"[연습 3] K(x,z) 가 가장 작은 reps = {int(np.argmin(ex_Kxz)) + 1},  reps 가 늘 때 K(x,z) 가 단조 감소? {bool(np.all(np.diff(ex_Kxz) < 0))}")
assert np.allclose(ex_Kxy, [0.955, 0.962, 0.981], atol=5e-4), "8절 문제 3: K(x, y) = 0.955 → 0.962 → 0.981 (reps = 1, 2, 3)"
assert np.allclose(ex_Kxz, [0.071, 0.725, 0.804], atol=5e-4), "8절 문제 3: K(x, z) = 0.071 → 0.725 → 0.804 (7절의 0.07, 0.72)"
assert int(np.argmin(ex_Kxz)) + 1 == 1 and not bool(np.all(np.diff(ex_Kxz) < 0)), "8절 문제 3: K(x, z) 는 reps = 1 에서 가장 작고 단조 감소가 아니다"
