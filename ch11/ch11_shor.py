"""11장 Shor 알고리즘 — 위수 찾기(모듈러 곱셈 유니터리 + QPE + 연분수)로 N=15, 21 을 인수분해하는 전 과정을 numpy 로 검증한다.
실행: python ch11_shor.py
환경: Python 3.13, numpy 2.5 (SDK 없이 동작). 큐비트 순서: 빅 엔디안(3장 규약). 계수 레지스터 t비트, 작업 레지스터 n = ⌈log2 N⌉ 비트.
"""
import numpy as np
from math import gcd
from fractions import Fraction

np.set_printoptions(precision=4, suppress=True, linewidth=120)
rng = np.random.default_rng(seed=7)

I = np.eye(2, dtype=complex); H = np.array([[1, 1], [1, -1]], dtype=complex)/np.sqrt(2)
P0 = np.diag([1, 0]).astype(complex); P1 = np.diag([0, 1]).astype(complex)
def Pgate(lam): return np.diag([1, np.exp(1j*lam)])
def kron(*ops):
    out = np.array([[1]], dtype=complex) if ops[0].ndim == 2 else np.array([1], dtype=complex)
    for o in ops: out = np.kron(out, o)
    return out
def single(g, q, n): ops = [I]*n; ops[q] = g; return kron(*ops)
def cphase(c, t_, lam, n):
    a = [I]*n; a[c] = P0; b = [I]*n; b[c] = P1; b[t_] = Pgate(lam); return kron(*a) + kron(*b)
def cnot(c, t_, n):
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    a = [I]*n; a[c] = P0; b = [I]*n; b[c] = P1; b[t_] = X; return kron(*a) + kron(*b)
def swap(a_, b_, n): return cnot(a_, b_, n) @ cnot(b_, a_, n) @ cnot(a_, b_, n)
def Hn(n): return kron(*([H]*n))
def seq(gates):
    M = np.eye(gates[0].shape[0], dtype=complex)
    for g in gates: M = g @ M
    return M
def qft(n):
    g = []
    for j in range(n):
        g.append(single(H, j, n))
        for k in range(j+1, n): g.append(cphase(k, j, 2*np.pi/2**(k-j+1), n))
    for j in range(n//2): g.append(swap(j, n-1-j, n))
    return seq(g)
def iqft(n): return qft(n).conj().T
def bits(i, n): return format(i, f'0{n}b')

# ── 1. 고전 환원: 인수분해 → 위수 찾기 ─────────────────────────────────
def order_classical(a, N):
    r = 1; x = a % N
    while x != 1: x = (x*a) % N; r += 1
    return r
def factor_from_order(a, r, N):
    if r % 2: return None, "r 홀수"
    y = pow(a, r//2, N)
    if y == N - 1: return None, "a^{r/2} ≡ -1"
    f = sorted({gcd(y-1, N), gcd(y+1, N)} - {1, N})
    return (f[0], N//f[0]) if f else (None, "자명")
N = 15
print(f"1) N = {N}: 모든 a 에 대해 위수와 인수분해 성공 여부")
good = 0; total = 0
for a in range(2, N):
    if gcd(a, N) != 1: print(f"   a={a:>2d}: gcd(a,N)={gcd(a,N)} → 바로 인수 (운)"); continue
    r = order_classical(a, N); f, why = factor_from_order(a, r, N); total += 1
    if f: good += 1; print(f"   a={a:>2d}: r={r}, a^(r/2)={pow(a, r//2, N):>2d} → gcd 로 인수 {f}")
    else: print(f"   a={a:>2d}: r={r} → 실패 ({why})")
print(f"   → 서로소 a 중 성공 {good}/{total} (일반 하한 1/2 = 1 − 1/2^(m−1): Shor 1997, N&C 정리 5.3 정오표)")
# 1b) 하한 검산: Z_N^* 전체(a = 1 포함)에서 "r 짝수이고 a^(r/2) ≢ −1" 인 a 의 비율. N = pq 의 일반 하한은 1/2 이고 실제로 닿는다
print("1b) 서로소 a 전체(a = 1 포함) 중 'r 짝수이고 a^(r/2) ≢ −1' 인 비율 (하한 1 − 1/2^(m−1) = 1/2):")
for N_ in (15, 21, 33, 57, 77):
    units = [a_ for a_ in range(1, N_) if gcd(a_, N_) == 1]
    ok = [a_ for a_ in units if factor_from_order(a_, order_classical(a_, N_), N_)[0]]
    bad = [a_ for a_ in units if a_ not in ok]
    print(f"   N = {N_}: {len(ok)}/{len(units)} = {len(ok)/len(units):.2f}" + (f"   실패 a = {bad}" if N_ in (15, 21) else ""))

# ── 2. 모듈러 곱셈 유니터리 U_a 와 고유상태 ─────────────────────────────
def U_mod(a, N, n):
    """U_a|y> = |a·y mod N>  (y < N),  y ≥ N 은 항등 (유니터리 유지)"""
    U = np.zeros((2**n, 2**n), dtype=complex)
    for y in range(2**n):
        U[(a*y) % N if y < N else y, y] = 1
    return U
n = 4; a = 7; r = order_classical(a, N)
Ua = U_mod(a, N, n)
print(f"\n2) U_{a} (mod {N}, {n}큐비트): 유니터리? {np.allclose(Ua.conj().T @ Ua, np.eye(16))},  U^{r} = I? {np.allclose(np.linalg.matrix_power(Ua, r), np.eye(16))}")
# 고유상태 |u_s> = (1/√r) Σ_k e^{-2πi s k / r} |a^k mod N>
def u_s(s, a, r, N, n):
    v = np.zeros(2**n, dtype=complex)
    for k in range(r): v[pow(a, k, N)] += np.exp(-2j*np.pi*s*k/r)
    return v/np.sqrt(r)
ok = all(np.allclose(Ua @ u_s(s, a, r, N, n), np.exp(2j*np.pi*s/r)*u_s(s, a, r, N, n)) for s in range(r))
tot = sum(u_s(s, a, r, N, n) for s in range(r))/np.sqrt(r)
print(f"   |u_s> 가 고유값 e^(2πi s/{r}) 의 고유상태? {ok};  (1/√r)Σ_s|u_s> = |1>? {np.allclose(tot, np.eye(16)[1])}  (연습문제 1)")

# ── 3. 위수 찾기 회로: QPE(U_a, |1>) → s/r ─────────────────────────────
def dft_matrix(T):                                           # 6장: QFT 게이트 조립 = DFT 행렬 (검증 완료) → 큰 t 에서는 직접 사용
    w = np.exp(2j*np.pi/T); j = np.arange(T)
    return w**np.outer(j, j) / np.sqrt(T)
def order_finding_probs(a, N, n, t):
    """상태를 (2^t, 2^n) 배열로 두고 게이트를 직접 적용 (전체 유니터리 행렬을 만들지 않는다 — 16큐비트면 64 GB)"""
    T, M = 2**t, 2**n
    psi = np.zeros((T, M), dtype=complex); psi[:, 1] = 1/np.sqrt(T)       # H^t|0..0> ⊗ |1>
    idx = np.arange(T)
    for j in range(t):                                                      # 계수 큐비트 j (q[0] 최상위) 가 U^{2^{t-1-j}} 제어
        ap = pow(a, 2**(t-1-j), N)                                          # '사전 제곱': a^{2^j} mod N 을 고전으로 — Shor 가 다항식인 이유
        perm = np.array([(ap*y) % N if y < N else y for y in range(M)])
        inv = np.empty(M, dtype=int); inv[perm] = np.arange(M)              # U|y> = |perm[y]>  ⇒  new[:, perm[y]] = old[:, y]
        mask = ((idx >> (t-1-j)) & 1) == 1
        psi[mask] = psi[mask][:, inv]
    psi = dft_matrix(T).conj().T @ psi                                      # 역 QFT (계수 레지스터)
    return (abs(psi)**2).sum(axis=1)
def continued_fraction_r(k, t, N):
    """k/2^t 의 연분수 수렴분수 중 분모 ≤ N 인 마지막 것의 분모를 r 후보로"""
    return Fraction(k, 2**t).limit_denominator(N).denominator
t = 2*n + 1                                                 # = 9 (N&C 권장 t = 2n+1)
p = order_finding_probs(7, 15, 4, t)
peaks = [int(k) for k in np.argsort(p)[::-1][:4]]
print(f"\n3) 위수 찾기 a=7, N=15, t={t}: 상위 4 피크 k = {sorted(peaks)} → k/2^t = {[round(k/2**t, 4) for k in sorted(peaks)]} = s/4 (r=4)")
print(f"   피크 확률 합 = {p[peaks].sum():.4f} (r | 2^t 이므로 정확히 1)")
# 측정 샘플링 → 연분수 → r → 인수
def shor_trial(a, N, n, t, p=None, multiples=1):
    """측정 → 연분수 r̂ → (선택) r̂ 의 작은 배수도 시도 → 인수.  multiples=1 이면 교과서 그대로."""
    if p is None: p = order_finding_probs(a, N, n, t)
    k = rng.choice(len(p), p=p/p.sum())
    r_hat = continued_fraction_r(k, t, N)
    for m in range(1, multiples + 1):                                   # 연분수가 r 의 약수 r/gcd(s,r) 를 주는 경우를 배수로 구제
        rr = m * r_hat
        if rr and pow(a, rr, N) == 1:
            return k, rr, factor_from_order(a, rr, N)[0]
    return k, r_hat, None
print("   샘플 10회 (k → 연분수 r̂ → 인수):")
hits = 0
for _ in range(10):
    k, r_hat, f = shor_trial(7, 15, 4, t, p); hits += f is not None
    print(f"     k={k:>3d} ({k/2**t:.4f}) → r̂={r_hat} → {f if f else '실패(재시도)'}")
trials = 500; succ = sum(shor_trial(7, 15, 4, t, p)[2] is not None for _ in range(trials))
print(f"   {trials}회: 인수분해 성공 {succ/trials:.3f}  (이론: k=0 과 s=2 (r̂=2) 가 실패 → 1/2)")

# ── 4. N = 21: r = 6 은 2^t 를 나누지 않으므로 위상이 비정확 ────────────────
N2 = 21; n2 = 5; a2 = 2; r2 = order_classical(a2, N2); t2 = 2*n2 + 1    # t = 11, 2^11 = 2048
print(f"\n4) N = 21, a = 2: r = {r2}, t = {t2} (계수 {2**t2}, 전체 {t2+n2} = 16 큐비트 상태벡터)")
p21 = order_finding_probs(a2, N2, n2, t2)
peaks = [int(k) for k in np.argsort(p21)[::-1][:6]]
print(f"   상위 6 피크 k/2^t = {sorted(round(k/2**t2, 4) for k in peaks)}  vs s/6 = {[round(s/6, 4) for s in range(6)]}")
print(f"   상위 6 피크 확률 합 = {p21[peaks].sum():.4f} (비정확 위상 → 1 보다 작다, 10장 Dirichlet 꼬리)")
succ = sum(shor_trial(a2, N2, n2, t2, p21)[2] is not None for _ in range(trials))
succ_m = sum(shor_trial(a2, N2, n2, t2, p21, multiples=6)[2] is not None for _ in range(trials))
print(f"   {trials}회: 인수분해 성공 {succ/trials:.3f}  (실패: s=0, gcd(s,6)>1 인 s=2,3,4 → 유효 s=1,5 → 1/3 규모 + 꼬리)")
print(f"   r̂ 의 배수(≤6배)도 시도하면 성공 {succ_m/trials:.3f}  → s=2,3,4 의 r̂=3,2,3 을 배수로 구제; k=0(r̂=1)조차 6배로 구제 — r 이 작을 때만 통하는 트릭")

# ── 5. t 에 따른 성공률 (N=21): 연분수가 s/r 를 복원하려면 2^t ≥ N² 필요 ──────
print(f"\n5) N=21, a=2: 계수 비트 t 에 따른 '연분수 r̂ = 6 복원률' (측정 500회)")
tres = {}
for tt in (5, 7, 9, 11):
    pp = order_finding_probs(a2, N2, n2, tt)
    ks = rng.choice(len(pp), size=trials, p=pp/pp.sum())
    rate = np.mean([continued_fraction_r(k, tt, N2) == 6 for k in ks]); tres[tt] = rate
    print(f"   t={tt:>2d} (2^t = {2**tt:>5d}, N² = {N2**2}): r̂=6 복원률 {rate:.3f}")

# ── 6. 그림 ───────────────────────────────────────────────────────
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Noto Sans CJK JP"; plt.rcParams["axes.unicode_minus"] = False
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 3.4))
xs = np.arange(16); fx = [pow(7, int(x), 15) for x in xs]
ax1.stem(xs, fx, basefmt=" "); ax1.set_xlabel("x"); ax1.set_ylabel("7^x mod 15"); ax1.set_title("주기 함수: r = 4")
ax2.plot(np.arange(2**t)/2**t, p, color="#4C6EF5"); ax2.set_xlabel("측정값 k/2^t"); ax2.set_ylabel("확률"); ax2.set_title("a=7, N=15: 피크 = s/4")
ax3.plot(np.arange(2**t2)/2**t2, p21, color="#F03E3E"); ax3.set_xlabel("측정값 k/2^t"); ax3.set_title("a=2, N=21: 피크 = s/6 (비정확)")
for s in range(6): ax3.axvline(s/6, ls=":", color="#868E96", lw=0.8)
plt.tight_layout(); plt.savefig("ch11_results.png", dpi=150); print("6) saved ch11_results.png")

# ── 게이트 비용 기록 (2큐비트 게이트 / 깊이 / T-count) ──
# 위 시뮬레이션은 제어 U_a 를 치환으로 직접 적용했다(게이트 없음). 여기서는 N = 15, a = 7 에 한해 '컴파일된' 게이트 회로를
# 실제로 짜서 치환 행렬·order_finding_probs 와 대조하고 센다. 난수 미사용.
# 규약(tools/gate_cost.py): 기저 CX + 1큐비트 게이트, 전결합. SWAP = 3 CX, Toffoli = 6 CX + 7 T, Fredkin = Toffoli + 2 CX,
# CP(θ) = 2 CX + 위상 회전 3 (CP(π) = CZ). 깊이 = ASAP 층 수. π/4 의 배수가 아닌 회전은 R 로 세고 T ≈ 3·log2(1/ε)·R (ε = 1e-6).
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools"))
import gate_cost as gc
_g = gc.G
def mul7_list(w):                                 # ×7 = ×8·(−1) mod 15: 오른쪽 순환(SWAP 3) 뒤 전부 반전(X 4). w = 작업 큐비트 4개(w[0] 최상위)
    return [_g("swap", w[2], w[3]), _g("swap", w[1], w[2]), _g("swap", w[0], w[1])] + [_g("x", q) for q in w]
def c_mul7_list(c, w):                            # 제어판: Fredkin 3 + CX 4
    return [_g("cswap", c, w[2], w[3]), _g("cswap", c, w[1], w[2]), _g("cswap", c, w[0], w[1])] + [_g("cx", c, q) for q in w]
def mul4_list(w): return [_g("swap", w[0], w[2]), _g("swap", w[1], w[3])]          # ×4 = 두 자리 순환: SWAP 2
def c_mul4_list(c, w): return [_g("cswap", c, w[0], w[2]), _g("cswap", c, w[1], w[3])]
def iqft_list(t_):                                # qft(t_) 의 게이트 순서를 뒤집은 것
    g = []
    for j in range(t_):
        g.append(_g("h", j))
        for k in range(j+1, t_): g.append(_g("cp", k, j, angle=2*np.pi/2**(k-j+1)))
    g += [_g("swap", j, t_-1-j) for j in range(t_//2)]
    return gc.inverse(g)
def order_finding_list(t_):                       # N = 15, a = 7: 계수 q[0:t_] (q[0] 최상위), 작업 q[t_:t_+4]
    w = list(range(t_, t_ + 4)); g = [_g("x", w[3])] + [_g("h", j) for j in range(t_)]
    for j in range(t_):
        ap = pow(7, 2**(t_-1-j), 15)                # 사전 제곱: 7, 4, 1, 1, …
        if ap == 7: g += c_mul7_list(j, w)
        elif ap == 4: g += c_mul4_list(j, w)
        else: assert ap == 1                        # ×1 = 항등 → 게이트 없음
    return g + iqft_list(t_)
def ctrl_block(M): return np.kron(P0, np.eye(len(M))) + np.kron(P1, M)
print("\n7) 게이트 비용 기록 — N = 15, a = 7 의 '컴파일된' 회로(이 a·N 전용). 기저 CX + 1큐비트, 전결합; SWAP = 3 CX, Fredkin = Toffoli(6 CX, 7 T) + 2 CX, CP = 2 CX + 위상 회전 3; 깊이 = ASAP 층 수; T ≈ 3·log2(1/ε)·R, ε = 1e-6")
W = [0, 1, 2, 3]; M7 = gc.unitary(mul7_list(W), 4); M4 = gc.unitary(mul4_list(W), 4); U7 = U_mod(7, 15, 4); U4 = U_mod(4, 15, 4)
ok7 = np.allclose(M7[:, 1:15], U7[:, 1:15]); swap015 = np.isclose(M7[15, 0], 1) and np.isclose(M7[0, 15], 1); ok4 = np.allclose(M4, U4); assert ok7 and swap015 and ok4
print(f"   ×7 mod 15 = SWAP 3 + X 4: {gc.fmt(mul7_list(W))}  [y = 1…14 에서 U_7 과 같은가: {ok7}; y = 0 ↔ 15 는 맞바꿈(U_7 은 둘 다 고정): {bool(swap015)} — |1> 의 궤도 {{1, 7, 4, 13}} 밖]")
print(f"   ×4 mod 15 = SWAP 2: {gc.fmt(mul4_list(W))}  [= U_4 (16개 전부): {ok4}]")
C7 = gc.unitary(c_mul7_list(0, [1, 2, 3, 4]), 5); C4 = gc.unitary(c_mul4_list(0, [1, 2, 3, 4]), 5)
okc7 = gc.same_up_to_phase(C7, ctrl_block(M7)) and np.allclose(C7[:, 17:31], ctrl_block(U7)[:, 17:31]); okc4 = gc.same_up_to_phase(C4, ctrl_block(U4)); assert okc7 and okc4
print(f"   제어 ×7 = Fredkin 3 + CX 4: {gc.fmt(c_mul7_list(0, [1, 2, 3, 4]))}  [= 제어(SWAP 3 + X 4), 제어 U_7 과 y = 1…14 에서 일치: {okc7}]")
print(f"   제어 ×4 = Fredkin 2: {gc.fmt(c_mul4_list(0, [1, 2, 3, 4]))}  [= 제어 U_4: {okc4}]")
# 작은 t 에서 전체 유니터리 경로를 행렬로 대조: |0…0>|0001> 에 본문식 행렬(H^t, 제어 U_7^{2^{t-1-j}}, QFT†)을 건 상태 = 게이트 리스트의 상태
t_s = 4; ref = np.kron(np.eye(2**t_s)[0], np.eye(16)[1]).astype(complex); ref = kron(Hn(t_s), np.eye(16)) @ ref
for j in range(t_s):
    a_ = [I]*t_s; a_[j] = P0; b_ = [I]*t_s; b_[j] = P1
    ref = (kron(*a_, np.eye(16)) + kron(*b_, U_mod(pow(7, 2**(t_s-1-j), 15), 15, 4))) @ ref
ref = kron(iqft(t_s), np.eye(16)) @ ref
ok_s = gc.same_up_to_phase(gc.apply(order_finding_list(t_s), np.eye(2**(t_s+4), dtype=complex)[0], t_s+4), ref); assert ok_s
for t_ in (9, 8):
    og = order_finding_list(t_); psi_g = gc.apply(og, np.eye(2**(t_+4), dtype=complex)[0], t_+4).reshape(2**t_, 16)
    ok_p = np.allclose((abs(psi_g)**2).sum(axis=1), order_finding_probs(7, 15, 4, t_)); assert ok_p
    assert gc.two_qubit_count(og) == 28 + 16 + t_*(t_-1) + 3*(t_//2)
    print(f"   위수 찾기 N = 15, a = 7, t = {t_} ({t_+4}큐비트{', 본문 3)' if t_ == 9 else ', 4절 Qiskit 스니펫의 t'}): {gc.fmt(og)}  [계수 분포 = order_finding_probs: {ok_p}; t = {t_s} 에서 상태 = 본문식 행렬: {ok_s}]")
print("     내역: 제어 ×7 (28) + 제어 ×4 (16) + 역 QFT (CP t(t−1)/2 개 = CX t(t−1), SWAP ⌊t/2⌋ 개 = 3⌊t/2⌋ CX); 제어 ×7^{2^j} (j ≥ 2) 는 ×1 이라 게이트 없음")
print("   N = 21 (a = 2, 16큐비트)과 일반 N: 이 스크립트는 제어 곱셈을 치환으로 직접 적용했고 N = 15 처럼 손으로 줄인 게이트 회로를 만들지 않았다 → 게이트 리스트가 없어 세 숫자를 정의하지 않는다")

# ── 게이트 비용 기록 (2): 컴파일하지 않은 회로 — Beauregard 2003 (quant-ph/0205095) 의 제어 모듈러 곱셈, N = 15 ──
# 위 7) 은 a = 7, N = 15 의 답(위수 4)을 알고 줄인 회로다. 여기서는 a 와 N 을 숫자로만 받는 일반 구성을 그대로 짜서 센다.
# 구성(논문 그림 3~7): Fourier 덧셈 φADD(a) → 이중 제어 모듈러 덧셈 φADD(a)MOD(N) → CMULT(a)MOD(N) → 제어 U_a = CMULT(a)·제어 SWAP·CMULT(a^{-1})†.
# 큐비트 2n + 3 = 11: 제어 c = 0 | x = 1..4 (x[0] 최상위) | b = 5..9 (n + 1 개, b[0] 최상위 = 넘침 비트) | 보조 = 10.
# 규약에 하나 추가: 이중 제어 위상 CCP(θ) = CP(θ/2)(c2,t)·CX(c1,c2)·CP(−θ/2)(c2,t)·CX(c1,c2)·CP(θ/2)(c1,t) (CP 3 + CX 2). QFT 는 SWAP 없이.
_n, _N = 4, 15; _C, _X, _B, _A, _NQ = 0, [1, 2, 3, 4], [5, 6, 7, 8, 9], 10, 11; _m = _n + 1
def _ccp(c1, c2, t, th):
    if gc._is_multiple(th, 2*np.pi): return []
    return [_g("cp", c2, t, angle=th/2), _g("cx", c1, c2), _g("cp", c2, t, angle=-th/2), _g("cx", c1, c2), _g("cp", c1, t, angle=th/2)]
def _qft_b():                                     # b 레지스터의 QFT (SWAP 없음)
    g = []
    for j in range(_m):
        g.append(_g("h", _B[j]))
        for k in range(j+1, _m): g.append(_g("cp", _B[k], _B[j], angle=2*np.pi/2**(k-j+1)))
    return g
def _phi_add(a, ctrls=()):                        # Fourier 공간에서 b ← b + a: b[j] 에 위상 2π·a/2^{m−j}
    g = []
    for j in range(_m):
        th = 2*np.pi*a/2**(_m-j)
        if len(ctrls) == 0: g.append(_g("p", _B[j], angle=th))
        elif len(ctrls) == 1: g.append(_g("cp", ctrls[0], _B[j], angle=th))
        else: g += _ccp(ctrls[0], ctrls[1], _B[j], th)
    return g
def _mod_add(a, c1, c2):                          # 이중 제어 φADD(a)MOD(N) (논문 그림 5)
    Q = _qft_b(); Qd = gc.inverse(Q)
    g = _phi_add(a, (c1, c2)) + gc.inverse(_phi_add(_N))
    g += Qd + [_g("cx", _B[0], _A)] + Q + _phi_add(_N, (_A,)) + gc.inverse(_phi_add(a, (c1, c2)))
    g += Qd + [_g("x", _B[0]), _g("cx", _B[0], _A), _g("x", _B[0])] + Q + _phi_add(a, (c1, c2))
    return g
def _cmult(a, c):                                 # b ← b + a·x mod N (c = 1 일 때)
    g = _qft_b()
    for i in range(_n): g += _mod_add((2**i * a) % _N, c, _X[_n-1-i])
    return g + gc.inverse(_qft_b())
def beauregard_cU(a, c=_C):                       # 제어 U_a: x ← a·x mod N
    g = _cmult(a, c) + [_g("cswap", c, _X[i], _B[i+1]) for i in range(_n)]
    return g + gc.inverse(_cmult(pow(a, -1, _N), c))
def _check_cU(a):                                 # b = 0, 보조 = 0 인 입력 2×15 개 전부: 제어 치환과 같고 열 사이 상대 위상이 없는가
    g = beauregard_cU(a); cols = [(cb << 10) | (xv << 6) for cb in (0, 1) for xv in range(_N)]
    tgts = [(cb << 10) | (((a*xv) % _N if cb else xv) << 6) for cb in (0, 1) for xv in range(_N)]
    out = gc.apply(g, np.eye(2**_NQ, dtype=complex)[:, cols], _NQ); ph = np.array([out[tg, k] for k, tg in enumerate(tgts)])
    return g, bool(abs(abs(ph[0]) - 1) < 1e-8 and np.allclose(ph, ph[0], rtol=0, atol=1e-8))
print("\n8) 게이트 비용 기록 (2) — 컴파일하지 않은 Beauregard 2003 회로, N = 15 (2n + 3 = 11큐비트). 추가 규약: 이중 제어 위상 = CP 3 + CX 2, QFT 는 SWAP 없이")
_bg = {}
for a_ in (7, 4, 1):
    _bg[a_], ok_b = _check_cU(a_); assert ok_b
    print(f"   제어 ×{a_} mod 15 한 번: {gc.fmt(_bg[a_])}  [입력 2×15 개 전부에서 제어 치환과 일치, 상대 위상 없음: {ok_b}]")
_seq = []
for j in range(2*_n): _seq += _bg[pow(7, 2**j, 15)]                 # 7^(2^j) mod 15 = 7, 4, 1, 1, … 을 숫자로만 넣는다(×1 도 그대로)
_cmp = gc.two_qubit_count(order_finding_list(2*_n)) - (2*_n)*(2*_n - 1) - 3*_n      # 컴파일된 회로의 제어 곱셈 부분(역 QFT 제외) = 28 + 16
print(f"   제어 곱셈 2n = 8번(a = 7; 역 QFT 제외, 제어 큐비트 하나를 되쓰는 구성): {gc.fmt(_seq)}")
print(f"   ×1 을 빼고 ×7·×4 만 넣으면(7⁴ ≡ 1 을 알아야 한다): {gc.fmt(_bg[7] + _bg[4])};  컴파일된 회로의 같은 부분은 2큐비트 {_cmp} → {gc.two_qubit_count(_seq)/_cmp:.0f}배")
_nadd = (2*_n)*(2*_n); print(f"   구조 검산: 모듈러 덧셈기 2n·2n = {_nadd}개, 덧셈기마다 {_m}큐비트 QFT 4번(CP {_m*(_m-1)//2}개씩) → QFT 의 CP 만 {4*_nadd*_m*(_m-1)//2}개 = CX {4*_nadd*_m*(_m-1)}")
