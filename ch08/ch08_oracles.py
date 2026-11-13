"""8장 오라클과 위상 되돌림 — Deutsch–Jozsa, Bernstein–Vazirani, Simon 을 numpy 로 조립해 검증한다.
실행: python ch08_oracles.py
환경: Python 3.13, numpy 2.5 (SDK 없이 동작). 큐비트 순서: 빅 엔디안(3장 규약). 입력 레지스터 q[0:n], 보조 q[n:].
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

def kron(*ops):
    out = np.array([[1]], dtype=complex) if ops[0].ndim == 2 else np.array([1], dtype=complex)
    for o in ops: out = np.kron(out, o)
    return out
def single(g, q, n): ops = [I]*n; ops[q] = g; return kron(*ops)
def cnot(c, t, n):
    a = [I]*n; a[c] = P0; b = [I]*n; b[c] = P1; b[t] = X; return kron(*a) + kron(*b)
def ket(bits):
    v = np.array([1], dtype=complex)
    for b in bits: v = np.kron(v, np.array([0, 1]) if b == '1' else np.array([1, 0]))
    return v
def Hn(n): return kron(*([H]*n))
def probs(psi): return abs(psi)**2
def bits(i, n): return format(i, f'0{n}b')

# ── 1. 오라클 만들기: 진리표 → 치환 행렬 ────────────────────────────
def oracle_bitflip(f, n, m=1):
    """U_f |x>|y> = |x>|y ⊕ f(x)>.  f: int → int (m비트).  (n+m)큐비트 치환 행렬."""
    N = 2**(n+m); U = np.zeros((N, N), dtype=complex)
    for x in range(2**n):
        for y in range(2**m):
            U[(x << m) | (y ^ f(x)), (x << m) | y] = 1
    return U
def oracle_phase(f, n):
    """U_f |x> = (-1)^{f(x)} |x>.  보조 큐비트 없는 버전 (연습문제 1)."""
    return np.diag([(-1)**f(x) for x in range(2**n)]).astype(complex)

# 위상 되돌림: 비트 뒤집기 오라클 + 보조 |-> == 위상 오라클 ⊗ |->
n = 3
f_rand = lambda x: int(rng.integers(0, 2))      # 아무 함수나 — 고정하려고 표로 만든다
table = [int(rng.integers(0, 2)) for _ in range(2**n)]
f = lambda x: table[x]
Uf = oracle_bitflip(f, n); Uph = oracle_phase(f, n)
minus = H @ np.array([0, 1], dtype=complex)
lhs = Uf @ np.kron(np.eye(2**n), minus.reshape(2, 1))            # |x>|-> 에 U_f
rhs = np.kron(Uph, minus.reshape(2, 1))                        # ((-1)^f |x>) |->
print("1) 위상 되돌림: U_f (I⊗|->) == (U_phase ⊗ |->) ?", np.allclose(lhs, rhs))
print("   진리표 f =", table, "→ 위상 오라클 대각 =", np.real(np.diag(Uph)).astype(int))
assert np.allclose(lhs, rhs), "5절 위상 되돌림: U_f·(I⊗|->) = U_φ⊗|->"
assert table == [1, 1, 1, 1, 1, 1, 1, 0], "5절 위상 되돌림: 난수 진리표 f = [1,1,1,1,1,1,1,0]"
assert np.real(np.diag(Uph)).astype(int).tolist() == [-1, -1, -1, -1, -1, -1, -1, 1], "5절 위상 되돌림: 위상 오라클의 대각"

# ── 2. Deutsch–Jozsa ───────────────────────────────────────────────
def deutsch_jozsa(f, n):
    """반환: 측정 확률 중 |0...0> 의 확률. 1 이면 상수, 0 이면 균형."""
    N = n + 1
    psi = ket('0'*n + '1')                       # 보조는 |1> → H 로 |->
    psi = Hn(N) @ psi
    psi = oracle_bitflip(f, n) @ psi
    psi = kron(Hn(n), I) @ psi
    p = probs(psi)
    p_in = p.reshape(2**n, 2).sum(axis=1)        # 보조 큐비트 추적
    return p_in[0], p_in
# 상수 함수 2개, 균형 함수 전부 (3비트: C(8,4)=70개)
n = 3
ok = True; cnt_c = cnt_b = 0
for c in (0, 1):
    p0, _ = deutsch_jozsa(lambda x, c=c: c, n); ok &= np.isclose(p0, 1.0); cnt_c += 1
for ones in __import__('itertools').combinations(range(2**n), 2**(n-1)):
    s = set(ones); p0, _ = deutsch_jozsa(lambda x, s=s: int(x in s), n); ok &= np.isclose(p0, 0.0); cnt_b += 1
print(f"\n2) Deutsch–Jozsa n={n}: 상수 {cnt_c}개 → P(0..0)=1, 균형 {cnt_b}개 → P(0..0)=0, 전부 통과? {ok}")
# 균형도 상수도 아닌 f (연습문제 2)
f_odd = lambda x: int(x == 0)                    # 8개 중 1개만 1
p0, p_in = deutsch_jozsa(f_odd, n)
print(f"   균형도 상수도 아닌 f (f(x)=1 은 x=0 뿐): P(0..0) = {p0:.4f}  (전체 분포: {np.round(p_in, 4)})")
assert ok and cnt_c == 2 and cnt_b == 70, "5절 Deutsch–Jozsa: 상수 2개 P(000) = 1, 균형 70개 P(000) = 0"
assert abs(p0 - 0.5625) < 5e-5 and abs(p0 - (6/8)**2) < 1e-12, "5절·8절 문제 2: 약속을 어긴 f 의 P(000) = 0.5625 = (6/8)²"
assert np.allclose(p_in[1:], 0.0625), "5절·8절 문제 2: 나머지 7개 y 에 0.0625 씩"

# ── 3. Bernstein–Vazirani ──────────────────────────────────────────
def bernstein_vazirani(s, n):
    f = lambda x: bin(x & s).count('1') % 2      # f(x) = s·x mod 2
    N = n + 1
    psi = Hn(N) @ ket('0'*n + '1')
    psi = oracle_bitflip(f, n) @ psi
    psi = kron(Hn(n), I) @ psi
    p_in = probs(psi).reshape(2**n, 2).sum(axis=1)
    return int(np.argmax(p_in)), p_in[np.argmax(p_in)]
n = 6; ok = True
for s in range(2**n):
    got, p = bernstein_vazirani(s, n); ok &= (got == s and np.isclose(p, 1.0))
print(f"\n3) Bernstein–Vazirani n={n}: 비밀 s 64가지 전부 1쿼리에 확률 1 로 복원? {ok}")
print(f"   예: s = {bits(0b101101, n)} → 측정 {bits(bernstein_vazirani(0b101101, n)[0], n)}  (고전은 n={n}번 쿼리 필요)")
assert ok and n == 6, "5절 Bernstein–Vazirani: n = 6, 비밀 s 64가지 전부 1쿼리·확률 1 로 복원"
assert bits(bernstein_vazirani(0b101101, n)[0], n) == "101101", "5절 Bernstein–Vazirani: s = 101101 → 측정 101101"

# ── 4. Simon ───────────────────────────────────────────────────────
def simon_oracle_fn(s, n):
    """f(x) = f(x⊕s) 인 2대1 함수를 난수로 만든다 (s≠0). 반환 f: int→int (n비트)."""
    vals = {}; next_val = list(rng.permutation(2**n))
    for x in range(2**n):
        if x in vals: continue
        v = next_val.pop(); vals[x] = v; vals[x ^ s] = v
    return lambda x: vals[x]
def simon_sample(f, n):
    """회로 1회 실행 → y (y·s = 0 인 n비트 문자열) 를 Born 규칙으로 샘플링."""
    N = 2*n
    psi = kron(Hn(n), np.eye(2**n)) @ ket('0'*N)
    psi = oracle_bitflip(f, n, m=n) @ psi
    psi = kron(Hn(n), np.eye(2**n)) @ psi
    p = probs(psi)
    idx = rng.choice(len(p), p=p/p.sum())
    return idx >> n                               # 입력 레지스터 값 y
def gf2_rank(rows, n):
    """GF(2) 가우스 소거. rows: int 리스트. 반환 (rank, 기약행 리스트)"""
    basis = []
    for r in rows:
        for b in basis:
            r = min(r, r ^ b)
        if r: basis.append(r); basis.sort(reverse=True)
    return len(basis), basis
def gf2_nullspace_1d(basis, n):
    """rank = n-1 인 basis 의 영공간 생성자 s 를 전수 탐색으로 (n ≤ 12 가정)"""
    for s in range(1, 2**n):
        if all(bin(s & b).count('1') % 2 == 0 for b in basis): return s
    return None
def simon_solve(f, n, max_queries=None):
    """y 를 모아 rank 가 n-1 이 될 때까지 쿼리. 반환 (s 추정, 쿼리 수)"""
    rows = []; q = 0
    while True:
        rows.append(simon_sample(f, n)); q += 1
        rank, basis = gf2_rank(rows, n)
        if rank == n - 1: return gf2_nullspace_1d(basis, n), q
        if max_queries and q >= max_queries: return None, q
n = 5
s_true = 0b10110
f = simon_oracle_fn(s_true, n)
ys = [simon_sample(f, n) for _ in range(8)]
print(f"\n4) Simon n={n}, s = {bits(s_true, n)}: 샘플 y 8개 = {[bits(y, n) for y in ys]}")
print("   모두 y·s = 0 mod 2 ?", all(bin(y & s_true).count('1') % 2 == 0 for y in ys))
assert [bits(y, n) for y in ys] == ['01111', '11100', '10100', '10010', '01111', '00111', '00000', '00111'], "5절 Simon: 샘플 y 8개"
assert all(bin(y & s_true).count('1') % 2 == 0 for y in ys), "5절 Simon: 모든 샘플이 y·s = 0 (mod 2)"
trials = 300; qs = []; succ = 0
for _ in range(trials):
    f = simon_oracle_fn(s_true, n)
    s_hat, q = simon_solve(f, n); qs.append(q); succ += (s_hat == s_true)
qs = np.array(qs)
print(f"   {trials}회 반복: 성공 {succ}/{trials}, 필요 쿼리 평균 {qs.mean():.2f}, 중앙값 {np.median(qs):.0f}, 최대 {qs.max()}  (n-1 = {n-1} 이 하한)")
hist = {k: int((qs == k).sum()) for k in range(n-1, qs.max()+1)}
print("   쿼리 수 분포:", hist)
# n 개의 쿼리 안에 끝날 확률 (연습문제 3): 이론은 ∏_{k=1}^{n-1}(1-2^{-k}) 이상
p_theory = np.prod([1 - 2.0**(-k) for k in range(1, n)])
print(f"   n-1 = {n-1} 번 만에 끝날 확률: 표본값 {(qs == n-1).mean():.3f}, 이론 ∏(1-2^-k) = {p_theory:.3f};  n+10 번 안에: 실측 {(qs <= n+10).mean():.4f}")
assert trials == 300 and succ == 300, "5절 Simon 표: s 복원 성공 300/300"
assert int(qs.sum()) == 1677 and abs(qs.mean() - 5.59) < 5e-3, "5절 Simon 표: 필요 쿼리 수 평균 5.59 (합 1677/300)"
assert np.median(qs) == 5 and qs.max() == 13, "5절 Simon 표: 중앙값 5, 최대 13"
assert hist == {4: 93, 5: 84, 6: 57, 7: 29, 8: 15, 9: 10, 10: 5, 11: 5, 12: 1, 13: 1}, "5절 Simon: 쿼리 수 분포(4번 93, 5번 84, 6번 57, 7번 29, 8번 15, … 13번 1)"
assert abs((qs == n-1).mean() - 0.310) < 5e-4 and abs(p_theory - 0.308) < 5e-4, "5절 Simon 표: n−1 = 4번에 끝난 비율 0.310 (이론 0.308)"
assert (qs <= n+10).mean() == 1.0, "5절 Simon 표: n+10 = 15번 안에 끝난 비율 1.000"
assert abs(sum(1/(1 - 2.0**(-j)) for j in range(1, n)) - 5.54) < 5e-3, "8절 문제 3 (b): 기대 쿼리 수 2 + 1.333 + 1.143 + 1.067 = 5.54"

# ── 5. 노이즈: 오라클 뒤 각 큐비트에 위상 감쇠(7장) → BV 성공률 ──────────
def phase_damp_all(rho, lam, N):
    K0 = np.array([[1, 0], [0, np.sqrt(1-lam)]]); K1 = np.array([[0, 0], [0, np.sqrt(lam)]])
    for q in range(N):
        A = single(K0, q, N); B = single(K1, q, N)
        rho = A @ rho @ A.conj().T + B @ rho @ B.conj().T
    return rho
def bv_noisy(s, n, lam):
    f = lambda x: bin(x & s).count('1') % 2
    N = n + 1
    psi = Hn(N) @ ket('0'*n + '1'); psi = oracle_bitflip(f, n) @ psi
    rho = np.outer(psi, psi.conj())
    rho = phase_damp_all(rho, lam, N)
    U = kron(Hn(n), I); rho = U @ rho @ U.conj().T
    p = np.real(np.diag(rho)).reshape(2**n, 2).sum(axis=1)
    return p[s]
n = 4; s = 0b1011
lams = [0, 0.02, 0.05, 0.1, 0.2, 0.4]
p_succ = [bv_noisy(s, n, l) for l in lams]
print(f"\n5) BV n={n} 노이즈: 오라클 뒤 각 큐비트 위상 감쇠 λ → 정답 확률")
for l, p in zip(lams, p_succ): print(f"   λ={l:.2f}: P(s)={p:.4f}   ((1+√(1-λ))/2)^n = {((1+np.sqrt(1-l))/2)**n:.4f}")
for l, p in zip(lams, p_succ):
    assert abs(p - ((1+np.sqrt(1-l))/2)**n) < 1e-9, "5절 노이즈 표: P(s) = ((1+√(1−λ))/2)ⁿ 이 정확히 맞는다"
    if l in (0, 0.05, 0.1, 0.2, 0.4):
        assert abs(p - {0: 1.0000, 0.05: 0.9503, 0.1: 0.9012, 0.2: 0.8050, 0.4: 0.6198}[l]) < 5e-5, "5절 노이즈 표: BV(n = 4, s = 1011) 정답 확률"

# ── 6. 그림 ───────────────────────────────────────────────────────
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Noto Sans CJK JP"; plt.rcParams["axes.unicode_minus"] = False
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 3.4))
_, p_c = deutsch_jozsa(lambda x: 1, 3); _, p_b = deutsch_jozsa(lambda x: bin(x).count('1') % 2, 3); _, p_o = deutsch_jozsa(f_odd, 3)
xs = np.arange(8); w = 0.27
ax1.bar(xs - w, p_c, w, label="상수 f", color="#4C6EF5"); ax1.bar(xs, p_b, w, label="균형 f (패리티)", color="#F03E3E"); ax1.bar(xs + w, p_o, w, label="둘 다 아님", color="#868E96")
ax1.set_xticks(xs); ax1.set_xticklabels([bits(i, 3) for i in xs], fontsize=7); ax1.set_ylabel("측정 확률"); ax1.set_title("Deutsch–Jozsa 출력 분포 (n=3)"); ax1.legend(frameon=False, fontsize=8)
ks = sorted(hist); ax2.bar(ks, [hist[k] for k in ks], color="#4C6EF5"); ax2.set_xlabel("필요 쿼리 수"); ax2.set_ylabel("횟수 (300회)"); ax2.set_title(f"Simon n={n+1}: rank n−1 까지 쿼리 수")
ax3.plot(lams, p_succ, 'o-', color="#4C6EF5", label="BV 정답 확률"); ax3.plot(lams, [((1+np.sqrt(1-l))/2)**4 for l in lams], '--', color="#868E96", label="((1+√(1−λ))/2)^n")
ax3.set_xlabel("위상 감쇠 λ (큐비트당)"); ax3.set_title("노이즈가 위상 되돌림을 지운다"); ax3.legend(frameon=False, fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "ch08_results.png"), dpi=150); print("6) saved ch08_results.png")
