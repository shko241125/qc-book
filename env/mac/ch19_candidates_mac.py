"""맥(Apple silicon, numpy 의 BLAS/LAPACK 이 Accelerate) 전용 실행본: ch19/ch19_candidates.py 를 그대로 돌리되
주성분(SVD)의 부호를 본문 수치를 낸 기준 환경(리눅스 OpenBLAS·WSL)에 맞춘다. 원본은 고치지 않는다.

왜 필요한가(2026-11-13 확인). 특이벡터의 부호는 수학적으로 정해지지 않아 LAPACK 구현마다 다르게 나온다.
같은 numpy 2.5.4 에서 리눅스 arm64·x86_64 는 기준 출력(ch19_stdout.txt)과 같고, 맥은 (100×192) 의 2·3번째 성분과
(200×192) 의 2·4·6번째 성분 부호가 반대다. 특이값은 끝자리까지 같다. ZZ 특성 맵은 [0, π] 로 옮긴 특성 x 와 π − x 를
구별하므로 부호가 바뀌면 5절 (2) 화자 임베딩 표의 커널 값·1-NN·정렬·샷 수가 달라진다(이 표에는 assert 가 없어 원본도 "통과"로 끝난다).

사용(저장소 루트에서): python env/mac/ch19_candidates_mac.py > runs/regress/ch19_candidates_mac.out
그다음 diff ch19/ch19_stdout.txt runs/regress/ch19_candidates_mac.out — 실행 시간 줄만 달라야 한다.
그림(ch19/ch19_results.png)은 원본과 같이 덮어쓴다.
"""
import os
import runpy
import sys

import numpy as np

# 기준 환경에서 잰 부호: 행렬 모양별로, 앞쪽 주성분마다 "절댓값이 가장 큰 원소"의 부호(스크립트는 앞 10개까지만 쓴다).
REF_SIGN = {(100, 192): "-+++++++++", (200, 192): "+++-+--++-", (275, 192): "---++++++-"}
_svd = np.linalg.svd
FLIPS = []


def svd(a, *args, **kw):
    r = _svd(a, *args, **kw)
    pat = REF_SIGN.get(np.shape(a))
    if pat is None or not kw.get("compute_uv", True) or kw.get("full_matrices", True):
        return r
    U, S, Vt = r
    U, Vt = U.copy(), Vt.copy()
    for k, want in enumerate(pat):
        v = Vt[k]
        if ("+" if v[np.argmax(np.abs(v))] > 0 else "-") != want:
            Vt[k] *= -1; U[:, k] *= -1; FLIPS.append((np.shape(a), k))
    return U, S, Vt


np.linalg.svd = svd
here = os.path.dirname(os.path.abspath(__file__))
target = os.path.normpath(os.path.join(here, "..", "..", "ch19", "ch19_candidates.py"))
os.chdir(os.path.dirname(target)); sys.argv = [target]
os.environ.setdefault("MPLBACKEND", "Agg")
runpy.run_path(target, run_name="__main__")
print(f"[mac] 부호를 뒤집은 주성분: {FLIPS or '없음'}", file=sys.stderr)
