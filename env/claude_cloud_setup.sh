#!/bin/bash
# Claude Code 클라우드 환경의 "Setup script" 칸에 붙여 넣는 내용(환경 대화상자 → Setup script).
# 세션이 시작되기 전에 root 로 한 번 돌고, 설치 결과는 환경 캐시에 남는다. 0 으로 끝나야 세션이 뜨므로 실패해도 넘어가게 했다.
# 작성됨, 미실행 — 저장소 파일에 기대지 않도록 패키지 이름을 직접 적었다. 범위는 requirements/sdk.txt 와 같게 유지한다.
# 주의: 대략 5분 안에 끝나야 환경 캐시가 만들어진다. 기본 네트워크 단계(Trusted)는 PyPI 는 열어 주지만 IBM Quantum 주소는 없으므로,
#       장치에 제출하려면 환경의 허용 도메인에 IBM Quantum 주소를 더한다.
PIP="python3 -m pip install --break-system-packages"
$PIP "numpy>=2.0" "scipy>=1.14" "matplotlib>=3.9" "markdown>=3.6" || true
ALL='qiskit>=2.5,<3 qiskit-aer>=0.17,<0.19 qiskit-ibm-runtime>=0.42,<0.53 pennylane>=0.45,<0.47 cirq-core>=1.7,<2 qutip>=5.2,<6 qutip-qip>=0.4,<0.5 pylatexenc qualtran>=0.7,<0.8'
if ! $PIP $ALL; then                       # 한꺼번에 안 풀리면 묶음별로 넣고, 실패한 묶음을 적어 둔다
  for grp in "qiskit>=2.5,<3 qiskit-aer>=0.17,<0.19 qiskit-ibm-runtime>=0.42,<0.53" "pennylane>=0.45,<0.47" \
             "cirq-core>=1.7,<2" "qutip>=5.2,<6 qutip-qip>=0.4,<0.5" "pylatexenc" "qualtran>=0.7,<0.8"; do
    $PIP $grp || echo "설치 실패: $grp" >> /root/qc-book-setup-failures.txt
  done
fi
exit 0
