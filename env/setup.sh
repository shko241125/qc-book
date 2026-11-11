#!/usr/bin/env bash
# 책 실행 환경 설치. 저장소 루트에서:  bash env/setup.sh [--hardware] [--gpu] [--no-venv]
#   --hardware  AWS Braket SDK 도 설치 (env/run_braket.py 용)
#   --gpu       NVIDIA CUDA-Q 설치 (Linux 또는 Apple silicon macOS)
#   --no-venv   가상환경을 만들지 않고 "지금 쓰는 python3" 에 설치 (Codespaces·Docker·Claude Code 클라우드 세션용)
set -euo pipefail
cd "$(dirname "$0")/.."
HARDWARE=0; GPU=0; VENV=1
for a in "$@"; do case "$a" in
  --hardware) HARDWARE=1;; --gpu) GPU=1;; --no-venv) VENV=0;;
  -h|--help) sed -n '2,5p' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
  *) echo "모르는 옵션: $a  (bash env/setup.sh --help)"; exit 2;;
esac; done

# 파이썬 고르기. venv 를 만들 때는 3.12 → 3.13 → 3.11 → python3 순(cirq·cudaq·qiskit-aer 가 3.11 이상, Braket SDK 는 3.13 까지 표시).
# --no-venv 는 "지금 쓰는 파이썬"에 설치해야 하므로 python3 → python 을 먼저 본다(나머지 스크립트가 python3 를 쓴다).
if [ "$VENV" = 1 ]; then CANDS="python3.12 python3.13 python3.11 python3"; else CANDS="python3 python python3.12 python3.13 python3.11"; fi
PY=""
for c in $CANDS; do
  if command -v "$c" >/dev/null 2>&1 && "$c" -c 'import sys; sys.exit(0 if (3,11) <= sys.version_info[:2] else 1)'; then PY="$c"; break; fi
done
[ -n "$PY" ] || { echo "Python 3.11 이상이 필요하다. (예: uv python install 3.12, brew install python@3.12, apt install python3.12 python3.12-venv)"; exit 1; }
echo "파이썬: $("$PY" -V) ($(command -v "$PY"))"
"$PY" -c 'import sys; sys.exit(0 if sys.version_info[:2] <= (3,13) else 1)' || echo "주의: 3.14 이상이다. 일부 SDK 가 아직 지원하지 않을 수 있다 — 설치가 실패하면 3.12 를 쓴다."

if [ "$VENV" = 1 ]; then
  [ -x .venv/bin/python ] || { rm -rf .venv; "$PY" -m venv .venv; }     # 반쯤 만들어진 .venv 는 다시 만든다
  # shellcheck disable=SC1091
  source .venv/bin/activate
  PY=python
  "$PY" -m pip install --upgrade pip
elif "$PY" -c 'import os,sys,sysconfig; sys.exit(0 if os.path.exists(os.path.join(sysconfig.get_path("stdlib"),"EXTERNALLY-MANAGED")) else 1)'; then
  # PEP 668 은 표식 파일로 판정한다(네트워크·pip 버전에 기대지 않는다). 환경 변수는 install 에만 적용되고 freeze 는 그대로 돈다.
  export PIP_BREAK_SYSTEM_PACKAGES=1
  echo "외부 관리(PEP 668) 파이썬이다 — PIP_BREAK_SYSTEM_PACKAGES=1 로 설치한다(시스템 pip 는 올리지 않는다)."
fi
"$PY" -m pip install -r requirements/sdk.txt
TAG=""
if [ "$HARDWARE" = 1 ]; then "$PY" -m pip install -r requirements/hardware.txt; TAG="$TAG-hardware"; fi
if [ "$GPU" = 1 ]; then
  if "$PY" -m pip install -r requirements/gpu.txt; then TAG="$TAG-gpu"; else echo "cudaq 설치 실패 — 이 플랫폼을 지원하지 않을 수 있다(요구 사항은 requirements/gpu.txt 주석)"; fi
fi

LOCK="requirements/lock-$(date +%Y%m%d)$TAG.txt"
"$PY" -m pip freeze > "$LOCK"
echo "설치된 버전을 $LOCK 에 적었다."
"$PY" env/smoke_test.py --strict || echo "점검에서 빠지거나 실패한 SDK 가 있다(위 표)."
if [ "$VENV" = 1 ]; then echo "다음부터는  source .venv/bin/activate  로 환경을 켠다."; fi
exit 0
