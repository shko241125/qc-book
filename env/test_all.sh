#!/usr/bin/env bash
# 책의 numpy 검증 스크립트(chNN/ 와 부록 appX/)를 전부 다시 돌린다(SDK 불필요). 저장소 루트에서: bash env/test_all.sh [--sdk] [ch10 appA ...]
# SDK를 불러오는 스크립트(예: ch14/ch14_submit.py — 인자 없이 돌면 fake 모드)는 --sdk 를 줄 때만 돌린다.
# 출력은 runs/regress/<스크립트 이름>.out 에 남고, 그림 파일은 덮어쓴다. 하나라도 실패하면 종료 코드 1.
set -u -o pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
SDK=0; args=()
for a in "$@"; do if [ "$a" = "--sdk" ]; then SDK=1; else args+=("$a"); fi; done
if [ -n "${PYTHON:-}" ]; then PY="$(command -v "$PYTHON")"; elif [ -x "$ROOT/.venv/bin/python" ]; then PY="$ROOT/.venv/bin/python"; else PY="$(command -v python3)"; fi
[ -n "$PY" ] || { echo "파이썬을 찾지 못했다"; exit 2; }
mkdir -p runs/regress; fail=0
targets=(${args[@]+"${args[@]}"}); [ ${#targets[@]} -eq 0 ] && targets=(ch*/ app*/)
for d in "${targets[@]}"; do
  d="${d%/}"
  [ -d "$d" ] || { echo "없는 폴더: $d"; fail=1; continue; }
  for f in "$d"/ch*.py "$d"/app*.py; do
    [ -f "$f" ] || continue
    name="$(basename "$f" .py)"
    sdk="$(grep -m1 -oE '^(from|import) +(qiskit|cirq|pennylane|qutip|cudaq|braket|qualtran)' "$f" | awk '{print $2}')"
    if [ -n "$sdk" ]; then
      if [ "$SDK" = 0 ]; then printf "%-34s %s\n" "$f" "건너뜀(SDK 스크립트 — --sdk 로 실행)"; continue; fi
      if ! "$PY" -c "import $sdk" >/dev/null 2>&1; then printf "%-34s %s\n" "$f" "건너뜀($sdk 미설치)"; continue; fi
    fi
    s=$(date +%s)
    if (cd "$d" && MPLBACKEND=Agg "$PY" "$name.py" > "$ROOT/runs/regress/$name.out" 2>&1); then st="통과"; else st="실패"; fail=1; fi
    printf "%-34s %s  %4ds\n" "$f" "$st" $(( $(date +%s) - s ))
  done
done
"$PY" tools/gate_cost.py | tail -1 || fail=1
exit $fail
