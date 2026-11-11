"""run_snippets.py — env/extract_snippets.py 가 뽑은 SDK 스니펫을 하나씩 돌려 통과·실패를 표로 남긴다.
실행: python env/run_snippets.py [--only ch10] [--sdk qiskit] [--hardware] [--gpu] [--timeout 300]
결과: runs/snippets-YYYYMMDD.md 와 .json (커밋 대상 — "미실행 → 실행 완료"로 바꾸는 근거가 된다)

- 필요한 SDK가 설치돼 있지 않은 스니펫은 "건너뜀(미설치)"로 적는다.
- 실기 호출(QiskitRuntimeService 등)이 들어 있는 스니펫은 --hardware 를 줄 때만 돌린다(무료 분량을 쓴다).
- GPU 대상(nvidia, tensornet*)을 쓰는 CUDA-Q 스니펫은 --gpu 를 줄 때만 돌린다.
- NameError 는 둘 중 하나다: 본문용 조각이라 앞 블록의 변수를 쓰는 경우, 또는 오타. 출처 줄을 보고 가린다.
- 폐기 예정 경고(DeprecationWarning)는 통과여도 따로 적는다 — "API 가 바뀌는 중"과 "스니펫 오류"를 가르는 단서다.
"""
import argparse
import importlib.util
import json
import os
import pathlib
import re
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parent.parent
MODULE = {"qiskit": "qiskit", "pennylane": "pennylane", "cirq": "cirq", "qutip": "qutip", "cudaq": "cudaq",
          "braket": "braket", "qualtran": "qualtran"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", help="장 이름으로 거르기 (예: ch10)")
    ap.add_argument("--sdk", help="SDK 이름으로 거르기 (예: qiskit)")
    ap.add_argument("--hardware", action="store_true", help="실기 호출이 든 스니펫도 돌린다")
    ap.add_argument("--gpu", action="store_true", help="GPU 대상(nvidia, tensornet*)을 쓰는 스니펫도 돌린다")
    ap.add_argument("--timeout", type=int, default=300)
    a = ap.parse_args()
    idx_path = ROOT / "snippets" / "index.json"
    if not idx_path.exists():
        sys.exit("snippets/index.json 이 없다. 먼저 python env/extract_snippets.py 를 돌린다.")
    index = json.loads(idx_path.read_text(encoding="utf-8"))
    env = dict(os.environ, MPLBACKEND="Agg", PYTHONWARNINGS="default")
    rows = []
    for e in index:
        if a.only and e["chapter"] != a.only: continue
        if a.sdk and a.sdk not in e["sdks"]: continue
        missing = [s for s in e["sdks"] if importlib.util.find_spec(MODULE[s]) is None]
        if missing:
            rows.append({**e, "status": "건너뜀(미설치)", "detail": ", ".join(missing), "seconds": 0}); continue
        if e["hardware"] and not a.hardware:
            rows.append({**e, "status": "건너뜀(실기)", "detail": "--hardware 로 실행", "seconds": 0}); continue
        src = (ROOT / e["file"]).read_text(encoding="utf-8")
        if "cudaq" in e["sdks"] and re.search(r"set_target\(\s*[\"'](nvidia|tensornet)", src) and not a.gpu:
            rows.append({**e, "status": "건너뜀(GPU)", "detail": "--gpu 로 실행", "seconds": 0}); continue
        t0 = time.time()
        try:
            r = subprocess.run([sys.executable, str(ROOT / e["file"])], cwd=ROOT / e["chapter"], env=env,
                               capture_output=True, text=True, timeout=a.timeout)
            ok = r.returncode == 0
            err_lines = r.stderr.strip().split("\n")
            last = err_lines[-1] if not ok else r.stdout.strip().split("\n")[-1]
            dep = sorted({ln.split(":", 2)[2].strip()[:100] for ln in err_lines if "DeprecationWarning" in ln and ln.count(":") >= 2})
            if ok:
                status = "통과(폐기 경고)" if dep else "통과"
            elif last.startswith("ModuleNotFoundError"):
                status = "건너뜀(미설치)"
            elif last.startswith("NameError"):
                status = "NameError(조각 또는 오타)"
            else:
                status = "실패"
            if dep:
                last = (last + " ‖ 경고: " + "; ".join(dep))[:300]
        except subprocess.TimeoutExpired:
            status, last = "시간 초과", f"{a.timeout}초"
        rows.append({**e, "status": status, "detail": last[:300], "seconds": round(time.time() - t0, 1)})
        print(f"{e['file']:<44} {status:<16} {rows[-1]['seconds']:>6}s  {last[:80]}")
    out = ROOT / "runs"; out.mkdir(exist_ok=True)
    stem = f"snippets-{time.strftime('%Y%m%d')}"
    (out / f"{stem}.json").write_text(json.dumps(rows, ensure_ascii=False, indent=1), encoding="utf-8")
    tally = {}
    for r in rows: tally[r["status"]] = tally.get(r["status"], 0) + 1
    lines = ["# SDK 스니펫 실행 보고", "", f"- 파이썬 {sys.version.split()[0]}, 실행 시각 {time.strftime('%Y-%m-%d %H:%M')}",
             "- 집계: " + ", ".join(f"{k} {v}" for k, v in sorted(tally.items())), "",
             "| 파일 | 출처 | SDK | 상태 | 초 | 마지막 줄 |", "| --- | --- | --- | --- | --- | --- |"]
    for r in rows:
        detail = r["detail"].replace("|", "/")
        lines.append(f"| {r['file']} | {r['md']}:{r['line']} | {', '.join(r['sdks'])} | {r['status']} | {r['seconds']} | {detail} |")
    (out / f"{stem}.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n집계:", ", ".join(f"{k} {v}" for k, v in sorted(tally.items())), f"→ runs/{stem}.md")


if __name__ == "__main__":
    main()
