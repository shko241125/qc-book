# R10 — Review of the execution-environment kit (`env/`, `requirements/`, Docker, devcontainer, docs)

Reviewer: R10 (independent). Date: 2026-11-10 UTC. Repo: `/home/claude/qc-book`. Nothing under the repo was edited by me
(`python3 env/extract_snippets.py` regenerated `snippets/`, as permitted).

**File state reviewed.** Four files were modified by another process while I was reviewing (mtime 23:21 UTC):
`env/test_all.sh`, `.gitignore`, `CLAUDE.md`, `README.md`. All line numbers below refer to the versions on disk at 23:28 UTC
(md5: test_all.sh `48417308…`, .gitignore `93c4ca8f…`, README.md `8e409a2b…`, CLAUDE.md `9367ad83…`, setup.sh `33d26067…`,
run_ibm.py `2417a5e3…`, run_braket.py `e3f8e7eb…`, smoke_test.py `d140f6b2…`). Three defects I found in the *earlier*
test_all.sh/.gitignore (log overwritten by the `"../$PY" || "$PY"` fallback, one `.out` per directory so `ch14_submit.py`
overwrote `ch14_hardware.py`'s log, `runs/*.json` ignored although CLAUDE.md calls it the evidence) are already fixed in
the 23:21 versions and are not repeated.

**What could and could not be run here.** No SDK is installed and PyPI is unreachable from the shell, so every SDK call
was checked against documentation only (web). Ran: `py_compile` on 5 scripts (OK), `bash -n` on 3 scripts (OK),
`smoke_test.py` (numpy row passes, 10 rows "미설치", exit 0; `--strict` exit 1), `extract_snippets.py` (37 snippets:
qiskit 16, pennylane 7, cirq 7, qutip 5, cudaq 1, qualtran 1; 1 hardware-flagged), `setup.sh --no-venv` in a scratch copy
(fails, see A1), `test_all.sh` only in a scratch fake repo (not the real one), unit checks of `load_dotenv`, `tvd`, the
Braket price lookup, and `git check-ignore` against a scratch copy of `.gitignore`. All proposed shell replacements were
syntax-checked and exercised in the scratch copy; proposed Python that needs an SDK was not executable.

Counts: **A = 4, B = 16, C = 12.**

---

## Grade A — will fail or mislead on the first real run

### A1. `env/setup.sh:27-31,37` — the `--no-venv` / PEP 668 branch always dies with "no such option"

```bash
27    PIP="$PY -m pip"
29    if ! $PIP install --dry-run numpy >/dev/null 2>&1; then PIP="$PIP --break-system-packages"; else PIP="$PY -m pip"; fi
31  $PIP install --upgrade pip        # -> python -m pip --break-system-packages install --upgrade pip
37  $PIP freeze > "$LOCK"             # -> python -m pip --break-system-packages freeze
```

`--break-system-packages` is an option of the `install` subcommand, not a global pip option, and `freeze` does not accept
it in any position. Reproduced here with both pip 24.0 (Ubuntu system python3.12) and pip 26.2.1:

```
$ python3 -m pip --break-system-packages install --dry-run numpy   ->  no such option: --break-system-packages
$ python3 -m pip freeze --break-system-packages                    ->  no such option: --break-system-packages
$ bash env/setup.sh --no-venv   (scratch copy, this sandbox)
파이썬: Python 3.12.3 (/usr/bin/python3.12)
no such option: --break-system-packages          exit=2
```

Two further problems in the same lines: (a) the probe on line 29 treats *any* failure of `pip install --dry-run numpy`
(no network, pip < 22.2 without `--dry-run`, proxy) as "externally managed"; (b) `pip install --upgrade pip` against a
distro-owned pip is unnecessary in this mode. The environment variable form works for `install` and is ignored by
`freeze` (verified: `PIP_BREAK_SYSTEM_PACKAGES=1 python3.12 -m pip freeze` runs on pip 24.0 and 26.2.1).

### A2. `env/setup.sh:13` — `--no-venv` installs into `python3.12`, not into the interpreter the rest of the kit uses

The chooser order `python3.12 python3.13 python3.11 python3` is right for creating a venv but wrong for
"현재 파이썬에 설치" (line 5). In this sandbox — Ubuntu 24.04, the image Claude Code cloud sessions run on, which line 5
names as a target — `python3.12` is `/usr/bin/python3.12` (system, EXTERNALLY-MANAGED, no numpy) while `python3`/`python`
is `/usr/local/bin/python3` 3.13.17 (has numpy 2.5.3). So even with A1 fixed, the SDKs go into 3.12, the script's own
smoke test (`"$PY" env/smoke_test.py`, line 39) passes, and then every documented command (`python env/smoke_test.py`,
`python env/run_snippets.py`, `env/test_all.sh`, and `env/claude_cloud_setup.sh` which uses `python3 -m pip`) runs on
3.13 and reports everything "미설치".

**Replacement for `env/setup.sh` lines 11-41** (fixes A1 + A2; syntax-checked and dry-run in a scratch copy, where it now
selects `/usr/local/bin/python3` 3.13.17):

```bash
# 파이썬 고르기. venv 를 만들 때는 3.12 → 3.13 → 3.11 → python3 순.
# --no-venv 는 "지금 쓰는 파이썬"에 설치해야 하므로 python3 → python 을 먼저 본다(나머지 스크립트가 python3 를 쓴다).
if [ "$VENV" = 1 ]; then CANDS="python3.12 python3.13 python3.11 python3"; else CANDS="python3 python python3.12 python3.13 python3.11"; fi
PY=""
for c in $CANDS; do
  if command -v "$c" >/dev/null 2>&1 && "$c" -c 'import sys; sys.exit(0 if (3,11) <= sys.version_info[:2] else 1)'; then PY="$c"; break; fi
done
[ -n "$PY" ] || { echo "Python 3.11 이상이 필요하다. (예: uv python install 3.12, brew install python@3.12, apt install python3.12)"; exit 1; }
echo "파이썬: $("$PY" -V) ($(command -v "$PY"))"

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
if [ "$HARDWARE" = 1 ]; then "$PY" -m pip install -r requirements/hardware.txt; fi
if [ "$GPU" = 1 ]; then "$PY" -m pip install -r requirements/gpu.txt || echo "cudaq 설치 실패 — 이 플랫폼을 지원하지 않을 수 있다(요구 사항은 requirements/gpu.txt 주석)"; fi

LOCK="requirements/lock-$(date +%Y%m%d).txt"
"$PY" -m pip freeze > "$LOCK"
echo "설치된 버전을 $LOCK 에 적었다."
"$PY" env/smoke_test.py --strict || echo "점검에서 빠지거나 실패한 SDK 가 있다(위 표)."
if [ "$VENV" = 1 ]; then echo "다음부터는  source .venv/bin/activate  로 환경을 켠다."; fi
exit 0
```

### A3. `env/run_ibm.py:96` and `:175` — the default Bell circuit reports 2 two-qubit gates (the `measure_all()` barrier is counted)

`build()` calls `qc.measure_all()` (line 75), which inserts a barrier across all qubits before the measurements. For
`n = 2` that barrier has `num_qubits == 2`, so `sum(1 for inst in isa.data if inst.operation.num_qubits == 2)` returns
2 (cx + barrier) for Bell and the correct n−1 for GHZ with n ≥ 3. The very first documented command
(`python env/run_ibm.py --local`, default `--circuit bell`) therefore prints "2큐비트 게이트 2", and the wrong value is
stored as `two_qubit_gates` in the `runs/*.json` that CLAUDE.md:38 designates as the source for "실측" entries. Barriers
survive transpilation, so the hardware path is affected too.
Basis: Qiskit's `measure_all` implementation (`circ.barrier()` then `circ.measure(...)`); the API pages I could fetch
do not mention the barrier, so confirm on the first machine with
`python -c "from qiskit import QuantumCircuit as Q; c=Q(2); c.measure_all(); print(c.count_ops())"` (expect `barrier: 1`).

Fix — add once and use in both places (logic checked with mock instructions: old count 2, new count 1):

```python
def two_qubit_gates(circ):
    """2큐비트 게이트 수. measure_all() 이 넣는 barrier(지시문)는 게이트가 아니므로 세지 않는다."""
    return sum(1 for inst in circ.data if inst.operation.num_qubits == 2 and inst.operation.name != "barrier")
```
```python
# run_ibm.py:96
    twoq = two_qubit_gates(isa)
# run_ibm.py:175
    twoq = two_qubit_gates(isa)
```

Related (C-level): `isa.depth()` includes the measurement layer, whereas CLAUDE.md:35 defines depth by the
`tools/gate_cost.py` convention (gates only). Either label it "깊이(측정 포함)" or use
`isa.depth(lambda inst: inst.operation.name not in ("barrier", "measure"))`.

### A4. `ch14/ch14_submit.py:103-107` (advertised at README.md:34, CLAUDE.md:14/26/45/53) — `real` mode loses the record after spending QPU time

```python
103    if mode != "fake":
105            m = j.metrics(); record.setdefault("quantum_seconds", []).append(m["usage"]["quantum_seconds"])
107    json.dump(record, open(f"ch14_record_{mode}.json", "w"), ...)
```

The qiskit-ibm-runtime 0.51.0 release notes state that jobs run on IBM Quantum Platform return `usage` containing
`qpu_charge_time_seconds` and `status` "instead of `quantum_seconds` and `seconds`" (local mode was changed to match).
So line 105 raises `KeyError` after both real jobs have finished, *before* line 107 writes the record, and the job IDs
were never printed (they exist only inside `record`). The counts are then recoverable only by finding the jobs in the
dashboard. Same stale key in `ch14/ch14_hardware.md:45, 101, 194, 199` and in the docstring at `ch14_submit.py:5`.
This file is outside the listed kit, but the 23:21 README/CLAUDE.md now route the reader to `ch14_submit.py real`.

Fix (replace lines 87-88 and 103-107):

```python
    record = {"backend": backend.name, "shots": shots, "t": time.strftime("%Y-%m-%d %H:%M"), "jobs": {}}
    job = sampler.run(list(isa.values()), shots=shots)   # 회로 4개를 한 작업(PUB 4개)으로 → 오버헤드 1회
    print("sampler job id:", getattr(job, "job_id", lambda: "fake")(), flush=True)   # 끊겨도 service.job(id) 로 되찾는다
```
```python
    def dump():
        with open(f"ch14_record_{mode}.json", "w") as fh:
            json.dump(record, fh, indent=1, ensure_ascii=False)

    dump()                                               # 사용 시간 조회가 실패해도 측정 결과는 남는다
    if mode != "fake":                                   # (3) 사용 시간 기록 — 14장 표의 '사용 분' 열
        for j in (job, ejob):
            try:
                u = j.metrics().get("usage", {})
                secs = u.get("qpu_charge_time_seconds", u.get("quantum_seconds"))   # 0.51 릴리스 노트: 플랫폼 키는 qpu_charge_time_seconds
                if secs is None:
                    secs = j.usage()                     # RuntimeJobV2.usage() -> float(초)
            except Exception as e:                       # 사용량 때문에 기록을 잃지 않는다
                secs = f"조회 실패: {type(e).__name__}: {e}"
            record.setdefault("qpu_seconds", []).append(secs)
        print("qpu_seconds:", record["qpu_seconds"])
        dump()
    return record
```
Also print the estimator job id after line 97 (`print("estimator job id:", getattr(ejob, "job_id", lambda: "fake")(), flush=True)`).
In the md, replace `job.metrics()['usage']['quantum_seconds']` by `job.usage()`.

---

## Grade B — wrong in an edge case, or misleading

### B1. `env/run_ibm.py:136, 160-171` — `--fetch` evaluates the job against the *command line's* circuit and re-transpiles
`qc, ideal = build(a.circuit, a.n)` runs before the fetch branch, so `--fetch <id>` of a GHZ job without repeating
`--circuit ghz --n 5` compares against the Bell ideal (I measured `tvd({'000':500,'111':500}, bell_ideal) == 1.0`),
names the file `…_bell2.json`, and line 170 re-runs the pass manager (no seed), so `two_qubit_gates`/`depth` describe a
new transpilation, not the circuit that ran. The docstring example on line 10 shows `--fetch <job id>` alone.

Fix: write a sidecar at submission and read it at fetch; refuse on width mismatch.
```python
# after line 180 (job submitted)
    pend = ROOT / "runs" / f"pending_{job.job_id()}.json"
    pend.parent.mkdir(exist_ok=True)
    pend.write_text(json.dumps({"job_id": job.job_id(), "backend": backend.name, "circuit": a.circuit, "n": qc.num_qubits,
                                "shots": a.shots, "two_qubit_gates": twoq, "depth": isa.depth(),
                                "ops": {k: int(v) for k, v in isa.count_ops().items()}}, indent=1), encoding="utf-8")
```
```python
# replace lines 160-171
    if a.fetch:
        pend = ROOT / "runs" / f"pending_{a.fetch}.json"
        meta = json.loads(pend.read_text(encoding="utf-8")) if pend.exists() else None
        if meta:
            qc, ideal = build(meta["circuit"], meta["n"])
        else:
            print("경고: runs/pending_<id>.json 이 없다. --circuit/--n 을 제출 때와 똑같이 줘야 TVD 가 맞는다.")
        job = service.job(a.fetch)
        print("상태:", job.status())
        counts = counts_of(job.result()[0])
        if len(next(iter(counts))) != qc.num_qubits:
            sys.exit(f"결과 비트 수({len(next(iter(counts)))})가 회로({qc.num_qubits}큐비트)와 다르다. 제출 때의 --circuit/--n 을 준다.")
        usage = final_usage(job)                                  # B4
        be = job.backend()
        isa = generate_preset_pass_manager(optimization_level=1, backend=be).run(qc)   # meta 가 없을 때만 참고값
        return report(job.job_id(), be.name if be else "unknown", qc, isa, counts, ideal, usage, sum(counts.values()), False, stats=meta)
```
and in `report(..., local, stats=None)` use `stats["two_qubit_gates"]`, `stats["depth"]`, `stats["ops"]` when `stats`
is given (else compute from `isa` and add `"transpile_note": "re-transpiled at fetch"` to the JSON).

### B2. `env/run_braket.py:61, 73-77` — same problem for `--fetch`, plus `None` result
`build(a.circuit, a.n)` supplies the ideal; `t.result()` returns `None` for a FAILED/CANCELLED task (AttributeError on
`.measurement_counts`); the device is saved as `"?"`. Fix:
```python
    if a.fetch:
        t = AwsQuantumTask(arn=a.fetch)
        state = t.state()
        print("상태:", state)
        if state in ("FAILED", "CANCELLED"):
            sys.exit(f"작업이 {state} 로 끝났다(요금 없음). 콘솔에서 사유를 확인한다.")
        r = t.result()
        if r is None:
            sys.exit("결과를 받지 못했다(대기 시간 초과). 나중에 다시 --fetch 한다.")
        counts = dict(r.measurement_counts)
        n_bits = len(next(iter(counts)))
        ideal = {"0" * n_bits: 0.5, "1" * n_bits: 0.5}            # Bell·GHZ 전용 스크립트이므로 비트 수로 정한다
        label = "bell2" if n_bits == 2 else f"ghz{n_bits}"
        device = t.metadata().get("deviceArn", "?")
        return save(a.fetch, device, label, sum(counts.values()), counts, ideal, None)
```
(`AwsQuantumTask.metadata()` returning the GetQuantumTask response with `deviceArn` is from my knowledge of the SDK, not
re-verified today.)

### B3. `SamplerV2` is deprecated since qiskit-ibm-runtime 0.50.0 and the range is open-ended
`run_ibm.py:135,142,180`, `smoke_test.py:66`, `requirements/sdk.txt:8`, `claude_cloud_setup.sh:8`.
Release notes 0.50.0 (2026-09-24): "`EstimatorV2` and `SamplerV2` are deprecated"; replacements are
`qiskit_ibm_runtime.executor_sampler.Sampler` / `executor_estimator.Estimator`, and "at the end of the deprecation period"
the top-level names will point to the client-side implementations. The current local-testing-mode guide already imports
`from qiskit_ibm_runtime.executor_sampler import Sampler` and labels the old import "Legacy server-side primitive before
qiskit-ibm-runtime 0.50.0". PyPI latest is **0.52.0 (2026-11-09)**; its release notes are not yet on the docs site, so I
could not check whether 0.52 changed anything further. `>=0.42` with no cap means the kit can silently switch
implementation on a future install. Today's code should still work (deprecated ≠ removed) — not verified by execution.
Fix in `requirements/sdk.txt:8` and `claude_cloud_setup.sh:8`:
```
qiskit-ibm-runtime>=0.42,<0.53     # SamplerV2/EstimatorV2 는 0.50 부터 폐기 예정 — 상한을 올리기 전에 run_ibm.py 를 --local 로 다시 확인
```
Do not switch to `executor_sampler.Sampler` blind: I could confirm `Sampler(mode=AerSimulator())` and `run([circuit])`
for it, but not `run(..., shots=)` nor `join_data()` on its result.

### B4. `env/run_ibm.py:165-169, 183-188` — `job.usage()` right after `result()` can be 0
`RuntimeJobV2.usage(partial=False)` → float seconds; docs: "By default, the job usage returned is `0` until the usage
calculation is completed." The script would print "QPU 사용량: 0.0 초" and store `usage_seconds: 0.0`. Fix:
```python
def final_usage(job, tries=6, wait=5):
    """job.usage() 는 집계가 끝나기 전에는 0 을 준다. 잠깐 기다렸다가 그래도 0 이면 None(미확정)으로 적는다."""
    import time
    for _ in range(tries):
        try:
            u = job.usage()
        except Exception:
            return None
        if u:
            return u
        time.sleep(wait)
    return None
```
use `usage = final_usage(job)` in both places and print "미확정 — 대시보드에서 확인" when it is `None`.

### B5. `env/smoke_test.py` — the bit-order check never fails
Docstring lines 7-8 and CLAUDE.md:22 promise a bit-order check ("'01'로 … '10'으로 보여야 한다"), but `order` is only
interpolated into the note; the status depends on `bell_ok` alone (lines 53, 85, 94, 104, 126, 135). A wrong order would
still print "통과". Fix each test's return, e.g.:
```python
# t_qiskit (line 53)
    return version("qiskit"), bell_ok(p, 1e-9) and order == "01", f"{order} (리틀 엔디안: 01 이어야 한다)"
# t_pennylane (85)
    return version("pennylane"), bool(np.allclose(p, REF_BELL)) and order == "10", f"{order} (빅 엔디안: 10 이어야 한다)"
# t_cirq (94)
    return version("cirq-core"), bool(np.allclose(abs(sv) ** 2, REF_BELL, atol=1e-6)) and order == "10", f"{order} (빅 엔디안: 10 이어야 한다)"
# t_qutip (104)
    return qutip.__version__, bool(np.allclose(p, REF_BELL)) and order == "10", f"{order} (tensor 순서: 10 이어야 한다)"
# t_cudaq (126)
    return version("cudaq"), bell_ok(probs_from_counts(counts, 4000), 0.05) and order == "10", f"{order} (q0 가 왼쪽: 10 이어야 한다); 대상 {cudaq.get_target().name}"
# t_braket (135)
    return version("amazon-braket-sdk"), bell_ok(probs_from_counts(dict(r.measurement_counts), 4000), 0.05) and order == "10", f"{order} (빅 엔디안: 10 이어야 한다)"
```
The expected strings themselves are right (see "Verified correct").

### B6. `env/run_braket.py` never reads `.env`; `.env.example:5-6` implies it does
`.env.example` says "AWS Braket: … (여기에는 리전만)  AWS_DEFAULT_REGION=us-east-1", but only `run_ibm.py` has
`load_dotenv()`. Outside Docker (`--env-file`) the value is ignored; without a region in `~/.aws/config`,
`AwsDevice.get_devices`/`AwsDevice(arn)` fail in boto3 with `NoRegionError` (boto3 behaviour from my knowledge; the
Braket README does say "Make sure that your AWS region is set"). Fix: add the same `load_dotenv()` (B16 version) to
`run_braket.py` and call it at the top of `main()`, or at minimum:
```python
    import os
    os.environ.setdefault("AWS_DEFAULT_REGION", "us-east-1")     # --list 는 어차피 모든 Braket 리전을 훑는다
```

### B7. `env/test_all.sh:18-23` runs an unexecuted SDK script inside the "numpy regression"; docs still say numpy-only
With qiskit installed (i.e. after the documented `setup.sh`), `ch14/ch14_submit.py` is executed with no argument →
`fake` mode (no hardware submission — good), but it is itself "작성됨, 미실행", simulates a 156-qubit fake backend's noise
model, and writes `ch14/ch14_record_fake.json` (tracked). A failure there turns the whole regression red (exit 1) although
all numpy scripts passed. README.md:20 ("numpy 스크립트(SDK 불필요, 10~15분)") and CLAUDE.md:23 contradict the script's
own header (lines 3-4). Also line 27: `"$PY" tools/gate_cost.py | tail -1 || fail=1` can never set `fail` (pipeline
status is `tail`'s; no `pipefail`), and `.venv` silently overrides an explicit `PYTHON=`.

Replacement for the whole file (tested in a scratch fake repo: SDK script skipped by default, missing folder and a
failing `gate_cost.py` both give exit 1):
```bash
#!/usr/bin/env bash
# 책의 numpy 검증 스크립트를 전부 다시 돌린다(SDK 불필요). 저장소 루트에서: bash env/test_all.sh [--sdk] [ch10 ch11 ...]
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
targets=(${args[@]+"${args[@]}"}); [ ${#targets[@]} -eq 0 ] && targets=(ch*/)
for d in "${targets[@]}"; do
  d="${d%/}"
  [ -d "$d" ] || { echo "없는 폴더: $d"; fail=1; continue; }
  for f in "$d"/ch*.py; do
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
```

### B8. `Dockerfile:1`, `README.md:41` — `docker run … --env-file .env` on a fresh clone
`.env` does not exist until the user copies `.env.example`; `docker run --env-file` with a missing file aborts
("open .env: no such file or directory" — Docker CLI behaviour from my knowledge, not run here). Fix both lines:
```
docker build -t qc-book . && docker run --rm -it -v "$PWD":/work $([ -f .env ] && echo --env-file .env) qc-book
```
Also note in `.env.example` that values must be **unquoted**: Docker's `--env-file` passes quotes literally, whereas
`load_dotenv()` strips them, so a quoted token works natively and fails in Docker.

### B9. `env/claude_cloud_setup.sh:5-10` — one resolver failure installs nothing, silently
All packages are in one `pip install … || true`. If any single requirement cannot be resolved (e.g. a qualtran
dependency conflict) pip installs none of them and the script still exits 0, so the session starts with no SDKs and no
message. The header claims were checked against the Claude Code docs and are right (Bash script, runs as root on
Ubuntu 24.04 before Claude Code starts, must exit 0, result cached as a filesystem snapshot) but two documented
constraints are missing: the snapshot is only cached if the script finishes in "roughly five minutes", and the default
"Trusted" network level allows PyPI and `*.amazonaws.com` but lists no IBM Quantum host, so `run_ibm.py` against real
devices needs the IBM hosts added to the environment's allowed domains. Replacement for lines 5-10:
```bash
PIP="python3 -m pip install --break-system-packages"
$PIP "numpy>=2.0" "scipy>=1.14" "matplotlib>=3.9" "markdown>=3.6" || true
ALL='qiskit>=2.5,<3 qiskit-aer>=0.17,<0.19 qiskit-ibm-runtime>=0.42,<0.53 pennylane>=0.45,<0.47 cirq-core>=1.7,<2 qutip>=5.2,<6 qutip-qip>=0.4,<0.5 pylatexenc qualtran'
if ! $PIP $ALL; then                       # 한꺼번에 안 풀리면 묶음별로 넣고, 실패한 묶음을 적어 둔다
  for grp in "qiskit>=2.5,<3 qiskit-aer>=0.17,<0.19 qiskit-ibm-runtime>=0.42,<0.53" "pennylane>=0.45,<0.47" \
             "cirq-core>=1.7,<2" "qutip>=5.2,<6 qutip-qip>=0.4,<0.5" "pylatexenc" "qualtran"; do
    $PIP $grp || echo "설치 실패: $grp" >> /root/qc-book-setup-failures.txt
  done
fi
exit 0
# 주의: 대략 5분 안에 끝나야 환경 캐시(스냅숏)가 만들어진다. IBM 장치에 제출하려면 환경의 허용 도메인에 IBM Quantum 호스트를 더한다.
```
(`$ALL`/`$grp` are intentionally unquoted so they split into arguments; no globbing characters are present.)

### B10. Version ranges: open upper bounds on fast-moving packages, stale "current" notes, all-or-nothing hardware install
Checked on PyPI today:

| requirement | latest on PyPI (date) | Requires-Python | note |
|---|---|---|---|
| `qiskit>=2.5,<3` | 2.6.2 (2026-11-05) | >=3.10 | matches sdk.txt:2 |
| `qiskit-aer>=0.17` | **0.18.2 (2026-10-26)** | **>=3.11** | sdk.txt:3 and the snippets say 0.17.x; 0.18.0 only dropped Py 3.9/3.10 wheels and i686; no API removals listed |
| `qiskit-ibm-runtime>=0.42` | **0.52.0 (2026-11-09)** | >=3.10 | see B3; ≥0.48 needs qiskit ≥2.3 |
| `pennylane>=0.45,<0.47` | 0.45.1 (2026-06-26) | >=3.11 | matches |
| `cirq-core>=1.7,<2` | 1.7.0 (2026-06-29) | >=3.11 | matches |
| `qutip>=5.2,<6` | **5.3.1** (date not shown) | not seen | sdk.txt:3 says 5.2.x |
| `qutip-qip>=0.4` | 0.4.2 (2026-06-23) | not seen | 0.5.0.dev docs: "Passing Gate as a string input has been deprecated"; `index`/`arg_value`/`control_value` deprecated |
| `qualtran` | 0.7.0 (2026-01-22) | badge: 3.10+ (main: 3.12+) | README: "experimental preview release" with no backwards-compatibility guarantees |
| `pylatexenc` | 2.11 (2026-07-25), universal wheel | — | ok |
| `cudaq>=0.16` | 0.16.0.post1 (2026-09-17) | >=3.11, pip >= 24 | gpu.txt:1-2 matches (Linux x86_64/arm64, macOS arm64 CPU-only) |
| `amazon-braket-sdk` | 1.126.0 (2026-08-24) | >=3.11 (classifiers ≤3.13) | exists |
| `qiskit-braket-provider` | 0.27.3 (2026-10-29) | >=3.11 (classifiers ≤3.13) | exists; badge "Qiskit >=2.2.0" |
| `pennylane-qiskit` | 0.45.0 (2026-05-12) | >=3.11 | exists; **Qiskit 2.x compatibility not stated anywhere I could read — unverified** |
| `qiskit-ionq` | 1.1.1 (2026-06-05) | (0.5.6 page: >=3.9) | exists; IonQ docs: "versions 1.0 and above are compatible with Qiskit v2.0" |

All package names exist. `import cirq` with only `cirq-core`: PyPI describes cirq-core as "all the code you need" and
`cirq` as "Cirq with all the optional modules, not just cirq-core"; that `cirq-core` ships the `cirq` import package is
from my knowledge (the page does not literally say it).
Recommended lines:
```
# requirements/sdk.txt
qiskit-aer>=0.17,<0.19
qiskit-ibm-runtime>=0.42,<0.53
qutip-qip>=0.4,<0.5
qualtran>=0.7,<0.8
# requirements/hardware.txt
amazon-braket-sdk>=1.126,<2
qiskit-braket-provider>=0.27,<0.28
pennylane-qiskit>=0.45,<0.47
qiskit-ionq>=1.1,<2
```
`setup.sh:33` installs hardware.txt as one resolver run under `set -e`; if pennylane-qiskit (or any one extra) conflicts
with `qiskit>=2.5`, the script exits before writing the lock file or running the smoke test. Only
`amazon-braket-sdk` is used by a kit script — consider moving the other three behind their own flag (C8).
Not verified: qualtran's dependency pins (it pulls a large tree incl. PennyLane; whether it forces numpy below the
2.5.3 the numpy scripts were run with could not be checked — compare `requirements/lock-*.txt` after the first install).

### B11. `.devcontainer/devcontainer.json:3` — tag `1-3.12` is not among the documented tags
The image's Docker Hub page documents `3`, `3.10`…`3.14`, `3.12-bookworm`, `3.12-trixie`, `3-3.14`, `3.2-3.14`, …
(image major version is now 3). `1-3.12` is a frozen v1 line if it still resolves at all (pullability unverified).
```json
  "image": "mcr.microsoft.com/devcontainers/python:3.12-bookworm",
```
The image runs as the non-root `vscode` user; pip then falls back to a user install, which works with
`setup.sh --no-venv`.

### B12. `env/run_snippets.py:35, 50-51` — classification hides what the task in CLAUDE.md:51-52 needs
* line 35 `PYTHONWARNINGS="ignore"` suppresses the DeprecationWarnings that distinguish "API changed between 2.5 and
  2.6 / runtime 0.50" from "snippet bug". Use `PYTHONWARNINGS="default"` and record them.
* line 51: any `NameError` on the last stderr line becomes "조각(NameError)" — a typo in a self-contained snippet is
  indistinguishable from a fragment; a `ModuleNotFoundError` for a *secondary* module (`qutip_qip`, `qiskit_aer`,
  `qiskit_ibm_runtime` — `MODULE` only checks the primary) becomes "실패"; `ch18_02_cudaq.py` (`set_target("tensornet-mps")`,
  GPU only) becomes "실패" on a CPU machine.
  Of the snippets I read (ch14_01, ch07_02, ch07_03, ch15_02, ch12_02, ch17_02, ch18_02, ch14_02): seven are
  self-contained; ch14_01 is a true fragment (`CIRCUITS`, `hva_tfim`, `H_TFIM` undefined) and is correctly
  hardware-flagged; with `--hardware` it authenticates and calls `least_busy`, then stops at the NameError on its line 11
  before submitting anything.

Replacement for lines 35 and 45-54:
```python
    env = dict(os.environ, MPLBACKEND="Agg", PYTHONWARNINGS="default")
```
```python
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
                status = "NameError(조각이거나 오타 — 출처 확인)"
            else:
                status = "실패"
            if dep:
                last = (last + " ‖ 경고: " + "; ".join(dep))[:300]
        except subprocess.TimeoutExpired:
            status, last = "시간 초과", f"{a.timeout}초"
        rows.append({**e, "status": status, "detail": last[:300], "seconds": round(time.time() - t0, 1)})
```
plus `import re` and `ap.add_argument("--gpu", action="store_true", help="GPU 대상(nvidia, tensornet*)을 쓰는 스니펫도 돌린다")`.
(Not executed against real SDK output.)

### B13. `env/run_ibm.py:149` prints the instance CRN; CLAUDE.md:44 forbids CRNs in logs
"API 키와 CRN 은 … 에만 둔다. 키를 명령줄 인자·코드·로그·본문에 쓰지 않는다" vs `print("활성 인스턴스:", service.active_instance())`
(returns the CRN string per the API docs). Mask it:
```python
            crn = service.active_instance()
            print("활성 인스턴스:", (crn[:24] + "…" + crn[-6:]) if crn and len(crn) > 34 else crn)
```

### B14. `.gitignore` — secret variants and a multi-GB artifact are tracked
`git check-ignore` on a scratch copy: `.env` ignored; **`.env.local`, `.env.bak` tracked**;
**`qc-book-image.tar.gz` tracked** (Dockerfile:2 tells the user to create exactly that file in the repo root);
`ch14/ch14_record_fake.json` tracked (regenerated by every `test_all.sh` run with qiskit present, see B7).
`runs/*.json`, `ch14/ch14_record_real.json`, `requirements/lock-*.txt` tracked — correct, they are the evidence.
```
.env
.env.*
!.env.example
.venv/
__pycache__/
*.pyc
.ipynb_checkpoints/
snippets/
runs/regress/
runs/pending_*.json
ch14/ch14_record_fake.json
*.tar.gz
*.b64
```

### B15. `env/setup.sh:22, 39` — re-run trap and always-green post-create (both fixed by the A1/A2 replacement)
`[ -d .venv ] || "$PY" -m venv .venv`: on Debian/Ubuntu without `python3.x-venv`, the first run fails after creating a
partial `.venv`; the second run skips creation and dies at `source .venv/bin/activate` or `python -m pip`.
Line 39 runs the smoke test without `--strict`, so `postCreateCommand` (devcontainer.json:4) reports success with
missing SDKs.

### B16. `env/run_ibm.py:33-40` — `.env` parser: `export` and inline comments
Measured: `export QISKIT_IBM_INSTANCE=crn:…` is stored under the key `"export QISKIT_IBM_INSTANCE"` (instance silently
ignored → auto-selection of another instance); `AWS_DEFAULT_REGION=us-east-1 # region` becomes `'us-east-1 # region'`.
Replacement for lines 33-40 (tested):
```python
    for line in p.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        if line.startswith("export "):
            line = line[len("export "):].lstrip()
        k, v = line.split("=", 1)
        k, v = k.strip(), v.strip()
        if len(v) >= 2 and v[0] == v[-1] and v[0] in "\"'":
            v = v[1:-1]                      # 따옴표로 싼 값은 그대로(안의 # 도 값)
        else:
            v = v.split(" #", 1)[0].strip()  # 따옴표 없는 값 뒤의 주석은 버린다
        if v and not os.environ.get(k):
            os.environ[k] = v
```

---

## Grade C — polish

- **C1. `requirements/base.txt:1`** says "책을 쓴 환경의 버전을 하한으로 둔다" but the book environment is numpy 2.5.3 /
  scipy 1.18.1 / matplotlib 3.11.2 (ch14_hardware.md:207, this sandbox) while the bounds are 2.0 / 1.14 / 3.9. Reword to
  "하한은 넉넉히 잡았다; 책을 쓴 환경은 numpy 2.5.3 · scipy 1.18.1 · matplotlib 3.11.2".
- **C2. Stale "current version" notes**: `requirements/sdk.txt:2-3` (add: qiskit-aer 0.18.2, qiskit-ibm-runtime 0.52.0,
  qutip 5.3.1, qutip-qip 0.4.2, qualtran 0.7.0 as of 2026-11-10); `ch14/ch14_submit.py:2` ("최신 0.50.x").
- **C3. `env/smoke_test.py:165`** `version(mod.replace("_", "-"))` gives "?" for `cirq` (dist `cirq-core`) and `braket`
  (dist `amazon-braket-sdk`) in error rows; add a `DIST = {"cirq": "cirq-core", "braket": "amazon-braket-sdk"}` lookup.
  `ljust` misaligns columns with Hangul (double-width) — cosmetic.
- **C4. `env/run_ibm.py`**: line 102 shot-noise hint — for a two-outcome 50/50 ideal the *mean* TVD is
  √(2/π)·0.5/√N ≈ 0.40/√N (0.5/√N is the standard deviation); line 171 `job.backend()` may be `None`; line 52 a missing
  account gives a raw traceback (catch and print "먼저 --save 또는 .env"); `--n 0` builds an invalid circuit
  (`ap.add_argument("--n", type=int, default=3)` → validate `n >= 2`); line 106 uses UTC while the file name (117) uses
  local time; line 58 accepts an empty token.
- **C5. `env/run_braket.py`**: `PER_SHOT` contains `Aquila` (QuEra analog device — a gate `Circuit` cannot run there;
  drop it so the lookup refuses); `--shots 0` or negative passes the guard (`cost` ≤ $0.30) — add
  `if a.shots < 1: sys.exit("--shots 는 1 이상")`; `versions` lacks the SDK version
  (`importlib.metadata.version("amazon-braket-sdk")`); device status/availability window is not shown before `--yes`.
  The price table itself matches aws.amazon.com/braket/pricing today (Cepheus 0.000425, Garnet 0.00145, Emerald 0.0016,
  Aquila 0.01, IBEX-Q1 0.0235, Forte 0.08, $0.30/task; free tier is simulators only).
- **C6. `env/extract_snippets.py`**: stale files from an earlier extraction are never removed (add
  `shutil.rmtree(OUT, ignore_errors=True)` at the top of `main()`); `HARDWARE` (line 23) misses other remote paths
  (`AwsQuantumTask(`, `qp.device("qiskit.remote"`/`"braket.aws.qubit"`, `cudaq.set_target("ionq"|"braket"|…)`) — none
  occur in today's 37 snippets. Extraction itself is sound: 56 python blocks, 19 numpy/stdlib blocks skipped, 37 SDK
  blocks = the 37 "미실행" markers; no indented, attribute-carrying or unbalanced fences in any chapter.
- **C7. `env/setup.sh`**: the lock file name has only the date (a `--hardware` run overwrites the plain run of the same
  day; in `--no-venv` mode it lists every system package); `-h/--help` prints "모르는 옵션"; no upper Python bound
  (falls through to a `python3` that may be 3.14/3.15 — amazon-braket-sdk and qiskit-braket-provider list classifiers
  only up to 3.13).
- **C8. `requirements/hardware.txt:4-6`, `.env.example:7-8`**: `qiskit-braket-provider`, `pennylane-qiskit`,
  `qiskit-ionq` and `IONQ_API_KEY` are not used by any kit script (`run_braket.py` needs only `amazon-braket-sdk`).
  They widen the resolver's conflict surface (B10). `IONQ_API_KEY` is the right variable name per IonQ's docs.
- **C9. Docker**: no `.dockerignore` (the build context ships `.venv/`, `.env`, `snippets/`, `runs/` to the daemon; the
  image itself only `COPY`s `requirements/`); the container runs as root, so files written to the bind mount
  (`runs/*.json`, `snippets/`) are root-owned on Linux hosts. Suggested `.dockerignore`: `.venv` `.env*` `snippets`
  `runs` `**/__pycache__` `*.tar.gz` `.git`.
- **C10. README.md:40** assumes a GitHub-hosted repository; `/home/claude/qc-book` is not a git repository yet
  (`git rev-parse` fails), so `.gitignore` has never been exercised in place.
- **C11. `snippets/report.md`** lives in an ignored directory, yet it is the evidence for flipping "미실행 → 실행 완료"
  (CLAUDE.md:36-37). Write the report to `runs/snippets-YYYYMMDD.md` or un-ignore `snippets/report.*`.
- **C12. `env/claude_cloud_setup.sh`** duplicates the ranges of `requirements/sdk.txt` by hand (already different from
  what B10 recommends); add a comment in sdk.txt pointing at it, or generate it.

---

## Statements in CLAUDE.md / README.md that contradict the scripts or the repo

| where | statement | reality |
|---|---|---|
| README.md:20, CLAUDE.md:23 | test_all.sh = "numpy 스크립트(SDK 불필요)" | runs `ch14_submit.py` (qiskit, Aer, runtime fake provider) whenever qiskit imports (B7) |
| CLAUDE.md:22, smoke_test.py:7-8 | "비트 순서 확인" | printed, never asserted (B5) |
| README.md:15, setup.sh:5 | `--no-venv` for Codespaces/Docker/Claude Code cloud | fails on any PEP 668 Python and targets the wrong interpreter in the cloud image (A1, A2) |
| CLAUDE.md:44 | CRN never in logs | `run_ibm.py --list` prints it (B13) |
| .env.example:5-6 | region for Braket is set here | `run_braket.py` does not read `.env` (B6) |
| README.md:41, Dockerfile:1 | `--env-file .env` | file absent on a fresh clone (B8) |
| CLAUDE.md:38 | measured values: device, date, job ID, shots "`runs/` 의 JSON 에 다 들어 있다" | `ch14_submit.py` writes to `ch14/ch14_record_<mode>.json`, and loses it in `real` mode (A4); Braket `--fetch` records device "?" (B2) |
| run_ibm.py:10 | `--fetch <job id>` alone | needs the original `--circuit/--n` to be right (B1) |
| requirements/base.txt:1 | lower bounds = book environment | they are not (C1) |

Consistent and confirmed: 37 snippets with the stated per-SDK split (CLAUDE.md:51); "IBM Open 플랜은 28일에 QPU 10분"
(IBM plans page: "up to 10 minutes per 28-day rolling window", us-east only); Braket "$0.30 + 샷당 요금";
`claude_cloud_setup.sh` header (root, exit 0, cached); timing "약 10~15분; ch18 이 6분" is consistent with the other
process's log (ch01-ch17 ≈ 256 s, ch18 ≈ 5 min).

---

## Verified correct against primary documentation (no change needed)

- **qiskit-ibm-runtime** (docs for 0.51): `channel="ibm_quantum_platform"` is the default/recommended value (`ibm_cloud`
  is the legacy alias, `local` exists); `save_account(token, instance, channel, overwrite, set_as_default, …)`;
  `least_busy(min_num_qubits=None, instance=None, filters=None, use_fractional_gates=False, **kwargs)` with
  `operational=True` in the docs' own example; `backends(..., **kwargs)` with the docs example
  `backends(simulator=False, operational=True)`; `backend(name)`; `job(job_id)`; `usage()` → dict for the active
  instance; `active_instance()` → CRN str; `RuntimeJobV2.usage(partial=False)` → float seconds; `status()` → one of
  INITIALIZING/QUEUED/RUNNING/CANCELLED/DONE/ERROR; `backend(timeout=None)` → Backend | None; `job_id()` → str.
  `SamplerV2(mode=backend).run([(isa,)], shots=…)` and `result[0].data.meas.get_counts()` appear in a current IBM guide.
- **Qiskit 2.6**: `from qiskit.transpiler import generate_preset_pass_manager` is the documented path;
  `optimization_level`, `backend`, `seed_transpiler` are valid; `SamplerPubResult.join_data(names=None)` exists
  ("bitstring concatenation").
- **qiskit-aer 0.18.x**: release notes list no API removals relative to the calls used.
- **PennyLane 0.45**: release notes — "The `qml` alias as in `import pennylane as qml` has been updated to `qp` in our
  source code and documentation"; no breaking change to `device("default.qubit", wires=)`, `qnode`, `probs`.
- **Cirq**: state-vector ordering is `numpy.kron` order with the first qubit most significant (flip example gives
  `[0 0 1 0]`) → "10"; `simulate(...).final_state_vector` is current.
- **qutip-qip 0.4.2**: `add_gate(gate, targets=None, controls=None, …)`, `gate` "string or Gate", targets/controls
  "int or list"; the gate table names Hadamard "SNOT".
- **CUDA-Q (latest docs, 0.16)**: `cudaq.make_kernel()` builder with `qalloc(2)`, `h`, `cx`, `mz(qubits)` (spec
  example identical to smoke_test.py:121); `cudaq.sample(kernel, *args, shots_count=1000, …)`;
  `SampleResult.items()`; `Target.name`; "`SampleResult` bitstring keys … place qubit 0 as the left-most character"
  → "10". The builder API is not marked deprecated.
- **Amazon Braket SDK**: `AwsDeviceType` (StrEnum, `QPU`/`SIMULATOR`) is defined in `braket.aws.aws_device`;
  `AwsDevice.get_devices(arns, names, types, statuses, provider_names, order_by="name", aws_session)`; properties
  `arn`, `provider_name`, `type`; `AwsDevice(arn)`; `run(task_specification, …, shots=None, …)`; `measurement_counts`
  is a `Counter` whose "Key is the measurements in a big endian binary string"; README example
  `Circuit().h(0).cnot(0, 1)` / `task.result().measurement_counts`.
- **Expected bit strings in smoke_test.py** (X on qubit 0 of 2): Qiskit "01"; PennyLane, Cirq, QuTiP tensor, CUDA-Q,
  Braket "10" — all correct. `REF_BELL`/`REF_X0` and `bell_ok` are correct; TVD is computed correctly and the
  Bell/GHZ ideal is invariant under bit reversal, so the Qiskit-order comparison in run_ibm.py is valid.
- **Braket cost guard** (run_braket.py:80-86): unknown device → exit; cost > `--max-usd` → exit before `--yes` is
  considered; substring match is case-insensitive (`Ibex-Q1` → 0.0235, `Forte-Enterprise-1` → 0.08, `Ankaa-3`/`sv1` →
  refused). 100 shots on IBEX = $2.65 and on Forte = $8.30 are both blocked by the default `--max-usd 2`.

## Could not verify from a primary source (stated rather than guessed)

1. qiskit-ibm-runtime **0.52.0** release notes (released 2026-11-09; docs site still ends at 0.51.0).
2. Whether constructing the deprecated `SamplerV2` emits a warning, and whether legacy `SamplerV2(mode=AerSimulator())`
   local mode is unchanged in ≥0.50 (the guide now shows only the client-side `Sampler`).
3. `measure_all()` inserting a barrier (A3) — from the implementation as I know it; API pages are silent.
4. `Statevector.probabilities_dict`, `CircuitInstruction.operation.num_qubits`, `count_ops`, `AerSimulator().run(qc,
   shots=, seed_simulator=)`, QuTiP `basis`/`tensor`, `qp.probs` ordering, `LocalSimulator().run(…).result()` — standard
   APIs, not re-fetched today.
5. `from braket.aws import AwsQuantumTask`, `AwsQuantumTask(arn=…)`, `.id`, `.state()`, `.metadata()`, `result()`
   returning `None` on failure; `AwsDevice.name` (inherited from `Device`) — the reference page for
   `aws_quantum_task` was not fetchable.
6. `import cirq` working with `cirq-core` alone (strongly implied, not literally stated on PyPI).
7. qualtran's dependency pins and pennylane-qiskit's Qiskit 2.x support (B10).
8. Pullability of `mcr.microsoft.com/devcontainers/python:1-3.12`; Docker's exact error for a missing `--env-file`;
   boto3 `NoRegionError` — tool behaviour from knowledge, not run.
9. Python 3.14 support of qutip, qutip-qip, qualtran, qiskit-ionq (aer, cirq-core, pennylane, cudaq, pennylane-qiskit
   list 3.14 classifiers; qiskit ships `Python :: 3 :: Only`).

## Sources

- https://quantum.cloud.ibm.com/docs/en/api/qiskit-ibm-runtime/qiskit-runtime-service.md
- https://quantum.cloud.ibm.com/docs/api/qiskit-ibm-runtime/release-notes
- https://quantum.cloud.ibm.com/docs/en/api/qiskit-ibm-runtime/release-notes.md
- https://quantum.cloud.ibm.com/docs/en/api/qiskit-ibm-runtime/runtime-job-v2.md
- https://quantum.cloud.ibm.com/docs/en/guides/migrate-from-sampler-to-executor.md
- https://quantum.cloud.ibm.com/docs/en/guides/local-testing-mode.md
- https://quantum.cloud.ibm.com/docs/guides/plans-overview
- https://quantum.cloud.ibm.com/docs/guides/measure-qubits
- https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.transpiler.generate_preset_pass_manager
- https://eu-de.quantum.cloud.ibm.com/docs/en/api/qiskit/qiskit.primitives.SamplerPubResult
- https://eu-de.quantum.cloud.ibm.com/docs/en/api/qiskit/1.0/qiskit.circuit.QuantumCircuit.md
- https://qiskit.github.io/qiskit-aer/release_notes.html
- https://docs.pennylane.ai/en/stable/development/release_notes.html
- https://quantumai.google/cirq/simulation
- https://quantumai.google/reference/python/cirq/final_state_vector
- https://qutip-qip.readthedocs.io/en/stable/apidoc/qutip_qip.circuit.html
- https://qutip-qip.readthedocs.io/en/stable/qip-basics.html
- https://qutip-qip.readthedocs.io/en/latest/_modules/qutip_qip/circuit/circuit.html
- https://nvidia.github.io/cuda-quantum/latest/api/languages/python_api.html
- https://nvidia.github.io/cuda-quantum/latest/specification/cudaq/dynamic_kernels.html
- https://nvidia.github.io/cuda-quantum/latest/using/examples/executing_kernels.html
- https://amazon-braket-sdk-python.readthedocs.io/en/latest/_modules/braket/aws/aws_device.html
- https://amazon-braket-sdk-python.readthedocs.io/en/latest/_apidoc/braket.tasks.gate_model_quantum_task_result.html
- https://amazon-braket-sdk-python.readthedocs.io/en/stable/_apidoc/braket.aws.html
- https://aws.amazon.com/braket/pricing/
- https://docs.ionq.com/sdks/qiskit.md
- https://code.claude.com/docs/en/cloud-environments.md
- https://hub.docker.com/r/microsoft/devcontainers-python
- https://github.com/quantumlib/Qualtran
- https://github.com/quantumlib/Qualtran/releases
- https://newreleases.io/project/github/PennyLaneAI/pennylane-qiskit/release/v0.45.0
- https://pypi.org/project/qiskit/
- https://pypi.org/project/qiskit-aer/
- https://pypi.org/project/qiskit-ibm-runtime
- https://pypi.org/project/pennylane/
- https://pypi.org/project/cirq-core
- https://pypi.org/project/qutip/4.7.6/
- https://pypi.org/project/qutip-qip/
- https://pypi.org/project/qualtran/
- https://pypi.org/project/pylatexenc/
- https://pypi.org/project/cudaq/
- https://pypi.org/project/amazon-braket-sdk/
- https://pypi.org/project/qiskit-braket-provider/0.26.1/
- https://pypi.org/project/PennyLane-qiskit/
- https://pypi.org/project/qiskit-ionq/0.5.6
