# 양자컴퓨팅 입문·심화서 — 원고와 실행 환경

장마다 본문(`chNN/*.md`), numpy 검증 스크립트(`*.py`), 그림(`*.png`), html 이 들어 있다. 실행 환경은 `env/` 의 키트로 만든다.
자세한 안내는 「양자 실행 환경 구축 가이드」 문서에 있고, 여기에는 명령만 모았다.

## 5분 설치

```bash
# Python 3.11 이상(3.12 권장)이 있는 Linux · macOS · WSL2 에서, 저장소 루트에서
bash env/setup.sh                 # .venv 생성 → SDK 설치 → 버전 기록 → 점검
source .venv/bin/activate
python env/smoke_test.py          # 표에 "통과"가 줄줄이 나오면 끝
```

옵션: `--hardware`(AWS Braket SDK), `--gpu`(NVIDIA CUDA-Q), `--no-venv`(Codespaces·Docker 처럼 이미 격리된 곳 — 지금 쓰는 `python3` 에 설치한다).

## 책 내용 다시 돌리기

```bash
bash env/test_all.sh                                   # 모든 장·부록의 numpy 스크립트(SDK 불필요, 10~15분). SDK 스크립트까지는 --sdk
python env/extract_snippets.py && python env/run_snippets.py     # 본문의 SDK 스니펫 41개 → runs/snippets-날짜.md
```

## 실제 장치에서 돌리기 (IBM Quantum, 무료 Open 플랜)

```bash
python env/run_ibm.py --local                # 계정 없이 흐름 확인
python env/run_ibm.py --save                 # API 키와 인스턴스 CRN 저장(한 번)
python env/run_ibm.py --list                 # 쓸 수 있는 장치
python env/run_ibm.py --circuit bell         # 제출 계획만 본다
python env/run_ibm.py --circuit bell --yes   # 제출. 결과는 runs/ 에 JSON 으로 남는다
```

14장의 회로 다섯 개(GHZ, BV, Grover, QPE, HVA)는 `cd ch14 && python ch14_submit.py fake` 로 먼저 로컬에서 돌려 보고, `python ch14_submit.py real` 로 계획을 본 뒤 `real --yes` 로 제출한다.
끝난 작업의 실제 회로와 그 시각의 보정값은 `python env/fetch_ibm_job.py <job id> --out runs/<이름>.json` 으로 받는다(QPU 사용 없음). 첫 실행(2026-11-11)의 기록과 해석은 부록 A(`appA/`).

AWS Braket 은 `bash env/setup.sh --hardware` 뒤 `python env/run_braket.py --local` 부터. QPU 는 유료다(작업당 $0.30 + 샷당 요금).

IBM 경로(`run_ibm.py`, `ch14_submit.py`)는 2026-11-11 qiskit-ibm-runtime 0.52.0 으로 실기까지 실행했다. 그 뒤 `ch14_submit.py` 의 기록 순서와 `fetch_ibm_job.py` 의 추가 항목을 고쳤고(검토 R11), 이 변경분은 아직 돌리지 않았다 — `ch14_submit.py fake` 로 먼저 확인한다. Braket 경로는 아직 실행하지 않았다.

## 다른 곳에서 열기

- **GitHub Codespaces**: 저장소를 GitHub 에 올린 뒤 Code → Codespaces → Create. `.devcontainer` 가 설치까지 한다.
- **Docker**: `docker build -t qc-book . && docker run --rm -it -v "$PWD":/work $([ -f .env ] && echo --env-file .env) qc-book`
- **Claude Code**: 이 폴더에서 `claude` 를 실행하면 `CLAUDE.md` 의 규칙을 읽고 시작한다.

## 비밀 관리

API 키는 `.env`(`.env.example` 을 복사) 또는 `--save` 로만 둔다. `.env` 는 `.gitignore` 에 들어 있다. 키가 노출됐으면 대시보드에서 바로 지우고 새로 만든다.
