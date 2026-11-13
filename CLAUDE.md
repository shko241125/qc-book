# 양자컴퓨팅 입문·심화서 — Claude Code 작업 규칙

한국어로 쓰는 양자컴퓨팅 책의 원고 저장소다. 장마다 본문(md), numpy 검증 스크립트(py), 그림(png), 배포용 html 이 한 폴더에 있다.
작업 지시는 한국어로 오고, 결과물(본문·보고)은 한국어로 쓴다. 코드, 명령, 식별자, 오류 문구는 원문 그대로 둔다.

## 구조

- `chNN/chNN_<slug>.md` — 본문. **이 파일이 원본이다.** html 과 Claude Docs 사본은 여기서 만든다.
- `chNN/chNN_<slug>.py` — numpy·scipy 만 쓰는 검증 스크립트. 본문의 모든 수치는 이 스크립트의 출력이다.
- `chNN/chNN_results.png` — 스크립트가 만드는 그림. `chNN/chNN_<slug>.html` — md 에서 만든 단일 파일 html.
- `ch00/` — 책 앞의 「기호와 비용 단위」. 기호가 겹치는 곳과 비용을 세는 규약이 여기 있다.
- `tools/md_to_html.py` — md → html. `tools/gate_cost.py` — 게이트 리스트에서 2큐비트 게이트 수·깊이·T-count 를 세는 도구.
- `env/` — 실행 환경 키트(아래 "명령"). `requirements/` — 버전 범위와 설치 기록(`lock-날짜.txt`). `runs/` — 장치 실행 결과(JSON), 실행 기록의 고정 사본(`*_ch14_record.json`), 실행 회로·보정값(`*_isa_calib.json`), 스니펫 실행 보고(`snippets-날짜.md`). 셋 다 커밋한다. `runs/` 의 기록은 고치지 않는다.
- `ch14/ch14_submit.py` — 14장 회로 다섯 개를 제출하는 스크립트. `fake`(가짜 Heron 장치의 보정 데이터로 로컬 실행), `real`(계획만), `real --yes`(IBM 장치 제출). 결과는 `ch14/ch14_record_<모드>.json`(다음 실행이 덮어쓴다; real 은 `runs/` 에 고정 사본도 남긴다).
- `appA/` — 부록 A 「첫 QPU 실행 기록」. `appA_qpu_runs.py` 는 `runs/` 의 고정 기록만 읽는 numpy 스크립트(당일 보정값 재현, BV 결함 열거). `env/fetch_ibm_job.py` — 끝난 작업의 ISA 회로와 그 시각 보정값을 받아 `runs/` 에 남긴다(QPU 사용 없음).
- `reviews/` — 독립 검토 보고서. `snippets/` — md 에서 뽑은 SDK 스니펫(생성물, 커밋하지 않는다).

## 명령

```bash
bash env/setup.sh                 # .venv 를 만들고 SDK 설치, 버전 기록, 점검까지. --hardware / --gpu / --no-venv
source .venv/bin/activate
python env/smoke_test.py          # SDK 별 Bell 상태 점검과 비트 순서 확인(기대와 다르면 "값·순서 불일치")
bash env/test_all.sh              # 모든 장·부록의 numpy 스크립트 재실행(약 10~15분; ch18 이 6~7분). 장 지정: bash env/test_all.sh ch10 appA. SDK 스크립트까지: --sdk(ch14_submit.py fake 가 돌아 ch14_record_fake.json 을 덮어쓴다)
python env/extract_snippets.py && python env/run_snippets.py     # 본문의 SDK 스니펫을 뽑아 돌리고 runs/snippets-날짜.md 에 표로 남긴다
python env/run_ibm.py --local     # 장치 실행 흐름을 Aer 로 확인. 실기는 --yes 를 붙일 때만 제출된다
cd ch14 && python ch14_submit.py fake   # 14장 회로를 가짜 장치로. real 은 계획만, real --yes 는 사용자가 승인했을 때만(QPU 약 10초)
python env/fetch_ibm_job.py <job id> ... --out runs/<이름>.json   # 끝난 작업의 ISA 회로·보정값(QPU 사용 없음)
cd chNN && python ../tools/md_to_html.py <slug> "<html 의 title>" chNN_results.png     # html 다시 만들기(그림이 없으면 - )
```

## 규약 (어기면 본문과 숫자가 어긋난다)

- **큐비트 순서는 빅 엔디안.** |q0 q1 …⟩ 에서 q0 가 왼쪽·상위 비트다. Qiskit 은 리틀 엔디안이므로 Qiskit 결과를 본문에 옮길 때 비트열을 뒤집는다(3장, `ch00`).
- **PennyLane 은 `import pennylane as qp`** 로 쓴다.
- **숫자는 스크립트에서만 나온다.** 본문에 수치를 새로 쓰거나 고치려면 먼저 스크립트가 그 값을 출력하게 하고, 출력과 본문을 대조한다. 기억이나 어림으로 수치를 쓰지 않는다.
- **게이트 비용**은 `tools/gate_cost.py` 의 규약(CX + 1큐비트 게이트, 전결합, ASAP 깊이)으로 세고, 게이트 리스트가 본문 회로와 같은지 수치로 확인한 뒤에만 적는다.
- **표시.** numpy 결과는 "numpy 구현(실행 완료)". SDK 코드 블록의 첫 줄 주석은 버전과 실행 여부를 적는다: `# qiskit==2.5.x — 작성됨, 미실행`.
  SDK 를 실제로 돌려 통과한 스니펫은 `# qiskit==<설치된 버전> — 실행 완료(YYYY-MM-DD)` 로 바꾸고, 출력이 본문 수치와 다르면 본문이 아니라 원인을 먼저 찾는다.
- **"예측"과 "실측".** 14장의 하드웨어 수치는 잡음 모형의 예측이다. 실제 장치에서 잰 값은 "실측"이라고 따로 표시하고, 장치 이름·날짜·작업 ID·샷 수를 함께 적는다(`runs/*.json` 과 `ch14/ch14_record_real.json` 에 들어 있다).
- md 본문에서 `~ _ [ * |` 는 `\~ \_ \[ \* \|` 로 쓴다. 표 칸 안에 맨 `|` 를 쓰면 표가 깨진다. 한글은 `\uXXXX` 로 쓰지 않는다.
- 장 구성은 10개 절(요약 / 직관 / 정의·수식 / 회로 구성 / 결과 / 다른 SDK / 한계와 반례 / 연습문제와 해답 / 학습 노트 / 출처와 환경)이다. 다른 SDK 병기는 장마다 하나, 고른 이유 한 문장.

## 하지 말 것

- **비밀을 커밋하지 않는다.** API 키와 CRN 은 `.env`(무시됨) 나 OS 환경 변수, `QiskitRuntimeService.save_account` 에만 둔다. 키와 CRN 을 명령줄 인자·코드·로그·본문에 쓰지 않는다(스크립트는 CRN 을 가려서 찍는다).
- **묻지 않고 장치에 제출하지 않는다.** `env/run_ibm.py`·`env/run_braket.py`·`ch14_submit.py real` 의 `--yes` 는 사용자가 그 실행을 승인했을 때만 쓴다. 유료 장치(Braket)는 예상 비용을 먼저 보여 준다. IBM Open 플랜은 28일에 QPU 10분이다.
- 통과시키려고 검증(assert, 대조)을 지우거나 허용 오차를 넓히지 않는다. 안 맞으면 원인을 보고한다.
- 본문을 고쳤으면 html 을 다시 만든다. 스크립트를 고쳤으면 그 장을 다시 돌려 본문 수치와 대조한다.

## 지금 남아 있는 일

- SDK 스니펫 37개가 전부 "작성됨, 미실행"이다(1~18장, Qiskit 16 · PennyLane 7 · Cirq 7 · QuTiP 5 · CUDA-Q 1 · Qualtran 1). `env/run_snippets.py` 로 돌려 통과·실패를 가리고 표시를 바꾼다.
  스니펫은 Qiskit 2.5 기준으로 썼고 지금 배포판은 2.6 이다. 실패하면 API 변경인지 스니펫 오류인지 구분해서 보고한다.
- (완료 2026-11-11) IBM 실기 경로: Bell(`ibm_kingston`) → `ch14_submit.py fake` → `real`(`ibm_fez`) → 14장 표의 "QPU 실측"·"사용 시간" 열 → 부록 A.
  Runtime 은 최상위 `Sampler`·`Estimator`(0.52 에서 클라이언트 쪽 executor 프리미티브)로 옮겼고 `requirements/sdk.txt` 는 `>=0.52,<0.54`.
- BV 의 오답 `0000`(5.2%)은 보조 큐비트(ibm_fez 물리 3번)가 |1⟩ 에서 시작한 것으로 설명되지만 크기가 보정값 상한의 2.65배라 확정하지 못했다(부록 A 5절).
  확정 실험(보조만 다른 큐비트에 두고 재실행, 3번 큐비트의 준비-측정 오류 직접 측정)은 QPU 수 초 — 사용자 승인 후.
- Braket 경로(`env/run_braket.py`)는 아직 실행하지 않았다.
- 다음 원고는 20장 「사례 1: 선형계와 블록 인코딩」.
