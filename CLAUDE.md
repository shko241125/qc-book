# 양자컴퓨팅 입문·심화서 — Claude Code 작업 규칙

한국어로 쓰는 양자컴퓨팅 책의 원고 저장소다. 장마다 본문(md), numpy 검증 스크립트(py), 그림(png), 배포용 html 이 한 폴더에 있다.
작업 지시는 한국어로 오고, 결과물(본문·보고)은 한국어로 쓴다. 코드, 명령, 식별자, 오류 문구는 원문 그대로 둔다.

## 구조

- `chNN/chNN_<slug>.md` — 본문. **이 파일이 원본이다.** html 과 Claude Docs 사본은 여기서 만든다.
- `chNN/chNN_<slug>.py` — numpy·scipy 만 쓰는 검증 스크립트. 본문의 모든 수치는 이 스크립트의 출력이다.
- `chNN/chNN_results.png` — 스크립트가 만드는 그림. `chNN/chNN_<slug>.html` — md 에서 만든 단일 파일 html.
  스크립트가 그림 말고 남기는 파일도 생성물이다(고치지 말고 스크립트를 다시 돌린다): `ch14/ch14_pred.json`(14장 예측과 판독 오류를 끈 예측, 가정한 CZ 수 — 부록 A 스크립트가 읽는다), `ch20/ch20_case.json`·`ch21/ch21_case.json`(20·21장 SDK 스니펫용 사례), `appA/appA_posthoc.json`(예측·실측·사후 확인값 병기 표 — 14장 5절과 13장 5절이 이 출력에서 옮긴다), `appA/appA_bv_plan.json`(BV 확정 실험의 회로·판정 규칙), `chNN_stdout.txt`.
- `ch00/` — 책 앞의 「기호와 비용 단위」. 기호가 겹치는 곳과 비용을 세는 규약이 여기 있다.
- `tools/md_to_html.py` — md → html. `tools/gate_cost.py` — 게이트 리스트에서 2큐비트 게이트 수·깊이·T-count 를 세는 도구.
- `env/` — 실행 환경 키트(아래 "명령"). `requirements/` — 버전 범위와 설치 기록(`lock-날짜.txt`). `runs/` — 장치 실행 결과(JSON), 실행 기록의 고정 사본(`*_ch14_record.json`), 실행 회로·보정값(`*_isa_calib.json`), 스니펫 실행 보고(`snippets-날짜.md`). 셋 다 커밋한다. `runs/` 의 기록은 고치지 않는다.
- `ch14/ch14_submit.py` — 14장 회로 다섯 개를 제출하는 스크립트. `fake`(가짜 Heron 장치의 보정 데이터로 로컬 실행), `real`(계획만), `real --yes`(IBM 장치 제출). 결과는 `ch14/ch14_record_<모드>.json`(다음 실행이 덮어쓴다). real 은 제출 직후부터 `runs/` 의 고정 사본에 쓰고, 본문과 `ch14_hardware.py` 는 고정 사본만 읽는다.
- `appA/` — 부록 A 「첫 QPU 실행 기록」. `appA_qpu_runs.py` 는 `runs/` 의 고정 기록만 읽는 numpy 스크립트(당일 보정값 재현, BV 결함 열거). `env/fetch_ibm_job.py` — 끝난 작업의 ISA 회로와 그 시각 보정값을 받아 `runs/` 에 남긴다(QPU 사용 없음).
  `appA_bv_confirm.py`(numpy) — BV '0000' 확정 실험의 변형 회로·실행 전 예측·판정 규칙을 만들고, 기록 파일을 주면 같은 규칙으로 판정한다. `appA_bv_submit.py`(Qiskit) — 첫 작업(회로 둘)을 제출한다: 인자 없이 Aer 확인, `real` 계획만, `real --yes` 제출, `real --resume <기록>` 결과만 다시 받기.
- `reviews/날짜/RN_findings.md` — 독립 검토 보고서(R11: 부록 A·14장, R12: 20장, R13: 21장 첫 판, R14: 2026-11-12 결정 반영분, R15: 21장 고쳐 쓴 판, R16: R15 반영분). `snippets/` — md 에서 뽑은 SDK 스니펫(생성물, 커밋하지 않는다).

## 명령

```bash
bash env/setup.sh                 # .venv 를 만들고 SDK 설치, 버전 기록, 점검까지. --hardware / --gpu / --no-venv
source .venv/bin/activate
python env/smoke_test.py          # SDK 별 Bell 상태 점검과 비트 순서 확인(기대와 다르면 "값·순서 불일치")
bash env/test_all.sh              # 모든 장·부록의 numpy 스크립트 재실행(약 35분; ch21 이 약 20분, ch18 이 6~7분). 장 지정: bash env/test_all.sh ch10 appA. SDK 스크립트까지: --sdk(ch14_submit.py fake 가 돌아 ch14_record_fake.json 을 덮어쓴다)
python env/extract_snippets.py && python env/run_snippets.py     # 본문의 SDK 스니펫을 뽑아 돌리고 runs/snippets-날짜.md 에 표로 남긴다
python env/run_ibm.py --local     # 장치 실행 흐름을 Aer 로 확인. 실기는 --yes 를 붙일 때만 제출된다
cd ch14 && python ch14_submit.py fake   # 14장 회로를 가짜 장치로. real 은 계획만, real --yes 는 사용자가 승인했을 때만(QPU 약 10초)
python env/fetch_ibm_job.py <job id> ... --out runs/<이름>.json   # 끝난 작업의 ISA 회로·보정값(QPU 사용 없음)
python appA/appA_bv_confirm.py && python appA/appA_bv_submit.py   # BV 확정 실험: 설계·예측(numpy) → Aer 로 회로 확인. real 은 계획만, real --yes 는 사용자가 승인했을 때만(QPU 약 4초)
cd chNN && python ../tools/md_to_html.py <slug> "<html 의 title>" chNN_results.png     # html 다시 만들기(그림이 없으면 - )
```

## 규약 (어기면 본문과 숫자가 어긋난다)

- **큐비트 순서는 빅 엔디안.** |q0 q1 …⟩ 에서 q0 가 왼쪽·상위 비트다. Qiskit 은 리틀 엔디안이므로 Qiskit 결과를 본문에 옮길 때 비트열을 뒤집는다(3장, `ch00`).
- **PennyLane 은 `import pennylane as qp`** 로 쓴다.
- **사전 등록의 판정은 실행 전에 정한 방식으로 적는다.** 결과나 검토를 본 뒤에 더한 분석은 "사후"로 나란히 적고 판정을 바꾸는 데 쓰지 않는다(21장 5절 (6)).
- **숫자는 스크립트에서만 나온다.** 본문에 수치를 새로 쓰거나 고치려면 먼저 스크립트가 그 값을 출력하게 하고, 출력과 본문을 대조한다. 기억이나 어림으로 수치를 쓰지 않는다.
- **게이트 비용**은 `tools/gate_cost.py` 의 규약(CX + 1큐비트 게이트, 전결합, ASAP 깊이)으로 세고, 게이트 리스트가 본문 회로와 같은지 수치로 확인한 뒤에만 적는다.
- **표시.** numpy 결과는 "numpy 구현(실행 완료)". SDK 코드 블록의 첫 줄 주석은 버전과 실행 여부를 적는다: `# qiskit==2.5.x — 작성됨, 미실행`.
  SDK 를 실제로 돌려 통과한 스니펫은 `# qiskit==<설치된 버전> — 실행 완료(YYYY-MM-DD)` 로 바꾸고, 출력이 본문 수치와 다르면 본문이 아니라 원인을 먼저 찾는다.
- **예측은 고치지 않고 병기한다(2026-11-12 결정).** 실측을 본 뒤에 한 계산은 "사후"라고 이름 붙여 예측·실측 옆에 놓는다(14장 5절, 부록 A 4.1절). 사후 값을 "예측"이라 부르지 않는다.
- **"예측"·"실측"·"표본값".** numpy 로 뽑은 표본의 값은 "표본값"이나 "계산값"이라 쓴다("실측"은 QPU 값에만). 14장의 하드웨어 수치는 잡음 모형의 예측이다. 실제 장치에서 잰 값은 "실측"이라고 따로 표시하고, 장치 이름·날짜·작업 ID·샷 수를 함께 적는다(`runs/*.json` 과 `ch14/ch14_record_real.json` 에 들어 있다).
- md 본문에서 `~ _ [ *` 는 `\~ \_ \[ \*` 로 쓴다. `|` 는 표 칸 안에서 `\|` 로 쓰고(맨 `|` 를 쓰면 표가 깨진다), 문단을 `|` 로 시작하지 않는다(`|x| ≥ …` 로 시작하면 표로 읽힐 수 있다). 한글은 `\uXXXX` 로 쓰지 않는다.
- 장 구성은 10개 절(요약 / 직관 / 정의·수식 / 회로 구성 / 결과 / 다른 SDK / 한계와 반례 / 연습문제와 해답 / 학습 노트 / 출처와 환경)이다. 다른 SDK 병기는 장마다 하나, 고른 이유 한 문장.

## 하지 말 것

- **비밀을 커밋하지 않는다.** API 키와 CRN 은 `.env`(무시됨) 나 OS 환경 변수, `QiskitRuntimeService.save_account` 에만 둔다. 키와 CRN 을 명령줄 인자·코드·로그·본문에 쓰지 않는다(스크립트는 CRN 을 가려서 찍는다).
- **묻지 않고 장치에 제출하지 않는다.** `env/run_ibm.py`·`env/run_braket.py`·`ch14_submit.py real`·`appA_bv_submit.py real` 의 `--yes` 는 사용자가 그 실행을 승인했을 때만 쓴다. 유료 장치(Braket)는 예상 비용을 먼저 보여 준다. IBM Open 플랜은 28일에 QPU 10분이다.
- 통과시키려고 검증(assert, 대조)을 지우거나 허용 오차를 넓히지 않는다. 안 맞으면 원인을 보고한다.
- 본문을 고쳤으면 html 을 다시 만든다. 스크립트를 고쳤으면 그 장을 다시 돌려 본문 수치와 대조한다.

## 지금 남아 있는 일

- SDK 스니펫 41개가 전부 "작성됨, 미실행"이다(1~18장 37개 + 20장의 Qiskit·PennyLane 2개 + 21장의 PennyLane·CUDA-Q 2개). `env/run_snippets.py` 로 돌려 통과·실패를 가리고 표시를 바꾼다.
  스니펫은 Qiskit 2.5 기준으로 썼고 지금 배포판은 2.6 이다. 실패하면 API 변경인지 스니펫 오류인지 구분해서 보고한다.
  20·21장 스니펫은 `ch20/ch20_case.json`·`ch21/ch21_case.json` 을 읽으므로 그 장 폴더 안에서 돌린다. 기대 출력은 스니펫의 주석(20장: 성공 확률 0.1249, 해의 오차 약 6.0e-04. 21장: ⟨C⟩ 0.2904, p* 0.0313).
- (완료 2026-11-11) IBM 실기 경로: Bell(`ibm_kingston`) → `ch14_submit.py fake` → `real`(`ibm_fez`) → 14장 표의 "QPU 실측"·"사용 시간" 열 → 부록 A.
  Runtime 은 최상위 `Sampler`·`Estimator`(설치된 0.52.0 에서 클라이언트 쪽 executor 프리미티브)로 옮겼고 `requirements/sdk.txt` 는 `>=0.52,<0.54`.
- **실행해 보지 않은 것(SDK 없는 환경에서 썼다).** `ch14/ch14_submit.py` 의 기록 순서(제출 직후부터 `runs/` 고정 사본, `primitives`·`bit_order_note` 항목),
  `env/fetch_ibm_job.py` 의 추가 항목(항목별 측정 시각, `init_error`, 게이트 길이, Sampler counts, 끝나지 않은 작업 건너뛰기, 이미 있는 파일 거절), `appA/appA_bv_submit.py` 전체.
  먼저 `cd ch14 && python ch14_submit.py fake`, `python env/fetch_ibm_job.py dbq0qmoqtr0c738f53u0 dbq0qpk4148s739vbffg --out runs/<새 이름>.json`(QPU 사용 없음), `python appA/appA_bv_submit.py`(Aer, 장치에 닿지 않는다)로 확인한다.
- **BV 확정 실험(2026-11-12 결정: 새 설계로 진행. 제출은 사용자가 그 실행을 승인한 뒤).** BV 의 오답 `0000`(5.15%)은 "첫 CZ 때 보조 큐비트(ibm_fez 물리 3번)가 |1⟩ 이었다"로 설명되고, 왜 |1⟩ 이었는지는 가르지 못했다.
  한 작업의 회로들이 어떤 순서로 도는지(같은 회로가 연달아인지, 샷마다 번갈아인지 — IBM 문서는 뒤쪽) 기록에 없으므로, 첫 작업은 순서에 기대지 않는 회로 둘(P0, P1)만 낸다(검토 R14).
  순서: ① `python env/fetch_ibm_job.py dbq0qmoqtr0c738f53u0 --out runs/<새 이름>.json --bitstrings`(`init_error` 와 샷별 비트열, QPU 불필요) → ② `python appA/appA_bv_submit.py`(Aer) → `real`(계획) → ③ 승인 후 `real --yes`(PUB 2개 × 4000샷, QPU 약 4초)
  → ④ `python appA/appA_bv_confirm.py runs/<기록>.json` 의 판정(S·X·I·T·약함·혼합)을 부록 A 5절에 "실측"으로 적는다. 경우별 예측과 판정 규칙은 실행 전에 고정돼 있다(부록 A 5절, `appA_bv_stdout.txt`) — 결과를 본 뒤 규칙을 바꾸지 않는다.
  둘째 작업(P2·P3·P4·PM·Qm 가운데)은 첫 작업의 결과를 본 뒤에 고르고, 읽는 법을 실행 전에 다시 고정한다.
- (결정 2026-11-12) 14장 5절의 예측값은 고치지 않는다. 게이트 중 열완화의 이중 계산을 뺀 값은 "사후 ②"로 병기했다. 다음 실기 실행에서는 두 모형(넣은 것·뺀 것)의 예측을 실행 전에 함께 적는다.
- Braket 경로(`env/run_braket.py`)는 아직 실행하지 않았다.
- (초안 완료 2026-11-11) 20장 「사례 1: 선형계와 블록 인코딩」 — 사전 등록 기준 둘 다 통과.
- (초안 완료 2026-11-12) 21장 「사례 2: 단백질 접힘 HP 격자」 — 판정이 훈련 방식에 달려 있다. 사전 등록의 판정은 설정을 검토 전에 정한 방식(잔기 12개에서 샷 훈련: 시드 20개의 중앙값 169,302, 담금질 15,601)으로 불통과. 작은 크기의 각도를 그대로 쓰는 방식(17,114 / 5,963)은 사후 분석으로 병기했다. 첫 판은 QAOA 출발점의 부호가 틀렸고 검토 R13 이 잡았다 — 변분 최적화의 값은 출발점 여럿·격자 대조 없이 "최솟값"이라 부르지 않는다. 둘째 판은 탐색을 80회에서 끊은 값과 시드 5개의 중앙값(90,328)을 적었고 검토 R15 가 잡았다 — 확률적 최적화는 범위를 함께 적고, 탐색을 횟수로 끊었으면 마지막 구간의 개선폭을 찍는다. 스크립트 끝의 "회귀 검사: 모두 통과" 줄이 `ch21_stdout.txt` 에 있어야 통과한 실행이다.
  `CH21_QUICK=1` 은 평가 횟수·시드를 줄인 개발용 점검이다(회귀 검사를 건너뛴다. 본문 수치의 출처로 쓰지 않는다). 다음 원고는 22장 「사례 3: 화자 임베딩 양자 커널」(19장 5절의 사전 등록 표). 등록 행에 "허용하는 절차와 그 비용을 세는 법"이 빠져 있지 않은지 실행 전에 먼저 본다.
- 17~19장의 Claude Docs 사본은 2026-11-11 의 용어 정리("실측" → "표본값", "6′ 같은 예산 비교") 전 판이다. 원본은 md.
