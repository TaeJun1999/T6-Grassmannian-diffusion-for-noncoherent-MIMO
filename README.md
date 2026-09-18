# T6 — Grassmannian diffusion prior for noncoherent MIMO: code & results archive

이 저장소는 **진행 중인 탐색 과제(T6)의 코드·수치 결과 보관소**다. 연구 노트·의사결정 기록·판정 근거는 Claude 프로젝트에 두고, 여기에는 재현에 필요한 스크립트와 결과 파일만 둔다.

> ⚠️ **공개 범위 주의.** 2026-09-18 확인 시 이 저장소는 **Public**이다. `records/`, `handoffs/`, `code/`는 아직 발표하지 않은 문제 설정·게이트 기준·베이스라인·상한 계산을 그대로 담고 있으므로, 논문 제출 전까지는 저장소를 **Private으로 바꾸는 것**을 권한다.
> Public을 유지해야 한다면 `records/`와 `handoffs/`를 빼고, 대신 `PREREGISTRATION_HASH.txt`만 올린다 — 내용을 드러내지 않고 "게이트 기준을 실험 전에 고정했다"는 시점만 커밋 타임스탬프로 증명한다.

## 1. 무엇이 어디에 있나

| 위치 | 내용 | 정본 |
|---|---|---|
| Claude 프로젝트 | 노트 10종(`T6_*.md`), 기록 4종, 논문 PDF 3, 세션에서 import하는 코드 일부, 요약 txt | **세션 중 작업본**(노트의 정본) |
| 이 저장소 | 모든 실행 코드, 결과 JSON·요약·그림, ray-tracing 채널 데이터, arXiv 원시 덤프, 세션 핸드오프, `records/`의 기록 4종 | **코드·데이터의 정본, 기록의 버전 이력** |

기록 4종은 양쪽에 둔다: 세션 중에는 프로젝트 파일을 고치고, 세션이 끝나면 여기에 커밋한다(규칙은 `records/README.md`). 커밋 이력이 곧 날짜별 diff이고, 프로젝트에서 `__1_` 접미사가 쌓이는 문제도 이걸로 해결된다.

노트에 인용된 수치는 전부 `results/**/*summary*.txt`에 있다. JSON은 그 요약을 다시 만들거나 곡선을 그릴 때만 필요하다.

## 2. 디렉터리

```
records/     기록 4종 (gate_evidence / decision_log / open_questions / scooping_log) + SHA256SUMS
code/        실행 스크립트 (아래 재현 절차 참조)
  ceiling/     게이트 (0) 상한: prior-free vs genie 비율·SER, 폐형 Bingham 적분 + 단위시험
  diff_ustm/   차동 USTM tall-block 확장, 2T-블록 참조
  energy/      energy-based 검출(MCG 2016)의 저랭크 분석
  q12_geom/    Q-12 (a) 기하(one-ring형) 모델: 이득 상관 vs 부분공간 지속성
  q12_rt/      Q-12 (b) Sionna RT: 궤적 생성·지표·요약
  competitor/  예측-CSI coherent 검출기의 alpha 의존 상한
  scooping/    월례 arXiv 확인 스크립트
results/     위 코드의 출력 (JSON) + 노트에 실린 요약 (txt/csv)
data/rt_channels/
  t6_q12_rt_H.npz        22개 RT 궤적의 채널 H_b (complex64, 키 = run 이름)
  raw/                   RT 실행 원본 JSON (지표 + H, 사람이 읽을 수 있는 형태)
figures/     노트에 실린 그림 3장
raw/arxiv/   2026-09-16 scooping 실행의 원시 응답
handoffs/    세션 핸드오프 (연구 상태 요약)
```

## 3. 환경

```
python 3.12, numpy, scipy, matplotlib
sionna-rt 2.1.0, mitsuba 3.9.1, drjit 1.5.0   # code/q12_rt 만 필요
pip install -r requirements.txt
```
CPU 1코어 기준: 상한 계산 ~5분/설정, RT 궤적 1개(80블록) 30–90초.

## 4. 재현 절차 (seed 고정)

```bash
# 게이트 (0) 상한 — 단위시험 먼저
python code/ceiling/test_t6_ceiling.py
python code/ceiling/t6_ceiling.py --configs "16,2,32,2;8,2,16,2" --seed 20260917 --out t6_ceiling_results.json
python code/ceiling/t6_summarize.py                     # -> t6_ceiling_summary.txt
# 장기 부분공간 genie 상한 (축소 문제; SNR + 10log10(N/r_eff) 로 환산)
python code/ceiling/t6_ceiling.py --configs "16,2,4,2;16,2,6,2;16,2,8,2;8,2,4,2" --snr="-4:14:2" \
       --n_ustm 8000 --n_fin 2000 --Ks 64,1024 --seed 20260917 --out t6_ceiling_LT.json

# 차동 USTM
python code/diff_ustm/t6_diff_ustm.py --selftest
python code/diff_ustm/t6_diff_ustm.py ; python code/diff_ustm/t6_diff_summarize.py

# energy-based (MCG 2016) 분석
python code/energy/t6_energy_det.py --selftest
python code/energy/t6_energy_det.py --out t6_energy_det.json

# Q-12 (a) 기하 모델
python code/q12_geom/t6_q12_geom.py --selftest
python code/q12_geom/t6_q12_geom.py --out t6_q12_geom.json
python code/q12_geom/t6_q12_longterm.py
python code/q12_geom/t6_q12_summarize.py                # -> t6_q12_summary.txt

# Q-12 (b) ray-tracing  (t6_q12_rt_cars.py / _nlos.py 는 t6_q12_rt.py 를 import)
python code/q12_rt/t6_q12_metrics.py                    # 단위시험
python code/q12_rt/t6_q12_rt.py                         # LOS 4점 × 산란계수 × 블록간격
python code/q12_rt/t6_q12_rt_cars.py                    # 자동차 클러터
python code/q12_rt/t6_q12_rt_nlos.py                    # 검증 NLOS 점
python code/q12_rt/t6_q12_rt_summarize.py               # -> t6_q12_rt_summary.txt

# 경쟁자 상한
python code/competitor/t6_alpha_competitor.py --n 3000

# 월례 scooping (companion §7 키워드)
python code/scooping/t6_arxiv_oai_check.py --since 2026-09-16
python code/scooping/t6_arxiv_listing_check.py
```

스크립트는 같은 디렉터리에서 실행하도록 쓰였다. 다른 폴더에서 돌릴 때는 `code/q12_rt`처럼 서로 import하는 경우(`t6_q12_rt_cars.py` → `t6_q12_rt.py`, `t6_q12_metrics.py`; `t6_alpha_competitor.py` → `t6_ceiling.py`)만 `PYTHONPATH`를 맞추면 된다.

## 5. RT 채널 데이터 쓰는 법

```python
import numpy as np
H = np.load('data/rt_channels/t6_q12_rt_H.npz')
print(list(H.keys())[:3])          # 예: 'rt_00_LOS_A_D10ms_S0.0_nocars'
Hb = H['nlos_00_NLOS_true_D10ms_S0.4_nocars']   # (B, M, N) = (80, 2, 32) complex64
```
지표만 필요하면 `results/q12_rt/t6_q12_rt_metrics.json`(40 kB)에 22개 run의 ρ_G, ρ_S, r_eff, PEF, LT 중첩과 설정이 전부 들어 있다. `data/rt_channels/raw/*.json`은 같은 내용의 원본(지표 + H 실수/허수 리스트)이다.

주의(노트에도 기록): 첫 스캔의 LOS 판정이 `paths.cir()`의 기본 `normalize_delays=True` 때문에 틀려, run 이름의 `NLOS_A`/`NLOS_B`는 실제로 **LOS 점**이다(P2, P3). 검증된 NLOS는 `NLOS_true`(P5)뿐이다. 요약 스크립트가 이 정정 라벨을 붙여 준다.

## 6. 현재 상태 (2026-09-18)

- 게이트 (0): 통과 확정. 장기 부분공간 genie 기준 간극 주 레짐 2.7 dB (r_eff=4) / 1.7–2.3 dB (r_eff=8).
- 게이트 (i)(ii)(iii): 근거 없음. 다음 실험은 Q-15 — `data/rt_channels`의 채널 위에서 검출기 6종 직접 비교.
- D-11(이득 탈상관 + 순간 부분공간 유지) 철회, D-15로 대체. D-14로 prior의 대상이 장기 수신 부분공간으로 바뀜.
- 자세한 근거·판정 규칙·미해결 질문은 Claude 프로젝트의 기록 4종에 있다.
