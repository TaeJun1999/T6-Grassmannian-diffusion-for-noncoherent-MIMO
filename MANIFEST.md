# 파일 목록과 출처

| 파일 | 만든 스크립트 | 노트 |
|---|---|---|
| results/ceiling/t6_ceiling_results.json, _summary.txt, _summary.csv | code/ceiling/t6_ceiling.py, t6_summarize.py | T6_ceiling_derivation.md |
| results/ceiling/t6_ceiling_LT.json | code/ceiling/t6_ceiling.py (축소 문제 (T,2,r_eff,2)) | T6_q12_geom_note.md §5 |
| results/diff_ustm/t6_diff_results.json, t6_twoblock_ref.json | code/diff_ustm/* | T6_diff_ustm_note.md |
| results/energy/t6_energy_det.json, _summary.txt | code/energy/t6_energy_det.py | T6_energy_based_note.md |
| results/q12_geom/t6_q12_geom.json, t6_q12_longterm.json, t6_q12_summary.txt | code/q12_geom/* | T6_q12_geom_note.md |
| results/q12_rt/t6_q12_rt_metrics.json, t6_q12_rt_summary.txt | code/q12_rt/* | T6_q12_rt_note.md |
| results/q12_rt/t6_q12_rt_scan2.json | code/q12_rt/t6_q12_rt_scan2.py | LOS/NLOS 점 탐색(정정판) |
| results/competitor/t6_alpha_competitor.json, _summary.txt | code/competitor/t6_alpha_competitor.py | T6_q12_rt_note.md §7 |
| data/rt_channels/t6_q12_rt_H.npz | code/q12_rt/* 결과에서 추출 (complex64) | Q-15 입력 |
| data/rt_channels/raw/*.json | code/q12_rt/* 원본 출력 | 위 npz와 동일 내용 |
| figures/t6_ceiling_plots.png | code/ceiling/t6_summarize.py | 게이트 (0) |
| figures/t6_diff_results_plots.png | code/diff_ustm/t6_diff_summarize.py | 차동 USTM |
| figures/t6_q12_plots.png | (세션 내 matplotlib, 코드 미보관) | Q-12 (a) 3패널 |
| raw/arxiv/*.json | code/scooping/* (2026-09-16 실행) | scooping_log_T6.md |
| handoffs/T6_session_handoff_2026-09-16b.md | 수기 | 세션 상태 요약 |
| records/*.md | 세션마다 갱신 | 기록 정본 사본, 규칙은 records/README.md |
| records/SHA256SUMS.txt, PREREGISTRATION_HASH.txt | sha256sum | 프로젝트 사본과의 동일성·기준 고정 시점 증명 |

프로젝트에만 두는 파일(여기 없음): 노트 10종 `T6_*.md`, 논문 PDF 3종(저작권), 세션에서 직접 import하는 코드 사본.
