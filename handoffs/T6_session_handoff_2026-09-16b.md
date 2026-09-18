# T6 세션 핸드오프 — 2026-09-16 (2차 세션: Phase 0 항목 1·3·4 완료) → 다음 챗

> **용도.** 챗 간 상태 인계. 상위 문서는 `T6_grassmannian_diffusion_noncoherent_handoff (1).md`(주제 핸드오프)와 `T6_research_companion.md`(동반 노트). 직전 문서는 `T6_session_handoff_2026-09-16.md`(1차 세션: 이론 상한 = 게이트 (0) 사전 점검)이며, **이 문서가 그 다음 상태**다. 이 세션의 산출물은 사용자가 모두 프로젝트 파일로 저장했다(§2).
> **언어 규약.** 대화는 한국어, 기술 용어는 영어, 수식은 LaTeX. 한 턴에 질문 하나 + 기본값 제안. go/no-go는 사용자의 명시적 결정 없이 넘어가지 않는다.
> **태그.** **[정확]** 가정 하 엄밀, **[근사]** 입력·수신기·동역학 고정 등, **[추측]** 수치·직관 근거만. **[V]** 프로젝트 파일 또는 명시된 사본에서 확인, **[M]** 기억 기반(인용 전 재확인), **[VERIFY]** 원문 없음.

---

## 0. 받는 Claude가 먼저 할 일

1. 프로젝트 파일에서 `gate_evidence_T6.md`, `decision_log_T6.md`, `open_questions_T6.md`, `scooping_log_T6.md`를 열어 현재 상태를 확인한다(모두 2026-09-16, 3차 항목까지 기록됨).
2. 사용자가 "이어서"/"다음 단계"라고 하면 §6의 순서대로 진행한다. **첫 항목: 동반 노트 Phase 0 항목 5 (Manolakos–Chowdhury–Goldsmith 2016과의 관계).**
3. `gate_evidence_T6.md` 상단의 **PRE-REGISTRATION(D-12)** 을 먼저 읽는다. 게이트 통과선·레짐·베이스라인·복잡도 예산이 실험 전에 고정되어 있으며, 사후 변경은 D-번호 + 재실행을 요구한다. 교수님 확인 대기 항목은 없다(모두 확정).
4. 매 세션 끝: 네 개의 기록 파일에 날짜가 붙은 항목을 추가하거나 "근거 추가 없음"을 적는다.
5. **파일명 주의**: 프로젝트에 `_T2`, `__1_` 접미사가 붙은 구버전 사본이 남아 있다(내용은 T6). **정본은 `*_T6.md`** 이며 이번 세션 갱신본이 최신이다. 다음 세션 산출물도 `_T6`로 통일할 것.

---

## 1. 이 세션에서 한 일 (Phase 0 항목 1, 3, 4)

### 1.1 항목 1 — arXiv 재확인 (월례 #1, 창: 2026-08-25 이후) **완료**
- 방법: arXiv search API가 샌드박스 공유 IP에서 **429 "Rate exceeded"** 로 차단 → **OAI-PMH 벌크 수확**으로 전환(D-09). eess+cs set에서 수정된 레코드 21,731건 수확, 그중 cs.IT/eess.SP/cs.LG/math.IT 6,569건의 제목+초록을 계층 키워드(CORE/TOPIC/GEN/COMM)로 채점. 보조: 월별 listing 스캔(cs.IT·eess.SP·cs.LG 2026-08/09, 7,512건 제목 + 후보 83건 초록), web search 3건, Semantic Scholar 피인용.
- **결과: 경쟁 논문(RED = 핵심 주제 + 생성모델 + 통신) 0건.** Grassmann/Stiefel 위 diffusion·flow prior를 noncoherent MIMO subspace prior 또는 constellation 설계에 적용한 사례 여전히 없음 → 주제 열림. 감시 논문(2510.15070 v3, 2605.04545 v1) 피인용 0.
- 새로 잡힌 도구·인접 논문:
  - **arXiv:2605.03588 v2 (2026-08-27) "Cartan flow matching"** (구제목 "Flow Matching on Symmetric Spaces"), Ruscelli–Zanchetta–Fioresi, cs.LG [V, 초록·버전 이력]. 대칭공간(Grassmannian 포함) 위 FM을 isometry group의 Lie algebra 부분공간 위 FM으로 환원 → geodesic interpolation 경로 구성 불필요. **실** Grassmannian만 시연. v2 comments: "Major revision. Fixed a mistake in section 4." → **Q-10**(생성 경로 선택).
  - arXiv:2603.24829 (2026-03-25) "Flow matching on homogeneous spaces", Ruscelli [V, 초록]. 코드 github.com/fresh999/HomogeneousFM. (창 밖, web search로 발견.)
  - arXiv:2608.28073 (2026-08-28) Jensen–Zimmermann, second-order Stiefel retraction의 조건수·보간 오차 한계 [V, 초록] → exp/log 대체 저비용 retraction 근거(기준 (iii)).
  - arXiv:2604.19904 v2 (2026-08-30) Khirwadkar–Rajamäki–Pal, "Grassmannian-Coded Beamforming for mmWave Channel Sensing with Unknown Complex Path Gain" [V, 초록]. 미지 복소 이득 하 DoA↔beamspace 부분공간의 Grassmann 기하 → 비경쟁이지만 우리 $\mathbf H_b=\mathbf G_b\mathbf U_b^H$와 같은 기하, **Q-02** 참고문헌.
- 한계: 어휘 기반 스캔이라 "spatial signature"류 표현만 쓴 논문은 놓칠 수 있음. 다음 재확인: **2026-10-16**, `--since 2026-09-16`.

### 1.2 항목 3 — 차동 USTM 사전 추정 (Q-06) **완료**
상세: `T6_diff_ustm_note.md`. 코드 `t6_diff_ustm.py`(+`--selftest`), `t6_diff_summarize.py`, `t6_twoblock_ref.py`. seed 20260916.

- **원문 확인 [V]**: Hochwald–Sweldens, IEEE TCOM 48(12):2041–2052, 2000. 사본 `http://rywei.ce.ncu.edu.tw/course/94/codedmodulation/dstm.pdf`(강의 자료). 핵심: 블록 길이 $=M$, $\Phi_\ell=\tfrac1{\sqrt2}[\mathbf I_M;\mathbf V_\ell]$(식 (17)), 송신 $\mathbf S_\tau=\mathbf V_{z_\tau}\mathbf S_{\tau-1}$(식 (18)), 수신 기본식 $\mathbf X_\tau=\mathbf V_{z_\tau}\mathbf X_{\tau-1}+\sqrt2\mathbf W'$(식 (22), 잡음 2배 ⇒ ≈3 dB), ML 복조 = $\arg\min_\ell\|\mathbf X_\tau-\mathbf V_\ell\mathbf X_{\tau-1}\|_F$(식 (21)). **프로젝트에 PDF 추가 권장([M]→[V] 확정용).**
- **우리 확장 [근사, 우리 것]**: $T>M$이므로 시간축 유니터리로 일반화 — $\{\mathbf V_z\}_{z=1}^K\subset U(T)$ Haar, $\mathbf X_b=\mathbf V_{z_b}\mathbf X_{b-1}$, rate $\log_2K/T$(per-block과 동일), 수신 $\hat z=\arg\max_z\mathrm{Re\,tr}(\mathbf V_z^H\mathbf Y_b\mathbf Y_{b-1}^H)$ [정확 항등식, 단위시험].
- **동역학 [근사]**: $\mathbf U$ 고정(부분공간 유지), $\mathbf G_b=\alpha\mathbf G_{b-1}+\sqrt{1-\alpha^2}\mathbf E_b$. Jakes 1-lag 대응 $\alpha=J_0(2\pi f_DTT_s)$ [M]. $T=16$: $\alpha=0.984\leftrightarrow f_DT_s=0.0025$(HS 시뮬값), $0.9\leftrightarrow0.0064$, $0.5\leftrightarrow0.0151$, $0\leftrightarrow0.0239$. ($T=8$은 2배.)
- **결과 — SER $10^{-2}$ 동작점 (dB; MC ±0.3–0.5, 꼬리 ±1)**

| 구성/$K$ | struct-ML | genie-U | diff α=1 | diff 0.984 | diff 0.9 | diff ≤0.5 | diffU α=1 | diffS2 α=1 |
|---|---|---|---|---|---|---|---|---|
| (8,2,16,2)/64 | −0.04 | −2.19 | **+2.17** | +1.58 | +4.61 | floor | −1.35 | +0.96 |
| (8,2,16,2)/1024 | +1.56 | 0.00 | **+3.41** | +3.39 | floor | floor | +0.55 | +2.74 |
| (16,2,32,2)/64 | −5.87 | −9.63 | **−1.95** | −2.95 | −0.68 | floor | −7.05 | −4.29 |
| (16,2,32,2)/1024 | −4.56 | −7.40 | **−0.80** | −0.71 | floor | floor | −5.46 | −1.48 |

(diffU = $\mathbf U$ 기지 + 차동, diffS2 = 두 블록 SVD 사영 + 차동. "floor" = $10^{-2}$ 위 오류 바닥.)

- **2T 참조**(채널이 두 블록 상수일 때 prior-free 상한: 블록 $2T$, $K^2$점, 같은 rate): $(16,2,16,2)$ $K{=}4096$ → struct **−2.40**, genie −3.97 (per-block 대비 **+2.4 dB**, genie 이득 +2.1과 동급). $(32,2,32,2)$ $K{=}4096$ → struct **−7.50**, genie −10.86 (**+1.6 dB**, genie 이득 +3.8의 절반 이하).
- **결론 [근사]**: 차동 USTM은 T6 레짐에서 약한 베이스라인 — 채널이 완전 상수여도 per-block struct-ML보다 +1.9~2.2 dB($N/r{=}8$), +3.8~3.9 dB($N/r{=}16$) 열세. 이유: $T\ge4M$에서 per-block 손실이 이미 작고, 차동은 $TN$차원 잡음 reference(2배 잡음)를 쓰며 저랭크 구조를 버림(HS는 $T=M$·$N=1$ 설계). 이동성은 악화만 시킴.
- **부분공간 지식의 값**: diffU − diff = 3.5/2.9 dB($N/r{=}8$), 5.1/4.7 dB($N/r{=}16$) → array gain 서사 유지. prior-free 두 블록 SVD는 0.7–2.3 dB만 회수.
- **레짐 단서 [근사→추측]**: $\alpha\approx1$이면 "블록을 $2T$로" 만으로 $N/r{=}8$에서 genie와 동급 이득 → **T6의 가치가 최대인 곳은 "이득 $\mathbf G_b$는 블록 간 탈상관, 부분공간 $\mathcal S_b$는 유지"인 레짐**(D-11). 동반 노트 Phase 1의 "산발 패킷이면 Framing A 근거가 사라진다"는 **정정 대상**(산발 패킷은 이득 상관을 죽이지만 AoA 상관은 남길 수 있음 → Q-12에서 검증).
- 부수: $\alpha=0$에서도 차동 SER은 우연($1-1/K$)보다 낮음 — 부분공간 지속성이 두 블록 수신기에 보인다는 확인 [정확 대수].

### 1.3 항목 4 — subspace tracking 동역학 가정 확인 **완료**
상세: `T6_tracker_assumptions_note.md`.
- **Saad-Falcon–Ancelin–Romberg, arXiv:2402.10352 v1 (2024, IEEE SAM 투고) [V, arXiv HTML 열람]**: 정규화 최소제곱 + Grassmann 매끄러움 벌칙. 동역학은 (a) static $\mathbf Y_{t+1}=\mathbf Y_t$, 벌칙 $\sum_t d(\mathbf Y_t,\mathbf Y_{t+1})^2$(식 (3)), (b) constant velocity $\mathbf Y_{t+1}=\exp_{\mathbf Y_t}(\mathbf H_t)$, $\mathbf H_{t+1}=\Gamma\mathbf H_t$(식 (4)–(5)) 둘뿐. **실수** Grassmannian. Gaussian·Kalman·process noise·jump·abrupt 언급 **0회**(전수 검색). 시연은 narrowband beamforming(검출 되먹임 없음).
- GROUSE, PETRELS [M, 원문 미확보]: 느린 변화 가정의 적응 필터, 명시적 동역학·잡음 모델 없음.
- **함의**: T6의 차별화 주장(점프·경로 생성/소멸·복소 $\mathcal G(N,r)$은 모델 밖) 성립 [V]. 단 공정성을 위해 tracker와 게이트 (ii)의 접공간 Gaussian AR prior 양쪽에 점프 인지(innovation 검정 + reset) 변형을 넣어야 함 → **Q-13**. D-11 레짐에서도 살아남는 유일한 prior-free/모델 기반 경쟁자가 tracker + 조건부 GLRT이므로 **이것이 (i)의 사실상 기준 베이스라인**이다.

---

## 2. 프로젝트 파일 맵 (다음 챗에서는 `/mnt/project/`에서 읽는다)

1차 세션까지: `T6_grassmannian_diffusion_noncoherent_handoff (1).md`, `T6_research_companion.md`, `research_report_2_diffusion_wireless_IT.md`, `Multiple-antenna_capacity_in_correlated_Rayleigh_fading_with_channel_covariance_information.pdf`(JG 2005), `T6_session_handoff_2026-09-16.md`, `T6_ceiling_derivation.md`, `t6_ceiling.py`, `test_t6_ceiling.py`, `t6_summarize.py`, `t6_ceiling_results.json`, `t6_ceiling_summary.txt/.csv`, `t6_ceiling_plots.png`.

이번 세션 추가:

| 파일 | 역할 |
|---|---|
| `T6_session_handoff_2026-09-16b.md` | 이 문서 |
| `T6_diff_ustm_note.md` | Q-06 유도·결과 노트(HS 재정리, 확장, 표, 해석, 한계) |
| `T6_tracker_assumptions_note.md` | Phase 0 항목 4 확인 노트 |
| `t6_arxiv_oai_check.py` | 월례 scooping 정본 도구(OAI-PMH 벌크 + 계층 채점) |
| `t6_arxiv_listing_check.py` | 보조 scooping(월별 listing 스크래핑) |
| `t6_oai_2026-09-16.json`, `t6_listing_2026-09-16.json` | scooping 결과 |
| `t6_diff_ustm.py` | 차동 USTM 시뮬레이션(`t6_ceiling.py` import, `--selftest`) |
| `t6_diff_summarize.py` | 표·그림 |
| `t6_twoblock_ref.py` | 2T 블록 참조 |
| `t6_diff_results.json/_summary.txt/_plots.png`, `t6_twoblock_ref.json` | 결과 |
| `gate_evidence_T6.md` / `decision_log_T6.md` / `open_questions_T6.md` / `scooping_log_T6.md` | 기록 파일(정본, 3차 항목까지) |

---

## 3. 결정 요약 (`decision_log_T6.md`)
1차 세션: D-01 prior 대상 정의 · D-02 $\gamma=N/r$ · D-03 게이트 (0) 기준 · D-04 입력 USTM + random-Haar 유한 $K$ · D-05 SNR 격자 · D-06 $F_D$ 폐형 · D-07 상한·정확 간극 병기 · D-08 게이트 (0) 조건부 통과(잠정).

이번 세션:
- **D-12** 게이트 pre-registration(§5 참조).
- **D-09** 월례 scooping 도구를 `t6_arxiv_oai_check.py`(OAI-PMH)로 고정. search API는 429로 차단됨.
- **D-10** 차동 USTM은 **2순위 베이스라인**(보고는 하되 (i)의 기준 아님). 기준 prior-free 베이스라인 = per-block struct-ML(1순위) + tracker/조건부 GLRT(Phase 3에서 순위 확정).
- **D-11**(확정, 단 Q-12 실패 시 자동 철회) Q-05 레짐에 동역학 조건 추가: "$\mathbf G_b$는 블록 간 탈상관($\alpha\lesssim0.5$), $\mathcal S_b$는 수십 블록 유지". 주 레짐 $(16,2,32,2)$류·−4~−8 dB·짧은 블록은 유지.

## 4. 열린 질문 (`open_questions_T6.md`)
- ~~Q-05~~ **닫힘**(2026-09-16 확정: 주 레짐 $(16,2,32,2)$류, 보조 $(8,2,16,2)$, D-11 동역학 조건 포함).
- **Q-12(최우선, Phase 1/2)** 그 레짐의 물리적 타당성(ray-tracing으로 블록 간 $|G$ 상관$|$ vs subspace chordal distance 측정). D-11의 존폐를 결정하는 관문 — **안 나오면 D-11 자동 철회**하고 $\alpha\approx1$ 레짐으로 복귀(Q-11의 경쟁자가 (i) 베이스라인에 진입).
- **Q-13(신설, Phase 3)** tracker·AR prior에 점프 인지(reset) 변형 포함 여부. 기본값: 포함(무변형/reset 2종씩).
- **Q-11(신설, Phase 1/3)** $\alpha\approx1$ 보조 레짐의 경쟁자: 2T-블록 USTM, full-channel Kalman + coherent 검출.
- **Q-10(신설, Phase 2)** 생성 경로: geodesic RFM(1순위 유지) vs Cartan FM(2605.03588) vs homogeneous-space FM(2603.24829).
- Q-06 **잠정 닫힘**(Phase 3에서 HS group constellation·구조 인식 pair-ML·$\mathbf U$ 드리프트 재확인).
- 계속 열림: Q-02(array-manifold-aware per-block 수신기; 참고 2604.19904), Q-09($\gamma=1$ 감도), Q-04(저 SNR 입력), Q-03($r$ 고정), Q-01(steering $\mathbf G_b$), Q-07(일반 $(M,r)$ 폐형), Q-08(연속 입력 prior-free 정확 MI).

## 5. 게이트 상태와 미결 (`gate_evidence_T6.md`)
- **(0) 조건부 통과(D-08, 잠정)** — 변동 없음.
- **(i)** 직접 근거 없음. 보조 근거만 추가: 베이스라인 순위(차동 제외), 부분공간 지식의 차동 대비 값, 2T 참조.
- **(ii)(iii)** 근거 없음.
- **미결 없음 (2026-09-16 4차 갱신).** 교수님은 게이트에 개입하지 않고 사용자 자율로 진행하므로 세 건 모두 확정: D-08 게이트 (0) 통과 확정, Q-05 주 레짐 확정(닫힘), D-11 채택(단 **Q-12 실패 시 자동 철회**).
- **D-12 PRE-REGISTRATION**: 게이트 통과선을 실험 전에 고정했다 — (i) ≥1.0 dB @ BLER $10^{-2}$(주 레짐, 최고 베이스라인 대비), (ii) 동일 예산 튜닝한 AR prior와 차이 ≥0.5 dB(미만이면 논문 축 변경), (iii) 블록당 GLRT 대비 ≤50×(샘플링 단계를 늘려 (i)를 사는 것은 위반), no-go 시 2쪽 노트 후 중단. **기준을 바꾸려면 D-번호 + 영향 실험 전면 재실행.** 전문은 `gate_evidence_T6.md` 상단.

---

## 6. 다음 챗의 진행 순서

**Phase 0 체크리스트 (동반 노트 §2)**: 1 arXiv ✔ · 2 JG ✔ · 3 차동 USTM ✔ · 4 tracker 가정 ✔ · **5 Manolakos–Chowdhury–Goldsmith 2016 [M] 관계 정리 — 미완(다음 첫 항목)** · 6 이론 상한 ✔.

1. **항목 5**: 각도 sparse 채널의 energy-based noncoherent/massive SIMO 검출과 우리 설정의 관계. 저 SNR·많은 수신 안테나 레짐이 D-11과 겹치므로 중요. 원문이 프로젝트에 없으므로 web으로 초록·모델 확인 + [M] 정리, 원문 PDF는 요청 목록에. 산출: "게이트 (i) 베이스라인에 energy-based 검출을 넣을 것인가" 판단.
2. **교수님 보고 3줄 갱신**(§8) — 위 3건 확인 요청 포함.
3. **Phase 1 착수**: Q-05 레짐 고정(D-11 포함) → $r$ 고정(Q-03, 기본 $r=2$) → 이득 정의(SER $10^{-2}$ dB 간극) → 정확/근사/추측 표.

선택적 후속(요청 시): (a) $N=64$ 구성, (b) $\gamma=1$ 감도, (c) Q-02용 DFT-beamspace 수신기, (d) 차동 구조 인식 pair-ML($\mathbf U$ 기지 시 우리 $F_T$ 폐형으로 계산 가능; §`T6_diff_ustm_note.md` §2), (e) Q-12용 ray-tracing 파일럿.

---

## 7. 코드 사용법과 주의

```bash
# scooping (월례)
python3 t6_arxiv_oai_check.py --since 2026-09-16 --sets eess,cs \
    --cats cs.IT,eess.SP,cs.LG,math.IT --out t6_oai_2026-10-16.json     # 약 5분

# 차동 USTM (재현)
python3 t6_diff_ustm.py --selftest                                       # 4개 시험
python3 t6_diff_ustm.py --snr="-8:6:1;-14:0:1" --R 200 --B 21 \
    --n_pb 4000 --seed 20260916 --out t6_diff_results.json               # 약 9분
python3 t6_diff_summarize.py t6_diff_results.json
python3 t6_twoblock_ref.py                                               # 약 15분
```
- `t6_diff_ustm.py`는 `/mnt/project`의 `t6_ceiling.py`를 import한다(`--ceiling_dir`).
- 음수 SNR 인자는 반드시 `--snr="-8:..."` 형태(등호 필수).
- 긴 실행은 `setsid nohup ... &`로 분리하고 로그를 폴링한다. **다만 이번 세션에서 분리 실행이 중간에 죽은 사례 있음**(`t6_twoblock_ref.py` 2번째 case가 −9 dB 이후 종료) → 로그의 마지막 SNR과 JSON의 case 수를 항상 대조하고, 빠진 점은 인라인으로 보충한다. (현재 `t6_twoblock_ref.json`의 case 2는 −14…−9 dB가 $n=1000$, −8…−6 dB가 별도 RNG·$n=600$ — 노트에 명시됨.)
- 함정: ① arXiv search API 429(도메인 전체 차단, 백오프 무효) → OAI-PMH 사용. ② 키워드 매칭은 짧은 토큰(≤5자)만 양쪽 word boundary, 어간(`grassmann`)은 접미사 허용 — 수정 전에는 `ustm`이 `adjustment`에, `prior`가 `priority`에 걸렸다. ③ 샌드박스 셸에 `time` 없음. ④ `brentq`의 Bessel 근 bracket은 $[10^{-9}, 2.41]$.
- 제약(1차 세션에서 이월): 폐형은 $M=r=2$ 전용(Q-07), 유한 $K$의 SER $10^{-2}$ 정밀도 ±0.3–0.7 dB.

---

## 8. 교수님 보고용 3줄 (갱신)
1. 차동 USTM(Hochwald–Sweldens 2000)을 우리 구조 채널에 맞게 구현해 재보니, 채널이 두 블록에서 완전 상수여도 per-block 구조 인식 ML보다 2 dB($N/r{=}8$)–4 dB($N/r{=}16$) 나쁘다. 따라서 차동은 우리 레짐의 기준 베이스라인이 아니며, 기준은 per-block 구조 인식 ML과 subspace tracking + 조건부 GLRT다.
2. 다만 채널 전체가 두 블록 동안 상수이면 "블록 길이를 2배로" 늘리는 것만으로 $N/r{=}8$에서는 부분공간을 아는 이득과 같은 +2.4 dB를 얻는다($N/r{=}16$에서는 +1.6 dB로 절반 이하). 그래서 T6의 가치가 가장 큰 레짐은 **채널 이득은 블록 간 탈상관하되 수신 부분공간(AoA)은 유지되는 경우**이며, 레짐 정의에 이 조건을 넣자고 제안드린다.
3. 모델 기반 subspace tracker(Saad-Falcon 2024)의 동역학은 정지·등속 geodesic 두 가지뿐이고 점프·경로 생성/소멸을 다루지 않음을 원문에서 확인했다 — 학습 prior의 차별화 지점이 여기이며, 공정 비교를 위해 tracker에도 점프 감지·재초기화 변형을 붙여 비교할 계획이다. 이번 달 arXiv 재확인에서 경쟁 논문은 없었다.
