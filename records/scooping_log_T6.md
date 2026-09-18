# scooping_log.md — T6

월 1회 arXiv(cs.IT, eess.SP, cs.LG) 재확인. 키워드(동반 노트 §7 / 핸드오프 §6):
`Grassmannian diffusion`, `Grassmannian diffusion constellation`, `Riemannian generative noncoherent`, `Riemannian generative noncoherent MIMO`, `subspace sequence generative prior`, `noncoherent detection learned prior`, `partially coherent MIMO learned`, `Grassmann flow matching`, `subspace tracking diffusion`, `generative constellation design`, (2026-09-16 추가) `Cartan flow matching`, `flow matching symmetric spaces`, `flow matching homogeneous spaces`, `beamspace Grassmannian`.
감시 그룹: Grassmannian constellation 설계 그룹(arXiv:2605.04545 "Z-Opt", arXiv:2510.15070 geodesic mapping 저자), 다양체 생성모델 그룹(Chen–Lipman, De Bortoli), Romberg(subspace tracking), Guillaud/Ngo(noncoherent 설계), (2026-09-16 추가) Ruscelli–Zanchetta–Fioresi(대칭공간/등질공간 flow matching), Khirwadkar–Rajamäki–Pal(beamspace Grassmannian sensing).

도구(D-09): `python3 t6_arxiv_oai_check.py --since <직전 점검일> --sets eess,cs --cats cs.IT,eess.SP,cs.LG,math.IT --out t6_oai_<날짜>.json` (OAI-PMH 벌크 수확, 초록 포함, ~5분). arXiv search API는 샌드박스에서 429로 차단됨. 보조: `t6_arxiv_listing_check.py`(월별 listing 페이지), web search.

```
[2026-08-25] 기준 조사(research_report_2): Grassmann/Stiefel/유니터리 군 위 score-based diffusion을 noncoherent MIMO subspace prior 또는 codebook 설계에 적용한 사례 미발견 → 열림.
[2026-09-16] 이번 세션은 이론 상한 사전 점검에 집중. arXiv 재확인 미실시(다음 세션 첫 항목, 동반 노트 Phase 0 체크 1).
[2026-09-16] 월례 재확인 #1 — 창: 2026-08-25 이후 신규·개정 레코드. 실행: t6_arxiv_oai_check.py(결과 t6_oai_2026-09-16.json) + t6_arxiv_listing_check.py + web search + Semantic Scholar.
방법: arXiv search API가 IP 단위 429("Rate exceeded")로 차단 → (a) OAI-PMH 벌크 수확: eess+cs set에서 2026-08-25 이후 수정된 레코드 21,731건, 그중 cs.IT/eess.SP/cs.LG/math.IT 6,569건의 제목+초록을 계층 키워드(CORE/TOPIC/GEN/COMM)로 채점; (b) 월별 listing 페이지(cs.IT, eess.SP, cs.LG의 2026-08/09, 7,512건) 제목 스캔 + 후보 83건 초록 확인; (c) web search 핵심 구문 3건; (d) Semantic Scholar 피인용(2510.15070, 2605.04545, 2605.03588 모두 0).
결과 — 경쟁(RED = 핵심 주제 + 생성모델 + 통신): 0건. Grassmann/Stiefel 위 diffusion·flow prior를 noncoherent MIMO의 subspace prior 또는 constellation 설계에 적용한 사례 여전히 미발견 → 열림 유지.
방법 쪽 도구 논문(프로젝트 파일에 없던 것, 인용 후보):
 - arXiv:2605.03588 v2 (2026-08-27) "Cartan flow matching"(구제목 "Flow Matching on Symmetric Spaces"), Ruscelli–Zanchetta–Fioresi, cs.LG [V, 초록·버전 이력 확인]. 대칭공간(Grassmannian 포함) 위 flow matching을 isometry group의 Lie algebra 부분공간 위 FM으로 재정식화하여 geodesic interpolation 경로 구성을 피함; 실 Grassmannian SO(n)/S(O(k)×O(n−k))에서 시연. v2 comments: "Major revision. Fixed a mistake in section 4." → 우리 생성기(동반 노트 결정: Chen–Lipman geodesic RFM)의 대안 후보. 복소 G(N,r)=U(N)/(U(r)×U(N−r))도 대칭공간이므로 원리상 적용 가능하나 논문은 실수 버전만 다룸. → Q-10 신설.
 - arXiv:2603.24829 v1 (2026-03-25) "Flow matching on homogeneous spaces", Ruscelli, cs.LG [V, 초록만]. G/H 위 FM을 G로 들어올려 Lie algebra 위 유클리드 FM으로 환원, premetric/geodesic 불필요; 코드 github.com/fresh999/HomogeneousFM. (창 밖, web search로 발견.)
 - arXiv:2608.28073 (2026-08-28) R. Jensen, R. Zimmermann, "Conditioning and interpolation error bounds for second-order Stiefel retractions with closed-form inverses", math.NA [V, 초록만]. Cayley·polar-light retraction의 보간 오차 한계 → exp/log 대신 쓸 저비용 retraction의 근거(기준 (iii), tracker/생성기 구현).
인접(비경쟁):
 - arXiv:2604.19904 v2 (2026-08-30) Khirwadkar–Rajamäki–Pal, "Grassmannian-Coded Beamforming for mmWave Channel Sensing with Unknown Complex Path Gain", eess.SP/cs.IT [V, 초록·버전 확인]. 미지 복소 이득 하에서 DoA 후보가 beamspace 응답을 통해 부분공간으로 사상되는 Grassmann 기하; Grassmannian packing을 beamspace code로 사용(sensing, 단일 RF chain). → 우리 H_b = G_b U_b^H의 "미지 이득 × 구조 부분공간" 관점과 같은 기하; Q-02(array-manifold-aware 수신기) 참고문헌.
감시 논문: 2510.15070 v3(2026-08-08, 창 밖) 외 변동 없음; 2605.04545 v1 그대로.
무관 매칭(기록용): 2609.08312 non-coherent AirFL; 2608.30104 noise-modulation FAMA; 2608.24880·2608.27041 D-MIMO noncoherent processing(localization/ISAC); 2609.15488 kernel-PCA subspace tracking; 2608.26076 Lie-group MeanFlow(robotics); Grassmann code 코딩이론 3건.
한계: 어휘 기반 스캔이라 핵심어 없이 "spatial signature"류 표현만 쓴 논문은 놓칠 수 있음(TOPIC 계층과 web search로 보완); 저자 감시는 성(姓) 스캔이라 동명이인 다수, 관련 히트는 Ruscelli·Zimmermann뿐.
다음 재확인: 2026-10-16 (--since 2026-09-16).
[2026-09-16 3차] scooping 추가 없음. 참고: Hochwald–Sweldens 2000 원문 사본 위치 http://rywei.ce.ncu.edu.tw/course/94/codedmodulation/dstm.pdf (강의 자료; 프로젝트에 PDF 추가 권장, 태그 [M]→[V]).
[2026-09-16 3차, 추가] Phase 0 항목 4 확인용으로 arXiv:2402.10352 v1 본문 열람(후속 버전 없음, v1만). 키워드에 `subspace tracking dynamical model` 추가 권장.
[2026-09-17] scooping 추가 없음(다음 월례 재확인 2026-10-16, --since 2026-09-16). 원문 PDF 2건이 프로젝트에 추가됨: Hochwald–Sweldens 2000 (Differential_unitary_space-time_modulation.pdf; 식 (17)(18)(19)–(20)(21)(22)와 "twice the variance / approximately a 3-dB performance loss" 서술 위치 재확인, 텍스트 레이어에 수식 기호 없음 → 기호 단위 재확인 시 래스터화), Manolakos–Chowdhury–Goldsmith 2016 (Energy-Based_Modulation_for_Noncoherent_Massive_SIMO_Systems.pdf; 식 (1)(2)(3)(8)(9)(21) 확인). 태그 [M]/[web 사본] → [V, 프로젝트 PDF]. 미확보 원문: Chowdhury–Manolakos–Goldsmith, IEEE TIT 62(4) 2016 [M](에너지 검출 = noncoherent ML 증명; 필요 시 요청 목록), GROUSE/PETRELS [M].
[2026-09-17 2차] scooping 추가 없음. Q-12 (a)에서 쓴 [M] 출처(인용 전 확인 목록): one-ring 국소 산란 모델 — D.-S. Shiu, G. Foschini, M. Gans, J. Kahn, IEEE TCOM 48(3), 2000 [M]; Clarke/Jakes J0 상관 [M]; 3GPP TR 38.901 클러스터 내 각퍼짐(BS 쪽) [M, VERIFY]. Q-12 (b) 도구: Sionna RT / DeepMIMO(설치 가능성 미확인).
[2026-09-17 3차] scooping 추가 없음. Q-12 (b) 도구 기록: sionna-rt 2.1.0 / mitsuba 3.9.1 / drjit 1.5.0 (pip, 샌드박스 1 CPU에서 동작). [VERIFY] 필요 문헌: 28 GHz 이동 환경 시변 채널 측정(보행자·차량 이동 시 coherence time, Doppler 퍼짐) — D-11을 "동적 클러터 환경"으로 되살리려면 필수; 후보 키워드 "28 GHz channel temporal variation measurement pedestrian", "mmWave Doppler spread measurement urban".
```
