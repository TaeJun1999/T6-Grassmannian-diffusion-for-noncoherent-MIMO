# decision_log.md — T6 (Riemannian diffusion prior on the receive-side subspace, noncoherent MIMO)

형식: `[날짜] D-번호 결정 / 이유 / 기각한 대안`. 기술 용어는 영어, 수식은 LaTeX.

```
[2026-09-16] D-01 결정: prior의 대상은 채널의 receive-side subspace S_b = span(H_b^H) ∈ G(N, r). 송신 부분공간 span(X_b) ∈ G(T, M)은 데이터이며 prior를 두지 않는다.
이유: 핸드오프 §1 Framing A·§2, 동반 노트 §1.2와 일치. 어긋나는 문장 목록은 T6_ceiling_derivation.md §0의 표에 기록(핸드오프 §4 M1의 G(4,1)/G(8,2) 표기, §0 규약의 G(N,r) 부재, 동반 노트 §1.2의 "partially coherent"·§1.3의 Zheng–Tse 점근식, research_report_2의 정정 전 표현).
대안: 없음(정의 확인).

[2026-09-16] D-02 결정: 채널 정규화는 gamma = N/r (H_b = G_b U_b^H, G_b iid CN(0, gamma), E||H_b||_F^2 = MN, 안테나당 평균 수신 SNR = rho). gamma = 1은 코드 옵션(--gamma_mode 1)으로만 남긴다.
이유: JG의 단위대각 상관행렬 관례(tr R^r = N)와 물리적 array gain에 부합. genie 채널은 SNR gamma·rho의 iid noncoherent (T, M, r) 채널이 된다(+9 dB, +12 dB).
대안: gamma = 1 (extra antennas give no array gain) — 비물리적이라 기각.

[2026-09-16] D-03 결정: 게이트 (0)의 기준 비교는 "genie(U_b 기지) vs 같은 구조 채널의 prior-free per-block 수신기"로 둔다. 1순위 structure-aware per-block ML(등방 prior 주변화, 정확 pdf), 2순위 SVD-projection(실용). i.i.d. Rayleigh 채널의 rate는 맥락용으로만 보고한다.
이유: 시간 prior의 가치는 정확히 I(X; U | Y) = I(X;Y|U) − I(X;Y)이고, 이는 같은 채널에서 per-block 수신기와 genie 사이의 간극이다. 동반 노트 §1.3의 "i.i.d. noncoherent rate와의 차이"는 다른 채널과의 비교라 게이트에 맞지 않는다.
대안: Zheng–Tse 점근식 대비 — 유한 SNR·다른 채널이라 기각.

[2026-09-16] D-04 결정: 입력은 isotropic USTM(연속)과 random-Haar 유한 constellation(K = 64, 1024 전 SNR; K = 8192 저 SNR)로 고정하고, 결과를 "입력 고정 achievable rate"로 표기한다(용량이 아님 → [근사]).
이유: 용량달성 입력(특히 저 SNR)은 미지. 연속 입력의 prior-free 정확 MI는 중첩 Haar 적분이라 계산 불가 → GMI 하한과 유한 K 정확 MI로 대체.
대안: 중첩 MC — 사후 집중으로 SNR ≥ 0 dB에서 실패, 기각. importance sampling — 추후 필요 시(Q-08).

[2026-09-16] D-05 결정: SNR 격자는 −10…20 dB(2 dB 간격). 사용자 요청 0–20 dB에 −10…−2 dB를 추가.
이유: 간극이 저 SNR에 집중될 것이라는 예상을 확인하기 위해(실제 결과가 그러함 — gate_evidence.md).

[2026-09-16] D-06 결정: rank-2 Bingham 적분 F_D의 자체 폐형(T6_ceiling_derivation.md §3.2)을 채택. 근거 밀도(truncated Haar unitary, det(I − Z Z^H)^{D−4})는 [M]이지만 단위 시험 T1/T2/T2b/T4/T6로 수치 검증(8/8 통과).
이유: 유한 K의 structure-aware ML을 정확 pdf로 계산할 수 있고, genie MI를 단일 수준 MC로 정확히 계산할 수 있다.
대안: 매 블록 U에 대한 MC 주변화 — 비용·분산 때문에 기각.

[2026-09-16] D-07 결정: 게이트 (0) 판정은 "상한(연속 USTM, R_genie − GMI_svd)"과 "정확 간극(유한 K, 저 SNR)"을 함께 보고하고, 판정 단위는 (i) 고정 SNR에서 rate 간극(%), (ii) SER 10^{-2}에서 SNR 간극(dB)으로 한다(동반 노트 §8과 같은 단위).
이유: 연속 입력에서는 상한만, 유한 K에서는 저 SNR에서만 정확값이 있으므로 둘을 병기해야 판정이 정직하다.

[2026-09-16] D-08 결정(잠정, 교수님 확인 전): 게이트 (0)은 "조건부 통과"로 기록. 조건 = 목표 레짐을 저 SNR 동작점(≤ 0 dB)·큰 N/r(≥ 8)·짧은 블록으로 고정. 기본 레짐 후보는 (16,2,32,2)류, (8,2,16,2)는 경계선.
이유: gate_evidence.md — 시간 prior의 최대 가치가 (i)의 목표를 넘는 레짐이 저 SNR·큰 N/r에만 있음. 10 dB 이상에서는 rate 간극 < 5%.
대안: 무조건 통과 — 근거 없음. no-go — 큰 N/r 저 SNR 레짐에서 3–4 dB 여유가 있으므로 이르다.

[2026-09-16] D-09 결정: 월례 scooping 도구를 arXiv OAI-PMH 벌크 수확 스크립트(t6_arxiv_oai_check.py; --since <직전 점검일>, sets eess,cs, cats cs.IT,eess.SP,cs.LG,math.IT)로 고정하고, 결과 JSON을 날짜를 붙여 프로젝트에 보관한다. web search는 보조.
이유: arXiv search API는 샌드박스 공유 IP에서 429("Rate exceeded")로 차단됨. OAI-PMH는 별도 제한이며 초록 포함 전수 수확이 5분 내 가능(2026-09-16: 21,731건). CORE/TOPIC/GEN/COMM 계층의 어휘 채점은 재현 가능하고 오탐을 눈으로 걸러낼 수 있는 규모(10여 건)로 줄여 준다.
대안: search API(차단), 월별 listing 페이지 스크래핑(제목만 전수, 초록은 후보만; t6_arxiv_listing_check.py로 보조 유지), web search 단독(재현 불가, 기각).

[2026-09-16] D-10 결정: 차동 USTM은 Phase 3에서 보고하되 게이트 (i)의 기준 베이스라인으로 삼지 않는다. 기준 prior-free 베이스라인 = per-block structure-aware ML(1순위), subspace tracker + 조건부 GLRT(Phase 3 측정 후 순위 확정). 차동은 HS 2000 [V, web 사본]의 tall-block 확장(Haar U(T))으로 구현하고, Phase 3에서 HS group/diagonal constellation과 구조 인식 pair-ML로 재확인한다.
이유: T6_diff_ustm_note.md §4-5 — 채널이 완전 상수여도 per-block struct-ML보다 2-4 dB 열세(T>=4M, 저 SNR, 큰 N/r에서 reference 잡음 2배·저랭크 구조 미사용). 이동성은 악화만 시킴.
대안: 차동을 최고 베이스라인으로 두고 우리 이득을 그 대비로 정의 — 기각(쉬운 상대와 비교하는 셈이 되어 (i)의 신뢰를 떨어뜨림).

[2026-09-16] D-11 결정(잠정, 교수님 확인 전): Q-05 레짐 기본값에 동역학 조건을 추가한다 — "블록 간격이 이득 G_b의 coherence time보다 길거나(산발 짧은 패킷) Doppler가 커서 G_b는 블록 간 탈상관(alpha <~ 0.5)하되, receive-side subspace S_b(AoA)는 수십 블록 동안 유지". 주 레짐 (16,2,32,2)류·동작점 -4~-8 dB·짧은 블록은 유지.
이유: 2T 참조(gate_evidence 3차 (c)) — alpha≈1이면 N/r=8에서 긴 블록/full-channel tracker가 T6와 대등. T6의 가치가 최대인 곳은 이득 탈상관·부분공간 유지 레짐이며, 여기서 차동·긴 블록·full-channel tracker는 모두 무력하고 per-block struct-ML과 subspace tracker만 남는다. 동반 노트 Phase 1의 "산발 패킷이면 Framing A 근거가 사라진다"는 이득 상관에만 해당하므로 정정 대상.
대안: 연속 블록·느린 페이딩(alpha≈1) 레짐 — N/r=16이면 여전히 T6가 유리(+3.8 vs +1.6 dB)하나 경쟁자가 늘어 (i) 입증 부담 큼; 보조 레짐으로만 유지.

[2026-09-16] D-12 결정: 게이트 기준을 실험 전에 고정한다(pre-registration, gate_evidence_T6.md 상단 블록). 통과선 (i) >= 1.0 dB @ BLER 1e-2, 주 레짐 (16,2,32,2)류, 최고 베이스라인 = per-block struct-ML과 tracker+조건부 GLRT 중 우수한 쪽; (ii) 동일 예산으로 튜닝한 접공간 Gaussian AR prior와의 차이 >= 0.5 dB, 미만이면 논문 축 변경; (iii) 블록당 GLRT 대비 <= 50x, 샘플링 단계 수를 늘려 (i)를 사는 것은 위반. no-go 시 2쪽 내부 노트 후 중단.
이유: 교수님이 게이트에 개입하지 않고 사용자 자율로 진행 → 결정자와 실행자가 같아 사후에 기준이 이동할 위험이 큼. (i)의 수치가 아직 하나도 없는 시점에 고정해야 "유리한 방향을 모르는 상태의 기준"이 된다. 변경 시 D-번호 + 영향 실험 전면 재실행.
대안: 실험 후 기준 조정(관행) — 기각. 자율 진행에서는 사실상 기준 소멸과 같음.

[2026-09-16] D-08 확정(잠정 해제): 게이트 (0) 조건부 통과를 사용자 권한으로 확정. 조건(저 SNR 동작점 <= 0 dB, N/r >= 8, 짧은 블록)은 그대로 유지.
[2026-09-16] Q-05 확정: 주 레짐 (16,2,32,2)류(수신 32, N/r = 16, 동작점 -4 ~ -8 dB, 짧은 블록), 보조 (8,2,16,2). Q-05 닫힘.
[2026-09-16] D-11 확정(조건부): 동역학 조건 "G_b 블록 간 탈상관(alpha <~ 0.5), S_b 수십 블록 유지"를 주 레짐 정의에 포함. **Q-12(ray-tracing 검증)에서 확인되지 않으면 자동 철회**하고 alpha ~ 1 레짐으로 되돌린다(이때 Q-11의 경쟁자 2종이 (i)의 베이스라인에 들어옴).
이유: 세 건 모두 교수님 검토를 전제로 "잠정"이었으나, 검토가 없을 예정이므로 상태를 명확히 한다. 미결로 남겨 두면 다음 세션이 "확인 대기"로 읽고 멈춘다.

[2026-09-17] D-13 결정: energy-based noncoherent 검출(Manolakos–Chowdhury–Goldsmith 2016 [V, 프로젝트 PDF])은 게이트 (i)의 베이스라인에 넣지 않는다. D-12(pre-registration)의 베이스라인 집합(per-block struct-ML, tracker+조건부 GLRT 중 우수한 쪽)은 그대로이며, D-13은 그 집합에 항목을 추가하지 않는 결정이므로 재실행 대상 실험은 없다. 관련 연구 지도에서 "statistics-only noncoherent, channel-hardening 기반, r≈N" 대조군으로만 위치시키고, OOK genie 간극(2.95/4.15/5.41 dB, N=16/32/64)은 Q-04(peaky 입력)의 근거로 쓴다. 리뷰어가 저복잡도 envelope 수신기 비교를 요구하면 MCG 블록 평균이 아니라 MF-energy prior-free ML(T6_energy_based_note.md §2 폐형)로 보고한다.
이유: (1) 입력 부류가 다름(진폭 L점 vs Grassmann K점; T=16에서 energy는 OOK 1/16 b/cu뿐). (2) 자기 입력에 대해서도 prior-free 블록 ML에 5.4 dB 열세(N=32,T=16) [정확]. (3) 저랭크에서 hardening이 대수적으로 없음(Var 신호항 p^2/r) → L>=3은 N,T,SNR 무관 오류 바닥 [정확]. MCG는 우리 모델의 r=N, M=1, T=1 모서리이며 T6 주 레짐(r=2, N=32, −4~−8 dB)과 겹치지 않음.
대안: OOK energy 검출을 (i)의 보조 베이스라인으로 보고 — 기각(rate 불일치, 이미 struct-ML이 지배). Rician/LOS 채널에서 재검토 — Q-12 결과로 레짐의 LOS 유무가 정해진 뒤 필요 시(Q-01 메모).
정정 항목(문서): 동반 노트 §2 Phase 0 항목 5와 §6 15번의 "각도 sparse 채널의 noncoherent/massive SIMO 검출(energy-based, MCG 2016)"은 부정확 — MCG 2016은 안테나 간 i.i.d. 채널을 가정하고 antenna correlation을 향후 과제로 명시(§II, §VIII). "안테나 간 i.i.d. 채널의 energy-based noncoherent massive SIMO 검출"로 고칠 것. 태그 갱신: HS 2000, MCG 2016 [M]/[web 사본] → [V, 프로젝트 PDF].

[2026-09-17] D-14 제안(결정 대기, 사용자 확정 필요): prior의 대상을 순간 부분공간 S_b = span(H_b^H) ∈ G(N,2)(D-01)에서 장기 수신 부분공간 S_b^LT ∈ G(N, r_eff)(창 공분산 sum_b H_b^H H_b의 상위 r_eff 고유공간; 등가로 저랭크 수신 공분산)로 바꾼다. 검출기는 tracked/sampled U^LT로 사영한 뒤 (T,M,r_eff,r) noncoherent 구조 인식 ML/GLRT(1차 세션 prior-free 수신기를 N'=r_eff에 적용; JG 2005 공분산 정보 수신기 구조). r_eff 기본값 주 4·보조 3(폐형 제약상 4로 계산), Q-12 (b) 실측으로 확정. D-11 문구를 "G_b 탈상관, S_b^LT 수십 블록 유지, 순간 2-평면은 그 안에서 블록마다 무작위"로 수정.
이유: T6_q12_geom_note.md — 국소 산란 기하 모델에서 순간 2-평면은 이득과 함께 탈상관(rho_S(1) ≈ 0.6 at N=32)하고 장기 부분공간만 수십 블록 유지. exact-U genie 상한(3.6–3.9 dB)은 어떤 시간 prior도 달성 불가; LT genie 상한 2.7 dB(r_eff=4)는 통과선 위. 생성 prior 서사(드리프트·blockage 점프·경로 생성/소멸)는 장기 부분공간에 더 자연스러움.
대안: (A) D-01 유지 + rank-2 tracked 부분공간(사영 손실 ≈ 1.1 dB, 상한 ≈ 2.5 dB) — 단순하나 손실이 (i) 여유를 먹음. (B) D-01 유지 + 레짐을 정반사 지배(PEF_2 ≥ 0.9)로 좁힘 — Q-12 (b)가 그런 채널을 보여줘야 하고 레짐이 좁아짐. 기본값: D-14.
pre-registration 정합: 게이트 기준·통과선·베이스라인 집합은 불변. 바뀌는 것은 prior 객체와 D-11 문구이며, 영향 실험(게이트 (0) 상한)은 t6_ceiling_LT.json으로 재실행 완료. 채택 시 T6_ceiling_derivation.md §0–2(정의·JG 대응)와 동반 노트 §1.2 서술을 갱신한다.

[2026-09-17] D-14 채택(기본값 승인): prior의 대상 = 장기 수신 부분공간 S_b^LT ∈ G(N, r_eff). r_eff 기본값은 RT 실측(NLOS 8–11, LOS 2–4)을 반영해 "주 레짐 4–8, 실험 격자에 {4, 8} 포함"으로 둔다. T6_ceiling_derivation.md §0–2·동반 노트 §1.2 서술 갱신은 다음 세션 문서 작업.
[2026-09-17] D-15 채택(기본값 승인, 사전 고정 규칙의 결과): D-11 철회. 주 레짐(Q-05)의 동역학 조건을 "블록 간격 2–10 ms(≤2.8λ)에서 이득 복소 상관 |alpha| ≈ 0.6–0.95(정지·정반사 지배 환경), 장기 수신 부분공간 r_eff ≈ 3–10, 순간 2-평면 비지속"으로 바꾼다. 게이트 (i) 베이스라인 집합에 Q-11 경쟁자 2종 — (a) 2T-블록 USTM(채널 2블록 상수 가정), (b) full-channel 예측(Kalman/AR on H_b) + coherent 또는 오차 인지 검출 — 을 추가한다(D-12의 "최고 베이스라인 = 그 설정에서 더 좋은 쪽" 정의는 유지; 집합 확장을 여기 기록). 통과선·(ii)(iii) 기준 불변.
이유: T6_q12_rt_note.md — Sionna RT 22개 궤적 중 rho_G(1) ≤ 0.5인 것 0개; D-11의 근거였던 기하 모델은 등방·동적 국소 산란을 가정한 것이고 RT 장면에는 그 요소가 없음. 규칙(2026-09-17 2차: "아니면 철회·alpha≈1 레짐 복귀·Q-11 진입")대로 처리.
대안: D-11을 "동적 클러터 환경" 조건부로 유지 — 외부 측정 근거 없이는 기각(규칙 위반). 사람·차량 이동을 RT에 추가해 재검증 — 가능하나(Sionna RT SceneObject.velocity) 이번 세션 범위 밖; Q-12 (c)로 남김.
영향 실험(재실행 필요 없음): (0) LT-genie 상한은 동역학과 무관. (i)는 아직 실험 없음. 차동 USTM 노트의 alpha 표는 그대로 유효(이제 주 레짐 alpha 0.6–0.95가 그 표의 "diff 0.9"·"diff ≤0.5" 열에 해당 — 차동은 여전히 열세).
[2026-09-18] 정리(결정 아님): 프로젝트 파일은 핵심(문서·기록 4종·PDF·import되는 코드·요약 txt)만 남기고, 결과 JSON·그림·보조 스크립트·arXiv 원시 덤프·세션 핸드오프는 GitHub로 이관. 기록 파일 정본은 이 4개(gate_evidence_T6 / decision_log_T6 / open_questions_T6 / scooping_log_T6)이며 구판(_T2, 무접미, __1_~__5_)은 삭제. 이관 목록은 2026-09-18 채팅 참조.
```
