# open_questions.md — T6

형식: `Q-번호 (상태, Phase): 질문. 담당/근거. 기본값.` 닫힐 때 날짜와 결론을 붙인다.

```
Q-01 (열림, Phase 1): 물리 채널에서 G_b는 iid Gaussian이 아니라 steering 구조(경로 이득 × 송신 응답)를 가진다. A3의 iid G_b 가정이 간극을 과대/과소평가하는가?
근거: 이번 계산은 A3 기준. 기본값: ray-tracing 채널로 Phase 2에서 재계산.
추가 2026-09-17 (메모): Rician/LOS 성분이 크면 ||G_b||^2 요동이 줄어 energy-based 검출(MCG 2016)이 일부 회복될 수 있음 [추측]. Q-12의 ray-tracing 궤적이 LOS를 포함하면 T6_energy_based_note.md §2 분해식에 LOS 항을 넣어 재확인(D-13에는 영향 없음: 베이스라인 집합은 D-12로 고정).

Q-02 (열림, Phase 1): prior-free 수신기가 U_b의 등방(Haar) 분포만 안다는 A4는 array-manifold 구조(AoA 부분다양체)를 무시한다. array-manifold-aware per-block 수신기(예: 빔공간 sparse 검출)를 넣으면 간극이 얼마나 줄어드는가?
근거: 계산된 간극은 상한 쪽 추정. 기본값: Phase 1에서 DFT-beamspace per-block 수신기를 베이스라인에 추가.
추가 2026-09-17: D-14 채택 시 prior-free의 "등방"은 G(N, r_eff) 위 등방이 되고, array-manifold-aware 수신기는 장기 부분공간이 steering 다양체 근방에 있다는 지식에 해당 — LT genie(부분공간 기지)와 prior-free 사이에 놓임. Q-12 (b) 채널로 DFT-beamspace 수신기의 PEF를 함께 잴 것.
추가 2026-09-16: 참고문헌 arXiv:2604.19904 v2(Khirwadkar–Rajamäki–Pal) — 미지 복소 이득 하 DoA↔beamspace 부분공간의 Grassmann 기하, beamspace Grassmannian code. 우리 A4(등방 U_b)와 array-manifold 구조 사이의 간극을 정량화할 때 인용 가능.

Q-03 (열림, Phase 1): r < M 또는 가변 r. 이번 계산은 r = M = 2 고정.
추가 2026-09-17: 재개. Q-12 (a)에서 순간 rank는 2(=M)이지만 지속 객체의 차원은 r_eff ≈ 3–4(N=32)로 드러남. D-14 채택 시 "r" = r_eff(prior 대상의 차원)와 "r = 2"(블록당 채널 rank)를 구분해 표기해야 함. 기본값: r_eff = 4 (주), Q-12 (b) 실측으로 확정. r < M이면 genie 채널의 pre-log가 min(M, r, ⌊T/2⌋)로 낮아지지만 prior-free도 같이 낮아지므로 간극 성격은 유지될 것으로 추측.
기본값: r = 2 고정 유지(동반 노트 Q-01과 동일 결론).

Q-04 (열림, Phase 0/1): 저 SNR 최적 입력(non-USTM, peaky)에서 간극이 커지는가 작아지는가? USTM 고정 결과는 용량 간극과 다를 수 있음.
근거: [M] 저 SNR noncoherent 최적 입력은 non-isotropic. 기본값: 게이트 (0)은 USTM 기준으로 두고, 판정이 경계선이면 on–off USTM으로 재계산.
추가 2026-09-17: OOK(M=1) 입력에서 부분공간 지식의 값은 정확히 2.95/4.15/5.41 dB (N=16/32/64, r=2; T6_energy_based_note.md 표 B) — USTM genie 간극(1.3–2.1/3.3–3.9 dB)보다 약간 크다. peaky 입력에서 간극이 줄지 않는다는 첫 정량 근거. on–off USTM(K점 중 하나를 0 codeword) 재계산은 (0)이 확정(D-08)된 지금 급하지 않으므로 기본값 유지: 판정에 필요할 때만.

Q-05 (열림, Phase 1, 최우선): 게이트 (0)/(i)의 레짐 고정. 전체 실행 결과 간극(상한)은 SNR ≤ 6 dB(T=8,N=16) / ≤ 2–4 dB(T=16,N=32)에 집중되고, SER 10^{-2} 동작점의 dB 간극은 (8,2,16,2)에서 1.3–2.1 dB, (16,2,32,2)에서 3.3–3.9 dB. 어느 (T, N/r, 동작 SNR) 레짐을 목표로 고정할 것인가?
근거: gate_evidence.md (2026-09-16). 기본값: (16,2,32,2)류(큰 N/r, 동작점 ≤ 0 dB, 짧은 블록)를 주 레짐으로, (8,2,16,2)는 보조. 교수님 확인 필요(동반 노트 §4 "관심 레짐").
추가 2026-09-16 (D-11): 동역학 조건 "G_b 블록 간 탈상관(alpha <~ 0.5), S_b 유지"를 기본값에 포함. 근거: 2T 참조 — alpha≈1이면 N/r=8에서 긴 블록이 T6와 대등. 물리적 타당성은 Q-12.

Q-06 (잠정 닫힘 2026-09-16, Phase 3 재확인): 차동 USTM(Hochwald–Sweldens 2000 [V, web 사본])이 같은 구조 채널에서 얼마나 얻는가? → 사전 추정 결과 "얻지 못한다": 채널 완전 상수여도 per-block struct-ML보다 2-4 dB 열세(T6_diff_ustm_note.md). 최고 베이스라인 후보에서 제외(D-10). 남은 확인: HS group constellation, 구조 인식 pair-ML, U 드리프트 하 거동(Phase 3).

Q-07 (열림, 도구): F_D 폐형은 M = r = 2에만 있다. 일반 (M, r)은 truncated Haar의 det(I − Z Z^H)^{D−2m} 밀도로 같은 방식의 폐형이 가능해 보이나 유도 필요.
기본값: 필요해질 때까지 보류.

Q-08 (열림, 도구): 연속 입력의 prior-free 정확 MI. Φ에 대한 Haar 적분(차원 T)에 importance sampling(top-2 left singular subspace 주변 제안분포)을 쓰면 계산 가능할 수 있음.
근거: 현재는 GMI_svd 하한(간극의 상한)만 있음. 기본값: 상한만으로 판정이 갈리지 않을 때 착수.

Q-09 (열림, Phase 0): γ = N/r 정규화가 결과를 어떻게 바꾸는가? γ = 1이면 genie의 유효 SNR 이득(9/12 dB)이 사라져 간극이 줄어들 것으로 추측.
기본값: γ = N/r 유지(D-02); 검토 필요 시 --gamma_mode 1로 재실행.

Q-10 (열림, Phase 2): 생성 경로 선택 — Chen–Lipman geodesic Riemannian flow matching(동반 노트 결정) vs Cartan flow matching(arXiv:2605.03588 v2: 대칭공간 위 FM을 Lie algebra 부분공간 위 FM으로 환원, geodesic interpolation 불필요) vs homogeneous-space FM(arXiv:2603.24829: G/H → G → Lie algebra). 복소 G(N,r)에 대한 대칭쌍 (U(N), U(r)×U(N−r))의 Cartan 분해 u(N) = (u(r) ⊕ u(N−r)) ⊕ p, p ≅ C^{r×(N−r)}는 표준 [M]이므로 원리상 적용 가능하나 논문은 실 Grassmannian만 시연. 학습 비용과 샘플링 단계 수가 기준 (iii)에 직결.
근거: scooping_log 2026-09-16. 기본값: geodesic RFM을 1순위로 유지(복소 Grassmann exp/log 폐형과 round-trip 시험이 이미 계획됨), Cartan FM은 Phase 2 ablation 후보. 2605.03588은 v2가 "§4 오류 수정"을 명시하므로 인용 전 v2 원문 확인 [VERIFY].

Q-11 (열림, Phase 1/3): alpha≈1(느린 페이딩, 연속 블록)에서는 "블록 길이 2T의 prior-free USTM"(+2.4 dB at N/r=8, +1.6 dB at N/r=16)과 full-channel tracker(Kalman on H_b)가 경쟁자다. Phase 3 베이스라인에 (a) 2T-블록 USTM(채널 상수 가정), (b) full-channel Kalman + coherent 검출을 넣어야 하는가?
근거: t6_twoblock_ref.json. 기본값: 주 레짐(D-11: alpha 작음)에서는 둘 다 무력하므로 보조 레짐(alpha≈1)에서만 (a)를 보고; (b)는 Phase 3에서 tracker 계열과 함께 구현 여부 결정.

Q-12 (열림, Phase 1/2): "이득 G_b는 블록 간 탈상관, receive-side subspace는 수십 블록 유지"인 물리 레짐이 실제로 있는가? 후보: mmWave 소수 경로, 블록 간격 >= 이득 coherence time(예: 슬롯 간격 산발 짧은 패킷, f_D·간격 >~ 0.3), AoA 변화율 << 1/블록 수. 산발 패킷은 이득 상관은 죽이지만 AoA 상관은 남긴다는 주장의 검증.
근거: T6_diff_ustm_note.md §5.3 [추측]. 기본값: Phase 2 ray-tracing(DeepMIMO/Sionna RT) 궤적에서 블록 간 |G 상관|과 subspace chordal distance를 같은 간격으로 재서 확인; 안 나오면 D-11 철회.
추가 2026-09-17 — (a) 기하 모델 사전 점검 완료(T6_q12_geom_note.md): 이득은 Jakes로 탈상관 ✔; 순간 rank-2 부분공간은 함께 탈상관 ✘(rho_S(1) ≈ 0.6, N rho_loc/R로 결정); 장기 부분공간(r_eff ≈ 3–4)은 48블록 유지 ✔ → D-11은 장기 부분공간으로 수정 필요(D-14 제안). (b) ray-tracing 측정 항목과 판정 규칙을 고정: v Delta/lambda ∈ [0.3, 3]에서 (i) rho_G(1) ≤ 0.5, (ii) r_eff95(창 16–32블록), (iii) PEF_2(k)·PEF_reff(k), k ≤ 48, (iv) LT 중첩(k), (v) rho_S(1), (vi) LOS/정반사 비율. 판정: LT 중첩(32) ≥ 0.8 이고 PEF_reff(32) ≥ 0.85 이면 D-11(수정본) 채택, 아니면 철회·alpha≈1 레짐 복귀(Q-11 진입); r_eff > 8이면 (0) 재판정. 도구 후보: Sionna RT(이동 궤적, diffuse scattering 옵션 필수 — 정반사만 켜면 국소 산란의 Doppler 퍼짐이 사라져 이득이 탈상관하지 않음), DeepMIMO(격자 횡단). 기본값: Sionna RT 설치 가능성 시험 → 불가 시 DeepMIMO.
추가 2026-09-17 3차 — (b) 완료(T6_q12_rt_note.md): Sionna RT 2.1.0 설치·실행 성공(1 CPU, 해당 1–4 s). 규칙 적용 결과 이득 탈상관 조건 통과 0/22 → D-11 철회(D-15). 남은 것: (c) 동적 클러터(SceneObject.velocity로 움직이는 차량·보행자 대용 물체) 추가 재검증, 다른 장면(etoile, san_francisco)·BS 높이(가로등 높이 5–10 m, NLOS 비율↑). 기본값: Q-15가 먼저; (c)는 (i) 결과가 경계선일 때.

Q-13 (열림, Phase 3): 모델 기반 tracker 베이스라인(Saad-Falcon 2024 static/constant-velocity RLS, GROUSE, PETRELS)은 점프·경로 생성/소멸을 모델하지 않는다(T6_tracker_assumptions_note.md). 공정 비교를 위해 점프 인지 변형(innovation 검정 + 재초기화)을 tracker와 게이트 (ii)의 접공간 Gaussian AR prior 양쪽에 넣어야 하는가?
근거: 점프에서 쉽게 무너지는 tracker와 비교하면 (i)·(ii)가 쉬운 상대와의 비교가 됨. 기본값: 넣는다 — tracker 2종(무변형/reset), AR prior 2종(무변형/reset)으로 Phase 4 실험 격자에 포함. 주 레짐(D-11)에서 살아남는 prior-free 경쟁자는 tracker+조건부 GLRT뿐이므로 이것이 (i)의 사실상 기준.

[2026-09-16] 상태 정리: 교수님 검토가 없을 예정이므로 "교수님 확인 필요/확인 전" 표기는 모두 무효. 해당 항목은 사용자 권한으로 확정됨(D-08, Q-05, D-11). 게이트 기준은 D-12(pre-registration)로 고정되었으므로, 아래 질문들의 기본값을 바꿀 때도 D-번호가 필요하다.

Q-12 승격(최우선, Phase 1/2): D-11이 조건부 확정이므로 Q-12는 이제 "선택적 검증"이 아니라 **D-11의 존폐를 결정하는 관문**이다. ray-tracing 궤적에서 블록 간 |G 상관|과 subspace chordal distance를 같은 간격으로 측정하여, "이득 탈상관 + 부분공간 유지" 구간이 실제로 존재하는지 확인한다. 없으면 D-11 자동 철회.
```
