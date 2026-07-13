# EXHALE 코드 검토 보고서 (한국어)

- 검토일: 2026-07-13
- 검토 기준 커밋: `cfc33b1fe36143e6192f2d854805e506ddcda966` (`update`)
- 범위: 주 실행 코드, Wind-AE 초기조건 코드, 정상상태 풀이기, 분자·금속 화학, 후처리 Python, 빌드·실행 스크립트
- 방법: 소스 정적 검토, 경고 강화 빌드, Python/셸 구문 검사, 정상 입력 스모크 테스트, 비정상 조합 재현

## 1. 요약

현재 코드는 경고 강화 디버그 빌드에 성공했고, `tutorial_nometals` 입력으로 474 step까지 bounds/FPE 검사 없이 진행되었다. 기존 검토에서 수정한 보간 NaN, 후처리 He/H 비율, `q_abs`/로그 보호도 현재 코드에 남아 있다. 특히 방사선 루프에서 과거에 쓰던 OpenMP `critical`은 제거되어 병렬 병목 하나가 이미 개선되었다.

다만 다음 항목은 우선 수정하는 것이 좋다.

| 우선순위 | 문제 | 성격 | 핵심 영향 |
|---|---|---|---|
| P0 | 분자 화학 호환성 검사가 옵션을 읽기 전에 실행됨 | 실행으로 재현 | 금속+분자 입력이 검증을 통과한 뒤 배열 경계 초과로 종료 |
| P0 | 동일 배열을 `intent(in)`과 `intent(out)` 인수에 동시에 전달 | 언어 규격상 오류 위험 | 최적화 수준·컴파일러에 따라 결과가 달라질 수 있음 |
| P1 | loaded SED의 하한이 low-IP 금속 임계값을 반영하지 않음 | 높은 확신의 물리/로직 오류 | Mg, Si, Ca, Na, K, Fe의 광이온화율을 과소평가할 수 있음 |
| P1 | 단색광 H-only 분기에서 He 질량과 EOS 처리가 불일치 | 높은 확신의 상태 불일치 | 초기 밀도·입자수·화학 상태가 서로 다른 조성을 가정 |
| P1 | transit 자동 파장창의 column 계산에 cm 변환 누락 | 확정 산술 오류 | column을 100배 작게 계산해 wing 파장창을 지나치게 좁힐 수 있음 |
| P1 | SED reader의 EOF/최소 행/정렬 검증 부족 | 높은 확신의 입력 오류 | 무한 반복 또는 `e_v(2)` 경계 초과 가능 |
| P1 | HLLC hot loop에서 비연속 row slice 임시 배열 생성 | 실행으로 확인 | 매 cell/interface마다 복사·임시 배열 비용 발생 |
| P2 | 문자열 선택값에 `case default`가 없음 | 높은 확신의 입력 오류 | 오타가 즉시 중단되지 않고 미정의 상태로 전파 |
| P2 | 빌드 stamp가 `FFLAGS` 변경을 추적하지 않음 | 명령으로 확인 | debug/release flag를 바꿔도 이전 object를 재사용 |
| P2 | 상수 literal 정밀도와 implicit interface 경고가 많음 | 빌드로 확인 | 수치 정밀도 저하 및 호출 불일치 탐지 실패 |

권장 순서는 P0 두 건, P1의 물리·단위 오류, SED 방어 로직, hot-loop profile/개선, 빌드·구조 정리 순이다.

## 2. 검토 및 검증 결과

### 2.1 경고 강화 빌드

다음 성격의 옵션으로 전체 Fortran 코드를 별도 `/tmp` object 디렉터리에 강제 재빌드했다.

```text
-O0 -g -fopenmp -Wall -Wextra -Wimplicit-interface
-Wconversion-extra -Wsurprising -fcheck=all -fbacktrace
-ffpe-trap=invalid,zero,overflow
```

빌드는 성공했다. 경고는 총 4,440건이며 주요 분류는 다음과 같다.

| 경고 | 개수 | 해석 |
|---|---:|---|
| 탭 문자 | 2,779 | 대부분 형식 문제이나 유의미한 경고를 가림 |
| 정밀도 변환 | 1,417 | 기본 실수 literal을 `real*8`에 넣는 경우가 다수 |
| implicit interface | 51 | MINPACK/LAPACK 호출의 형식·rank 검사가 불완전 |
| 실수 동등 비교 | 51 | 일부는 의도적 sentinel 비교, 일부는 정리 권장 |
| 미사용 변수 | 40 | 오래된 경로와 중복 구현 정리 후보 |
| do-subscript | 26 | 확인한 주요 항목은 endpoint 분기로 보호된 false positive |

Wind-AE standalone build도 성공했다. Python 파일은 모두 `py_compile`, 추적 중인 셸 스크립트는 `bash -n`을 통과했다.

### 2.2 실행 검증

`examples/tutorial_nometals/input.inp`를 디버그 executable로 실행하여 20초 동안 474 step까지 진행했다. 경계 검사나 FPE 오류는 없었다. 다만 런타임이 다음 위치에서 비연속 배열 section을 위한 임시 배열을 반복 보고했다.

- `src/modules/states/PLM_rec.f90:26`
- `src/modules/time_step/RK_rhs.f90:41,49,53`

분자 화학과 금속을 동시에 활성화한 입력은 parser가 거부하지 않았고, 시간 적분 시작 후 아래 오류로 종료되었다.

```text
Fortran runtime error: Index '9' of dimension 1 of array 'sys_x'
above upper bound of 8
at src/modules/radiation/ionization_equilibrium.f90:409
```

따라서 아래 3.1은 잠재적 위험이 아니라 재현된 결함이다.

## 3. 정확성 및 안정성 문제

### 3.1 [P0, 재현됨] 분자 화학 호환성 검사의 실행 위치가 잘못됨

근거:

- `src/modules/files_IO/input_read.f90:262-276`에서 분자 화학과 He, He diffusion, metals의 조합을 검사한다.
- 그러나 `Molecular chemistry`는 같은 파일 `:394-397`, `He diffusion`은 `:428-432`에서 뒤늦게 읽는다.
- 즉 검사 시점에는 두 옵션이 기본값인 `False`여서 제한 조건이 작동하지 않는다.
- 실제 분자+금속 입력은 `ionization_equilibrium.f90:409`에서 `sys_x(9)`를 접근하다 상한 8을 넘었다.

영향:

- 문서상 금지된 `molecular + no He`, `molecular + He diffusion`, `molecular + metals` 조합이 parser를 통과한다.
- 실패 지점이 입력 단계가 아니라 시간 적분 내부이므로 원인 파악도 어렵다.

권장 수정:

1. 모든 선택 옵션과 `base.inp`를 읽은 뒤, allocation/초기화 전에 단일 `validate_configuration()`을 호출한다.
2. 금지된 조합은 설명을 포함한 `error stop`으로 비영(非零) 종료한다.
3. 검사를 여러 parsing 지점에 흩뜨리지 말고 한 함수에 모은다.

필수 회귀 테스트:

- molecular+metals, molecular+He diffusion, molecular+HeH=0이 모두 입력 단계에서 실패해야 한다.
- 허용되는 molecular 예제는 기존과 동일하게 시작되어야 한다.

### 3.2 [P0] `intent(in)`/`intent(out)` 인수 aliasing

대표 사례:

- `Apply_BC(u,u)`가 `EXHALE_main.f90`, `steady_newton.f90`, `init.f90`에서 반복된다.
- `Apply_BC` 자체도 `Apply_BC_W(W,W)`를 호출한다.
- `ioniz_eq(T,rho,f_sp,rho,f_sp,...)`는 입력 `rho/f_sp`와 출력 `rho_out/f_sp_out`에 같은 actual argument를 준다.
- `post_process.f90`에서는 필요 없는 여러 `intent(out)` 값에 같은 scratch 변수 `dum_v`를 동시에 전달한다.

Fortran에서 한 procedure 호출 중 한 dummy가 정의되는 동안 다른 dummy를 통해 같은 actual object를 참조하는 것은 제한된다. 현재 구현이 우연히 입력을 local로 먼저 복사하더라도 interface 계약 자체가 비정상이며, compiler 최적화가 alias가 없다고 가정할 수 있다.

권장 수정:

- 실제 in-place 동작인 boundary routine은 명시적 `intent(inout)` API 하나로 만든다.
- 입력과 출력을 모두 지원해야 하면 별도 buffer를 쓰는 wrapper와 in-place wrapper를 분리한다.
- `ioniz_eq` 호출부는 `rho_new`, `f_sp_new`를 사용한 뒤 명시적으로 대입한다.
- 버리는 출력은 각각 다른 local을 쓰거나, optional output/dedicated rate API를 제공한다.
- gfortran `-O0/-O3`와 Intel ifx 중 가능한 두 compiler에서 golden regression을 수행한다.

### 3.3 [P1] loaded SED가 low-ionization-potential 금속의 광자를 잘라냄

`src/modules/radiation/set_energy_vectors.f90:31-53`의 power-law 경로는 `thereis_lowIP_metal`일 때 에너지 grid 하한을 활성 금속의 ionization threshold까지 낮춘다. 반면 `src/modules/radiation/sed_read.f90:32-33`의 loaded SED 경로는 He I triplet만 고려하고 low-IP 금속은 고려하지 않는다.

영향:

- SED 파일에 13.6 eV 아래 광자가 존재해도 선택 단계에서 제외될 수 있다.
- Mg, Si, Ca, Na, K, Fe 등의 광이온화율과 전자 밀도, cooling, transit ion fraction이 연쇄적으로 편향될 수 있다.

권장 수정:

- 활성 species metadata(`species_table`)에서 최소 ionization threshold를 계산해 power-law와 loaded SED가 같은 하한 정책을 쓰게 한다.
- SED 자체가 그 범위를 포함하지 않으면 경고 또는 명시적 오류를 낸다.
- low-IP metal 각각에 대해 threshold 양쪽에 line/continuum bin을 둔 rate regression을 추가한다.

### 3.4 [P1] 단색광 H-only 분기의 helium 상태 불일치

`input_read.f90:174-176`은 monochromatic energy가 24.6 eV보다 작으면 `thereis_He=.False.`로 바꾼다. 그러나 `HeH`, `mass_per_H=1+4*HeH`, boundary density에는 helium 질량이 남으며, `set_IC.f90:151-156`도 helium fraction을 기록한다. 반대로 `utilities.f90`의 `calc_rho`와 `calc_ntot`은 `thereis_He=False`일 때 helium을 완전히 제외한다.

이 상태는 “helium은 중성이지만 질량은 존재한다”와 “helium 자체가 없는 H-only 유체”를 혼합한다.

권장 수정:

- 물리적으로는 24.6 eV 미만이어도 neutral He의 질량·입자수는 유지하고, photoionization channel만 비활성화하는 편이 자연스럽다.
- 정말 H-only 계산을 의도했다면 `HeH=0`, `mass_per_H`, 초기 abundance와 EOS를 모두 같은 모드로 바꾼다.
- `composition_has_He`와 `He_photoactive`를 별도 논리값으로 나누는 것이 가장 명확하다.

### 3.5 [P1] SED reader의 EOF 및 입력 구조 검증 부족

`sed_read.f90:46`에서 `iostat`과 함께 읽은 뒤, `:49`에서 값을 사용하고 `:53`의 `cycle`을 거쳐서야 `:59`에서 `iostat`을 검사한다. 모든 행이 선택 상한 밖에 있으면 EOF에서 이전/미정 값으로 계속 `cycle`할 가능성이 있다. 또한 선택된 행이 0개 또는 1개여도 `:90`에서 `e_v(2)`를 참조한다.

추가로 wavelength의 양수성, 엄격한 정렬, 중복 행, malformed row에 대한 검증이 없다.

권장 수정:

- `read` 직후 `iostat`을 먼저 처리하고 EOF와 malformed row를 구분한다.
- 선택 결과가 최소 2개인지 검사한다.
- wavelength가 양수이고 엄격히 증가하는지 검증한다.
- 실패 메시지에 파일명과 행 번호, 필요한 에너지 범위를 포함한다.

### 3.6 [P1, 확정] transit 자동 파장창 column의 단위가 100배 틀림

`EXHALE_transit.py:406-409`에서 `_dr_cm = np.gradient(r) * Rp`로 계산한다. 이때 `Rp`는 meter 단위인데 density는 cm\(^{-3}\)이다. 이름과 주석은 cm를 가정하지만 `1e2` 변환이 빠졌다. 같은 파일의 실제 LOS 적분 경로 `:472-475`는 올바르게 `r * Rp * 1e2`를 사용한다.

영향:

- 자동 파장창용 column density가 100배 작다.
- damping-wing 추정 폭은 대략 \(\sqrt{N}\)에 비례하므로 cap에 걸리기 전 최대 약 10배 좁은 창을 선택할 수 있다.
- 실제 optical-depth 적분보다 주로 He/Lyα 자동 파장 범위 선택에 영향을 주며, 날개가 잘릴 수 있다.

수정은 `_dr_cm = np.abs(np.gradient(r)) * Rp * 1e2`로 하고, 균일한 analytic atmosphere로 단위 테스트한다.

### 3.7 [P2] enum/string 입력에 기본 오류 분기가 없음

다음 `select case`에는 `case default`가 없다.

- reconstruction: `src/modules/states/Reconstruction.f90`
- numerical flux: `src/modules/states/Num_Fluxes.f90`
- grid type: `src/modules/grid/define_grid.f90`
- spectrum type: `src/modules/files_IO/input_read.f90:151-178`

오타가 난 입력은 즉시 실패하지 않고 flag나 출력 배열을 정의하지 않은 채 진행할 수 있다.

권장 수정은 모든 user-facing 선택값을 읽는 즉시 허용 목록과 비교하고, 원래 값 및 허용 값을 보여주는 `case default; error stop`을 두는 것이다.

### 3.8 [P2] `base.inp` override 후 조성 상태를 다시 확정해야 함

lower-atmosphere 경로에서 `base.inp`가 `HeH`를 바꿀 수 있다. 최종 조성 flag와 `mass_per_H`, molecular compatibility validation은 override까지 끝난 값을 기준으로 한 번만 재계산해야 한다. 현재처럼 parsing 도중 상태를 부분 갱신하면 `HeH`, `thereis_He`, species count가 불일치하기 쉽다.

## 4. 성능 최적화 후보

### 4.1 [우선] HLLC/PLM hot loop의 비연속 row slice 제거

배열이 `(cell, component)` 순서인데 `W(j,:)`, `u(k,:)`를 procedure에 전달한다. Fortran column-major layout에서 이 section은 비연속이며, 디버그 런타임이 `PLM_rec.f90:26`, `RK_rhs.f90:41,49,53`마다 임시 배열 생성을 실제 보고했다.

권장 접근:

1. 먼저 production 크기에서 profiler로 임시 배열과 flux/reconstruction의 비중을 측정한다.
2. 작은 변경으로는 scalar 인수 또는 명시적인 길이-3 local을 사용해 숨은 allocation/copy를 제거한다.
3. 장기적으로는 primitive/conserved state를 `(component, cell)`로 전치하거나 SoA 구조로 바꾸고 전 구간 conversion을 한 번에 수행한다.
4. 변경 전후 mass/energy conservation과 wall time을 함께 비교한다.

세 번째 방법은 영향 범위가 크므로 P0/P1 수정 후 독립 PR로 하는 것이 안전하다.

### 4.2 `EXHALE_transit.py`의 중첩 loop와 중복 line 계산

현재 스크립트는 impact parameter × wavelength × LOS의 Python loop를 수행하며, 일부 doublet는 single-component 계산을 먼저 한 뒤 유사 계산을 다시 한다. 또한 `np.append`를 반복해 좌표를 만든다.

권장 개선:

- chord geometry, temperature, bulk velocity, abundance를 impact parameter당 한 번 계산한다.
- wavelength 축은 memory 상한을 둔 chunk vectorization을 사용한다.
- `np.append` loop는 사전 크기 할당 또는 vector 식으로 바꾼다.
- line 정보를 표 구조로 만들고 singlet/doublet을 하나의 solver로 처리한다.
- 최종 스펙트럼의 상대 오차와 peak memory를 함께 기록한다.

### 4.3 화학·cooling rate의 중복 평가

H/He, triplet, metal, molecular system이 같은 cell state에서 유사한 rate와 electron density를 반복 계산한다. 단, 이 부분은 먼저 profiler로 확인해야 한다. 개선 시에는 한 cell의 named rate cache/derived type을 만들고, Newton residual과 Jacobian이 같은 계산 결과를 공유하도록 한다. 캐시 수명과 온도·밀도 의존성을 명확히 하지 않으면 오히려 stale value 오류가 생길 수 있다.

### 4.4 OpenMP 상태

현재 radiation loop의 과거 `critical` 직렬화는 제거되어 있다. 다음 최적화는 무작정 directive를 추가하기보다 thread scaling(1, 2, 4, 8 threads), scheduling imbalance, chemistry solver 비용을 측정한 뒤 결정하는 것이 좋다. CI에는 1-thread와 multi-thread 결과 허용오차 비교를 추가한다.

## 5. 중복 및 간략화 후보

### 5.1 HLLC와 ROE wave-speed 코드 공통화

`speed_estimate_HLLC.f90`와 `speed_estimate_ROE.f90`에는 primitive 분해, PVRS 추정, TRRS/TSRS 분기, star pressure/velocity 계산이 상당 부분 반복된다.

공통 routine이 다음을 반환하도록 분리할 수 있다.

- left/right primitive state
- `p_star`, `u_star`
- shock/rarefaction 보정값

각 solver는 마지막 wave-speed 조합만 유지한다. 먼저 representative Riemann states에 대한 golden test를 만든 뒤 refactor해야 수치 분기 변화를 막을 수 있다.

### 5.2 H/He/TR/metal/molecular system의 중복과 `params(N)` 제거

`System_HeH*`, `System_HeH_TR*`, metal/implicit-advection 변형은 H/He 반응식을 반복하며, 위치 의존적인 `params(N)` 배열을 사용한다. 이 방식은 인덱스 하나만 어긋나도 컴파일러가 잡지 못한다.

권장 구조:

- `type(ion_cell_state)`와 `type(ion_rates)`에 이름 있는 field를 둔다.
- species별 기여를 작은 pure routine으로 분리한다.
- residual/Jacobian assembly는 활성 species metadata를 순회한다.
- 한 번에 모두 합치기보다 H/He 공통부 → metal → molecule 순으로 golden test와 함께 이동한다.

Triplet 경로에서 collisional ionization을 생략하는 현재 선택은 코드 주석에 의도라고 적혀 있다. 이를 즉시 오류로 분류하기보다, 적용 온도 범위에서 영향이 무시 가능한지를 수치 실험으로 확인하고 문서화할 필요가 있다.

### 5.3 composition 계산의 단일화

`calc_rho`, `calc_ntot`, 초기화, boundary, chemistry가 helium/metal/molecule 포함 여부를 각각 분기한다. 3.4와 같은 불일치는 이 분산된 로직에서 생긴다. `species_table`을 H/He/molecule까지 확장하고, 질량·입자수·전자수 계산을 metadata 기반 routine으로 모으는 것이 좋다.

### 5.4 `post_process` 출력 API 단순화

필요 없는 출력 때문에 같은 `dum_v`를 여러 `intent(out)` 인수에 전달한다. 계산 목적별 API(heating only, cooling components, diagnostic rates)를 제공하거나 optional output을 쓰면 aliasing도 제거되고 호출부도 읽기 쉬워진다.

### 5.5 `EXHALE_transit.py` 모듈화

약 1,400줄의 top-level script에 입력 parsing, 물리 계산, convolution, plotting, interactive input이 섞여 있다.

권장 분리:

- `main()`과 argparse/config
- EXHALE profile/header reader
- line metadata table
- pure Voigt/LOS/transit solver
- instrument convolution
- plotting

또한 metal column을 `17, 23, 25`처럼 고정하지 말고 출력의 `# columns` header를 읽어 이름으로 선택해야 한다. species 순서가 바뀌면 현재 코드는 조용히 잘못된 column을 읽을 수 있다.

### 5.6 입력 parser 통합

Fortran parser는 앞부분의 위치 기반 필수 행과 뒤쪽의 key 검색형 옵션을 혼합한다. Python utilities도 같은 파일을 별도로 해석한다. 한 개의 key-value schema와 validation 규칙을 정하고, legacy name alias를 명시적으로 지원하면 중복과 순서 의존성을 줄일 수 있다.

## 6. 빌드, 이식성, 유지보수

### 6.1 [확인됨] `FFLAGS` 변경이 재빌드를 유발하지 않음

Makefile의 compiler stamp는 compiler 이름은 기록하지만 flag 전체를 dependency로 추적하지 않는다. 디버그 flag로 object를 만든 뒤 같은 `OBJDIR`에서 `FFLAGS='-O3 -fopenmp' make -n`을 실행했더니 `Nothing to be done for 'all'`이 나왔다.

권장 수정:

- compiler, `FFLAGS`, module flag, preprocessor option을 내용으로 하는 configuration stamp를 생성한다.
- link library 변경도 relink dependency에 포함한다.
- `debug`, `release`, `check`가 서로 다른 object directory를 기본 사용하게 한다.

### 6.2 실수 상수 정밀도

`src/parameters.f90:230-265` 등은 `real*8` 변수에 기본 정밀도 literal을 대입한다. literal이 먼저 real32로 반올림된 뒤 real64로 승격될 수 있으며, 1,417개의 conversion warning 대부분과 연결된다.

`iso_fortran_env, only: real64`와 `real(real64)`, `1.0_real64` 형식을 사용한다. 결과의 마지막 몇 bit가 바뀔 수 있으므로 물리량 허용오차 기반 regression과 함께 점진적으로 바꾼다.

### 6.3 implicit interface와 선형대수 ABI

MINPACK 및 `dgbtrf/dgbtrs` 호출에서 51개의 implicit-interface 경고가 발생한다. 명시적 module/interface를 제공하면 rank, kind, intent 오류를 컴파일 시 잡을 수 있다.

현재 환경에서는 실행 파일이 `libgfortran.so.5`와 system LAPACK이 요구하는 `libgfortran.so.4`를 동시에 load한다는 linker 경고도 확인했다. 당장 빌드는 성공했지만 ABI/exception/runtime state 위험이 있으므로, production 환경은 주 compiler와 같은 runtime으로 빌드한 BLAS/LAPACK을 쓰는 것이 안전하다.

### 6.4 hard-coded 경로와 실행 위치

`input_read.f90:818-831`에는 lower-atmosphere helper의 개발자 절대 경로 fallback이 있다. 다른 clone에서는 `EXHALE_ROOT`를 설정하지 않으면 실패하며, shell command 인수가 quote되지 않아 공백이 있는 경로도 깨진다.

`run_EXHALE.sh`는 script 위치가 아니라 `pwd`를 project root로 간주한다. `${BASH_SOURCE[0]}`에서 script directory를 구하고 그곳으로 이동하는 편이 견고하다.

### 6.5 Python 의존성 문서 불일치

`EXHALE_transit.py:4`는 `astropy.convolution`을 무조건 import하지만 README requirements에는 numpy, scipy, matplotlib, tkinter만 있다. `astropy`를 requirements에 추가하거나, 사용 기능을 SciPy convolution으로 대체해 의존성을 제거해야 한다.

### 6.6 오류 종료 방식

일부 입력 오류는 일반 `stop`을 사용해 shell에서 성공 status로 보일 수 있다. 모든 fatal configuration/I/O 오류는 공통 error routine 또는 `error stop`으로 비영 종료하고, filename/key/value를 메시지에 포함하는 것이 batch workflow에 적합하다.

## 7. 테스트 전략 제안

현재 Wind-AE에 standalone test program이 일부 있지만 최상위 `make check`나 CI가 없고, benchmark input/reference가 자동 assertion으로 연결되어 있지 않다. 다음의 작은 test pyramid를 권장한다.

1. 단위 테스트
   - EOS/composition: H-only, H+He, triplet, molecule, metals
   - SED selection: low-IP threshold, EOF, malformed/1-row/non-monotonic file
   - reconstruction/flux: constant state, shock tube, positivity
   - transit: 단위가 알려진 uniform slab와 analytic column
2. 입력 검증 테스트
   - 금지 조합, 잘못된 enum, 누락 파일이 모두 빠르고 비영 status로 실패
3. 짧은 통합 테스트
   - tutorial_nometals, one metal, molecular, lower-atmosphere 각각 제한 step 실행
   - mass/energy와 abundance 범위, NaN/Inf 부재를 assertion
4. compiler matrix
   - gfortran debug: `-fcheck=all`, FPE trap
   - gfortran optimized
   - 가능하면 ifx 한 구성으로 aliasing/표준 위반 탐지
5. 성능 회귀
   - 고정 thread 수와 입력으로 step/s, radiation/chemistry/flux 시간을 기록
   - 결과 허용오차를 통과한 경우에만 최적화 성능을 비교

초기에는 4,440개 경고를 모두 CI 실패로 만들기보다 `strict-check` target에서 분류하고, 새 경고 금지 → 기존 conversion/interface 경고 감축 순으로 운영하는 편이 현실적이다.

## 8. 권장 실행 로드맵

### 1단계: 즉시 안정화

- 최종 `validate_configuration()`을 도입하고 분자 조합 crash 회귀 테스트 추가
- `Apply_BC`, `ioniz_eq`, post-process aliasing 제거
- transit cm 변환 수정 및 analytic test 추가
- SED EOF/최소 행/정렬 방어 로직 추가

### 2단계: 물리 일관성

- loaded SED와 power-law의 low-IP threshold 정책 통합
- composition 존재 여부와 photoactive 여부 분리
- 모든 enum에 fail-fast validation 추가

### 3단계: 빌드와 수치 위생

- configuration stamp 수정
- `real64` 상수로 점진적 전환
- MINPACK/LAPACK 명시적 interface와 일관된 runtime 구성
- `make check`와 CI 도입

### 4단계: 성능·구조 개선

- hot-loop array temporary를 profile하고 layout/API 개선
- wave-speed 및 chemistry 공통부를 golden test 아래에서 refactor
- transit script를 line-table 기반 모듈로 분리하고 vectorize

## 9. 범위와 한계

이 검토는 전체 production 규모의 수렴 결과나 관측 스펙트럼을 기준 자료와 장시간 비교한 물리 검증은 아니다. 정상 예제는 20초 스모크 실행까지만 확인했고, compiler warning 중 endpoint 분기로 보호된 `do-subscript` 항목은 false positive로 판정했다. 성능 항목은 runtime temporary가 확인된 부분을 제외하면 profiler로 우선순위를 재확인해야 한다.

그럼에도 3.1의 배열 경계 초과, 3.6의 단위 누락, 6.1의 flag 재빌드 실패는 직접 재현하거나 산술적으로 확정했으며, 최우선 수정 대상으로 볼 수 있다.
