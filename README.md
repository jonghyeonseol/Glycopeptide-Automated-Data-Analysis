# pGlyco Auto Combine

자동 glycoproteomics 데이터 통합 및 분석 파이프라인

## 기능

1. **데이터 통합**: Dataset 폴더의 모든 CSV 파일을 자동으로 통합
2. **Glycan Annotation**: Fucosylation 및 Sialylation 여부 자동 주석
3. **통계 분석**: PCA 및 그룹 간 통계 비교
4. **시각화**: PCA plot, Boxplot, Heatmap, 분포도 자동 생성

## 프로젝트 구조

```
pGlyco_auto_combine/
├── config.yaml              # 설정 파일
├── main.py                  # 메인 파이프라인 실행 스크립트
├── requirements.txt         # Python 패키지 의존성
├── src/
│   ├── data_loader.py       # CSV 통합 모듈
│   ├── annotator.py         # Glycan annotation 모듈
│   ├── analyzer.py          # 통계 분석 모듈
│   └── visualizer.py        # 시각화 모듈
├── Dataset/                 # 입력 CSV 파일 (C_01.csv ~ N_24.csv)
└── Results/                 # 출력 결과
    ├── integrated_example.csv       # 통합 및 주석된 데이터
    ├── analysis_summary.txt         # 분석 요약 리포트
    ├── glycan_type_statistics.csv   # Glycan type별 통계
    ├── pca_plot.png                 # PCA 시각화
    ├── pca_samples.png              # PCA 샘플 분포
    ├── boxplot_glycan_types.png     # Boxplot
    ├── glycan_type_distribution.png # Glycan type 분포
    └── heatmap_top_glycopeptides.png # Heatmap
```

## 설치

```bash
pip install -r requirements.txt
```

## 사용 방법

### 1. 데이터 준비
- `Dataset/` 폴더에 CSV 파일들을 배치 (파일명: `C_01.csv`, `C_02.csv`, ..., `N_01.csv`, ...)
- 각 CSV 파일은 다음 열을 포함해야 함:
  - `Peptide`: 펩타이드 서열
  - `GlycanComposition`: Glycan 조성 (예: H(5)N(4)A(2)F(1))
  - `IsotopeArea`: 정량 값

### 2. 파이프라인 실행

```bash
python3 main.py
```

### 3. 결과 확인
- `Results/` 폴더에 모든 출력 파일이 생성됨

## Annotation 규칙

### Sialylation
- GlycanComposition에 `A`가 포함된 경우 → **Sialylated**
- 예: `H(5)N(4)A(2)` → Sialylated (2개의 sialic acid)

### Fucosylation
- GlycanComposition에 `F`가 포함된 경우 → **Fucosylated**
- 예: `H(5)N(4)F(1)` → Fucosylated (1개의 fucose)

### Glycan Type
- `Non`: Sialylation 없음, Fucosylation 없음
- `Sialylated`: Sialylation만 있음
- `Fucosylated`: Fucosylation만 있음
- `Both`: Sialylation과 Fucosylation 둘 다 있음

## 출력 파일 설명

### integrated_example.csv
- 모든 샘플의 데이터를 통합한 파일
- 구조:
  ```
  Peptide | GlycanComposition | C1 | C2 | ... | N1 | N2 | ... | Sialylation | Fucosylation | GlycanType
  ```

### analysis_summary.txt
- 전체 분석 결과 요약
- 샘플 수, Glycan type 분포, PCA 결과, 통계 등

### 시각화 파일
- **pca_plot.png**: Cancer vs Normal 샘플의 PCA 분포
- **boxplot_glycan_types.png**: Glycan type별 강도 분포
- **glycan_type_distribution.png**: Glycan type 개수 분포
- **heatmap_top_glycopeptides.png**: 상위 50개 glycopeptide heatmap

## 설정 변경

`config.yaml` 파일을 수정하여 다양한 파라미터를 조정할 수 있습니다:

- 경로 설정
- QC 필터링 기준
- Annotation 마커
- PCA 파라미터
- 시각화 옵션

## 개발자

개발 날짜: 2025-10-01

## 라이센스

MIT License
