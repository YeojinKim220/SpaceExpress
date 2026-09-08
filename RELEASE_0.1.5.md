# SpaceExpress 0.1.5

English version: [RELEASE_0.1.5-en.md](RELEASE_0.1.5-en.md)

## 변경 사항

- 실패한 유전자-embedding 차원 적합(음수 또는 비유한 통계량)의 FDR은 1입니다. 유효한 empirical-null 적합이 불가능한 차원도 FDR=1로 반환합니다.
- 최종 FDR을 [0, 1]로 제한합니다. 기존 empirical-null 추정식을 다른 보정법으로 교체하지 않습니다.
- `select_hvg_after_outlier`는 정규화·log 변환된 공통 유전자에서 두 표본을 합친 평균 + 4 표본 표준편차 이상의 값을 0으로 바꾼 뒤 batch-aware HVG를 선택합니다. 입력 객체 대신 복사본을 수정합니다.
- 초기 전처리 이력을 `.uns['spaceexpress_preprocessing']`에 저장합니다. 모든 입력에 이 이력이 있으면 DSE에서 중복 평균+4SD 제거를 건너뜁니다. 이력이 없는 기존 입력은 기존 동작을 유지하고, 이력이 섞여 있으면 명시적 설정을 요구합니다.
- `SpaceExpress_DSE(..., remove_mean_sd_outliers=False)`로 후속 평균+4SD 제거를 명시적으로 끌 수도 있습니다. 학습 함수의 기존 95% clipping은 변경하지 않았습니다. 다중 반복 표본 DSE의 별도 99% 필터도 유지합니다.
- 반환 AnnData의 `.varm`에 `DSE-statistic`, `DSE-fit-failed`를 추가하여 실패를 0 예측값으로 추측하지 않고 확인할 수 있습니다.

## 설치

Python 3.11 이상과 R이 필요합니다. R 패키지 `lmtest`, `fitdistrplus`, `dplyr`, `lme4`를 먼저 설치하고 Python에서 R 공유 라이브러리를 찾을 수 있도록 설정하십시오. Python 의존성은 패키지 메타데이터에 선언했습니다. CUDA에 맞는 PyTorch 설치는 사용하는 GPU 환경을 따릅니다.

```bash
Rscript -e 'install.packages(c("lmtest", "fitdistrplus", "dplyr", "lme4"), repos="https://cloud.r-project.org")'
python -m pip install 'spaceexpress[notebooks]==0.1.5'
```

`examples/install_release.sh BASE_PYTHON NEW_ENV`는 기존 과학 계산·R 환경을 공유하는 새 venv에 PyPI 릴리스를 설치합니다. 패키지 다운로드 출처는 `pip_install_report.json`, 의존성은 `pip_freeze.txt`에 기록합니다. 전체 의존성을 처음부터 새로 설치하는 격리 환경과는 다릅니다.

## 데이터 실행

`examples/public_pair.py`는 config 기준 상대 경로를 사용하며 원본을 덮어쓰지 않습니다. 원본 count에서 최소 3개 검출 유전자·유한 spatial 좌표를 확인하고, 필요하면 공간 구획별 비례 표집합니다. 두 표본에서 각각 3개 이상 관측치에 검출된 공통 유전자를 정규화(총량 10,000)하고 log1p 후 위 전처리로 HVG 200개를 선택합니다. 원본 선택 데이터와 QC 수치도 저장합니다.

```bash
bash examples/submit_public_pair.sh /path/to/env/bin/python /path/to/config_v015.json /path/to/results_v015
```

Slurm 기본 자원은 전처리 12 CPU/128 GB, 학습 H100 80 GB 1장·8 CPU/256 GB RAM, DSE 작업당 12 CPU/128 GB, 노트북 4 CPU/64 GB입니다. 각 작업의 시간 한도는 4시간입니다. `SE_ACCOUNT`, `SE_QOS`, `SE_CPU_PARTITION`, `SE_GPU_PARTITION`, `SE_GPU_CONSTRAINT`로 클러스터 설정을 변경할 수 있습니다. CPU 기본 파티션은 `cpu-medium`입니다. A100 사용 시 `SE_GPU_PARTITION=gpu-a100 SE_GPU_CONSTRAINT=A100-80GB`를 설정하십시오. 이는 요청량이며 측정된 필요량이 아닙니다.

전처리, 전체 관측치 2-epoch GPU probe, 본 학습, k=30/50/100 DSE, 결과 노트북 순서로 실행합니다. 앞 단계 성공을 의존 조건으로 사용합니다. 본 학습은 config의 `epochs`, `patience`를 사용합니다. 실행 폴더 재사용은 거부하므로 새 실험에는 새 출력 경로를 지정하십시오.

결과에는 `sample_summary.csv`, QC 및 유전자 진단 CSV, 선택 원본·전처리·embedding H5AD, 모델 상태, FDR·통계량·실패 마스크 CSV, DSE pickle, 단계별 JSON·시간·메모리 로그, 그림이 포함된 실행 완료 `results.ipynb`가 포함됩니다. Notebook에서는 원본 대비 실제 관측치 수, 공간 분포, counts·검출 유전자 수, embedding, 실패 적합, k 민감도, 유전자 발현과 적합 결과를 확인합니다. Pickle은 신뢰하는 실행에서 생성한 파일만 여십시오.

## 해석과 한계

95% clipping은 희소 유전자의 상한이 0이 되면 그 표본에서 해당 유전자의 학습 신호 전체를 없앨 수 있습니다. 초기 4SD 제거 이후에도 발생할 수 있으며, `gene_qc_*.csv`와 notebook에 기록합니다. 이번 릴리스는 합의된 두 전처리를 유지하므로 이를 자동으로 다른 방법으로 바꾸지 않습니다.

Slide-seq는 bead, Stereo-seq는 spatial bin일 수 있으므로 관측치 수를 실제 단일 세포 수와 동일시하지 않습니다. Brain tissue annotation만 있는 입력에서 세부 세포 유형을 추정해 표시하지 않습니다. 미토콘드리아 비율이나 검출량에 추가 임계값을 자동 적용하지 않으며, 분포와 표본별 QC를 검토한 뒤 결정해야 합니다.

FDR 범위·실패 처리·관측치 보존 검사는 수치적 동작 검사이지 생물학적 검증이 아닙니다. k 안정성과 공간 이웃 보존율도 함께 검토하십시오. 조건당 생물학적 표본 1개인 비교는 반복 표본 기반의 집단 수준 추론을 대신하지 못합니다.
