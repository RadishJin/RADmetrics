#!/bin/bash

# 에러 발생 시 즉시 스크립트 중단 (-e) 및 미정의 변수 사용 금지 (-u), 뒤 옵션 지정 (-o)
set -euo pipefail

python_exe = $(which python3)

echo "[Step 1/4] Downloading raw dataset..."
bash src/gathering.sh

echo "[Step 2/4] Parsing..."
 src/parsing.py

# 3. Decoy 생성
echo "[Step 3/4] Generating decoys..."
python_exe src/preparing.py

# 4. Metric 계산 및 시각화 리포트 생성
echo "[Step 4/4] Calculating metrics and plotting results..."
python_exe src/eval.py

echo "=========================================="
echo " Pipeline Finished Successfully!"
echo "=========================================="
