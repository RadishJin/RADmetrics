#!/bin/bash

# 에러 발생 시 즉시 스크립트 중단 (-e) 및 미정의 변수 사용 금지 (-u), 뒤 옵션 지정 (-o)
set -euo pipefail

echo "[Step 1/4] Downloading raw dataset..."
bash src/gathering.sh

echo "[Step 2/4] Parsing files..."
python src/parsing.py

echo "[Step 3/4] Generating decoys..."
python src/preparing.py

echo "[Step 4/4] Calculating metrics and plotting results..."
python src/eval.py

echo "=========================================="
echo "Successfully finished"
echo "=========================================="
