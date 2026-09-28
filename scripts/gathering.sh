#!/bin/bash

# 데이터 넣을 디렉토리 정하고, 없으면 생성
TARGET_DIR="./raw_dataset"
# mkdir -p(Parent: 이미 있으면 그냥 넘어가)
mkdir -p "$TARGET_DIR"


# Bash에서는 리스트를 공백으로 구분
TARGETS=("1CRN" "1CLL" "5DK3")

# API로 다운로드
for ID in "${TARGETS[@]}"; do
    echo "Downloading ${ID}.cif.gz..."
    curl -fs -o "${TARGET_DIR}/${ID}.cif.gz" "https://files.rcsb.org/download/${ID}.cif.gz"
    if [[ $? -eq 0 ]]; then
        echo "  [SUCCESS] Saved to ${TARGET_DIR}/${ID}.cif.gz"
    else
        echo "  [ERROR] Failed to download ${ID}"
    fi
done
