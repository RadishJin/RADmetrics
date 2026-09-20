#!/bin/bash

# 데이터 넣을 디렉토리 정하고, 없으면 생성
TARGET_DIR="./raw_dataset"
# mkdir -p(Parent: 이미 있으면 그냥 넘어가)
mkdir -p "$TARGET_DIR"


# Target Protein
    # Standard Rigid Target - (1CRN - Crambin)
    # Domain-Flexible Target - (1CLL - Calmodulin)
    # Antibody - (5DK3 - Pembrolizumab Fab)

# Bash에서는 리스트를 공백으로 구분
TARGETS=("1CRN" "1CLL" "5DK3")


# curl(Client URL : 웹 서버와 통신하는 표준 CLI)
# curl -o(output : 원하는 형식으로 데이터 받기)
# curl -f(Fail Silently: 서버 에러 발생하면 즉시 실패 처리)
# curl -s(Silent: 다운로드 진행 바 안 띄움)
# echo (프린트랑 비슷함)

# ⚬	-eq : Equal (==)
# ⚬	-ne : Not Equal (!=)
# ⚬	-gt : Greater Than (>)
# ⚬	-ge : Greater or Equal (>=)
# ⚬	-lt : Less Than (<)
# ⚬	-le : Less or Equal (<=)
# $? : 저장된 종료 코드. 0이면 정상작동

# if 뒤의 대괄호 조건 앞위로 공백 필수, 뒤에 then 붙이고 fi로 마무리

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