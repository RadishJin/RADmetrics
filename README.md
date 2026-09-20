# RADmetrics

Benchmarking protein structural evaluation metrics (**RMSD**, **lDDT**, **TM-score**)

- **Author**: Mu Jin Kim (B.S. student @ Kyungpook National University, Department of Biotechnology)
- **Contact**: radishj24@gmail.com
- **Git Strategy**: Managed under GitFlow branching model

---

## Environment & Dependencies

| Tool / Library | Version |
| :--- | :--- |
| **Python** | `3.11` |
| **Biotite** | `1.6.0` |
| **NumPy** | `1.26.4` |
| **PyTorch** | `2.13.0` |
| **Pandas** | `3.0.5` |

---

## Directory Structure

```text
RADetector/
├── .vscode/
│   └── settings.json
├── raw_dataset/        # Raw dataset from PDB [Ignored by Git]
├── test_dataset/       # Noise-injected decoy datasets [Ignored by Git]
├── result/             # Final evaluation summary [Ignored by Git]
├── src/
│   ├── gathering.sh    # Bash script for fetching datasets
│   ├── parsing.py      # Data parsing and tensor conversion
│   ├── preparing.py    # Torsion angle 기반 decoy 생성
│   └── eval.py         # Metric calculation with decoy (TM-score, lDDT, RMSD)
├── .gitignore
├── LICENSE
├── README.md
└── run.sh

## Overview

### Pipeline Logic

1. **데이터 전처리 (정답 데이터셋 만들기)**
   1. Bash 이용 Protein Dataset 수집
   2. Python 이용 Data Parsing

2. **테스트 데이터셋 만들기**
   1. `AtomArray` torch 이용 텐서화
   2. 정답 데이터에 특정 강도 노이즈 주입
   3. 테스트 데이터로 저장

3. **결과 정리하기**
   1. 정답 데이터 및 테스트 데이터 로드
   2. TM-score, lDDT, RMSD 산출
   3. 결과 요약 출력 (`.csv`)

### Target Proteins

* **Standard Rigid Target**: `1CRN` (Crambin)
* **Domain-Flexible Target**: `1CLL` (Calmodulin)
* **Antibody**: `5DK3` (Pembrolizumab Fab)

---

## Result & Review

### Limitations

1. **타깃 분류 의미 퇴색**: 기존 계획은 Target 특성을 고려한 노이즈 주입이었으나, 모든 대상에 모든 노이즈 타입을 주어 Decoy 수를 늘림. Decoy 다양성 측면에서는 성공적이었으나 Rigid / Non-rigid 구분 의미가 옅어짐.
2. **Lever-Arm Effect**: Global Torsion Angle에 Gaussian Noise를 주입하면서 Lever-arm effect가 발생함. 이로 인해 RMSD, lDDT, TM-score 산출값이 모두 폭발적인 비정상치를 보임.
3. **물리학적 타당성 부족**: Decoy가 수학적 맥락만 고려되어 생성됨. 이에 따라 물리학적으로 형성이 불가능한 구조가 대부분을 차지함.
4. **확장성 및 성능 병목**: 일부 하드코딩된 로직과 무거운 반복문이 산재하여, 대용량 데이터로의 확장성이 떨어짐.
5. **구조 표현범위 축소**: NSAA(Non-Standard Amino Acid)까지 고려한 All-Atom 설계를 추진했으나, 확장성과 연산 편의성 간 타협으로 Backbone 및 $\text{C}\alpha$ 수준으로 축소됨.
6. **Multi-Chain Perturbation의 한계**: Multi-chain (`5DK3`) local perturbation decoy 생성 시 4개 체인 중 임의의 1개 체인 중심부에만 변형을 가함. 결과적으로 체인 간 Intermolecular force를 예측/평가하는 데이터셋으로 활용하기에 유의미한 Question data가 되지 못함.
7. **이식성 검증 미비**: 실행용 Bash script 내 Interpreter 설정이 다른 사용자 환경에서도 에러 없이 작동하는지 검증할 수 있는 테스트 체계가 부족함.
8. **단일 책임 원칙(SRP) 위배**: Python 파일 내 모듈화 미흡. (예: `eval.py`에서 Metric 계산과 Pandas DataFrame 가공 작업을 동시에 수행함)

### Potential Improvements

1. **Physical Decoy Set 구축**: Ramachandran plot 및 B-factor를 도입하여 물리·화학적으로 타당한 Decoy 생성을 구현.
2. **All-Atom 표현 복원**: 연산량 최적화 수행 후, 오픈소스를 활용하여 NSAA를 가장 가까운 Parental Standard Amino Acid로 근사(Approximation)하여 All-Atom 모델링 도전.
3. **Structure Selection 도입**: Target protein의 생체 내 구조적 특성을 반영하여 실제 형성 가능한 구조 위주로 선별하는 Structure Selection 프로세스를 구축. Decoy의 다양성, 규모, Physical/Chemical 타당성을 동시에 확보 가능.
