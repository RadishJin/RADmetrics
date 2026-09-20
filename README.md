# RADmetrics

---
RADmetrics is a lightweight pipeline designed for benchmarking protein structural evaluation metrics:
RMSD (Backbone, Ca), lDDT (Backbone, Ca), TM-score

Made by Mu Jin Kim (BS student @ Kyungpook National University, Department of Biotechnology)
Contact : radishj24@gmail.com
Managed under GitFlow branching model
---


 - Versions
    python     3.11
    biotite    1.6.0
    numpy      1.26.4
    pytorch    2.13.0
    pandas     3.0.5

 - Directories
    RADetector/

        .vscode/
            settings.json

        raw_dataset/        # Raw dataset from PDB [Ignored by Git]
        test_dataset/       # Noise-injected decoy datasets [Ignored by Git]

        result/             # Final evaluation summary [Ignored by Git]
        src/
            gathering.sh    # Bash script for fetching datasets
            parsing.py      # Data parsing and tensor conversion
            preparing.py    # Torsion angle 기반 decoy 생성
            eval.py         # Metric calculation with decoy (TM-score, lDDT, RMSD)

        .gitignore
        LICENSE
        README.md
        run.sh


 - Brief Explanation
    Algorithm
    1. 데이터 전처리 (정답 데이터셋 만들기)
        1. Bash, Protein Dataset 가져오기
        2. Python, Data Parsing

    2. 테스트 데이터셋 만들기
        1. AtomArray torch 이용 텐서화
        2. 정답 데이터에 특정 강도 노이즈 주입
        3. 테스트 데이터로 저장

    3. 결과 정리하기
        1. 정답 데이터, 테스트 데이터 가져오기
        2. TM-score, lDDT, RMSD 산출
        3. 결과 요약 출력 (.csv)

    Target Protein
	    Standard Rigid Target - (1CRN - Crambin)
	    Domain-Flexible Target - (1CLL - Calmodulin)
	    Antibody - (5DK3 - Pembrolizumab Fab)

------

Result & Review

Limitations
    1. 기존 계획은 target 특성 고려한 노이즈 주입이었으나, 모든 대상에 모든 노이즈 타입을 주어 decoy 수를 늘렸음.  
       decoy 다양성 측면에서는 성공적이었으나 타깃 프로틴을 rigid non-rigid 로 나눈 의미가 옅어짐  
    2. global torsion angle 에 gaussian noise 주입한 것이 lever-arm effect 가져옴. rmsd, lddt, tmscore 전부 값이 폭발함.  
    3. decoy가 수학적 맥락만 고려하여 만들어짐. 따라서 물리학적으로 대부분이 불가능한 형태임.  
    4. 일부 하드코딩된 부분과 무거운 반복문이 산재. 대용량 데이터로 확장성이 떨어지는 상태임.  
    5. NSAA 까지 고려하여 All-Atom으로 최초에 설계하였으나, 확장성과 연산 편의성 간의 저울질 이후 BackBone과 AlphaCarbon으로 타협함.  
    6. MultiChain(5DK3)의 경우 local perturbation decoy를 만들 때 4개의 체인 중 랜덤하게 한 체인의 중심부에만 변형을 가하였고,  
       따라서 체인 간의 intermolecular force를 예측하는 데이터셋으로 이 decoy를 제공한다고 했을 때 유의미한 Question data가 되지 못함.  
    7. 실행용 bash script 에서 interpreter 설정이 다른 사용자 pc에서도 결과물 산출까지 과정을 에러 없이 거치는지에 대해 검증할 방법이 현재로서는 없음.  
    8. python file이 몇몇 서로 다른 기능을 같이 가지고 있는 경우가 있음. 이를테면 eval.py는 pandas dataframe으로 만드는 작업과 metric 계산을 한 파일 안에서 같이 수행함.  
  
Potential Improvements  
    1. Ramachandran plot, B-factor 이용하여 조금 더 Physical 한 decoy set을 만드는 것에 이용 가능해보임.  
    2. 연산량을 최적화 한 이후, 오픈소스 이용하여 NSAA를 가장 가까운 parental standard AA로 근사하여 All-Atom 구현에 도전해볼 수 있어보임.  
    3. Target Protein의 생체 내에서의 구조적 특징을 따와서, 가장 가질 법한 구조들을 decoy로 만드는 일종의 structure selection을 decoy set 구축 이후 도입할 시  
       decoy 다양성, 규모와 physical, chemical 구조를 동시에 잡는 것이 가능해보임.  
    
