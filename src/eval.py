import biotite.structure.io.pdbx as pdbx    # 읽어오기    
from pathlib import Path
import biotite.structure as struc
import numpy as np
import pandas as pd

####################
# 함수들

def calc_metrics(ans: struc.AtomArray, decoy: struc.AtomArray) -> tuple[float, float, float]:

    rmsd = f'{struc.rmspd(ans, decoy):.3f}'
    lddt = f'{struc.lddt(ans, decoy):.3f}'

    sup_pos = struc.superimpose(ans, decoy)
    n_range = np.arange(len(ans))
    tm = f'{struc.tm_score(ans, sup_pos[0], n_range, n_range):.3f}'

    return (rmsd, lddt, tm)

#####################
# 실행 코드


# 정답 데이터 가져오기

# Answer - Backbone
id_list = ["1CRN", "1CLL", "5DK3"]
ans_dict = {}
for id in id_list:
    with open(f"raw_dataset/parsed_{id}.cif", "rt", encoding = "utf-8") as f:
        if id == "1CRN":
            n1CRN = pdbx.CIFFile.read(f)
            n1CRN = pdbx.get_structure(n1CRN, model = 1)
            ans_dict[f"{id}"] = n1CRN

        if id == "1CLL":
            n1CLL = pdbx.CIFFile.read(f)
            n1CLL = pdbx.get_structure(n1CLL, model = 1)
            ans_dict[f"{id}"] = n1CLL

        if id == "5DK3":
            n5DK3 = pdbx.CIFFile.read(f)
            n5DK3 = pdbx.get_structure(n5DK3, model = 1)
            ans_dict[f"{id}"] = n5DK3
# print(ans_list["5DK3"])
ans_dict_bb = ans_dict

# Answer - Alpha Carbon
ans_dict_ca = {}
for id, ans in ans_dict.items():
    is_ca = (ans.get_annotation("atom_name") == "CA")
    ca_ans = ans[is_ca]
    ans_dict_ca[f"{id}"] = ca_ans

# test
# print(ans_dict_ca["5DK3"])



# decoy 가져오기

# decoy - backbone, Ca
decoy_1CRN_bb = {}
decoy_1CLL_bb = {}
decoy_5DK3_bb = {}
decoy_1CRN_ca = {}
decoy_1CLL_ca = {}
decoy_5DK3_ca = {}
decoy_dict = {
    "1CRN_bb" : decoy_1CRN_bb,
    "1CLL_bb" : decoy_1CLL_bb,
    "5DK3_bb" : decoy_5DK3_bb,
    "1CRN_ca" : decoy_1CRN_ca,
    "1CLL_ca" : decoy_1CLL_ca,
    "5DK3_ca" : decoy_5DK3_ca
    }

directory = Path("test_dataset/")
for file in sorted(directory.iterdir()):

    if not file.is_file:
        continue

    for id in id_list:
        if id in file.name:
            with open(file, "rt", encoding = "utf-8") as f:
                pre0 = pdbx.CIFFile.read(f)
                pre0 = pdbx.get_structure(pre0, model = 1)
                decoy_dict[f"{id}_bb"][file.stem] = pre0
                is_ca = (pre0.get_annotation("atom_name") == "CA")
                pre0 = pre0[is_ca]
                decoy_dict[f"{id}_ca"][file.stem] = pre0
                decoy_dict[file.name] = pre0
# print(decoy_dict["1CLL_bb"].keys())



########
# 데이터 구조가..
# ans_dict_bb[id] = 정답데이터_backbone
# ans_dict_ca[id] = 정답데이터_alphacarbon
# decoy_dict[id_bb] = { (파일이름 : decoy데이터)  이렇게 id의 backbone decoy 파일들 저장 }
# decoy_dict[id_ca] = { (파일이름 : decoy데이터)  이렇게 id의 backbone decoy 파일들 저장 }
########



# Metric 계산하기
rmsd_bb_list = []
lddt_bb_list = []
tm_list = []
rmsd_ca_list = []
lddt_ca_list = []

# backbone
for id in id_list:
    decoy_list = list(decoy_dict[f"{id}_bb"].keys())

    for decoy in decoy_list:
        ans = ans_dict_bb[f"{id}"]
        decoy = decoy_dict[f"{id}_bb"][f"{decoy}"]

        rmsd, lddt, _ = calc_metrics(ans, decoy)

        rmsd_bb_list.append(rmsd)
        lddt_bb_list.append(lddt)

# alpha carbon
for id in id_list:
    decoy_list = list(decoy_dict[f"{id}_ca"].keys())

    for decoy in decoy_list:
        ans = ans_dict_ca[f"{id}"]
        decoy = decoy_dict[f"{id}_ca"][f"{decoy}"]

        rmsd, lddt, tm = calc_metrics(ans, decoy)

        rmsd_ca_list.append(rmsd)
        lddt_ca_list.append(lddt)
        tm_list.append(tm)



# .csv 파일로 뽑기

# 저장순서
# ['decoy_1CRN_0.5angstrom', 'decoy_1CRN_1.0angstrom', 'decoy_1CRN_10degree', 'decoy_1CRN_2.0angstrom', 'decoy_1CRN_5.0angstrom', 'decoy_1CRN_5degree', 'decoy_1CRN_ex30', 'decoy_1CRN_ex60']


#                      |rmsd-bb | rmsd-ca | tmscore | lddt-bb | lddt-ca
#  noise   0.5A
#          1.0A
#          2.0A
#          5.0A  
#  torsion 5degree
#  angle   10degree
#  local   30degree
# pertur.  60degree


pre_dataset = np.column_stack((rmsd_bb_list, rmsd_ca_list, tm_list, lddt_bb_list, lddt_ca_list))
# print(pre_dataset)
i = 0
data_dict = {}
for id in id_list:
    dataset = pre_dataset[[i, i + 1, i + 3, i + 4, i + 5, i + 2, i + 6, i + 7]]
    data_dict[f"{id}"] = dataset
    i += 8
    # print(dataset)
# print(data_dict)

columns = ["rmsd-bb", "rmsd-ca", "tmscore", "lddt-bb", "lddt-ca"]
indices = pd.MultiIndex.from_tuples([
    ('noise', '0.5Å'),
    ('noise', '1.0Å'),
    ('noise', '2.0Å'),
    ('noise', '5.0Å'),
    ('global torsion', "5°"),
    ('global torsion', "10°"),
    ('local torsion', "30°"),
    ('local torsion', "60°")
], names = ['Type', 'Perturbation'])

df_dict: dict[str, pd.DataFrame] = {}
for id, data in data_dict.items():
    df = pd.DataFrame(data, index= indices, columns = columns)
    df_dict[id] = df
# print(data_dict["1CRN"])
# print(type(df_dict["1CRN"]))

for id, df in df_dict.items():
    df.to_csv(f'result/{id}_result.csv', index= True, encoding='utf-8-sig', )


