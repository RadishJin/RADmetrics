import biotite.structure.io.pdbx as pdbx    # 읽어오기    
from pathlib import Path
import biotite.structure as struc
import numpy as np


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

# metric 결과값 저장소
decoy_rmsd_dict = {}
decoy_lddt_dict = {}
decoy_tmscore_dict = {}
for id in id_list:
    decoy_rmsd_dict.setdefault(f"{id}_bb", [])
    decoy_rmsd_dict.setdefault(f"{id}_ca", [])
# print(decoy_rmsd_dict.keys())
for id in id_list:
    decoy_lddt_dict.setdefault(f"{id}_bb", [])
    decoy_lddt_dict.setdefault(f"{id}_ca", [])
# print(decoy_lddt_dict.keys())
for id in id_list:
    decoy_tmscore_dict.setdefault(f"{id}", [])
# print(decoy_tmscore_dict.keys())


# test
# ['decoy_1CRN_0.5angstrom', 'decoy_1CRN_1.0angstrom', 'decoy_1CRN_10degree', 'decoy_1CRN_2.0angstrom', 'decoy_1CRN_5.0angstrom', 'decoy_1CRN_5degree', 'decoy_1CRN_ex30', 'decoy_1CRN_ex60']
# 0.7153824065272497
# 1.4128974324698986
# 58.361695975918884
# 2.8623942985583017
# 7.160457759030418
# 28.658485820659063
# 7.191420607572877
# 13.848104800380229

# rmsd - backbone
rmsd_bb_list = []
decoy_list = list(decoy_dict)
for id in id_list:
    decoy_list = list(decoy_dict[f"{id}_bb"].keys())
    # print(decoy_list)
    for decoy in decoy_list:
        rmspd = struc.rmspd(ans_dict_bb[f"{id}"], decoy_dict[f"{id}_bb"][f"{decoy}"])
        rmsd_bb_list.append(rmspd)

# rmsd - Ca
rmsd_ca_list = []
decoy_list = list(decoy_dict)
for id in id_list:
    decoy_list = list(decoy_dict[f"{id}_ca"].keys())
    # print(decoy_list)
    for decoy in decoy_list:
        rmspd = struc.rmspd(ans_dict_ca[f"{id}"], decoy_dict[f"{id}_ca"][f"{decoy}"])
        rmsd_ca_list.append(rmspd)

# lDDT - backbone
lddt_bb_list = []
decoy_list = list(decoy_dict)
for id in id_list:
    decoy_list = list(decoy_dict[f"{id}_bb"].keys())
    # print(decoy_list)
    for decoy in decoy_list:
        lddt = struc.lddt(ans_dict_bb[f"{id}"], decoy_dict[f"{id}_bb"][f"{decoy}"])
        lddt_bb_list.append(lddt)

# lDDT - Ca
lddt_ca_list = []
decoy_list = list(decoy_dict)
for id in id_list:
    decoy_list = list(decoy_dict[f"{id}_ca"].keys())
    # print(decoy_list)
    for decoy in decoy_list:
        lddt = struc.lddt(ans_dict_ca[f"{id}"], decoy_dict[f"{id}_ca"][f"{decoy}"])
        lddt_ca_list.append(lddt)

# TMscore
tm_list = []
decoy_list = list(decoy_dict)
for id in id_list:
    decoy_list = list(decoy_dict[f"{id}_ca"].keys())
    # print(decoy_list)
    for decoy in decoy_list:
        sup_pos = struc.superimpose_structural_homologs(ans_dict_ca[f"{id}"], decoy_dict[f"{id}_ca"][f"{decoy}"],)
        # print(type(sup_pos))
        n_range = np.arange(len(ans_dict_ca[f"{id}"]))
        # print(n_range)
        tm = struc.tm_score(ans_dict_ca[f"{id}"], sup_pos[0], n_range, n_range)
        tm_list.append(tm)
# print(tm_list)

# .csv 파일로 뽑기
# 이미지로 뽑기

