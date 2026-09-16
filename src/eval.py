import biotite.structure.io.pdbx as pdbx    # 읽어오기    
from pathlib import Path
import biotite.structure as struc

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



# Metric 계산하기


# rmsd - backbone

rmspd = struc.rmspd(ans_dict_bb["1CRN"], decoy_dict["1CRN_bb"]["decoy_1CRN_ex30"])
lddt = struc.lddt(ans_dict_bb["1CRN"], decoy_dict["1CRN_bb"]["decoy_1CRN_ex30"])
print(rmspd, lddt)


# rmsd - Ca

# lDDT - backbone

# lDDT - Ca

# TMscore


# + RADetector

