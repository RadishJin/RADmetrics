# 목표 : (Backbone, Ca) RMSD, (Backbone, Ca) lDDT, TM-score 를 위한 데이터 전처리

import gzip                                 # .gz 압축된 파일 읽기용
import biotite.structure.io.pdbx as pdbx    # .cif 파일 읽기용
import biotite.structure as struc           # 텍스트 파일 파싱용


id_list = ["1CRN", "1CLL", "5DK3"]

for id in id_list:
        
    # gzip으로 압축파일 간단하게 열고, biotite 라이브러리 이용해서 객체로 바로 읽어오기
    with gzip.open(f"raw_dataset/{id}.cif.gz", "rt", encoding = "utf-8") as f:
        raw = pdbx.CIFFile.read(f)
    # print(raw)


    # pdbx.get_structrue로 텍스트데이터를 AtomArray로 파싱
    atoms = pdbx.get_structure(raw, model= 1) # 1D로 파싱하기 위해 model = 1
    # print(atoms)

    # Biotite.structure 이용 residue와 관련된 원자만 남기기 (비표준 아미노산도 포함)
    residue_mask = struc.filter_amino_acids(atoms)
    residue_atom = atoms[residue_mask]


    # print(list(set(residue_atom.get_annotation("res_id"))))
    # print(list(set(residue_atom.get_annotation("res_name"))))
    # 앞뒤로 Residue 4개 missing. 앞으로 고려 필요함. 끝부분이라서 달랑거렸을테니 Crystallography의 한계인듯? 


    # Backbone atoms

    is_backbone = struc.filter_peptide_backbone(residue_atom)
    bb_atoms = residue_atom[is_backbone]
    # print(bb_atoms)

    # .cif 파일로 저장
    cif_file = pdbx.CIFFile()
    pdbx.set_structure(cif_file, bb_atoms, data_block=f"backbone_{id}")
    cif_file.write(f"raw_dataset/parsed_{id}.cif")



        # Zero - Centering
        # center = struc.centroid(bb_atoms)
        # print(center)
        # superposition 할거라면 의미없을듯









# Ca 만 따로 남기기
# is_ca = (bb_atoms.get_annotation("atom_name") == "CA")
# print(is_ca)
# ca_atoms = bb_atoms[is_ca]
# print(ca_atoms)


# # All-Atom trial

# # Calmodulin의 경우엔 필요없지만, residue 전부 canonical로 변환
# # 어짜피 RMSD나 lDDT나 위치가 중요하므로 바꿔도 문제가 적을듯
# # 실제로 모든 NSAA를 매핑하는건 말이 안되니까 이 하드코딩 방법이 무식한것도 맞지만 PDBfixer로 데이터 상 오류나 위치 상 오류를 감수하는것보단 나음
# nsaa_map = {
#     "MSE" : "MET",
#     "SEP" : "SER",
#     "TPO" : "THR",
#     "PTR" : "TYR",
#     "HYP" : "PRO",
#     "KCX" : "LYS",
#     "CSO" : "CYS",
#     "CSX" : "CYS",
#     "CME" : "CYS",
#     "MLY" : "LYS",
#     "PCA" : "GLU"
# }

# #        A         4                    LEU                     CA           C        -6.696   22.003   26.447
# # ['chain_id', 'res_id', 'ins_code', 'res_name', 'hetero', 'atom_name', 'element']

# # 전후비교용
# print(residue_atom.get_annotation("res_name"))

# # numpy 기반이라 가능한 속성 이용 불리언마스킹
# aa_list = residue_atom.get_annotation("res_name")
# for nsaa, saa in nsaa_map.items():
#     aa_list[aa_list == nsaa] = saa

# # get_annotation으로 가져오고 set_annotation으로 변경
# residue_atom.set_annotation("res_name", aa_list)
# print(residue_atom.get_annotation("res_name"))

# # 아직 이름만 바뀜.. 원자는 도대체 어떻게 다 바꿔야할까