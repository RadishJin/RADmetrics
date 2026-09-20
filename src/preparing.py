import biotite.structure.io.pdbx as pdbx    # 읽어오기     
import biotite.structure as struc
import torch                                # 텐서 변환, 노이즈 생성용 
from math import pi
import random


id_list = ["1CRN", "1CLL", "5DK3"]
data_dict = {}
for id in id_list:    
    with open(f"raw_dataset/parsed_{id}.cif", "rt", encoding = "utf-8") as f:
        data_dict[f"{id}"] = pdbx.get_structure(pdbx.CIFFile.read(f))

# test
# print(data_dict["1CLL"])



# 각 target의 특징에 따라 decoy 생성
id_single_chain = ["1CRN", "1CLL"]
multi_chain = data_dict["5DK3"][0]
decoy_dict = {}

# 0. 멀티체인 전처리 (체인별)

chain_start = struc.get_chain_starts(multi_chain)
# print(chain_start)
chain_mask = struc.get_chain_masks(multi_chain, chain_start)
chain_list = []
for i in range(len(chain_mask)):
    strand = multi_chain[chain_mask[i]]
    chain_list.append(strand)
# print(chain_list)

# 1. 1CRN Crambin (가장 기본형)


# (1) Gaussian Coordiante Noise (sigma = 0.5, 1.0, 2.0, 5.0 angstrom)

# single chain
for id in id_single_chain:
    bb_coord = struc.coord(data_dict[f"{id}"])
    bb_tensor = torch.tensor(bb_coord)
    sigmas = [0.5, 1.0, 2.0, 5.0]
    for sigma in sigmas:
        noise = torch.randn_like(bb_tensor) * sigma
        noise_coord = bb_tensor + noise
        noise_coord = noise_coord.detach().cpu().numpy()
        pre1 = data_dict[f"{id}"].copy()
        struc.coord(pre1)[:] = noise_coord
        decoy_dict[f"{id}_{sigma}angstrom"] = pre1

# Multi chain

# 체인별로 나눠서 노이즈 주입
chain_decoy = {}
n = 0
for chain in chain_list:
    bb_coord = struc.coord(chain)
    bb_tensor = torch.tensor(bb_coord)
    sigmas = [0.5, 1.0, 2.0, 5.0]
    for sigma in sigmas:
        noise = torch.randn_like(bb_tensor) * sigma
        noise_coord = bb_tensor + noise
        noise_coord = noise_coord.detach().cpu().numpy()
        pre1 = chain.copy()
        struc.coord(pre1)[:] = noise_coord
        chain_decoy[f"n{sigma}angstrom_chain{n}"] = pre1
    n += 1

# print(chain_decoy.keys())
# print(chain_list[1][:10])
# print(chain_decoy["chain1_2.0angstrom"][:10])

# 하나로 합치기
ang05, ang10, ang20, ang50 = [], [], [], []
for id, chain in chain_decoy.items():
    if id.startswith(f"n{sigmas[0]}"):
        ang05.append(chain)
    if id.startswith(f"n{sigmas[1]}"):
        ang10.append(chain)
    if id.startswith(f"n{sigmas[2]}"):
        ang20.append(chain)
    if id.startswith(f"n{sigmas[3]}"):
        ang50.append(chain)

n = 0
for i in ang05, ang10, ang20, ang50:
    combined_chain = sum(i, struc.AtomArray(0))
    decoy_dict[f"5DK3_{sigmas[n]}angstrom"] = combined_chain
    n += 1

# print(decoy_dict.keys())
# print(decoy_dict["5DK3_0.5angstrom"][:10])
# print(data_dict["5DK3"][0][:10])



    
# test
# print(noise_coord)
# print(decoy_dict["1CRN_5.0angstrom"])
# print(decoy_dict.keys())


# (2) Global Perturbation (Gaussian, sigma = 5, 10 degree)

# single chain
for id in id_single_chain:
    # 바꿀 토션앵글 값 구하기
    bb_phi, bb_psi, bb_omega = struc.dihedral_backbone(data_dict[f"{id}"])
    bb_phi_tensor = torch.tensor(bb_phi)
    bb_psi_tensor = torch.tensor(bb_psi)
    bb_omega_tensor = torch.tensor(bb_omega)
    tensor_list = [bb_phi_tensor, bb_psi_tensor, bb_omega_tensor]
    angle_noise_list = [pi/36, pi/18]
    pre_decoy = []
    for angle in angle_noise_list:
        for tensor in tensor_list:
            noise = torch.randn_like(tensor) * angle
            noise = noise.detach().cpu().numpy().tolist()
            pre_decoy.append(noise[0])

    # print(bb_omega)
    # print(pre_decoy[5])

    # 델타 토션앵글 리스트화
    five_degree = [i for j in zip(pre_decoy[0], pre_decoy[1], pre_decoy[2]) for i in j]
    ten_degree = [i for j in zip(pre_decoy[3], pre_decoy[4], pre_decoy[5]) for i in j]
    # print(ten_degree)

    # 쓰이지 않는, 의미없는 값 잘라내기
    five_degree = five_degree[1:-1]
    ten_degree = ten_degree[1:-1]
    # print(five_degree, ten_degree)


    # print(pre[0].coord)
    # N-term부터 하나씩 돌리면서 새 구조 만들기
    j = 0
    for angle in five_degree, ten_degree:
        pre2 = data_dict[f"{id}"][0].copy()
        for i in range(len(angle)):
            if i+2 > len(pre2):
                break
            axis = pre2[i+1].coord - pre2[i].coord
            support = pre2[i+1].coord
            downstream = pre2[i+2:]
            downstream = struc.rotate_about_axis(
                downstream,
                angle = angle[i],
                axis = axis,
                support = support
            )
            struc.coord(pre2[i+2:])[:] = struc.coord(downstream)
        j += 5
        decoy_dict[f"{id}_{j}degree"] = pre2


# test
# print(data_dict["1CRN"][0])
# print(pre1[0][-5:])
# print(five_degree)
# print(pre2)

# print(decoy_dict.keys())

# Multi Chain
chain_decoy = {}
n = 0
for chain in chain_list:
    bb_phi, bb_psi, bb_omega = struc.dihedral_backbone(chain)
    bb_phi_tensor = torch.tensor(bb_phi)
    bb_psi_tensor = torch.tensor(bb_psi)
    bb_omega_tensor = torch.tensor(bb_omega)
    tensor_list = [bb_phi_tensor, bb_psi_tensor, bb_omega_tensor]
    angle_noise_list = [pi/36, pi/18]
    pre_decoy = []
    for angle in angle_noise_list:
        for tensor in tensor_list:
            noise = torch.randn_like(tensor) * angle
            noise = noise.detach().cpu().numpy().tolist()
            pre_decoy.append(noise)

    five_degree = [i for j in zip(pre_decoy[0], pre_decoy[1], pre_decoy[2]) for i in j]
    ten_degree = [i for j in zip(pre_decoy[3], pre_decoy[4], pre_decoy[5]) for i in j]

    five_degree = five_degree[1:-1]
    ten_degree = ten_degree[1:-1]

    j = 0
    for angle in five_degree, ten_degree:
        pre2 = chain.copy()
        for i in range(len(angle)):
            if i+2 > len(pre2):
                break
            axis = pre2[i+1].coord - pre2[i].coord
            support = pre2[i+1].coord
            downstream = pre2[i+2:]
            downstream = struc.rotate_about_axis(
                downstream,
                angle = angle[i],
                axis = axis,
                support = support
            )
            struc.coord(pre2[i+2:])[:] = struc.coord(downstream)
        j += 5
        chain_decoy[f"n{j}degree_{n}"] = pre2
    n += 1
# print(chain_decoy.keys())

five, ten = [], []
angles = [5, 10]
for id, chain in chain_decoy.items():
    if id.startswith(f"n{angles[0]}"):
        five.append(chain)
    if id.startswith(f"n{angles[1]}"):
        ten.append(chain)

n = 0
for i in five, ten:
    combined_chain = sum(i, struc.AtomArray(0))
    decoy_dict[f"5DK3_{angles[n]}degree"] = combined_chain
    n += 1


# (3) Local Extreme Perturbation (대충 가운데 부근 한 곳에서 결합각 대폭 조정, 30도 60도.)

# single chain
for id in id_single_chain:
    # 원본 데이터셋 가져와서
    pre3 = data_dict[f"{id}"][0].copy()
    pre4 = data_dict[f"{id}"][0].copy()

    # 돌릴 각도 저장해놓고
    torsion = [pi/6, pi/3]

    # 대충 중간 지점 정하고
    mid = int(len(pre3)/2)
    # print(mid)

    # 변형 - 30degree
    downstream = pre3[mid:]
    axis = pre3[mid - 1].coord - pre3[mid - 2].coord
    support = pre3[mid-1].coord
    downstream = struc.rotate_about_axis(
        downstream,
        angle = torsion[0],
        axis = axis,
        support = support
    )
    struc.coord(pre3[mid:])[:] = struc.coord(downstream)
    decoy_dict[f"{id}_ex30"] = pre3

    # 변형 - 60degree
    downstream = pre4[mid:]
    axis = pre4[mid-1].coord - pre4[mid-2].coord
    support = pre4[mid-1].coord
    downstream = struc.rotate_about_axis(
        downstream,
        angle = torsion[1],
        axis = axis,
        support = support
    )
    struc.coord(pre4[mid:])[:] = struc.coord(downstream)
    decoy_dict[f"{id}_ex60"] = pre4




# Multi Chain
chain_decoy = {}

chain = random.choice(chain_list)
chain_name = struc.get_chains(chain)[0]
# print(chain_name)

pre3 = chain.copy()
pre4 = chain.copy()

mid = int(len(pre3)/2)

downstream = pre3[mid:]
axis = pre3[mid - 1].coord - pre3[mid - 2].coord
support = pre3[mid-1].coord
downstream = struc.rotate_about_axis(
    downstream,
    angle = torsion[0],
    axis = axis,
    support = support
)
struc.coord(pre3[mid:])[:] = struc.coord(downstream)
chain_decoy[f"{chain_name}_ex30"] = pre3

downstream = pre4[mid:]
axis = pre4[mid-1].coord - pre4[mid-2].coord
support = pre4[mid-1].coord
downstream = struc.rotate_about_axis(
    downstream,
    angle = torsion[1],
    axis = axis,
    support = support
)
struc.coord(pre4[mid:])[:] = struc.coord(downstream)
chain_decoy[f"{chain_name}_ex60"] = pre4
# print(chain_decoy)
# print(struc.get_chains(chain_list[0]))

chain_num_map = {
    "A" : 0,
    "B" : 1,
    "F" : 2,
    "G" : 3
}
alpha = list(chain_decoy.keys())[0]
alpha = alpha.strip().split("_")[0]
chain_num = chain_num_map[alpha]
# print(chain_num)

pre8 = chain_list.copy()
pre9 = chain_list.copy()
pre8[chain_num] = chain_decoy[f"{alpha}_ex30"]
pre9[chain_num] = chain_decoy[f"{alpha}_ex60"]

ex_list = [30, 60]
n = 0
for i in pre8, pre9:
    combined_chain = sum(i, struc.AtomArray(0))
    decoy_dict[f"5DK3_ex{ex_list[n]}"] = combined_chain
    n += 1


# test
# print(decoy_dict.keys())
# print(decoy_dict["5DK3_ex60"] == data_dict["5DK3"][0])
for name, array in decoy_dict.items():
    cif_file = pdbx.CIFFile()
    pdbx.set_structure(cif_file, array, data_block=f"decoy_{name}")
    cif_file.write(f"test_dataset/decoy_{name}.cif")