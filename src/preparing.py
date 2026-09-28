import biotite.structure.io.pdbx as pdbx    # 읽어오기     
import biotite.structure as struc
import torch                                # 텐서 변환, 노이즈 생성용 
from math import pi
import numpy as np
import random

###################
# 함수들



# Gaussian noise - Cartesian
def noise_gaussian(bb_atoms: struc.AtomArray, sigma: float) -> struc.AtomArray:

    pre = bb_atoms.copy()

    bb_coord = struc.coord(pre)
    bb_tensor = torch.tensor(bb_coord)

    noise = torch.randn_like(bb_tensor) * sigma
    noise_coord = bb_tensor + noise
    noise_coord = noise_coord.detach().cpu().numpy()

    struc.coord(pre)[:] = noise_coord

    return pre


# Global Perturbation - Torsion Angle
def noise_global_torsion(bb_atoms: struc.AtomArray, angle: float) -> struc.AtomArray:

    pre = bb_atoms.copy()

    bb_phi, bb_psi, bb_omega = struc.dihedral_backbone(pre)

    bb_phi_tensor = torch.tensor(bb_phi)
    bb_psi_tensor = torch.tensor(bb_psi)
    bb_omega_tensor = torch.tensor(bb_omega)

    tensor_list = [bb_phi_tensor, bb_psi_tensor, bb_omega_tensor]

    pre_decoy = []
    for tensor in tensor_list:
        noise = torch.randn_like(tensor) * angle
        noise = noise.detach().cpu().numpy().tolist()
        pre_decoy.append(noise)

    pertub = [i for j in zip(pre_decoy[0], pre_decoy[1], pre_decoy[2]) for i in j]
    pertub = pertub[1:-1]

    for i in range(len(pertub)):
        axis = pre[i+1].coord - pre[i].coord
        support = pre[i+1].coord
        downstream = pre[i+2:]
        downstream = struc.rotate_about_axis(
            downstream,
            angle = pertub[i],
            axis = axis,
            support = support
        )
        struc.coord(pre[i+2:])[:] = struc.coord(downstream)

    return pre


# Local Extreme Perturbation - Torsion Angle
def noise_local_torsion(bb_atoms: struc.AtomArray, torsion: float) -> struc.AtomArray:

    pre = bb_atoms.copy()

    mid = int(len(pre)/2)

    downstream = pre[mid:]
    axis = pre[mid - 1].coord - pre[mid - 2].coord
    support = pre[mid-1].coord
    downstream = struc.rotate_about_axis(
        downstream,
        angle = torsion,
        axis = axis,
        support = support
    )
    struc.coord(pre[mid:])[:] = struc.coord(downstream)

    return pre

#######################
# 실행 코드



# 체인 데이터 딕셔너리로 정리
id_list = ["1CRN", "1CLL", "5DK3"]
data_dict = {}
for id in id_list:    
    with open(f"raw_dataset/parsed_{id}.cif", "rt", encoding = "utf-8") as f:
        data_dict[f"{id}"] = pdbx.get_structure(pdbx.CIFFile.read(f), model = 1)


# 싱글체인 데이터 정리
id_single_chain = ["1CRN", "1CLL"]
single_dict = {}
for id in id_single_chain:
    single_dict[f"{id}"] = data_dict[f"{id}"]

# 멀티체인 데이터 정리
multi_chain = data_dict["5DK3"]
chain_list = list(struc.chain_iter(multi_chain))

# 그릇
decoy_dict = {}



# (1) Gaussian Coordiante Noise (sigma = 0.5, 1.0, 2.0, 5.0 angstrom)
sigmas = [0.5, 1.0, 2.0, 5.0]


# single chain
for id, bb_atoms in single_dict.items():
    for sigma in sigmas:
        decoy = noise_gaussian(bb_atoms, sigma)
        decoy_dict[f"{id}_{sigma}angstrom"] = decoy

# Multi chain
for sigma in sigmas:
    chain_decoy = []
    for chain in chain_list:
        decoy = noise_gaussian(chain, sigma)
        chain_decoy.append(decoy)
    pre = struc.concatenate(chain_decoy)
    decoy_dict[f"5DK3_{sigma}angstrom"] = pre



# (2) Global Perturbation (Gaussian, sigma = 5, 10 degree)
angles = [pi/36, pi/18]


# single chain
for id, bb_atoms in single_dict.items():
    for angle in angles:
        decoy = noise_global_torsion(bb_atoms, angle)
        deg = int(np.rad2deg(angle))
        decoy_dict[f"{id}_{deg}degree"] = decoy

# Multi Chain
for angle in angles:
    chain_decoy = []
    for chain in chain_list:
        decoy = noise_global_torsion(chain, angle)
        chain_decoy.append(decoy)
    pre = struc.concatenate(chain_decoy)
    deg = int(np.rad2deg(angle))
    decoy_dict[f"5DK3_{deg}degree"] = pre



# (3) Local Extreme Perturbation (대충 가운데 부근 한 곳에서 결합각 대폭 조정, 30도 60도.)
torsions = [pi/6, pi/3]


# single chain
for id, bb_atoms in single_dict.items():
    for torsion in torsions:
        decoy = noise_local_torsion(bb_atoms, torsion)
        deg = int(np.rad2deg(torsion))
        decoy_dict[f"{id}_ex{deg + 1}"] = decoy


# Multi Chain
for torsion in torsions:
    chain = random.choice(chain_list)
    decoy = noise_local_torsion(chain, torsion)

    chain_decoy = chain_list.copy()
    idx = chain_decoy.index(chain)
    chain_decoy[idx] = decoy

    pre = struc.concatenate(chain_decoy)
    deg = int(np.rad2deg(torsion))
    decoy_dict[f"5DK3_ex{deg + 1}"] = pre



# 출력
for name, array in decoy_dict.items():
    cif_file = pdbx.CIFFile()
    pdbx.set_structure(cif_file, array, data_block=f"decoy_{name}")
    cif_file.write(f"test_dataset/decoy_{name}.cif")

########################

# 추후 멀티체인과 싱글체인 동시에 다룰 수 있도록 2차 리팩토링 가능성