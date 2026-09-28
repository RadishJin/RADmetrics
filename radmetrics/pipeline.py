from .io import bb_parser, ca_parser, load_structure
from .decoy import noise_gaussian, noise_global_torsion, noise_local_torsion

import biotite.structure as struc
from math import pi
import numpy as np
import random

def run_pipeline():

    ###
    # 참조 데이터셋
    id_list = ["1CRN", "1CLL", "5DK3"]
    sigmas = [0.5, 1.0, 2.0, 5.0]
    angles = [pi/36, pi/18]
    torsions = [pi/6, pi/3]
    ###


    ###
    # cif.gz 파일 읽어와서, BackBone AtomArray, Ca AtomArray 정답 데이터셋 만들기
    ans_dict_bb = {}
    ans_dict_ca = {}
    for id in id_list:

        # cif.gz -> AtomArray
        atoms = load_structure(id)

        # AtomArray(All-Atom) -> AtomArray(Backbone Atoms)
        bb_atoms = bb_parser(atoms)

        # AtomArray(All-Atom) -> AtomArray(Backbone Atoms)
        ca_atoms = ca_parser(atoms)

        ans_dict_bb[f"{id}"] = bb_atoms
        ans_dict_ca[f"{id}"] = ca_atoms
    ###


    ###
    # 싱글체인, 멀티체인 분리

    # 싱글체인 모데이터
    multichain = dict(list(ans_dict_bb.items()[:2]))
    single_dict = dict(list(ans_dict_bb.items()[2:]))
    #이거검증필요오늘여기까지

    # 멀티체인
    chain_list = list(struc.chain_iter(multichain))
    ###


    ###
    # 그릇
    decoy_dict = {}
    ###


    ###
    # Global Gaussian Noise, Cartesian(0.5, 1, 2, 5 angstrom)

    # Single Chain
    for id, bb_atoms in single_dict.items():
        for sigma in sigmas:
            decoy = noise_gaussian(bb_atoms, sigma)
            decoy_dict[f"{id}_{sigma}angstrom"] = decoy

    # Multi Chain
    for sigma in sigmas:
        chain_decoy = []
        for chain in chain_list:
            decoy = noise_gaussian(chain, sigma)
            chain_decoy.append(decoy)
        pre = struc.concatenate(chain_decoy)
        decoy_dict[f"5DK3_{sigma}angstrom"] = pre
    ###


    ###
    # Global Gaussian Noise, Torsion Angle(5, 10 degree)

    # Single Chain
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
    ###


    ###
    # Local Extreme Perturbation, Torsion Angle(30, 60 degree)

    # Single Chain
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
    ###

    print(decoy_dict.keys())
    

    



