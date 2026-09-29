from .io import bb_parser, ca_parser, load_structure
from .decoy import noise_gaussian, noise_global_torsion, noise_local_torsion
from .metrics import calc_metrics

from collections import defaultdict
import biotite.structure as struc
from math import pi
import numpy as np
import random
import pandas as pd

def run_pipeline():

    ###
    # 참조 데이터셋
    id_list = ["1CRN", "1CLL", "5DK3"]
    sigmas = [0.5, 1.0, 2.0, 5.0]
    angles = [pi/36, pi/18]
    torsions = [pi/6, pi/3]

    ans_dict_bb = {}
    ans_dict_ca = {}
    decoy_dict = {}
    decoy_bb = {}
    decoy_ca = {}
    rmsd_bb_dict = defaultdict(list)
    lddt_bb_dict = defaultdict(list)
    tm_dict = defaultdict(list)
    rmsd_ca_dict = defaultdict(list)
    lddt_ca_dict = defaultdict(list)
    ###


    ###
    # cif.gz 파일 읽어와서, BackBone AtomArray, Ca AtomArray 정답 데이터셋 만들기
    for id in id_list:

        # cif.gz -> AtomArray
        atoms = load_structure(id)

        # AtomArray(All-Atom) -> AtomArray(Backbone Atoms)
        bb_atoms = bb_parser(atoms)

        # AtomArray(All-Atom) -> AtomArray(Backbone Atoms)
        ca_atoms = ca_parser(atoms)

        ans_dict_bb[id] = bb_atoms
        ans_dict_ca[id] = ca_atoms
    ###


    ###
    # 싱글체인, 멀티체인 분리

    # 싱글체인 모데이터
    multichain = dict(list(ans_dict_bb.items())[2:])
    single_dict = dict(list(ans_dict_bb.items())[:2])

    # 멀티체인
    multi_chain = multichain["5DK3"]
    chain_list = list(struc.chain_iter(multi_chain))
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


    ###
    # Decoy backbone, alphacarbon parsing

    # 백본은 그대로 가져오기
    decoy_bb = decoy_dict

    # 디코이 백본에서 알파카본만 남긴 새로운 딕셔너리
    for id, decoy in decoy_bb.items():
        ca_decoy = ca_parser(decoy)
        decoy_ca[id] = ca_decoy
    ###

    
    ###
    # Metric 계산하기

    # backbone
    for id in id_list:
        ans = ans_dict_bb[id]
        for key, decoy in decoy_bb.items():
            target = key.strip().split('_')[0]
            if target == id:
                rmsd, lddt, _ = calc_metrics(ans, decoy)
                rmsd_bb_dict[id].append(rmsd)
                lddt_bb_dict[id].append(lddt)

    # alpha carbon
    for id in id_list:
        ans = ans_dict_ca[id]
        for key, decoy in decoy_ca.items():
            target = key.strip().split('_')[0]
            if target == id:
                rmsd, lddt, tm = calc_metrics(ans, decoy)
                rmsd_ca_dict[id].append(rmsd)
                lddt_ca_dict[id].append(lddt)
                tm_dict[id].append(tm)
    ###

    
    ###
    # .csv 파일로 출력

    # pandas DataFrame 조건 설정
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

    # id 별 데이터 분리, 변환
    for id in id_list:
        data = np.column_stack((rmsd_bb_dict[id], rmsd_ca_dict[id], tm_dict[id], lddt_bb_dict[id], lddt_ca_dict[id]))
        df = pd.DataFrame(data, index= indices, columns = columns)
        df.to_csv(f'result/{id}_result.csv', index= True, encoding='utf-8-sig')
    ###




