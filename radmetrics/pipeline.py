from .io import bb_parser, load_structure

def run_pipeline():

    # 참조 데이터셋
    id_list = ["1CRN", "1CLL", "5DK3"]


    # cif.gz 파일 읽어와서, BackBone AtomArray만 남기기
    bb_dict = {}
    for id in id_list:
        # cif.gz -> AtomArray
        atoms = load_structure(id)
        # AtomArray(All-Atom) -> AtomArray(Backbone Atoms)
        bb_atoms = bb_parser(atoms)
        bb_dict[f"{id}"] = bb_atoms

    여기서부터
    # chain_list = {bb_dict.pop("5DK3")}
    # single_dict = bb_dict

    # print(chain_list, single_dict.keys())



