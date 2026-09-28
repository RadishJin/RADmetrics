import biotite.structure.io.pdbx as pdbx    # 읽어오기     
import biotite.structure as struc
import torch                                # 텐서 변환, 노이즈 생성용 


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
