import biotite.structure as struc
import numpy as np

# Metric 계산기 (rmsd, lddt, tmscore)
def calc_metrics(ans: struc.AtomArray, decoy: struc.AtomArray) -> tuple[float, float, float]:

    rmsd = f'{struc.rmspd(ans, decoy):.3f}'
    lddt = f'{struc.lddt(ans, decoy):.3f}'

    sup_pos = struc.superimpose(ans, decoy)
    n_range = np.arange(len(ans))
    tm = f'{struc.tm_score(ans, sup_pos[0], n_range, n_range):.3f}'

    return (rmsd, lddt, tm)
