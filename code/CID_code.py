import numpy as np
np.set_printoptions(precision=5, linewidth=200, suppress=True)
import psi4

from old_code.helper_CI import Determinant, HamiltonianGenerator
from itertools import combinations

def calculate_CID(mol, scf_e, wfn):

    # Memory for Psi4 in GB
    psi4.core.set_output_file('output.dat', False)

    # Memory for numpy in GB
    numpy_memory = 2

    # Grab data from wavfunction class
    C = wfn.Ca()
    ndocc = wfn.doccpi()[0]
    nmo = wfn.nmo()
    nvirt = nmo - ndocc
    nDet_S = ndocc * nvirt * 2

    # Integral generation from Psi4's MintsHelper
    mints = psi4.core.MintsHelper(wfn.basisset())
    H = np.asarray(mints.ao_kinetic()) + np.asarray(mints.ao_potential())

    #Make spin-orbital MO
    MO = np.asarray(mints.mo_spin_eri(C, C))

    # Update H, transform to MO basis and tile for alpha/beta spin
    H = np.einsum('uj,vi,uv', C, C, H)
    H = np.repeat(H, 2, axis=0)
    H = np.repeat(H, 2, axis=1)

    # Make H block diagonal
    spin_ind = np.arange(H.shape[0], dtype=np.int) % 2
    H *= (spin_ind.reshape(-1, 1) == spin_ind)


    occList = [i for i in range(ndocc)]
    det_ref = Determinant(alphaObtList=occList, betaObtList=occList)
    detList = det_ref.generateDoubleExcitationsOfDet(nmo)
    detList.append(det_ref)

    Hamiltonian_generator = HamiltonianGenerator(H, MO)
    Hamiltonian_matrix = Hamiltonian_generator.generateMatrix(detList)

    e_cid, wavefunctions = np.linalg.eigh(Hamiltonian_matrix)

    hartree2eV = 27.211

    # for i in range(1, len(e_cis)):
    #     excit_e = e_cid[i] + mol.nuclear_repulsion_energy() - scf_e

    return e_cid[0] + mol.nuclear_repulsion_energy() # scf_e