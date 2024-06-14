from qiskit import QuantumCircuit, Aer, BasicAer, execute, transpile
from qiskit.circuit import Parameter
from qiskit_nature.second_q.circuit.library import HartreeFock
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit.algorithms.optimizers import SLSQP
from qiskit.algorithms.minimum_eigensolvers import VQE

from qiskit.primitives import Estimator, BackendEstimator
from qiskit_aer.primitives import Estimator as AerEstimator

from qiskit_nature.second_q.algorithms import GroundStateEigensolver
from qiskit.circuit.library import RYGate
from qiskit_aer import noise, AerSimulator
from qiskit_aer.noise import NoiseModel
from qiskit.providers.fake_provider import FakeVigo
from qiskit.utils import algorithm_globals
# seed = 170
# algorithm_globals.random_seed = seed

from itertools import combinations, product

import numpy as np


class CustomAnsatz:

    def __init__(self, mol, singles=False, pair_doubles=False, doubles=False, k=1):
        
        self.mol = mol

        self.circuit = self.create_circuit()
        self.circuit = self.add_reference_state()

        if singles:
            for exc in self.generate_CCS_excitations()*k:
                self.add_excitation(exc[0], exc[1])

        if pair_doubles:
            for exc in self.generate_pCCD_excitations()*k:
                self.add_excitation(exc[0], exc[1])

        if doubles:
            for exc in self.generate_CCD_excitations()*k:
                self.add_excitation(exc[0], exc[1])        
        
        return

    def create_circuit(self):
        mol = self.mol
        return QuantumCircuit(mol.num_spin_orbitals)

    def generate_CCS_excitations(self):
        """
        function generating single excitations
        """
        mol = self.mol
        # list of occupations
        occ = list(mol.orbital_occupations) + list(mol.orbital_occupations_b)

        lst = [[[[i],[j]] for j,o in enumerate(occ) if o==0] for i,s in enumerate(occ) if s==1]

        full_singles = []
        for lst in lst:
            full_singles += lst

        return full_singles

    def generate_pCCD_excitations(self):
        """
        function generating pair double excitations
        """
        mol = self.mol
        # list of occupations
        occ = mol.orbital_occupations + mol.orbital_occupations_b

        # pair doubles, 
        pair_doubles = [[[[i,i+mol.num_spatial_orbitals],[j,j+mol.num_spatial_orbitals]] 
                    for j,o in enumerate(occ) if o==0] for i,s in enumerate(occ) if s==2]
        
        # concatenating the lists of spatial orbital excitations into one
        full_pair_doubles = []
        for lst in pair_doubles:
            full_pair_doubles += lst

        return full_pair_doubles
    
    def generate_CCD_excitations(self):
        """
        function generating all double excitations
        """
        mol = self.mol

        occ = list(mol.orbital_occupations) + list(mol.orbital_occupations_b)

        # all occupied
        occ_list = [ i for i,o in enumerate(occ) if o==1]

        # all virtual
        vir_list = [ i for i,o in enumerate(occ) if o==0]

        occ_pairs = list([list(c) for c in combinations(occ_list, 2)])
        vir_pairs = list([list(c) for c in combinations(vir_list, 2)])

        cartesian_product = list(product(occ_pairs, vir_pairs))

        return cartesian_product


    def add_reference_state(self):
        mol = self.mol
        circuit = self.circuit
        # function adding the HF reference state to a circuit
        mapper = JordanWignerMapper()
        initial_state = HartreeFock( mol.num_spatial_orbitals, mol.num_particles, mapper)

        for q, occ in enumerate(initial_state._bitstr):
            if occ:
                circuit.x(q)

        return circuit

    def add_compressed_single(self, from_ind, to_ind):
        # ! if the gates are malfunctioning a looping cnot ladder might
        # need to be added, going from excitations[1] +1 to excitations[0] -1 (DOI:https://doi.org/10.1103/PhysRevA.102.062612)
        circuit = self.circuit

        excitations = from_ind + to_ind
        excitations.sort()

        # rotation angle
        angle = Parameter(f'θ_{circuit.num_parameters}')

        circuit.cnot(excitations[0], excitations[1])
        circuit.cry(angle, excitations[1], excitations[0])
        circuit.cnot(excitations[0], excitations[1])
        return circuit

    # def add_compressed_double(self, from_ind, to_ind):
    #     # compressed double, modified with looping cnot ladders, currently depricated (DOI:https://doi.org/10.1103/PhysRevA.102.062612)

    #     circuit = self.circuit

    #     from_ind.sort()
    #     to_ind.sort()
    #     excitations = from_ind + to_ind

    #     # rotation angle
    #     angle = Parameter(f'θ_{circuit.num_parameters}')

    #     circuit.cx(excitations[0], excitations[1])
    #     circuit.cx(excitations[2], excitations[3])
    #     circuit.cx(excitations[0], excitations[2])

    #     control1 =  RYGate(angle).control(3, None, '010')
    #     circuit.append(control1, [excitations[1],excitations[2],excitations[3],excitations[0]])

    #     circuit.cx(excitations[0], excitations[2])
    #     circuit.cx(excitations[0], excitations[1])
    #     circuit.cx(excitations[2], excitations[3])

    #     return circuit

    # def add_compressed_double(self, from_ind, to_ind):
    #     # compressed double

    #     circuit = self.circuit

    #     from_ind.sort()
    #     to_ind.sort()
    #     excitations = from_ind + to_ind

    #     index_pairs = list(range(circuit.num_qubits))
    #     index_pairs += index_pairs
    #     index_pairs

    #     # pairs for looping cnot ladders
    #     loop_1 = index_pairs[(excitations[1]+1) : (int(len(index_pairs)/2)+excitations[0])]
    #     loop_2 = index_pairs[(excitations[3]+1) : (int(len(index_pairs)/2)+excitations[2])]

    #     loop_1_s =  loop_1[1:] + [loop_1[0]]
    #     loop_2_s =  loop_2[1:] + [loop_2[0]]

    #     loop_1 = loop_1[:-1]
    #     loop_1_s = loop_1_s[:-1]

    #     loop_2 = loop_2[:-1]
    #     loop_2_s = loop_2_s[:-1]

    #     # cnot ladder 1
    #     for o, n in zip(loop_1, loop_1_s):
    #         circuit.cx(o, n)

    #     # cnot ladder 2
    #     for o, n in zip(loop_2, loop_2_s):
    #         circuit.cx(o, n)

    #     # rotation angle
    #     angle = Parameter(f'θ_{circuit.num_parameters}')

    #     circuit.cx(excitations[0], excitations[1])
    #     circuit.cx(excitations[2], excitations[3])
    #     circuit.cx(excitations[0], excitations[2])

    #     control1 =  RYGate(angle).control(3, None, '010')
    #     circuit.append(control1, [excitations[1],excitations[2],excitations[3],excitations[0]])

    #     circuit.cx(excitations[0], excitations[2])
    #     circuit.cx(excitations[0], excitations[1])
    #     circuit.cx(excitations[2], excitations[3])

    #     # cnot ladder 1
    #     for o, n in zip(loop_1[::-1], loop_1_s[::-1]):
    #         circuit.cx(o, n)

    #     # cnot ladder 2
    #     for o, n in zip(loop_2[::-1], loop_2_s[::-1]):
    #         circuit.cx(o, n)

    #     return circuit

    def add_excitation(self, from_ind, to_ind, test_operator=False):

        circuit = self.circuit

        from_ind.sort()
        to_ind.sort()
        excitations = from_ind + to_ind

        if len(excitations) == 4:
            if test_operator:
                operator = ["+XXYX"]
            else:
                operator = ["+XYXX", "+XXXY", "+YYYX", "+YYXY", "-XXYX", "-YXXX", "-YXYY", "-YYYX"] # correct

        elif len(excitations) == 2:
            excitations.sort()
            operator = ["+YX", "-XY"]

        else:
            print("unsupported excitation")

        # rotation angle
        angle = Parameter(f'θ_{circuit.num_parameters}')


        for ladder in operator:
            # in single excitations the cnot ladders are between the target qubits

            # changing bases
            for q, op in zip(excitations, ladder[1:]):
                if op == "X":
                    circuit.h(q)
                elif op == "Y":
                    circuit.rx(np.pi/2, q)
            
            if len(operator) == 2:
                # cnot ladder
                for o in range(excitations[0], excitations[-1]):
                    circuit.cx(o, o+1)

                # rotation
                sign = -1 if ladder[0] == "-" else 1
                circuit.rz( sign * angle/len(operator), excitations[-1])

                # cnot ladder
                for o in range(excitations[-1]-1, excitations[0]-1, -1):
                    circuit.cx(o, o+1)

            if len(operator) == 8:
                # in double excitations the cnot ladders have to go through all qubits

                # cnot ladder
                for o in range(0, circuit.num_qubits-1):
                    circuit.cx(o, o+1)

                # rotation
                sign = -1 if ladder[0] == "-" else 1
                circuit.rz( sign * angle/len(operator), circuit.num_qubits-1)

                # cnot ladder
                for o in range(circuit.num_qubits-2, -1, -1):
                    circuit.cx(o, o+1)

            # changing back bases
            for q, op in zip(excitations, ladder[1:]):
                if op == "X":
                    circuit.h(q)
                elif op == "Y":
                    circuit.rx(-np.pi/2, q)

        return circuit
    


def solve(mol, ansatz):
    """
    function performing vqe using the provided ansatz
    """
    mapper = JordanWignerMapper()

    # noiseless_estimator = AerEstimator(run_options={"seed": seed, "shots": 100000}) # the AerEstimator does not seem to work with the VQE
    noiseless_estimator = Estimator()

    vqe = VQE( noiseless_estimator, ansatz, SLSQP())

    vqe.initial_point = np.zeros(ansatz.num_parameters)
    
    solver = GroundStateEigensolver(mapper, vqe)
    result = solver.solve(mol)

    return result.total_energies[0]

def solve_noisy(mol, ansatz):
    """
    function performing vqe on a noisy backend using the provided ansatz, the AerEstimator does not seem to work with the VQE
    """
    # using a noisy backend, based on data from IBM
    device = FakeVigo()
    coupling_map = device.configuration().coupling_map
    noise_model = NoiseModel.from_backend(device)

    # noisy_estimator = AerEstimator(
    # backend_options={
    #     "method": "density_matrix",
    #     "coupling_map": coupling_map,
    #     "noise_model": noise_model,
    # },
    # run_options={"seed": seed, "shots": 10240},
    # transpile_options={"seed_transpiler": seed},
    # )

    noisy_estimator = BackendEstimator(AerSimulator(noise_model=noise_model))
    # noisy_estimator = BackendEstimator(AerSimulator())

    mapper = JordanWignerMapper()

    vqe = VQE(noisy_estimator, ansatz, SLSQP())

    vqe.initial_point = np.zeros(ansatz.num_parameters)
    
    solver = GroundStateEigensolver(mapper, vqe)
    result = solver.solve(mol)

    return result.total_energies[0]


from qiskit.quantum_info import SparsePauliOp

def solve_noisy_2(r):
    import tequila as tq
    # define the active space
    # active_orbitals = {"A1":[1], "B1":[0]}
    samples = 1000000
    backend = "qiskit"
    device = "fake_rome"
    # define the molecule
    molecule = tq.chemistry.Molecule(geometry = f"H 0.0 0.0 0.0\nH 0.0 0.0 {r}",
                                    basis_set="sto-3g",
                                    #  active_orbitals=active_orbitals,
                                    transformation="jordan_wigner")

    fci = molecule.compute_energy("fci")

    H = molecule.make_hamiltonian()

    # Toy circuit (no deeper meaning)
    U = tq.gates.Ry(angle="a", target=0)
    U += tq.gates.X(target=1, control=0)
    # E = tq.ExpectationValue(H=H, U=U, optimize_measurements=True)

    # vqe = tq.minimize(method="cobyla", objective=E, initial_values=0.0)
    # noisy_vqe = tq.minimize(method="cobyla", objective=E, samples=samples, backend=backend, device=device, initial_values=0.0)

    # The same with UpCCGSD
    UpCCGSD = molecule.make_upccgsd_ansatz(name="UpCCGSD")
    E2 = tq.ExpectationValue(H=H, U=UpCCGSD, optimize_measurements=True)
    # ucc = tq.minimize(method="cobyla", objective=E2, initial_values=0.0)
    noisy_ucc = tq.minimize(method="cobyla", objective=E2, samples=samples,  backend=backend, device=device, initial_values=0.0)

    # print("VQE         = {:2.8f}".format(min(vqe.history.energies)))
    # print("VQE (noisy) = {:2.8f}".format(min(noisy_vqe.history.energies)))
    # print("UCC         = {:2.8f}".format(min(ucc.history.energies)))
    print("UCC (noisy) = {:2.8f}".format(min(noisy_ucc.history.energies)))

    # Retrieve the energies from the optimization history
    # energies = noisy_ucc.history.energies
    
    # # Compute the mean and standard deviation of the energies
    # mean_energy = np.mean(energies)
    # std_energy = np.std(energies)

    return noisy_ucc