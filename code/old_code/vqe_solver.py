from qiskit import QuantumCircuit, Aer, execute
from qiskit.circuit import Parameter
import numpy as np

from qiskit_nature.second_q.drivers import PySCFDriver

from qiskit_nature.second_q.mappers import JordanWignerMapper

from qiskit.algorithms.optimizers import SLSQP
from qiskit.algorithms.minimum_eigensolvers import VQE
from qiskit.primitives import Estimator
from qiskit_nature.second_q.algorithms import GroundStateEigensolver


def solve(mol, ansatz):

    mapper = JordanWignerMapper()

    vqe = VQE(Estimator(), ansatz, SLSQP())

    vqe.initial_point = np.zeros(ansatz.num_parameters)
    
    solver = GroundStateEigensolver(mapper, vqe)
    result = solver.solve(mol)

    return result.total_energies[0]