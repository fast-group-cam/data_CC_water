from ase import io
from ase.neighborlist import neighbor_list
from pymatgen.analysis.local_env import JmolNN
import itertools
import numpy as np
from typing import List, Dict, Tuple


def extract_water_cluster(
    filename: str,
    output_filename: str = "test_cluster.xyz",
    cutoff_radius: float = 5.0,
    supercell_min_length: float = 20.0,
    random_seed: int | None = None,
) -> None:
    """
    Extracts a molecular cluster centered around a random O atom
    from a periodic water system, including nearby O and H atoms.

    Parameters
    ----------
    filename : str
        Path to the input structure file (e.g., 'water_example.extxyz').
    output_filename : str, optional
        Path for the output cluster file. Default is 'test_cluster.xyz'.
    cutoff_radius : float, optional
        Cutoff distance (Å) to include neighboring O atoms. Default is 5.0 Å.
    supercell_min_length : float, optional
        Minimum length (Å) of the constructed supercell in each dimension. Default is 20.0 Å.
    random_seed : int | None, optional
        Optional random seed for reproducibility. Default is None.

    Returns
    -------
    None
        Writes the extracted cluster structure to `output_filename`.
    """
    if random_seed is not None:
        np.random.seed(random_seed)

    # Load structure and build supercell
    unit_cell = io.read(filename)
    unit_cell_lengths = unit_cell.cell.cellpar()[:3]
    supercell = unit_cell * tuple(
        int(supercell_min_length / l) + 1 for l in unit_cell_lengths
    )

    # Identify atomic species and generate bond pairs
    atomic_species = list(set(supercell.get_chemical_symbols()))
    bond_pairs = list(itertools.combinations(atomic_species, 2)) + [
        (x, x) for x in atomic_species
    ]

    # Determine cutoff distances for each pair using Jmol radii
    cutoff_dict: Dict[Tuple[str, str], float] = {
        pair: JmolNN().get_max_bond_distance(pair[0], pair[1]) for pair in bond_pairs
    }

    # Get nearest neighbor info (ASE neighbor_list)
    i_list, j_list, S_list = neighbor_list("ijS", supercell, cutoff=cutoff_dict)

    ase_nn_info: List[Dict[int, Dict[str, Tuple[int, Tuple[int, int, int]]]]] = [
        {} for _ in range(len(supercell))
    ]
    for i, j, S in zip(i_list, j_list, S_list):
        ase_nn_info[i][len(ase_nn_info[i])] = {"site_index": j, "image": tuple(S)}

    # Choose a random O atom as cluster center
    O_indices = [i for i, atom in enumerate(supercell) if atom.symbol == "O"]
    center_O_idx = np.random.choice(O_indices)

    # Translate center O atom to the geometric center of the cell
    cell = supercell.get_cell()
    cell_center = 0.5 * (cell[0] + cell[1] + cell[2])
    supercell.translate(cell_center - supercell[center_O_idx].position)
    supercell.wrap()

    # Collect O atoms within cutoff radius
    cluster_O_indices = [
        i
        for i in O_indices
        if np.linalg.norm(supercell[i].position - supercell[center_O_idx].position)
        <= cutoff_radius
    ]

    cluster_H_indices: List[int] = []
    for O_idx in cluster_O_indices:
        for nn in ase_nn_info[O_idx].values():
            neighbor_idx = nn["site_index"]
            if supercell[neighbor_idx].symbol == "H":
                cluster_H_indices.append(neighbor_idx)

    # Combine O and H atoms to form the final cluster
    cluster_indices = sorted(set(cluster_O_indices + cluster_H_indices))
    cluster = supercell[cluster_indices].copy()

    # Center cluster and remove PBC
    cluster.translate(-cell_center)
    cluster.pbc = False
    cluster.set_cell([0, 0, 0])

    # Write output
    io.write(output_filename, cluster)
    print(f"Cluster with {len(cluster)} atoms written to {output_filename}.")