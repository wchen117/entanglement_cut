import itertools
import multiprocessing as mp
import os
from collections import defaultdict
from warnings import filterwarnings

import numpy as np
from MDAnalysis import Universe
from numba import njit
from scipy.spatial.distance import pdist, squareform
from topoly import lasso_type

filterwarnings("ignore")


@njit(fastmath=True)
def helper_dot(Runit: np.ndarray, dR_cross: np.ndarray) -> list:

    return [np.dot(x, y) for x, y in zip(Runit, dR_cross)]


def point_60_rounding(num: float) -> float:

    """
    This function perform rounding to the nearest 0.6. 

    Ex:  0.61 rounds up to 1 but 0.59 rounds down to 0

    """
    if len(str(num).split("e")) == 2:
        num = 0.0

    if num % 1 >= 0.60:
        rounded_num = round(num)

    else:
        rounded_num = int(str(num).split(".")[0])

    return rounded_num


def get_entanglements(
    coor: np.ndarray,
    l: int,
    termini_threshold: list,
    chain: str,
    resids: np.ndarray,
    uniprotid: str,
) -> dict:

    """
    Find crossing residues that pierces the loop formed by native contact i and j

    """
    Nterm_thresh = termini_threshold[0]
    Cterm_thresh = termini_threshold[1]

    # make native contact contact map
    dist_matrix = squareform(pdist(coor))
    native_cmap = np.where(
        dist_matrix <= 8.0, 1, 0
    )  # if true then 1 will appear otherwise zero
    native_cmap = np.triu(
        native_cmap, k=4
    )  # element below the 4th diagonal starting from middle are all zeros; # protein contact map

    nc_indexs = np.stack(np.nonzero(native_cmap)).T  # stack indices based on rows

    # make R coordinate and gradient of R length N-1
    range_l = np.arange(0, l - 1)
    range_next_l = np.arange(1, l)

    coor = coor.astype(np.float32)
    R = 0.5 * (coor[range_l] + coor[range_next_l])
    dR = coor[range_next_l] - coor[range_l]

    # make dRcross matrix
    pair_array = np.asarray(
        list(itertools.product(dR, dR))
    )  # combination of elements within array

    x = pair_array[:, 0, :]
    y = pair_array[:, 1, :]

    dR_cross = np.cross(x, y)

    # make Rnorm matrix
    pair_array = np.asarray(list(itertools.product(R, R)))
    diff = pair_array[:, 0, :] - pair_array[:, 1, :]
    diff = diff.astype(np.float32)

    Runit = diff / np.linalg.norm(diff, axis=1)[:, None] ** 3
    Runit = Runit.astype(np.float32)

    # make final dot matrix
    dot_matrix = helper_dot(Runit, dR_cross)
    dot_matrix = np.asarray(dot_matrix)
    dot_matrix = dot_matrix.reshape((l - 1, l - 1))

    nc_gdict = {}

    for i, j in nc_indexs:

        loop_range = np.arange(i, j)
        nterm_range = np.arange(Nterm_thresh, i - 5)  # buffer of 6
        cterm_range = np.arange(j + 6, l - (Cterm_thresh + 1))

        gn_pairs_array = np.fromiter(
            itertools.chain(*itertools.product(nterm_range, loop_range)), int
        ).reshape(-1, 2)
        gc_pairs_array = np.fromiter(
            itertools.chain(*itertools.product(loop_range, cterm_range)), int
        ).reshape(-1, 2)

        if gn_pairs_array.size != 0:
            gn_vals = dot_matrix[gn_pairs_array[:, 0], gn_pairs_array[:, 1]]
            gn_vals = gn_vals[~np.isnan(gn_vals)]  # get rid of nan values
            gn_val = np.sum(gn_vals) / (4.0 * np.pi)
        else:
            gn_val = 0

        if gc_pairs_array.size != 0:
            gc_vals = dot_matrix[gc_pairs_array[:, 0], gc_pairs_array[:, 1]]
            gc_vals = gc_vals[~np.isnan(gc_vals)]  # get rid of nan values
            gc_val = np.sum(gc_vals) / (4.0 * np.pi)

        else:
            gc_val = 0

        rounded_gc_val = point_60_rounding(np.float64(abs(gc_val)))
        rounded_gn_val = point_60_rounding(np.float64(abs(gn_val)))

        total_link = rounded_gn_val + rounded_gc_val

        if np.abs(rounded_gn_val) >= 1 or np.abs(rounded_gc_val) >= 1:
            nc_gdict[(int(i), int(j))] = (
                gn_val,
                gc_val,
                rounded_gn_val,
                rounded_gc_val,
            )

    nc_data, missing_residues = eliminate_missing_residues(
        nc_gdict, resids, Nterm_thresh, Cterm_thresh
    )

    entangled_res = find_crossing(coor.tolist(), nc_data, chain, resids)

    return entangled_res, missing_residues


def eliminate_missing_residues(
    native_contacts: dict, resids: np.ndarray, Nterm_thresh: int, Cterm_thresh: int
):

    """
    Eliminate native contacts that has missing residues in loop and/or threading segment 
    based on the partial linking number

    Capture and output missing residues in protein
    
    """
    check_all_resids = np.arange(resids[0], resids[-1] + 1)

    missing_residues = np.setdiff1d(check_all_resids, resids)

    check = set()

    for ij, values in native_contacts.items():

        native_i = resids[ij[0]]

        native_j = resids[ij[1]]

        rounded_gn = values[-2]

        rounded_gc = values[-1]

        check_loop = np.arange(native_i, native_j + 1)
        missing_res_loop = np.setdiff1d(check_loop, resids[ij[0] : ij[1] + 1])

        if abs(rounded_gn) >= 1:

            check_resids_N = np.arange(
                resids[Nterm_thresh], native_i - 5
            )  # the N terminus to i - 6
            missing_res_N = np.setdiff1d(
                check_resids_N, resids[Nterm_thresh : ij[0] - 5]
            )

            if missing_res_N.size or missing_res_loop.size:

                native_contacts[ij] = None

        if abs(rounded_gc) >= 1:

            check_resids_C = np.arange(
                native_j + 6, resids[-Cterm_thresh]
            )  # the j + 6 to C terminus - 5
            missing_res_C = np.setdiff1d(
                check_resids_C, resids[ij[1] + 6 : -Cterm_thresh]
            )

            if missing_res_C.size or missing_res_loop.size:

                native_contacts[ij] = None

    return native_contacts, missing_residues


def find_crossing(coor: np.ndarray, nc_data: dict, chain: str, resids: np.ndarray):

    entangled_res = {}

    native_contacts = [[ij[0], ij[1]] for ij, values in nc_data.items() if values]

    # reduction:
    # 1. each crossing must be 10 residues apart [default]
    # 2. first crossing should be at least 6 residues from the loop
    # 3. first crossing should be at least 5 residues from the closest termini

    data = lasso_type(
        coor,
        loop_indices=native_contacts,
        pdb_chain=chain,
        more_info=True,
        precision=0,
        density=0,
        min_dist=[10, 6, 5],
    )
    # high precision, low denisty

    for native_contact in native_contacts:

        crossings = []

        native_contact = tuple(native_contact)

        if abs(nc_data[native_contact][-2]) >= 1:  # if rounded_gn >= 1

            crossingN = [
                resids[int(cr[1:])] for cr in data[native_contact]["crossingsN"]
            ]

            crossings += crossingN

        if abs(nc_data[native_contact][-1]) >= 1:  # if rounded_gc >= 1

            crossingC = [
                resids[int(cr[1:])] for cr in data[native_contact]["crossingsC"]
            ]

            crossings += crossingC

        gn = nc_data[native_contact][0]

        gc = nc_data[native_contact][1]

        ij_gN_gC = (resids[native_contact[0]], resids[native_contact[1]]) + (gn, gc)

        entangled_res[ij_gN_gC] = np.unique(crossings)

    return entangled_res


def gaussian_entanglement(ref_univ:Universe):

    pdb = "default"
    uniprotid = "default"

    ref_calphas_dup = ref_univ.select_atoms("name CA and protein")

    termini_threshold = [5, 5]

    entangled_info = ""
    for chain in set(ref_calphas_dup.chainIDs):
        print(chain)

        # get rid of duplicate resids (which will have duplicate coordinates)
        PDB_resids, i = np.unique(
            ref_calphas_dup.select_atoms(f"segid {chain}").resids, return_index=True
        )

        # no duplicate resids in this Universe obj
        ref_calphas = ref_calphas_dup.select_atoms(f"segid {chain}")[i]

        # x y z cooridnates of chain
        coor = ref_calphas.positions

        chain_res = PDB_resids.size

        if PDB_resids.size:

            ent_result, missing_residues = get_entanglements(
                coor, chain_res, termini_threshold, chain, PDB_resids, uniprotid
            )

            if ent_result:

                for ij_gN_gC, crossings in ent_result.items():
                    if crossings.size:
                        #  with open(f"unmapped_GE_alpha/{gene}_GE.txt", "a") as f:
                        #  f.write(f"Chain {chain} | ({ij_gN_gC[0]}, {ij_gN_gC[1]}, {crossings}) | {ij_gN_gC[2]} | {ij_gN_gC[3]} | No | No\n")
                        entangled_info += f"[{pdb}, {chain}, {ij_gN_gC[0]}, {ij_gN_gC[1]}, {crossings}, 99999]\n"

    return entangled_info


def get_pdb_files(directory):
    from os import listdir
    from os.path import isfile, join

    only_files = [
        join(directory, f) for f in listdir(directory) if isfile(join(directory, f))
    ]
    return only_files


def exit():
    import sys

    sys.exit()
    return


def date_time():
    from datetime import datetime

    now = datetime.now()
    return now.strftime("%d-%m-%Y %H:%M:%S")


output_folder = "../../Entanglement_Information_Without_TER_PDB/"


def write_entanglements(input_pdb_file):
    try:
        gene_name = input_pdb_file.split("/")[-1].split("_")[0]
        with open(f"{output_folder}{gene_name}.txt", "a") as output:

            log_file = open(f"{output_folder}logging.txt", "a")
            log_file.write(
                f"{date_time()}| Writing entanglements for {input_pdb_file}\n"
            )
            log_file.close()

            line = gaussian_entanglement(f"{input_pdb_file}")

            output.write(line)
    except Exception as e:
        error_file = open(f"{output_folder}error.txt", "a")
        error_file.write(
            f"Error entanglement calculation in {input_pdb_file}\n Error info:\n {e}\n\n"
        )
        error_file.close()
    return


if __name__ == "__main__":

    cores = 1  # change the core. I recommend at least 12
    result_obj = set()
    input_folder = "../../../Before_TER_PDBs/"
    file_list = get_pdb_files(input_folder)

    # Make sure you have pdbs in Before_TER_PDBs
    # in the name of gene.pdb. This will allow easy comparison between fold's data and my data.

    # gaussian_entanglement("P0A6L9.pdb") # if you want to run a single protein at a time

    with mp.get_context("spawn").Pool(cores) as p:
        chunk = len(file_list) // (cores ** 2) + 1
        results = p.map_async(write_entanglements, iterable=file_list, chunksize=chunk)
        #  append the mapresult obj
        result_obj.add(results)

        #  get the mapresult obj (only for async methods)
        for result in result_obj:
            result.get()


# chuck is variable that is used to feed each core.
# if we have 12 cores and have 4000 pdbs then each core will take 28 pdbs at time
# based on the formula: number of PDBs // (cores * cores) + 1
# once those 28 pdbs are completed then that core will take another 28 pdbs.
# This continues until mp is done
