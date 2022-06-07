from ast import parse
import MDAnalysis
import numpy as np
import gaussian_entanglements as ge
import clustering_updated as cluster
import subprocess
import sys
import re


def parse_stride_output(stride_output):
    """parse the stride output to find the secondary structure indices"""
    # use regex expression to find expression of "LOC NAME resname resid - resname resid -"
    exp = re.findall("LOC \D+ \d{1,5} \w \D+ \d{1,5} \w", stride_output)

    indices = np.zeros((len(exp), 2))
    #indices = []
    for idx, each_line in enumerate(exp):
        # for now, only alphahelix and strand is considered forbidden secondary
        # structure

        tmp = each_line.split()
        if (tmp[1] == 'AlphaHelix' or tmp[1] == 'Strand'):
            indices[idx, 0] = int(tmp[3])
            indices[idx, 1] = int(tmp[6])
            print(each_line)
   
    return indices
 
def stride_wrapper(pdb_file_name: str, pdb_length: int):
    """A python wrapper that calls the stride program 
        and obtain the secondary structure indices"""

    proc = subprocess.Popen("stride -o %s"%(pdb_file_name), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc.wait()
    (stdout, stderr) = proc.communicate()

    # if not normal return
    if proc.returncode != 0:
        # if error, print error message and quit
        sys.exit(stderr)
    else:
        #print("success")
        stride_output = stdout.decode('ascii')
        parsed_indices = parse_stride_output(stride_output)

        # hard code the size for now, TBD
        masked_array = np.ones(pdb_length)
        for each_pair in parsed_indices:
            # the resid starts from 1 
            ss_start_index = int(each_pair[0] - 1)
            ss_end_index = int(each_pair[1] - 1)
            masked_array[ss_start_index:ss_end_index] = 0
        return masked_array

def load_pdb(pdb_file_name: str):
    """
    load the pdb struct, input: file_name
    outputs: atoms, list of residue id, list of residue name
    """
    struct = MDAnalysis.Universe(pdb_file_name, format="PDB")

    return struct

def get_residue_resname_list(pdb_struct: MDAnalysis.Universe):

    """Given the pdb_struct loaded by load_pdb() function,
       return the list of atoms, residue id and residue name"""
    protein = pdb_struct.select_atoms("protein")
    residue_list = list(protein.residues.resids)
    resname_list = list(protein.residues.resnames)

    
    return protein, residue_list, resname_list

def ent_calc_wrapper(pdb_struct: MDAnalysis.Universe):
    """
    wrapper for entanglement wrapper calculation, returns only 
    the number of entanglement
    """
    # note the here input arguments for gaussian_entanglement have been modified
    # to by of type: str to type: Universe
    list_of_ents = ge.gaussian_entanglement(pdb_struct)

    if list_of_ents:
        ent_dict,ent=cluster.cluster_entanglements(list_of_ents, 55)
    return len(ent_dict)
    
def cut_protein(pdb_struct:MDAnalysis.Universe, site_index:int):
    """ cut the protein sturcture at a given index, returns a list of protein structure
    for input site_index = N, the cut is placed between N and N+1. Therefore the maximum site_index allowed is n_residue - 1"""

    if (site_index >= pdb_struct.residues.n_residues):
        exit("cut site index must be smaller than number of residues in protein ")
    else:
        # site_index is allowed
        protein, resid_list, resname_list = get_residue_resname_list(pdb_struct)

        # cut the protein at site_index

        rearanged_resid_list = []
        rearanged_resid_list.append(resid_list[:site_index])
        rearanged_resid_list.append(resid_list[site_index:])

        return rearanged_resid_list

def main():

    pdb_file_name = "./pdbs/native_chain_B.pdb"
    # the starting cut 
    n_cut = 1
    # load the structure into our system
    pdb_struct = load_pdb(pdb_file_name)
    protein, resid_list, resname_list = get_residue_resname_list(pdb_struct)
    n_residue = pdb_struct.residues.n_residues
    # the maximal number of cut sites available
    n_cut_sites = n_residue - 1

    # a generic python wrapper for stride program, note that 
    protein_mask = stride_wrapper(pdb_file_name, n_residue)
    print(protein_mask)

    # calculate the inital number of entanglement, start with this

    initial_ent_num = ent_calc_wrapper(pdb_struct)

    # we need a deep copy, not a reference
    tmp_cut_struct = pdb_struct.copy()
    # a place to keep track of the cutted sites

    already_cutted_sites = []

    for cut_site_index in range(n_cut_sites):
        # see if the cut_site_index cuts into the secondary structure identified by stride program
        if (protein_mask[cut_site_index] == 1 or cut_site_index not in already_cutted_sites):
            protein_pieces =  cut_protein(tmp_cut_struct, cut_site_index)




    


        
    


    
    return

if __name__ == '__main__':
    main();
    