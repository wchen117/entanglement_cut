##############################################################
# Cut.py, a simple program to iteratively cut protein to 
# eliminate entanglement
# Author: Weinan Chen
# May 2022
##############################################################

from ast import parse
import MDAnalysis
import numpy as np
import gaussian_entanglements as ge
import clustering_updated as cluster
import subprocess
import sys
import re
import itertools
from anytree import Node, RenderTree

class IndexNode(Node):
    """
       A derived class based on Node containing two additional quantities: 
       Index (resid) and Ent (num of entanglement/partial linking number)
    """
    def __init__(self, name, Ent=-1, Index=-1,  parent=None, children=None, **kwargs):
        super().__init__(name, parent, children, **kwargs)
        # set default value to -1 for error handling
        self.Ent = Ent
        self.Index = Index


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

def permute_connect(nested_list: list):
    """ take a nested list produced by cut_protein function, iteratve through all its permutations (A_N^N -1)
    """

    all_connections = []
    all_permutations = list(itertools.permutations(nested_list))
   
    # we should always discard the first one since it's unchanged
    for each_perm in all_permutations[1:]:
        left = []
        for each_piece in each_perm:
            left += each_piece
        all_connections.append(left)
            
    return all_connections  

def cut_protein(resid_list:list, site_indices:list, max_site: int):
    """ cut the protein sturcture at a given index.
        Note that the only input needed for cutting is resids/ list of resids, which is the first argument
        of this function.
        Please note that at this stage, the first input must be a complete resid_list, not a nested_list
        max_site must be smaller than n_residue - 1"""

    # check if input exceed max_site
    mask = np.array(site_indices) >= max_site
    # all elements must be 0 (site_index >= max_site is false for all input site indices)
    if (mask.prod() != 0):
        exit("input site_indices >= max_site is not allowed")
    else:
        # all indices are legit, sort the indices so that we start from small indices to large
        # essentially left to right
        site_indices = sorted(site_indices)
        # cutted list will be put into this nested list, with each element
        # being a piece of cutted resid_list
        nested_list = []
        # variable needed to reset the index of right array
        tmp = 0

        for each_site in site_indices:
            each_site = each_site - tmp
            left = resid_list[:each_site]
            nested_list.append(left)
            right = resid_list[each_site:]
            # not sure if we need a deep copy but just in case
            resid_list = right.copy()
            tmp = each_site + tmp
        del each_site, left, right, tmp
            
        nested_list.append(resid_list)


    return nested_list

def construct_pdb(permuated_list:list, original_list:list, new_pdb: MDAnalysis.Universe):
    """construct a pdb based on permutated resid_list, also need a copy of pdb_struct in the main func"""
    
    rewired_residue_positions=list()
    for resid in original_list:
        position_in_rewired_list = permuated_list.index(resid)
        resid_corresponding_to_rewired_list = original_list[position_in_rewired_list]
        rewired_residue_positions.append(resid_corresponding_to_rewired_list)

    rewired_array = np.array(rewired_residue_positions)
    new_pdb.residues.resids = rewired_array
    
    # return the new_pdb with modified resids
    return new_pdb

def main():

    pdb_file_name = "./pdbs/native_chain_B.pdb"
    # the starting number of cut sites to consider
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
    tmp_cut_list = resid_list
    # a place to keep track of the cutted sites

    # a place to hold cut site_indices, only hold 
    # optimal value for each n_cut value
    # in a tree implementation, this is replaced by one of the transversed path
    # for each n_cut iteration now, we need to iteratiev through all of the tranversed paths
    
    already_cutted_sites = []
    # the global tree structure that keep tracks of cutting indices

    root = Node("root")
 
    
    # for now limit the maximum concurrent cut to 5
    while (n_cut < 5):
        # need a place to hold temp indices
        local_index = -1
        # here is the place to tranverse the paths of the node, 

        for cut_site_index in range(n_cut_sites):
            # the first n elements of local_cut_sites are from previous 
            local_cut_sites = []
            # taking in account cut sites from the previous iterations 
            local_cut_sites = already_cutted_sites + [cut_site_index,]
            # see if the cut_site_index cuts into the secondary structure identified by stride program
            if (protein_mask[cut_site_index] == 1 and cut_site_index not in already_cutted_sites):
                nested_resid_pieces = cut_protein(tmp_cut_list, local_cut_sites, n_cut_sites)
                # each configuration is 
                all_configs = permute_connect(nested_resid_pieces)
                for each_config in all_configs:
                    # compute the number of entanglement (perhaps G. score) per configuration
                    # need a copy of the current pdb_struct to create new_pdb 
                    tmp_pdb = construct_pdb(each_config, resid_list, pdb_struct.copy())
                    tmp_ent_num = ent_calc_wrapper(tmp_pdb)
                    print("cut_site_index {cut_site} return {num_ent} entanglements".format(cut_site = cut_site_index, num_ent= tmp_ent_num))

                    # could be equal to what we have now
                    if tmp_ent_num <= initial_ent_num:
                        # pick one that gives rise to a smaller number of entanglement
                        initial_ent_num = tmp_ent_num.copy()
                        # stash its index in a tmp variable
                        local_index = cut_site_index.copy()
                    del tmp_pdb
                    
            else:
                # if site has been considered or in a prohibited secondary structure
                continue
            
        # we actually one or more index
        if (local_index != -1):
            already_cutted_sites.append(local_index)
            n_cut = n_cut + 1
        # well... shall we keep on searching or just quit?
        else:
            exit("no candidate found in n_cut = {n_cut} ".format(n_cut = n_cut) )
                        
                        
        
    
    return

if __name__ == '__main__':
    main()
    
