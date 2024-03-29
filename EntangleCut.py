##############################################################
# Cut.py, a simple program to iteratively cut protein to 
# eliminate entanglement
# Author: Weinan Chen
# May 2022
##############################################################

from ast import arg, parse
import MDAnalysis
import numpy as np
import gaussian_ent.gaussian_entanglements as ge
import gaussian_ent.clustering_updated as cluster
import subprocess
import sys
import re
import itertools
from IndexNode import IndexNode
from anytree import RenderTree
from gaussian_ent.rewiring_new_without_invert import write_pdb
import argparse


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
        masked_array = np.ones(pdb_length, dtype=int)
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
        ent_dict, ent=cluster.cluster_entanglements(list_of_ents, 55)

        if ent_dict is not None:
            return len(ent_dict)
    else:
        # no entanglement, yay!
        return -1

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

def compute_residue_distance(pdb_struct: MDAnalysis.Universe, resid_config: list, max_distance: int):
    """Given a reconstructed pdb, compute the max inter-residue distances
        so that the maximum inter-residue distance is <= max_distance (angstrom)"""
    residue_pos = pdb_struct.select_atoms("name CA").positions
    # reconnected cut sites are those whose indices changes abruptly
    cut_sites = np.where(np.diff(resid_config) !=1)
    
    for each_site in cut_sites[0]:
        #site_index = resid_config[each_site]
        left = resid_config[each_site] - 1
        right = resid_config[each_site+1] - 1
        tmp_dist = np.linalg.norm(residue_pos[left] - residue_pos[right])
        #import ipdb; ipdb.set_trace()
        print("dist = %.2f \AA"%(tmp_dist))
        #print("dist = {dist} \AA, sites = "tmp_dist, left, right)
        if tmp_dist > max_distance:
            #print(each_site)
            return -1

    return 0

def cut_protein(resid_list:list, site_indices:list, max_site: int):
    """ cut the protein sturcture at a given index.
        Note that the only input needed for cutting is resids/ list of resids, which is the first argument
        of this function.
        Please note that at this stage, the first input must be a complete resid_list, not a nested_list
        max_site must be smaller than n_residue - 1"""

    # check if input exceed max_site
    mask = np.array(site_indices) > max_site
    # all elements must be 0 (site_index >= max_site is false for all input site indices)
    if (mask.prod() != 0):
        exit("input site_indices > max_site is not allowed")
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

def cut_site_separation_constrain(previous_site_list: list, new_site: int, separation: int):
    """A constrain such that the new site is at least n (n = separation) 
       site away from sites in previous_site_list
       return 0 or 1"""

    mask = np.abs(new_site - np.array(previous_site_list)) >= separation

    # if true for all elements, return 1
    # otherwise return 0
    return np.prod(mask)


def iterative_cut(already_cutted_sites: list, protein_mask: np.ndarray, resid_list: list, \
                  n_cut_sites: int, pdb_struct: MDAnalysis.Universe, initial_ent_num: int, \
                  site_separation: int, max_distance: int):

    """A function wrapper for looping through all the cut sites in protein chain,
       taking into account already_cutted_sites from previous n_cut numbers, 
       and try to find new configurations that has smaller or equal number of entanglements,
       and return them as a list of new nodes"""
    
    list_of_new_nodes = []

    # can't cut after the last residue
    for cut_site_index in resid_list[:-1]:
    
        #print("working on cut_site_index = %d ..."%(cut_site_index))
        # need to unpack the tree structure, so that for each path
        # one could get a list of indices on the path, to form a local_cut_sites list

        # the first n elements of local_cut_sites are from previous 
        #local_cut_sites = []
        # taking in account cut sites from the previous iterations 
        local_cut_sites = already_cutted_sites + [cut_site_index,]

        site_boolean = cut_site_separation_constrain(already_cutted_sites, cut_site_index, site_separation)
        if site_boolean == 0:
            break
         
        #local_cut_sites = [13,8]
        # see if the cut_site_index cuts into the secondary structure identified by stride program
        if (protein_mask[cut_site_index] == 1 and cut_site_index not in already_cutted_sites):
            nested_resid_pieces = cut_protein(resid_list, local_cut_sites, n_cut_sites)
            # each configuration is 
            all_configs = permute_connect(nested_resid_pieces)
            for each_config in all_configs:
                # compute the number of entanglement (perhaps G. score) per configuration
                # need a copy of the current pdb_struct to create new_pdb 
                tmp_pdb = construct_pdb(each_config, resid_list, pdb_struct.copy())
                print("local_cut_sites = ", local_cut_sites)
                # we can try to enforce the constrain (<= 8 angstrom inter-residue distance) here
                distance_constrain = compute_residue_distance(tmp_pdb, each_config, max_distance)
                # the reconnected residues are further apart than 8 angstrom
                if distance_constrain < 0:
                    print("reconnect constrain not feasible")
                    continue
                tmp_ent_num = ent_calc_wrapper(tmp_pdb)
               
                print("config = ", each_config)

                #print("cut_site_index {cut_site} return {num_ent} entanglements".\
                #        format(cut_site = cut_site_index, num_ent= tmp_ent_num))

                # found an ent = 0 case
                if tmp_ent_num == -1:
                    # found a candidate, write it out
                    out_name = "_".join(str(x) for x in local_cut_sites)
                    write_pdb(tmp_pdb, "{name}.pdb".format(name = out_name))
                    # can we just .... exit here?
                    print("found a candidate, writing its pdb to current folder, exit!")
                    exit(0) 
                    #import ipdb; ipdb.set_trace()
                
                if tmp_ent_num <= initial_ent_num:
                    # pick one that gives rise to a smaller number of entanglement
                    #initial_ent_num = tmp_ent_num.copy()
                    # stash its index in a tmp variable
                    #local_index = cut_site_index.copy()
                    newNode = IndexNode("{index}".format(index=cut_site_index), Ent=tmp_ent_num, Index=cut_site_index)
                    list_of_new_nodes.append(newNode)
                del tmp_pdb
                    
        else:
            # if site has been considered or in a prohibited secondary structure
            continue
            
    # we actually have one or more index
    
    return list_of_new_nodes

def attach_nodes(node: IndexNode, list_of_nodes: list):

    """attach list_of_nodes (list) as new leaves of input node (IndexNode)"""

    for each_node in list_of_nodes:
         each_node.parent = node

    return


def parse_arg():
    """parse the input argument. for now, the code accepts: 1. input file name
       2. cut site index separations 3. max separation distance between two connected sites
       4. max number of iterations"""
    helper_message = "A simple python program to iteratively cut protein chains. Usage: Python EntangleCut.py"
    parser = argparse.ArgumentParser(description=helper_message)
    
    parser.add_argument("-f", type=str, default = "./pdbs/native_chain_B.pdb", help="required field for input pdb file")
    parser.add_argument("-s", type=int, default = 4, help="minimum cut site indices separations, default = 4")
    parser.add_argument("-m", type=float, default = 8, help="maximum geometric distance threshold between two re-connected sites in Angstrom, default = 8")
    parser.add_argument("-i", type=int, default = 5, help="maximum number of concurrent cut sites to consider, default = 3")

    return parser.parse_args()

def pack_all_nodes(root_node: IndexNode, resid_list: list, ent_num: int):
    """attach all resid in the resid_list to the root_node,
       while assigning the same ent_num to all attached node.
       This function is mainly used to bypass the situation where a single cut
       site doesn't produce a valid candidate
    """
    # can't cut after the last residue
    for cut_site_index in resid_list[:-1]:
        node = IndexNode("{index}".format(index=cut_site_index), Ent=ent_num, Index=cut_site_index, parent = root_node)

    return

def enforce_secondary_structure_constrain(resid_list: list, protein_mask: np.ndarray):
    """ used explicitely to enforce the secondary structure constrain 
        at the end of n_cut = 1, if no viable candidate is selected
    """
    local_resid = np.array(resid_list)
    mask = local_resid * protein_mask

    return list(local_resid[mask != 0])

def main():

    # parse command line argument
    parsed_arg = parse_arg()
   
    pdb_file_name = parsed_arg.f
    separation = parsed_arg.s
    max_distance = parsed_arg.m
    max_iteration = parsed_arg.i
    

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
    #print(protein_mask)

    # calculate the inital number of entanglement, start with this
    initial_ent_num = ent_calc_wrapper(pdb_struct)

    # the global tree structure that keep tracks of cutting indices
    # we use the IndexNode class for this
    # the root node has a default index of -1 
    tree_struct = IndexNode("root", Ent=initial_ent_num)
    
    
    # for now limit the maximum concurrent cut to 5
    while (n_cut <= max_iteration):

        # if we are starting, no need to check for leaves, directly continue to test
        # all 302 cut sites for n_cut = 1 scenarios
        if n_cut == 1:
            previous_cut_sites = []
            #list_of_new_nodes = iterative_cut(previous_cut_sites, protein_mask, resid_list, 6, pdb_struct, tree_struct.Ent)

            list_of_new_nodes = iterative_cut(previous_cut_sites, protein_mask, resid_list,\
                                              n_cut_sites, pdb_struct, initial_ent_num, separation, max_distance)
            #attach_nodes(tree_struct, list_of_new_nodes)
            # if no candidate from this round
            if (len(list_of_new_nodes) == 0):
                # 
                #import ipdb; ipdb.set_trace()
                allowed_resids = enforce_secondary_structure_constrain(resid_list, protein_mask)
                pack_all_nodes(tree_struct, allowed_resids, 99)
                #previous_cut_sites = resid_list.copy()
            else:
                
                attach_nodes(tree_struct, list_of_new_nodes)


        else: 
            
             # for each terminal nodes in this tree, for n = 2, it's layer 2, for n =3, it's layer 3 
            for leaf in tree_struct.leaves:
                # only care about the most extended layers, avoid recalculating the previous layers
                if leaf.depth + 1 == n_cut:

                    previous_cut_sites = leaf.getIndexPath()
                    list_of_new_nodes = iterative_cut(previous_cut_sites, protein_mask, resid_list,\
                                                      n_cut_sites, pdb_struct, leaf.Ent, separation, max_distance)
                    attach_nodes(leaf, list_of_new_nodes)
                else:
                    continue


        tree_struct.visualize()

        n_cut = n_cut + 1


                        
        
    
    return

if __name__ == '__main__':
    main()
    
