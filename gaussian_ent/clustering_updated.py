from collections import defaultdict
from scipy.spatial.distance import squareform, pdist, cdist
from geom_median.numpy import compute_geometric_median
import numpy as np
import itertools
import random
import copy

import sys


def group_using_pdb(file_name):
    """
    This function groups the list of the entangled list by the PDB id.
    
    argument:
        file_name: file containg the information about the entanglement 
    return:
        entangled_info: defaultdict where the key is the PDB id and the value are list of list of the information about the loop and cross_over
    """

    start_query = False
    entangled_info = defaultdict(str)
    uniprot_id = file_name.split("/")[-1]
    uniprot_id = uniprot_id.split(".")[0]

    with open(file_name, "r") as f:
        for line in f.readlines():
            pdb,*garbage=line.split(",")
            entangled_info[pdb] += line
    return entangled_info


def loop_distance(entangled_A: tuple, entangled_B: tuple):

    smaller_ent = min(len(entangled_A), len(entangled_B))
    bigger_ent = max(len(entangled_A), len(entangled_B))

    if len(entangled_A) == smaller_ent:

        cr_values = list(entangled_A[2:])

        med = np.median(cr_values)

        entangled_A += (med, ) * (bigger_ent - len(entangled_A))

    else: 

        cr_values = list(entangled_B[2:])

        med = np.median(cr_values)

        entangled_B += (med, ) * (bigger_ent - len(entangled_B))

    difference = np.asarray(entangled_A, dtype=int) - np.asarray(entangled_B, dtype=int)

    squared_diference = difference ** 2

    return np.sqrt(np.sum(squared_diference))


def cluster_entanglements(entanglement_info: str, cut_off: int):


    full_entanglement_data = defaultdict(list)

    rep_chain_ent = defaultdict(list)

    GE_data = entanglement_info.strip().split("\n")
    for line in GE_data:
       pdb,chain,native_contact_i,native_contact_j,*threading,garbage = line.replace("[","").replace("]","").replace(","," " ).split()

       crossing_res = list(map(lambda x: int(x), threading))
       crossing_res = np.array(crossing_res)
       full_entanglement_data[chain].append((int(native_contact_i), int(native_contact_j), *crossing_res))
    for chain in full_entanglement_data.keys():

        length_key = defaultdict(list)
        loop_dist = defaultdict(list)
        dups = []
        clusters = {} 
        cluster_count = 0

        pairwise_entanglements = list(itertools.combinations(full_entanglement_data[chain], 2))

        if pairwise_entanglements:

            for i, pairwise_ent in enumerate(pairwise_entanglements):

                dist = loop_distance(pairwise_ent[0], pairwise_ent[1])

                if dist <= cut_off and pairwise_ent[0] not in dups and pairwise_ent[1] not in dups:
                    # 1. pair must be <= =cut_off
                    # 2. the neighbor cannot be the next key and it cannot be captured by another key

                    loop_dist[pairwise_ent[0]].append(pairwise_ent[1])
                    dups.append(pairwise_ent[1])
                    
            key_list = list(loop_dist.keys())

            for key in key_list:

                length_key[len(loop_dist[key])].append(key)

            # create clusters

            while len(length_key.values()) > 0:

                max_neighbor = max(length_key.keys())

                selected_ent = random.choice(length_key[max_neighbor])

                cluster = copy.deepcopy(loop_dist[selected_ent])
                cluster.append(selected_ent)

                clusters[cluster_count] = cluster
                cluster_count += 1

                length_key[max_neighbor].remove(selected_ent)

                if len(length_key[max_neighbor]) == 0:
                    length_key.pop(max_neighbor)
            
        # create single clusters

        if clusters:
            clusters_ijr_values = list(itertools.chain.from_iterable(list(clusters.values())))
        else:
            clusters_ijr_values = []

        full_ent_values = np.asarray(full_entanglement_data[chain], dtype=object)

        difference_ent = np.zeros(len(full_ent_values), dtype=bool)

        for k, ijr in enumerate(full_ent_values):

            if tuple(ijr) in clusters_ijr_values:
                difference_ent[k] = True
            else:
                difference_ent[k] = False

        i = np.unique(np.where(difference_ent == False)[0])

        next_cluster_count = cluster_count

        for single_cluster in full_ent_values[i]:
            
            single_cluster_list = []
            single_cluster_list.append(tuple(single_cluster))

            clusters[next_cluster_count] = single_cluster_list

            next_cluster_count += 1

        # pick representative entanglement per cluster

        for counter, ijr_values in clusters.items():

            # clusters contain many entanglements
            if len(ijr_values) > 1:

                max_length = max([len(point) for point in ijr_values])

                ijr = []
                remove_extra_cr = {}

                for point in ijr_values:

                    if len(point) != max_length:
                    
                        med = np.median(point[2:])

                        extra_cr = [med] * (max_length - len(point))

                        point += (med, ) * (max_length - len(point))

                        remove_extra_cr[point] = extra_cr
                        ijr.append(point)
                    
                    else:
                        ijr.append(point)

                ijr = np.asarray(ijr)

                cr_values = ijr[:, 2:]

                median_cr = compute_geometric_median(cr_values).median

                distances = cdist(cr_values, [median_cr])

                minimum_distances_i = np.where(distances == min(distances))[0]

                possible_cand = ijr[minimum_distances_i]

                loop_lengths = np.abs(possible_cand[:, 0] - possible_cand[:, 1])

                smallest_loop_length = min(loop_lengths)

                rep_entanglement = possible_cand[random.choice(np.where(smallest_loop_length == loop_lengths)[0])]

                # remove extra padding if rep_entanglement has padding
                if tuple(rep_entanglement) in remove_extra_cr:

                    extra_padding = remove_extra_cr[tuple(rep_entanglement)]

                    remove_indices = -len(extra_padding)

                    rep_entanglement = rep_entanglement[:remove_indices]

                rep_chain_ent[f"{chain}_{counter}"].append(rep_entanglement)

            # clusters with a single entnalgement
            else:
                rep_chain_ent[f"{chain}_{counter}"].append(ijr_values[0])
    line=""
    for chain_counter, ijrs in rep_chain_ent.items():

        chain, counter = chain_counter.split("_")

        for ijr in ijrs:
            crossings = list(ijr[2:])
            line += f"[{pdb}, {chain}, {ijr[0]}, {ijr[1]}, {crossings}, 99999:{counter}\n"
    return rep_chain_ent, line

def clustering(file_name,cutoff):
    total_info = ""
    entanglement_pdb_info = group_using_pdb(file_name)
    for pdb, ent_info in entanglement_pdb_info.items():
        rep_chain_info, ent_string = cluster_entanglements(ent_info,cutoff)
        total_info += ent_string

    return total_info
            

def date_time():
    from datetime import datetime
    now = datetime.now()
    return now.strftime("%d-%m-%Y %H:%M:%S")

def get_files(directory):
    from os import listdir
    from os.path import isfile, join
    onlyfiles = [join(directory,f) for f in listdir(directory) if isfile(join(directory,f))]
    return onlyfiles

def write_clustering(info_dict):
    file_name = info_dict["ent_file"]
    out_file_name = info_dict["out_file"]
    cutoff = info_dict["cutoff"]

    try:
        clustered_string = clustering(file_name,cutoff)
    except Exception as e:
        print(f"Error: {e} in file {file_name}")

    with open(out_file_name,"w") as outfile:
        outfile.write(clustered_string)
        print(f"Success: Clustering for {file_name}")


if __name__ == "__main__":

    # only change that I made now is at line 133 and made it an object array to avoid a deprecation warning.
    # If this cause a problem (i.e code cannot run) just remove dtype=object and ignore the warning.  

    # Example how to run:
   #  string_entanglement=clustering("../../GENE_INFORMATION/Entanglement_Information/P0A794.txt", 55)
   #  print(string_entanglement)

    import argparse
    import multiprocessing as mp


    parser = argparse.ArgumentParser(
        prog='ent_clustering',
        formatter_class = argparse.RawDescriptionHelpFormatter,
        description = """Entanglement clustering""",
        epilog = """
        Example: python clustering_updated -i ~/Entanglement/ -o ~/Out_dir/ -n 16 
        """)

    parser.add_argument('-i', '--input', dest="input_directory", type = str, action='store', help="provide directory where entanglement files are located", required=True)
    parser.add_argument('-o','--output', dest="output_directory", type = str, action="store", help="provide directory where the images is to be saved", default = "./", required = False)
    parser.add_argument('-n','--process', dest="process", type = int,action="store", help="provide number of processors", default = 1, required= False)
       
    args = parser.parse_args() 
        
    input_dir = args.input_directory
    output_dir_for_data_base = args.output_directory

    process = args.process
    multiprocess = False

    if process > 1:
        multiprocess = True



    processed_ent = get_files("/gpfs/group/epo2/default/hkb5295/from_Viraj/GENE_INFORMATION/Entanglement_Information_Clustered_New/")
    get_uniport = lambda x: x.split("/")[-1].split(".")[0]
    processed_ent = list(map(get_uniport, processed_ent))
    


    entanglement_list = get_files(input_dir)

    input_arguments = list()
    for count, entanglements_file in enumerate(entanglement_list):
        gene_name = entanglements_file.split("/")[-1].split(".")[0]
        try:
            if gene_name in processed_ent:
                continue

            enteries_input_list = dict()
            enteries_input_list["ent_file"] = entanglements_file
            enteries_input_list["out_file"] = f"{output_dir_for_data_base}{gene_name}.clustered_txt"
            enteries_input_list["cutoff"] = 55
            input_arguments.append(enteries_input_list)

        except Exception as e:
            print(f"{entanglements_file}-> {e}")
            continue 

    if multiprocess:
        pool = mp.Pool(process)
        pool.map(write_clustering,input_arguments)
        pool.close()
    else:
        for entries in input_arguments:
            write_clustering(entries)
            #  if "P76318" in entries["ent_file"]:
                #  clustered_string = clustering(entries['ent_file'],55)
                #  print(clustered_string)
                #  sys.exit()


    
        
    
      
    
    

