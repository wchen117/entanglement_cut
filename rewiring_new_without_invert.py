# coding: utf-8
import sys

sys.path.append("../../Codes")


import gaussian_ent.gaussian_entanglements as ge
import gaussian_ent.clustering_updated as cluster


from warnings import filterwarnings

filterwarnings("ignore")


def rewiring(pdb, r1, r2, r3, r4):
    """[summary]

    Args:
        pdb ([type]): [description]
        r1 ([type]): [description]
        r2 ([type]): [description]
        r3 ([type]): [description]
        r4 ([type]): [description]

    Returns:
        [type]: [description]
    """

    import pickle

    import MDAnalysis as mda
    import numpy as np

    pdb_rewired = pdb.split(".pdb")[0] + "_rewiring_new.pdb"
    pdb_rewired_transform = pdb.split(".pdb")[0] + "_rewiring.transform_dict"
    u = mda.Universe(pdb, format="PDB")
    protein = u.select_atoms("all")
    residue_list = list(protein.residues.resids)
    resname_list = list(protein.residues.resnames)
    #  print(list(zip(resname_list,residue_list)))


    

    #  print(r1,r2,r3,r4)


    m1 = residue_list.index(r1)
    m2 = residue_list.index(r2)
    assert(m1<m2),"The first cut site is larger than second cut site"
    m3 = m1 + 1
    m4 = m2 + 1
    
    non_loop_chain_list = residue_list[:m1+1] + residue_list[m4:] 

    loop_chain_list = residue_list[m3:m2+1]




    c1 = non_loop_chain_list.index(r3)
    c2 = loop_chain_list.index(r4)
    list_rewired = non_loop_chain_list[:c1+1] + loop_chain_list[c2:] + loop_chain_list[:c2] + non_loop_chain_list[c1+1:]

    #  print(list_rewired)

    rewired_residue_positions=list()
    for resid in residue_list:
        position_in_rewired_list = list_rewired.index(resid)
        resid_corresponding_to_rewired_list = residue_list[position_in_rewired_list]
        rewired_residue_positions.append(resid_corresponding_to_rewired_list)


    rewired_array = np.array(rewired_residue_positions)
    u.residues.resids = rewired_array

    residue_list = list(protein.residues.resids)
    resname_list = list(protein.residues.resnames)
    #  print(list(zip(resname_list,residue_list)))

    #  u.select_atoms("protein").write(pdb_rewired)
    #  u.select_atoms("protein").write("test_rewiring.pdb")
    #  u.select_atoms("all").write("test_rewiring.pdb")

    return u, pdb_rewired


def get_a_pdb_line(atom, index):
    tag = "ATOM"
    name = atom.name
    resname = atom.resname
    chain = atom.segid
    resid = atom.resid
    x, y, z = atom.position
    occupancy = atom.occupancy
    tempfactor = atom.tempfactor
    type = atom.type
    pdb_line = f"{tag:<6}{index:>5}  {name:<3} {resname:<3} {chain:<1}{resid:>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f}{occupancy:>6.2f}{tempfactor:>6.2f}      {chain:<4}{type:>2}\n"
    return pdb_line, resname, chain, resid


def write_pdb(u, pdb_rewired):
    outfile = open(pdb_rewired, "w")
    chain_list = sorted(set(u.select_atoms("name CA").segids))

    index = 1
    for chain_id in chain_list:
        resids = sorted(u.select_atoms(f"name CA and segid {chain_id}").resids)

        for resid in resids:
            atoms = u.select_atoms(f"resid {resid} and segid {chain_id}")
            for atom in atoms:
                pdb_line, resname, chain, resid = get_a_pdb_line(atom, index)
                index += 1
                outfile.write(pdb_line)

    tag = "TER"
    name = " "
    pdb_line = f"{tag:<6}{index:>5}  {name:<3} {resname:<3} {chain:<1}{resid:>4}\n"
    outfile.write(pdb_line)
    outfile.write("END")
    outfile.close()


if __name__ == "__main__":
    pdb = sys.argv[1]
    ci1 = int(sys.argv[2])
    ci2 = int(sys.argv[3])

    chain_cut = int(sys.argv[4])
    loop_cut = int(sys.argv[5])
    


    universe, out_name = rewiring(pdb, ci1, ci2,chain_cut,loop_cut)
    if(sys.argv[6]):
        out_name = sys.argv[6]
    print(out_name)
    write_pdb(universe, out_name)
    ent_str=ge.gaussian_entanglement(out_name)
    ent="No entanglement"
    if ent_str:
        ent_dict,ent=cluster.cluster_entanglements(ent_str, 55)
    print(ent)
    ent_info_file=out_name.replace(".pdb",".ent_info")
    with open(ent_info_file,'w') as out:
        out.write(ent)
