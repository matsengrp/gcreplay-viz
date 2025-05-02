import os,sys
import argparse
import pandas as pd
import glob
from Bio.PDB import PDBParser,PDBIO
from Bio.PDB.DSSP import DSSP


def pdb_get_df(pdb_path, chainids=None):
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure(pdb_path, pdb_path)
    pdb_sites = []
    for model in structure:
        for chain in model:
            if (chainids is None) or (chain.id in chainids):
                min_num = None
                for site,residue in enumerate(chain):
                    # residue sequence number
                    res_num = int(residue.id[1])
                    if min_num is None:
                        min_num = res_num
                        min_num = int(min_num / 100) * 100
                    # residue insertion code
                    res_ins = (residue.id[2].strip() or "-")
                    # residue id
                    _res_ins = (res_ins if (res_ins != "-") else "")
                    res_id = f"{res_num}{_res_ins}"
                    # residue amino short code
                    aa_long = residue.resname
                    aa_short = Encoder.long2short(residue.resname)
                    pdb_sites.append((site+1, chain.id, res_id, res_num, res_ins, aa_long, aa_short))
    df = pd.DataFrame(pdb_sites, columns=['site', 'chainid', 'res_id', 'res_num', 'res_ins', 'aa_long', 'aa_short'])
    return df


### MAIN ###

bench_dir = os.getcwd()
pdb_paths = os.listdir(bench_dir)
print(f"{pdb_paths=}")

pdb_dfs = {}
for pdb_path in pdb_paths:
    pdb_dfs[pdb_path] = pd.read_csv(pdb_path)

