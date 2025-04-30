#!/usr/bin/env python3

import os,sys
import argparse
import shutil
import glob
import numpy as np
import pandas as pd
import json
import pprint
from pathlib import Path
import matplotlib.pyplot as plt

from Bio.PDB import PDBParser,PDBIO
from Bio.PDB.DSSP import DSSP
from utility import *

# palette for binding
DARK_PALETTE = ["#7570b3", "#808080"]  # original Dark2 pallete
VISIBLE_PALETTE = ["#675ed6", "#808080"]  # more visible pallete
ALT_PALETTE = ["#6A5ACD", "#B22222", "#2E8B57"]
COLOR_PALETTE = ALT_PALETTE
COLOR_MAP = 'brg'
# amino acid codes
AA_ALPHABET = sorted(list("RKHDEQNSTYWFAILMVGPC"))
AA_COUNT = len(AA_ALPHABET)


def mpl_rgba_to_hex(rgba):
    r, g, b, a = [int(x * 255) for x in rgba]
    return "#{:02X}{:02X}{:02X}".format(r, g, b)


def generate_color_palette(n_colors=10, colormap='gist_earth', as_hex=True):
    cmap = plt.get_cmap(colormap)
    colors = [cmap(i / max(n_colors - 1, 1)) for i in range(n_colors)]

    if as_hex:
        # Convert RGBA to hex
        return [mpl_rgba_to_hex(color) for color in colors]
    else:
        # Return RGB tuples without alpha
        return [tuple(color[:3]) for color in colors]


def pdb_get_chainids(pdb_path):
    chainids = []
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure(pdb_path, pdb_path)
    pdb_sites = []
    for model in structure:
        for chain in model:
            chainids.append(chain.id)
    return chainids


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


def pdb_get_flat_df(pdb_path):
    with open('yourfile.pdb', 'r') as f:
        lines = f.readlines()

    atom_lines = [line for line in lines if line.startswith(('ATOM', 'HETATM'))]
    records = []
    for line in atom_lines:
        record = {
            'record_name': line[0:6].strip(),
            'atom_serial': int(line[6:11]),
            'atom_name': line[12:16].strip(),
            'alt_loc': line[16],
            'res_name': line[17:20].strip(),
            'chain_id': line[21],
            'res_seq': int(line[22:26]),
            'i_code': line[26],
            'x': float(line[30:38]),
            'y': float(line[38:46]),
            'z': float(line[46:54]),
            'occupancy': float(line[54:60]),
            'temp_factor': float(line[60:66]),
            'element': line[76:78].strip(),
            'charge': line[78:80].strip()
        }
        records.append(record)

    df = pd.DataFrame(records)
    return df


def metric_get_binding_df(pdb_df, metric_path, chainids=None, metric_names=None):
    raw_metric_df = pd.read_csv(metric_path)
    raw_metric_df.loc[raw_metric_df["wildtype"] == raw_metric_df["mutant"], "mutant"] = "-"
    if chainids is not None:
        raw_metric_df = raw_metric_df[raw_metric_df["chain"].isin(chainids)]
        raw_metric_df["site"] = [x for x in range(1,len(set(raw_metric_df['position']))+1) for _ in range(AA_COUNT)]

    metric_df = pd.melt(
        raw_metric_df,
        id_vars=["site", "position", "position_IMGT", "chain", "wildtype", "mutant"],
        # value_vars=["bind_CGG", "expr"],
        value_vars=["single_nt",
                    "bind_CGG", "delta_bind_CGG", "n_bc_bind_CGG", "n_libs_bind_CGG",
                    "expr", "delta_expr", "n_bc_expr", "n_libs_expr"],
        var_name="condition",
        value_name="factor")
    metric_df["position_IMGT"] = metric_df["position_IMGT"].astype(int)
    if metric_names:
        metric_df = metric_df[metric_df["condition"].isin(metric_names)]
    return metric_df


def write_sitemap_csv(pdb_df, output_path, site_count=None):
    res_ins = [x if not x.startswith("-") else "" for x in pdb_df["res_ins"]]
    protein_sites = [f"{num}{ins}" for num, ins in zip(pdb_df["res_num"], res_ins)]
    sitemap_df = pd.DataFrame({
        'sequential_site': range(1, len(protein_sites) + 1),
        'reference_site': range(1, len(protein_sites) + 1),
        'protein_site': protein_sites,
    })
    if site_count:
        assert len(sitemap_df) >= site_count
        sitemap_df = sitemap_df[0:site_count].copy()

    if output_path:
        sitemap_df.to_csv(output_path, index=False)
    return sitemap_df


def write_metric_csv(pdb_df, metric_df, output_path, site_count=None, metric_cols=None):
    if site_count:
        assert len(metric_df) >= site_count
        metric_df = metric_df[metric_df["site"] > site_count]

    if metric_cols:
        metric_df = metric_df.copy()
        metric_df = metric_df[metric_df["condition"].isin(metric_cols)]

    if output_path:
        metric_df.to_csv(output_path, index=False)
    return metric_df


def configure_dms_viz(
    name,
    plot_colors,
    metric,
    input_metric_path,
    input_sitemap_path,
    output_path,
    included_chains=[],
    excluded_chains=[],
    add_options="",
    local_pdb_path=None,
):
    program = "configure-dms-viz format"
    plot_colors_str = ",".join(plot_colors)
    base_options = f'--name "{name}" \
                     --title "{name}" \
                     --description "GCReplay: {name}" \
                     --colors "{plot_colors_str}" \
                     --metric "{metric}" '
    if local_pdb_path:
        base_options += f'--structure "{local_pdb_path}" '
    if len(included_chains) != 0:
        base_options += f'--included-chains "{" ".join(included_chains)}" '
    if len(excluded_chains) != 0:
        base_options += f'--excluded-chains "{" ".join(excluded_chains)}" '

    cmd = f'{program} {base_options} \
            --input "{input_metric_path}" \
            --sitemap "{input_sitemap_path}" \
            --output "{output_path}" \
            {add_options} '

    cmd = " ".join(cmd.split())
    output = run_command(cmd)

    if output is None:
        raise Exception("ERROR: configure-dms-viz format run failed.")
    pass


def compare_seqs(aa_seqs):
    all_equal = True
    for i, seq_i in enumerate(aa_seqs):
        for j, seq_j in enumerate(aa_seqs):
            if i >= j:
                continue
            is_equal = (seq_i == seq_j)
            # print(f"seqs {i},{j}: {is_equal}")
            if not is_equal:
                all_equal = False
                print(seq_i)
                print(seq_j)
                for k in range(len(seq_i)):
                    if seq_i[k] == seq_j[k]:
                        print(seq_i[k], end="")
                    else:
                        print("-", end="")
                print()
    print(f"seqs::all_equal: {all_equal}")


def chainids_get_other_chainids(heavy_chainids=[], light_chainids=[], all_chainids=[]):
    other_chainids = []
    for chainid in all_chainids:
        if chainid not in (heavy_chainids + light_chainids):
            other_chainids.append(chainid)
    return other_chainids


### MAIN ###


def parse_args(args):
    arg_parser = argparse.ArgumentParser("gcreplay-viz pipeline")
    arg_parser.add_argument("--input-dir", type=Parser.parse_input_dir(), help="input directory for pdbs")
    arg_parser.add_argument("--output-dir", type=Parser.parse_output_dir(), help="output directory for dms-viz jsons")
    arg_parser.add_argument("--temp-dir", type=Parser.parse_output_dir(), help="temporary directory", default="_temp")
    arg_parser.add_argument("--chain-id", type=Parser.parse_list(str), help="heavy chain ids", default=["H"])
    arg_parser.add_argument("--light-chain-id", type=Parser.parse_list(str), help="light chain ids", default=["L"])
    parser = Parser(arg_parser=arg_parser)
    args = parser.parse_args()
    return args


def main(args=sys.argv):
    args = parse_args(args)
    input_dir = args['input_dir']
    output_dir = args['output_dir']
    temp_dir = args['temp_dir']
    # heavy_chainids = args["chain_id"]
    # light_chainids = args["light_chain_id"]
    heavy_chainids = ['H']
    light_chainids = ['L']
    focal_chainids = (heavy_chainids + light_chainids)
    # focal_chainids = (heavy_chainids)

    summary_data = {
        "dmsviz_filepath": [],
        "pdb_filepath": [],
        "pdbid": [],
        "chainid": [],
        "metricid": [],
        "metric_full_name": [],
        "description": [],
    }

    aa_seqs = {}
    all_chainids = []
    other_chainids = []

    # load pdb files
    input_pdb_paths = glob.glob(f"{input_dir}/*.pdb")
    print(input_pdb_paths)

    # get all chain ids
    for input_pdb_path in input_pdb_paths:
        all_chainids += pdb_get_chainids(pdb_path=input_pdb_path)
    # other chainids include chainids not in heavy or light chain
    all_chainids = list(set(all_chainids))
    print(f"all_chainids: {all_chainids}")
    for chainid in all_chainids:
        aa_seqs[chainid] = []

    # get all pdbs
    all_pdb_dfs = {}
    for pdb_path in input_pdb_paths:
        for chainid in all_chainids:
            pdb_df = pdb_get_df(
                pdb_path=pdb_path,
                chainids=chainid)
            if len(pdb_df) > 0:
                all_pdb_dfs[tuple([pdb_path, chainid])] = pdb_df

            aa_seq = ''.join(list(pdb_df['aa_short']))
            # print(f"aa_seq: {len(aa_seq)} {aa_seq}")
            if len(aa_seq) > 0:
                aa_seqs[chainid].append(aa_seq)

    # get first pdb_df from file
    first_key = next(iter(all_pdb_dfs))
    pdb_df = all_pdb_dfs[first_key]

    # parse metric files
    input_metric_paths = glob.glob(f"{input_dir}/*.csv")
    input_metric_path = input_metric_paths[0]

    metric_names = {
        # "all_metrics": ["bind_CGG", "expr", "delta_bind_CGG", "delta_expr", "n_bc_bind_CGG", "n_bc_expr", "n_libs_bind_CGG", "n_libs_expr", "single_nt"],
        "binding": ["bind_CGG"],
        "expression": ["expr"],
        "metric": ["bind_CGG", "expr"],
        "delta": ["delta_bind_CGG", "delta_expr"],
        "n_bc": ["n_bc_bind_CGG", "n_bc_expr"],
        "n_libs": ["n_libs_bind_CGG", "n_libs_expr"],
        "single_nt": ["single_nt"],
    }
    metric_full_names = {
        # "all_metrics": "All Binding/Expression metrics",
        "binding": "Binding",
        "expression": "Expression",
        "metric": "Binding/Expression",
        "delta": "Binding/Expression: Delta Change Relative to Wildtype",
        "n_bc": "Binding/Expression: Number of Barcodes",
        "n_libs": "Binding/Expression: Number of Libraries",
        "single_nt": "Binding/Expression: Mutation by Single Nucleotide Change",
    }

    all_metric_dfs = {}
    for (pdb_path, chainid), pdb_df in all_pdb_dfs.items():
        metric_df = metric_get_binding_df(
                    pdb_df=pdb_df,
                    metric_path=input_metric_path,
                    chainids=[chainid],
                    metric_names=None)
        if len(metric_df) > 0:
            all_metric_dfs[chainid] = metric_df

        tmp_metric_df = metric_df.drop_duplicates(subset=["position"])
        aa_seq = ''.join(list(tmp_metric_df['wildtype']))
        # print(f"aa_seq: {chainid} {len(aa_seq)} {aa_seq}")
        if len(aa_seq) > 0:
            aa_seqs[chainid].append(aa_seq)

    pprint.pp(aa_seqs)
    # compare_seqs(aa_seqs=aa_seqs)

    # build csvs and dmsviz jsons
    for (pdb_path, chainid), pdb_df in all_pdb_dfs.items():
        pdb_prefix = os.path.basename(pdb_path).split(".")[0]
        # pdb_prefix = "CGG_naive_DMS"
        chain_str = f"{''.join(chainid)}"
        print(f"pdb: {pdb_prefix=} {chainid=}")

        # skip if not in focal_chainids
        if chainid not in focal_chainids:
            continue
        other_chainids = chainids_get_other_chainids(
            heavy_chainids=heavy_chainids, light_chainids=light_chainids, all_chainids=all_chainids)

        for metric_name, metric_cols in metric_names.items():
            print(f"metric: {metric_name=} {metric_cols=}")
            metric_full_name = metric_full_names[metric_name]
            metric_df = all_metric_dfs[chainid]

            # get number of metrics
            metric_df = metric_df[metric_df["condition"].isin(metric_cols)]
            metric_types = set(metric_df["condition"])
            num_metrics = len(metric_types)

            # prune down to only common IMGT sites
            # if only_common_sites:
            pdb_sites = set(pdb_df.res_id.astype(str))
            metric_sites = set(metric_df.position_IMGT.astype(str))
            union_sites = pdb_sites & metric_sites
            xor_sites = pdb_sites ^ metric_sites
            metric_df = metric_df[metric_df.position_IMGT.astype(str).isin(union_sites)]
            print(f"omitted_sites: {len(xor_sites)} {sorted(list(xor_sites))}")
            metric_sites = sorted(list(set(metric_df.site)))
            metric_site_map = {x: y for x, y in zip(metric_sites, range(1, len(metric_sites)+1))}
            metric_df["site"] = [metric_site_map[x] for x in metric_df["site"]]

            # build sitemap csv
            sitemap_path = f"{temp_dir}/{pdb_prefix}.{chain_str}.sitemap.csv"
            sitemap_df = write_sitemap_csv(
                pdb_df=pdb_df,
                output_path=sitemap_path)

            # build metric csv
            metric_path = f"{temp_dir}/{pdb_prefix}.{chain_str}.{metric_name}.csv"
            metric_df = write_metric_csv(
                pdb_df=pdb_df,
                metric_df=metric_df,
                output_path=metric_path,
                metric_cols=metric_cols,)

            add_options = ""
            condition_options = '--condition "condition" '
            condition_options += '--condition-name "Factor" '
            add_options += condition_options

            try:
                # build dms-viz json
                full_description = f"{pdb_prefix} :: {chain_str} :: {metric_full_name}"
                dmsviz_path = f"{temp_dir}/{pdb_prefix}.{chain_str}.{metric_name}.dmsviz.json"
                # COLOR_PALETTE = generate_color_palette(
                #     n_colors=num_metrics, colormap=COLOR_MAP, as_hex=True)
                COLOR_PALETTE=ALT_PALETTE[:num_metrics]
                configure_dms_viz(
                    name=full_description,
                    plot_colors=COLOR_PALETTE,
                    metric="factor",
                    input_metric_path=metric_path,
                    input_sitemap_path=sitemap_path,
                    output_path=dmsviz_path,
                    included_chains=chainid,
                    excluded_chains=other_chainids,
                    add_options=add_options,
                    local_pdb_path=input_pdb_path)

                # add summary data entry
                summary_data["metric_full_name"].append(metric_full_name)
                summary_data["dmsviz_filepath"].append(os.path.basename(dmsviz_path))
                summary_data["pdb_filepath"].append(os.path.basename(pdb_path))
                summary_data["pdbid"].append(pdb_prefix)
                summary_data["chainid"].append(chain_str)
                summary_data["metricid"].append(metric_name)
                summary_data["description"].append(full_description)

            except Exception as e:
                cprint(f"[ERROR] {pdb_prefix} {chainid} {metric_name}", color=colors.RED)
                cprint(f"[ERROR] error occurred during configure-dms-viz: {e}", color=colors.RED)
            else:
                cprint(f"[SUCCESS] configure-dms-viz completed successfully!", color=colors.GREEN)

        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(f"{temp_dir}/summary.csv", index=False)
        summary_json = summary_df.to_json(orient='records')
        with open(f"{temp_dir}/summary.json", "w") as file:
            # json.dump(summary_json, file)
            file.write(f"{summary_json}\n")
        print(summary_df)

    if args['output_dir'] is not None:
        temp_jsons = glob.glob(f"{temp_dir}/*.dmsviz.json")
        for temp_json in temp_jsons:
            shutil.copy(temp_json, f"{output_dir}/dmsviz-jsons/")
        shutil.copy(f"{temp_dir}/summary.csv", f"{output_dir}/metadata/summary.csv")
        shutil.copy(f"{temp_dir}/summary.json", f"{output_dir}/metadata/summary.json")
    return


if __name__ == "__main__":
    print("[BEGIN] main")
    main()
    print("[END] main")
