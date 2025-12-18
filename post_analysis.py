import os
import subprocess
import shutil
import glob
import json
import csv
import gemmi
import numpy as np
import pandas as pd
from pathlib import Path

def setup_results_dir(source_dir, results_dir):
    """Collect PDB and JSON files into results directory."""
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    
    print(f"Collecting files from {source_dir} to {results_dir}...")
    
    # Copy PDBs
    for pdb_file in glob.glob(os.path.join(source_dir, "*.pdb")):
        shutil.copy(pdb_file, results_dir)
        
    # Copy JSONs
    for json_file in glob.glob(os.path.join(source_dir, "*.json")):
        shutil.copy(json_file, results_dir)

def get_ca_coords(structure, chain_id):
    """Extract CA coordinates for a specific chain."""
    coords = []
    model = structure[0]
    for chain in model:
        if chain.name == chain_id:
            for res in chain:
                ca = res.find_atom("CA", '*')
                if ca:
                    coords.extend([ca.pos.x, ca.pos.y, ca.pos.z])
            break
    return coords

def align_and_extract(results_dir, output_csv):
    """Align structures and extract coordinates/scores."""
    pdb_files = sorted(glob.glob(os.path.join(results_dir, "*.pdb")))
    if not pdb_files:
        print("No PDB files found.")
        return

    ref_pdb = pdb_files[0]
    print(f"Using {os.path.basename(ref_pdb)} as reference for alignment.")
    
    ref_st = gemmi.read_structure(ref_pdb)
    
    # Prepare CSV header
    # We need to know the number of CA atoms in Chain B to generate headers
    ref_b_coords = get_ca_coords(ref_st, "B")
    print(f"Reference Chain B coords length: {len(ref_b_coords)}")
    
    num_atoms = len(ref_b_coords) // 3
    
    headers = ["file", "cluster", "ipTM", "pair_chains_ipTM_CRBN_VAV1_avg"]
    for i in range(1, num_atoms + 1):
        headers.extend([f"{i}x", f"{i}y", f"{i}z"])
        
    data_rows = []
    
    for pdb_file in pdb_files:
        basename = os.path.basename(pdb_file)
        print(f"Processing {basename}...")
        
        # Load structure
        st = gemmi.read_structure(pdb_file)
        
        # Align to reference Chain A
        try:
            # Perform superposition using Chain objects
            # Use CaP selection (CA for proteins, P for nucleic acids)
            sup = gemmi.calculate_superposition(ref_st[0]["A"].whole(), st[0]["A"].whole(), gemmi.PolymerType.PeptideL, gemmi.SupSelect.CaP)
            st[0].transform_pos_and_adp(sup.transform)
        except Exception as e:
            print(f"Alignment failed for {basename}: {e}")
            continue
            
        # Extract Chain B CA coords
        b_coords = get_ca_coords(st, "B")
        
        if len(b_coords) != len(ref_b_coords):
            print(f"Warning: Atom count mismatch for {basename}. Ref B: {len(ref_b_coords)}, Curr B: {len(b_coords)}. Skipping.")
            continue
            
        # Load scores
        # The JSON file has the same basename but with .json extension and 'confidence_' prefix
        # e.g. CRBN_vav1_template1_model_0.pdb -> confidence_CRBN_vav1_template1_model_0.json
        json_basename = basename.replace(".pdb", ".json").replace("CRBN_", "confidence_CRBN_")
        json_file = os.path.join(results_dir, json_basename)
        # print(f"Checking JSON: {json_file}")
        
        iptm = 0.0
        pair_iptm = 0.0
        
        if os.path.exists(json_file):
            with open(json_file, 'r') as f:
                scores = json.load(f)
                iptm = scores.get("iptm", 0.0)
                # Calculate average of pair_chains_iptm 0->1 and 1->0
                pair_chains = scores.get("pair_chains_iptm", {})
                val_0_1 = pair_chains.get("0", {}).get("1", 0.0)
                val_1_0 = pair_chains.get("1", {}).get("0", 0.0)
                pair_iptm = (val_0_1 + val_1_0) / 2.0
        else:
            print(f"Warning: JSON score file not found for {basename} at {json_file}")

        try:
            cluster = basename.split("template")[1].split("_")[0]
        except:
            cluster = "0"

        row = [basename, cluster, iptm, pair_iptm] + b_coords
        data_rows.append(row)
        # print(f"Added row for {basename}")
        
    # Write CSV
    print(f"Writing {len(data_rows)} rows to {output_csv}")
    with open(output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(data_rows)
        
    print(f"Saved aligned coordinates and scores to {output_csv}")

def main():
    # Paths
    base_dir = os.getcwd()
    combined_dir = os.path.join(base_dir, "..", "CRBN_VAV1_template_noMSA_20runs", "combined")
    results_dir = os.path.join(base_dir, "results")
    output_csv = "atom.csv"
    
    # 1. Collect results
    setup_results_dir(combined_dir, results_dir)
    
    # 2. Align and extract
    align_and_extract(results_dir, output_csv)
    
    # 3. Run UMAP visualization
    print("Running UMAP visualization...")
    import sys
    try:
        subprocess.run([sys.executable, "umap_coord_energy.py"], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running UMAP script: {e}")
    
    print("Post-analysis complete.")

if __name__ == "__main__":
    main()
