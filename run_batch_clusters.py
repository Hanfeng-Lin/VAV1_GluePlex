import os
import glob
import subprocess
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
# ---------------------------------------------------------------
# ------------------------ Configuration ------------------------
# ---------------------------------------------------------------
CIF_FILES = [
    "cluster1_1.pdb.cif",
    "cluster2_1.pdb.cif",
    "cluster3_1.pdb.cif",
    "cluster4_1.pdb.cif",
    "cluster5_1.pdb.cif",
    "cluster6_1.pdb.cif",
]
# Set specific GPUs to use (e.g., [0, 1, 2, 6, 7]). 
# If empty, will attempt to detect available GPUs.
GPUS = [0, 1, 2, 6, 7] 

DIFFUSION_SAMPLES = 20
PROTEIN_A_NAME = "CRBN"
PROTEIN_A_SEQUENCE = "IINFDTSLPTSHTYLGADMEEFHGRTLHDDDSCQVIPVLPQVMMILIPGQTLPLQLFHPQEVSMVRNLIQKDRTFAVLAYSNVQEREAQFGTTAEIYAYREEQDFGIEIVKVKAIGRQRFKVLELRTQSDGIQQAKVQILPECVLPSTMSAVQLESLNKCQIFPSKPVSREDQCSYKWWQKYQKRKFHCANLTSWPRWLYSLYDAETLMDRIKKQLREWDENLKDDSLPSNPIDFSYRVAACLPIDDVLRIQLLKIGSAIQRLRCELDIMNKCTSLCCKQCQETEITTKNEIFSLSLCGPMAAYVNPHGYVHETLTVYKACNLNLIGRPSTEHSWFPGYAWTVAQCKICASHIGWKFTATKKDMSPQKFWGLTRSALLPTIPDTEDEISPDKVILCL"
PROTEIN_B_NAME = "VAV1"
PROTEIN_B_SEQUENCE = "KYFGTAKARYDFCARDRSELSLKEGDIIKILNKKGQQGWWRGEIYGRVGWFPANYVEEDYS"
LIGAND_SMILES = "O=C1CCC(c2cccc(-c3ccc4[nH]ccc4c3)c2)C(=O)N1"
# ---------------------------------------------------------------
# ---------------------------------------------------------------

# Template for the run script
SCRIPT_TEMPLATE = """import os
import subprocess
import sys
from pathlib import Path

YAML_CONTENT = \"\"\"
version: 1  # Optional, defaults to 1
sequences:
  - protein:
      id: A  # {protein_a_name}
      sequence: {protein_a_seq}
      msa: empty
  - protein:
      id: B  # {protein_b_name}
      sequence: {protein_b_seq}
      msa: empty
  - ligand:
      id: C
      smiles: '{ligand_smiles}'
templates:
    - cif: {cif_path} # if a pdb path is provided, Boltz will incrementally assign template chain ids based on the chain names in the PDB file (A1, A2, B1, etc)
      chain_id: [A, B]
      template_id: [A, B]
\"\"\"

def main(yaml_file, out_dir, output_format="pdb", diffusion_samples={diffusion_samples}, gpu_id=0):
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    with open(yaml_file, "w") as f:
        f.write(YAML_CONTENT)
    print(f"Created {{yaml_file}}")

    cmd = [
        "boltz",
        "predict",
        yaml_file,
        "--out_dir",
        out_dir,
        "--output_format",
        output_format,
        "--diffusion_samples",
        str(diffusion_samples)
    ]
    
    print(f"Running command: {{' '.join(cmd)}}")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running boltz: {{e}}")
        sys.exit(1)

if __name__ == "__main__":
    main("{yaml_filename}", "{out_dir}", 
        output_format="pdb", diffusion_samples={diffusion_samples}, gpu_id={gpu_id})
"""

def get_visible_devices():
    """Get list of visible CUDA devices."""
    if GPUS:
        return GPUS
        
    cuda_visible = os.environ.get("CUDA_VISIBLE_DEVICES")
    if cuda_visible:
        return [int(x.strip()) for x in cuda_visible.split(",")]
    else:
        # Fallback: try to count GPUs using nvidia-smi or just assume 0 if fails
        try:
            result = subprocess.run(["nvidia-smi", "-L"], capture_output=True, text=True)
            if result.returncode == 0:
                count = len(result.stdout.strip().split('\\n'))
                return list(range(count))
        except FileNotFoundError:
            pass
        return [0] # Default to device 0 if detection fails

def process_cif(cif_file, gpu_id):
    """Generate script and run it."""
    cif_path = Path(cif_file)
    if not cif_path.exists():
        print(f"Warning: {cif_file} does not exist.")
        return

    # Extract cluster number from filename (e.g., cluster1_1.pdb.cif -> 1)
    try:
        cluster_name = cif_path.name.split('_')[0] # cluster1
        cluster_num = cluster_name.replace("cluster", "")
    except Exception:
        cluster_name = cif_path.stem
        cluster_num = "X"

    script_name = f"run_boltz2_{cluster_name}.py"
    yaml_filename = f"{PROTEIN_A_NAME}_{PROTEIN_B_NAME}_template{cluster_num}.yaml"
    out_dir = f"./{PROTEIN_A_NAME}_{PROTEIN_B_NAME}_template_noMSA_20runs" # Using same output dir as requested

    script_content = SCRIPT_TEMPLATE.format(
        cif_path=f"./{cif_path.name}",
        yaml_filename=yaml_filename,
        out_dir=out_dir,
        gpu_id=gpu_id,
        protein_a_name=PROTEIN_A_NAME,
        protein_a_seq=PROTEIN_A_SEQUENCE,
        protein_b_name=PROTEIN_B_NAME,
        protein_b_seq=PROTEIN_B_SEQUENCE,
        ligand_smiles=LIGAND_SMILES,
        diffusion_samples=DIFFUSION_SAMPLES
    )

    with open(script_name, "w") as f:
        f.write(script_content)
    
    print(f"Generated {script_name} for {cif_file} on GPU {gpu_id}")
    
    # Run the script
    cmd = ["python", script_name]
    try:
        subprocess.run(cmd, check=True)
        print(f"Finished {script_name}")
    except subprocess.CalledProcessError as e:
        print(f"Error running {script_name}: {e}")

def main():
    cif_files = CIF_FILES
    if not cif_files:
        print("No CIF files defined in CIF_FILES list.")
        return

    gpus = get_visible_devices()
    print(f"Detected GPUs: {gpus}")
    
    with ThreadPoolExecutor(max_workers=len(gpus)) as executor:
        futures = []
        for i, cif_file in enumerate(cif_files):
            gpu_id = gpus[i % len(gpus)]
            futures.append(executor.submit(process_cif, cif_file, gpu_id))
        
        for future in futures:
            future.result()

if __name__ == "__main__":
    main()
