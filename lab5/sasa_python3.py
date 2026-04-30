import sys
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

# --- CONFIGURATION ---
PROBE_RADIUS = 1.1  # Standard water probe size in Angstroms
N_POINTS = 100      # Higher = more accurate, slower (e.g., 960)
# ---------------------

def calculate_sasa(pdb_file):
    # 1. Initialize the parser and structure
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("protein", pdb_file)
    except Exception as e:
        print(f"Error: Could not read file {pdb_file}. {e}")
        return

    # 2. Setup the ShrakeRupley algorithm
    sr = ShrakeRupley(probe_radius=PROBE_RADIUS, n_points=N_POINTS)

    # 3. Run the calculation at the Residue level
    sr.compute(structure, level="R")

    # 4. Print Header
    print(f"Analysis for: {pdb_file}")
    print(f"Probe Radius: {PROBE_RADIUS} Å")
    print(f"{'-'*35}")
    print(f"{'Residue':<12} {'Chain':<8} {'SASA (Å²)':<10}")
    print(f"{'-'*35}")

    total_protein_sasa = 0.0

    # 5. Iterate and output
    for model in structure:
        for chain in model:
            for residue in chain:
                # Filter: Only process actual amino acids (ignores HOH/Heteroatoms)
                if residue.id[0] == ' ':
                    res_name = residue.get_resname()
                    res_num = residue.id[1]
                    sasa_val = residue.sasa
                    
                    total_protein_sasa += sasa_val
                    
                    # Print residue row
                    print(f"{res_name}{res_num:<7} {chain.id:<8} {sasa_val:<10.2f}")

    # 6. Print Total
    print(f"{'-'*35}")
    print(f"TOTAL PROTEIN SASA: {total_protein_sasa:.2f} Å²")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 sasa_py3.py <file.pdb>")
    else:
        calculate_sasa(sys.argv[1])
