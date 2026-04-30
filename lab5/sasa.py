import sys
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

# --- DEFAULT CONFIGURATION ---
DEFAULT_PROBE = 1.4  
# -----------------------------

def calculate_sasa(pdb_file, probe_size):
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("protein", pdb_file)
    except Exception as e:
        print(f"Error: Could not read file {pdb_file}. {e}")
        return

    # Setup the algorithm with the specific probe size
    sr = ShrakeRupley(probe_radius=probe_size, n_points=100)
    sr.compute(structure, level="R")

    print(f"Analysis for: {pdb_file}")
    print(f"Probe Radius: {probe_size} Å")
    print(f"{'-'*35}")
    print(f"{'Residue':<12} {'Chain':<8} {'SASA (Å²)':<10}")
    print(f"{'-'*35}")

    total_protein_sasa = 0.0

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == ' ': # Standard amino acids only
                    res_name = residue.get_resname()
                    res_num = residue.id[1]
                    sasa_val = residue.sasa
                    total_protein_sasa += sasa_val
                    print(f"{res_name}{res_num:<7} {chain.id:<8} {sasa_val:<10.2f}")

    print(f"{'-'*35}")
    print(f"TOTAL RELATIVE SASA: {total_protein_sasa:.2f} Å²")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 sasa.py <file.pdb> [probe_size]")
    else:
        filename = sys.argv[1]
        # Use 2nd argument as probe size if provided, otherwise default to 1.4
        p_size = float(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_PROBE
        calculate_sasa(filename, p_size)
