import sys
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

# CONFIGURATION: Change your probe size here
MY_PROBE_SIZE = 1.4  # Value in Angstroms
MY_RESOLUTION = 100  # Number of points (higher = more accurate but slower)

if len(sys.argv) < 2:
    print "Usage: python code.py <your_file.pdb>"
    sys.exit(1)

pdb_filename = sys.argv[1]

try:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_filename)
except Exception as e:
    print "Error loading PDB file: {}".format(e)
    sys.exit(1)

# Initialize with your custom probe size
sr = ShrakeRupley(probe_radius=MY_PROBE_SIZE, n_points=MY_RESOLUTION)

# Compute
sr.compute(structure, level="R")

print "Using Probe Radius: {} A".format(MY_PROBE_SIZE)
print "{:<10} {:<6} {:<10}".format("Residue", "Chain", "SASA (A^2)")
print "-" * 30

total_sasa = 0.0
for model in structure:
    for chain in model:
        for residue in chain:
            if residue.id[0] != ' ':
                continue
            res_label = "{}{}".format(residue.get_resname(), residue.id[1])
            sasa = residue.sasa
            total_sasa += sasa
            print "{:<10} {:<6} {:<10.2f}".format(res_label, chain.id, sasa)

print "-" * 30
print "TOTAL PROTEIN SASA: {:.2f} A^2".format(total_sasa)
