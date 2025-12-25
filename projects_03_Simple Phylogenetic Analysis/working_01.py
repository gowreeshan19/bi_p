from Bio import SeqIO
from Bio import Entrez
from Bio.Align.Applications import ClustalOmegaCommandline
import os
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
import matplotlib.pyplot as plt
import time

Entrez.email = "gowreeshangowree@gmail.com"  # Always tell NCBI who you are
search_term=input("Enter the search term: ")
handle=Entrez.esearch(db="nucleotide", term=search_term, retmax=10)
record=Entrez.read(handle)
id_list=record["IdList"]



# Fetch the FASTA records and save them to a file

# Fetch the FASTA records
with Entrez.efetch(db="nuccore", id=id_list, rettype="fasta", retmode="text") as net_handle:
    # Open the output file using 'with'
    with open("coi_sequences.fasta", "w") as out_handle:
        out_handle.write(net_handle.read())

print("File 'coi_sequences.fasta' has been saved successfully.")



# Set the paths to your files
# If clustalo.exe is in the same folder, just use 'clustalo'
# Otherwise, provide the full path like: r"C:\path\to\clustalo.exe"
# Update this line to include 'clustalo.exe' at the end
clustal_exe = r"C:\Users\T.GOWREESHAN\Documents\Python\software\clustal-omega\clustalo.exe" 
in_file = "coi_sequences.fasta"
out_file = "aligned_coi.fasta"

# Create the command
# 'auto=True' lets Clustal choose the best settings automatically
clustal_cline = ClustalOmegaCommandline(cmd=clustal_exe, infile=in_file, outfile=out_file, verbose=True, auto=True)

# Run the command
print("Starting alignment... this may take a moment.")
stdout, stderr = clustal_cline()

print("Alignment complete! Results saved to", out_file)



# 1. Load the aligned sequences
# Use 'fasta' or 'clustal' depending on your Step 2 output format
alignment = AlignIO.read("aligned_coi.fasta", "fasta")

# 2. Calculate the Distance Matrix
# The 'identity' model calculates distance based on sequence mismatch
calculator = DistanceCalculator('identity')
distance_matrix = calculator.get_distance(alignment)

# 3. Construct the Tree
# We will use the UPGMA algorithm for this example
constructor = DistanceTreeConstructor(calculator)
tree = constructor.upgma(distance_matrix) # You can also use 'nj' for Neighbor-Joining

# 4. Save the tree to a Newick file
Phylo.write(tree, "coi_tree.nwk", "newick")

# 5. Visualize the Tree
# Draw a simple text version in the terminal
print("\n--- Phylogenetic Tree (ASCII) ---")
Phylo.draw_ascii(tree)

# Draw a nice graphical version (requires matplotlib)
print("\nGenerating graphical tree...")
Phylo.draw(tree)