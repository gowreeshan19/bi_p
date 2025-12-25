import os
import subprocess
import matplotlib.pyplot as plt
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# --- 1. LOCAL SETTINGS ---
local_fasta = r"C:\Users\T.GOWREESHAN\Documents\Python\bioinformatics_project\msa_sequence_02.fasta" 
clustal_exe = r"C:\Users\T.GOWREESHAN\Documents\Python\software\clustal-omega\clustalo.exe" 
output_alignment = "aligned_local.fasta"

# --- 2. ALIGNMENT (Local) ---
print(f"\nStep 1: Aligning local file: {local_fasta}...")
cmd = [clustal_exe, "-i", local_fasta, "-o", output_alignment, "--auto", "--force"]

try:
    # Run the alignment
    subprocess.run(cmd, check=True) 
    print(f"Alignment complete! Results saved to '{output_alignment}'.")
    
    # --- 3. TREE CONSTRUCTION ---
    print("\nStep 2: Building Phylogenetic Tree...")
    alignment = AlignIO.read(output_alignment, "fasta") #
    
    # FIX: Use 'identity' model to avoid "Model not supported" error
    calculator = DistanceCalculator('identity') 
    dm = calculator.get_distance(alignment) 

    # Construct tree using Neighbor-Joining (NJ)
    constructor = DistanceTreeConstructor() 
    tree = constructor.nj(dm) 

    # --- 4. AUTOMATED ROOTING ---
    # Automatically select the last sequence as outgroup
    outgroup_id = alignment[-1].id 
    print(f"Step 3: Rooting the tree with: {outgroup_id}")
    
    try:
        tree.root_with_outgroup(outgroup_id) #
    except Exception as root_error:
        print(f"Rooting failed, showing unrooted tree: {root_error}")

    # --- 5. VISUALIZATION ---
    print("\nGenerating graphical tree window...")
    # Use Matplotlib for professional output
    fig = plt.figure(figsize=(10, 8), dpi=100) 
    ax = fig.add_subplot(1, 1, 1)
    plt.title(f"Phylogeny: {os.path.basename(local_fasta)}")
    
    Phylo.draw(tree, axes=ax)
    plt.show() # This opens the pop-up window

except Exception as e:
    print(f"\nAn error occurred: {e}")