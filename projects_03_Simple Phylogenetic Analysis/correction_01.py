import os
import time
import subprocess
import matplotlib.pyplot as plt
from Bio import SeqIO, Entrez, AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# --- 1. SETTINGS & SEARCH ---
Entrez.email = "gowreeshangowree@gmail.com"
search_term = input("Enter the search term (e.g., COI[Gene] AND primates[Organism]): ")

# Define the path to your Clustal Omega executable
clustal_exe = r"C:\Users\T.GOWREESHAN\Documents\Python\software\clustal-omega\clustalo.exe" 

print(f"\nStep 1: Searching NCBI for: {search_term}...")
with Entrez.esearch(db="nucleotide", term=search_term, retmax=10) as handle:  #data retrieval limit set to 10 for demo purposes and handle NCBI server load
    search_record = Entrez.read(handle) #  handle result in to readeable format
    id_list = search_record["IdList"] # extract list of sequence IDs from search result

if not id_list: # check if any IDs were found 
    print("No sequences found. Please try a different search term.")
    exit() # exit program if no sequences found

# --- 2. FETCH GENBANK RECORDS & RENAME WITH SPECIES NAMES ---
# We fetch GenBank (.gb) instead of FASTA to access the 'source' organism field

records = [] # list to hold renamed SeqRecord objects
print(f"Found {len(id_list)} IDs. Fetching names and sequences...") # iterate over each sequence ID 

for seq_id in id_list: # fetch GenBank record for each ID
    with Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text") as handle:# fetch GenBank record 
        seq_record = SeqIO.read(handle, "genbank")# read record into SeqRecord object
        
        # Extract the scientific name and replace spaces with underscores
        species_name = seq_record.annotations.get("source", "Unknown_Species") # get species name from annotations
        clean_name = species_name.replace(" ", "_") # clean species name for file compatibility
        
        # Rename the ID so the tree labels show the species name
        seq_record.id = f"{clean_name}_{seq_record.id}" # rename SeqRecord ID with species name 
        seq_record.description = ""  # Clear description for clean FASTA headers
        records.append(seq_record)# add renamed record to list
        print(f"  Successfully identified: {species_name}")
    
    # Brief pause to respect NCBI server limits
    time.sleep(0.3) # pause between requests to avoid overloading NCBI servers

# Save renamed records to FASTA for alignment
SeqIO.write(records, "my_data.fasta", "fasta") # write records to FASTA file 
print("\nFile 'my_data.fasta' saved with species names.")

# --- 3. ALIGNMENT (Using Modern Subprocess) ---
# Subprocess is the modern way to run external tools in Biopython
print("\nStep 2: Starting Multiple Sequence Alignment (Clustal Omega)...")# prepare command for Clustal Omega
cmd = [clustal_exe, "-i", "my_data.fasta", "-o", "aligned_data.fasta", "--auto", "--force"]

try:
    subprocess.run(cmd, check=True)# run the command and check for errors
    print("Alignment complete! Results saved to 'aligned_data.fasta'.")
    
    # --- 4. SCIENTIFIC TREE CONSTRUCTION ---
    print("\nStep 3: Building Phylogenetic Tree...")
    alignment = AlignIO.read("aligned_data.fasta", "fasta") # read aligned sequences from file 
    
    # Calculate genetic distance using the 'identity' model
    calculator = DistanceCalculator('identity') # identity model for distance calculation
    dm = calculator.get_distance(alignment) # compute distance matrix from alignment


    # Construct tree using Neighbor-Joining (NJ)
    # NJ is more accurate for evolution than UPGMA
    constructor = DistanceTreeConstructor() # initialize tree constructor  
    tree = constructor.nj(dm) # build tree using NJ method

# --- 4.5 AUTOMATED ROOTING ---
    # Automatically pick the very last sequence downloaded as the outgroup
    outgroup_record = records[-1] # select last record as outgroup
    outgroup_id = outgroup_record.id # get the ID of the outgroup record

    print(f"\nStep 4: Automatically rooting the tree with: {outgroup_id}")
    
    try:
        # Perform the rooting
        tree.root_with_outgroup(outgroup_id) # root the tree using the outgroup ID
        print("Tree successfully rooted!") # confirm successful rooting
    except Exception as root_error:
        print(f"Automatic rooting failed: {root_error}")
        print("Displaying unrooted tree instead.")
    # Save the tree structure to a Newick file
    Phylo.write(tree, "coi_tree.nwk", "newick") # write tree to Newick file
    print("Phylogenetic tree constructed and saved to 'coi_tree.nwk'.")

    # --- 5. VISUALIZATION ---
    print("\n--- Phylogenetic Tree (ASCII) ---") # display ASCII representation of the tree
    Phylo.draw_ascii(tree) # print ASCII tree to console

    # This creates the graphical window (requires 'pip install matplotlib')
    print("\nGenerating graphical tree window...") # visualize tree using Matplotlib
    Phylo.draw(tree) # draw tree graphically
    plt.show() # show the plot window


except Exception as e: # catch any errors during alignment or tree construction
    print(f"\nAn error occurred: {e}") # print error message
    print("Please ensure Clustal Omega is correctly installed and the path is set.")