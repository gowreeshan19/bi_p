from Bio.Seq import Seq
from Bio import SeqIO

for Seq_record in SeqIO.parse("/Users/T.GOWREESHAN/Documents/Python/New folder/sequence_01.fasta", "fasta"):  
    my_coading_seq=Seq_record.seq
    print("Original DNA sequence:", my_coading_seq)
    print(my_coading_seq.translate(to_stop=False)) # Translate the entire coding sequence to protein
# List to store all potential protein sequences
    all_proteins =[]
#intract through all three reading frames (0,1,2)   
    for frame in range(3):
        remainder= len(my_coading_seq[frame:])%3
        if remainder !=0:
            coding_seq=my_coading_seq[frame:-remainder]
        else:
            coding_seq=my_coading_seq[frame:]
        protein_seq=coding_seq.translate(to_stop=False) # Translate the coding sequence to protein
        true_protein=protein_seq.split("*") # Split protein sequence at stop codons
        all_proteins.extend(true_protein) # Add proteins from this frame to the list


        if all_proteins:
         longest_protein = max(all_proteins, key=len)          
         print("Total potential ORFs found:",len(all_proteins))
         print("Longest Protein Sequence:",len(longest_protein))
         print(longest_protein)
        
        else:
         print("No potential ORFs found.")
