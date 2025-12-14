from Bio import SeqIO
from Bio import Entrez
from Bio.SeqUtils import gc_fraction

for Seq_record in SeqIO.parse("/Users/T.GOWREESHAN/Documents/Python/New folder/sequence_01.fasta", "fasta"):  
    my_seq=Seq_record.seq
    print("length of sequence:",len(my_seq)) # Calculate length of the sequence
    print("gc content:",gc_fraction(my_seq)*100,"%")# Calculate GC fraction and convert to percentage
    print("Original DNA sequence:", my_seq)
    print("reverse complement:",my_seq.reverse_complement()) # Get reverse complement of the DNA sequence

    print("transcription to mRNA:",my_seq.transcribe()) # Transcribe DNA to mRNA
    remainder= len(my_seq)%3
    if remainder !=0:
        protein=my_seq[:-remainder]
    else:
        protein=my_seq
    print("translation to Protein:",protein.translate()) # Translate mRNA to Protein

