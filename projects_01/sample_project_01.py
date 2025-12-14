from Bio import Entrez, SeqIO
from Bio.SeqUtils import gc_fraction

# 1. SETUP
Entrez.email = "gowreeshangowree@gmail.com"

# 2. ASK USER FOR SEARCH TERM
# Example: Type "Insulin Homo sapiens" or "Hemoglobin E. coli"
term = input("Enter gene name/organism to search for: ")

print("Searching for:",term)

# 3. SEARCH (Find the ID)
# retmax=1 means "just give me the top result"
search_handle = Entrez.esearch(db="nucleotide", term=term, retmax=20)
search_results = Entrez.read(search_handle)
search_handle.close()

# Check if we found anything
if not search_results["IdList"]:
    print("No results found. Try a different search term.")
else:
    # Get the first ID found
    best_id = search_results["IdList"][10]
    print("Found Best Match ID:",best_id)

    # 4. FETCH (Download the sequence using that ID)
    print("Downloading sequence...")
    fetch_handle = Entrez.efetch(db="nucleotide", id=best_id, rettype="fasta", retmode="text")
    record = SeqIO.read(fetch_handle, "fasta")
    fetch_handle.close()

    # 5. ANALYZE
    my_seq = record.seq
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