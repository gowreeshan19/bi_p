"""from rdkit import Chem
from rdkit.Chem import AllChem

# 1. Create a molecule for Aspirin using its SMILES string
smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
mol = Chem.MolFromSmiles(smiles)

# 2. Add Hydrogen atoms (crucial for docking!)
mol_with_h = Chem.AddHs(mol)

# 3. Generate a 3D shape (conformation)
AllChem.EmbedMolecule(mol_with_h, AllChem.ETKDG())

print(f"Successfully processed molecule with {mol_with_h.GetNumAtoms()} atoms.")

try:
    from openbabel import pybel
    print("Success! Open Babel is ready for docking.")
except ImportError:
    print("Installation failed. Ensure the .exe is installed and PATH is set.")"""
from rcsbsearchapi.search import AttributeQuery

# 1. Define your variables
ECNumber = "3.4.21.4"

# 2. Build queries manually using strings (Avoids the 'rcsb_attributes' crash)
# Note: We use the exact field names as strings
q1 = AttributeQuery("rcsb_polymer_entity.rcsb_ec_lineage.id", operator="exact_match", value=ECNumber)
q2 = AttributeQuery("chem_comp.formula_weight", operator="greater_or_equal", value=300)
q3 = AttributeQuery("chem_comp.formula_weight", operator="less_or_equal", value=800)

# 3. Combine them
query = q1 & q2 & q3

resultL = list(query())           # assign the results of the query to a list variable

print(resultL[0:10])              # list the first 10 results

print("There are", len(resultL), "trypsin structures that contain ligands in the RCSB PDB.")