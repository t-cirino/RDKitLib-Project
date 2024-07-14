import csv
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, Lipinski
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize

""" Pre-processing of data: Validate and clean SMILES by striping out salts, 
neutralising charged compounds, removing metals."""

#Function to read a list of SMILES strings and return it in a list
def readcsv(infile_name):
    smiles = []
    with open(infile_name,'r',newline='') as f:
        next(f)
        reader = csv.reader(f)

        for row in reader:
            smiles.append(row[0])

    return smiles

#Function to create the outfile and write down the data
def writecsv(outfile_name, data):
    with open(outfile_name,'w',newline='') as f:
        writer = csv.writer(f)
        for row in data:
            writer.writerow(row)

#Function to validate and "clean" the SMILES 
def processMolecule(smiles_str) -> str:
    # Create RDKit mol object
    mol = Chem.MolFromSmiles(smiles_str, sanitize=True)
    if mol is None:
        return None #SMILES string invalid

    # Remove salts and minor fragments
    remover = SaltRemover()
    mol = remover.StripMol(mol, dontRemoveEverything=True, sanitize=True)

    # Remove heavy metals (atomic number > 50)
    mol = Chem.DeleteSubstructs(mol, Chem.MolFromSmarts('[#51-]'))

    # Remove stereochemistry
    Chem.RemoveStereochemistry(mol)

    # Normalize charges
    uncharger = rdMolStandardize.Uncharger()
    mol = uncharger.uncharge(mol)

    clean_smiles = Chem.MolToSmiles(mol, canonical=True)

    return clean_smiles