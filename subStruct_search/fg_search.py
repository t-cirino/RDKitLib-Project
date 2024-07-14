import csv
from rdkit import Chem

"""This script analyzes SMILES structures for the presence of 
common functional groups. It reads SMILES strings from a CSV file, 
identifies specific functional groups using SMARTS patterns with 
RDKit library, and outputs the results to a new CSV file."""

# A dictionary with SMARTS patterns for some common functional groups
functional_groups = {
    'Alcohol': '[OX2H]',
    'Aldehyde': '[CX3H1](=O)[#6]',
    'Ketone': '[#6][CX3](=O)[#6]',
    'Carboxylic Acid': '[CX3](=O)[OX2H1]',
    'Amine': '[NX3;H2,H1;!$(NC=O)]',
    'Amide': '[NX3][CX3](=[OX1])[#6]',
    'Ether': '[OD2]([#6])[#6]',
    'Ester': '[#6][CX3](=O)[OX2H0][#6]',
    'Nitro': '[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]',
    'Sulfonamide': '[#16X4](=[OX1])(=[OX1])([NX3R0,NX4R0;!$([N][!#6])][#6])',
    'Aminoacid': '[NH2,NHR,NR2][C@H,CH]([*])[C](=O)[OH]'
}

# Function to search for substructures in a molecule (considers that SMILES is already validated)
def findFunctionalGroup(mol, substructure_dict):

    functg_dict = {}
    for functg_group, smarts in substructure_dict.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            print(f"Error: Invalid SMARTS pattern for {functg_group}: {smarts}")
            continue
        if mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
            functg_dict[functg_group] = len(matches)  # Count of matches
    return functg_dict

#Function to read a list of SMILES strings and return it in a list
def readcsv(infile_name):
    smiles = []
    with open(infile_name,'r',newline='') as f:
        next(f)
        reader = csv.reader(f)

        for row in reader:
            smiles.append(row[0:2])

    return smiles

#Function to create the outfile and write down the data
def writecsv(outfile_name, data):
    with open(outfile_name,'w',newline='') as f:
        writer = csv.writer(f)
        for row in data:
            writer.writerow(row)


# Creating the output nested list to write in the file
output = [["SMILES","name"]+[fg for fg in functional_groups]]
smiles = readcsv("smiles_names.csv")

for smi,name in smiles:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        output.append(["Invalid SMILES entry",name,"-","-","-","-","-","-","-","-","-","-"])
    else:   
        Chem.RemoveStereochemistry(mol)
        fg_dict = findFunctionalGroup(mol, functional_groups)
        output_row = [smi, name] + [fg_dict.get(fg, 0) for fg in functional_groups]
        output.append(output_row)

writecsv("fg_matches.csv", output)

