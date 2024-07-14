import csv
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, Lipinski
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from validSMILES import processMolecule

""" Calculate Lipinski's Rule of Five parameters (molecular weight, number of 
hydrogen bond acceptors and donors, and LogP), count Lipinski 
violations, and output the results to a new CSV file."""

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

#Function to calculate the RO5 and count the Lipinski violations
def lipinskiRules(smiles_str):
    mol = Chem.MolFromSmiles(smiles_str)
    molec_formula = CalcMolFormula(mol)
    counter = 0
    if mol is not None:

        mw = Descriptors.ExactMolWt(mol)
        if mw > 500:
            counter += 1

        hba = Descriptors.NOCount(mol)
        if hba > 10:
            counter += 1

        hbd = Descriptors.NHOHCount(mol)
        if hbd > 5:
            counter += 1
        
        logP = Descriptors.MolLogP(mol)
        if logP > 5:
            counter += 1
        
        return molec_formula, mw, hba, hbd, logP, counter



#Create the output file with "clean" SMILES and Lipinski RO5 violations count 
output_data = [["SMILES","QSAR_SMILES","molec_form","MW","HBA","HBD","logP","Lip_violation"]]
smiles = readcsv("example_smiles.csv")
for smi in smiles:
    clean_smi = processMolecule(smi)
    smi_lip_rules = [clean_smi]
    if clean_smi is None:
        output_data.append([smi,"-", "-","-","-","-","-"])
    else:
        smi_lip_rules.append(clean_smi)
        smi_lip_rules.extend(lipinskiRules(smi))
        output_data.extend([smi_lip_rules])

writecsv("outfile.csv", output_data)









