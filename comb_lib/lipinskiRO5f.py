import csv
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, Lipinski
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
#from validSMILES import processMolecule

""" Calculate Lipinski's Rule of Five parameters (molecular weight, number of 
hydrogen bond acceptors and donors, and LogP), count Lipinski violations"""

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
        
        return [molec_formula, mw, hba, hbd, logP, counter]









