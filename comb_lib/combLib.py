#Thalita Cirino on Monday, April 10th 2023 writes:

"""Develop a Python program that generates a combinatorial library of chemical 
structures based on a given molecular backbone and sets of substituents. Validate 
the final structures and write the valid ones in an external file with their RO5
properties and violation"""

import re
import csv
from typing import List, Dict
from rdkit import Chem
from lipinskiRO5f import lipinskiRules

# L-1/4: function to ask the user to input a SMILES backbone string.
# The backbone SMILES is input by the user. e.g. CC1C(X)CC(W)CC(Z)C1C
def import_backbone() -> str:
    while True:
        backbone = input('Insert your backbone structure with up to 4 insertion points: (W), (X), (Y) and (Z). > ').upper()
        insertion_points = {'(W)', '(X)', '(Y)', '(Z)'}
        if any(ip in backbone for ip in insertion_points):
            return backbone
        print('Structure invalid. At least one insertion point must be included.')

# L-2/4: function to import substituents list from the text file.
# A text file for each insertion point must be created and placed in the same folder as the code.
def import_substituents(subst_filename: str) -> List[str]:
    try:
        with open(subst_filename, 'r') as s:
            subst_list = [line.strip() for line in s if line.strip()]
        if not subst_list:
            raise ValueError("The text file cannot be empty.")
        return subst_list
    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}. Try again.")
        return import_substituents(input("Enter the correct filename: "))

# L-3/4: function to update numbers in case substituent and backbone have cycles.

def updateCycles(backbone, subs):

    bb_numbers = [int(num) for num in re.findall(r'\d+', backbone)] #Find numbers in the backbone
    max_bb_number = max(bb_numbers) if bb_numbers else 0 # Total number of cycles contained in the backbone.
    
    subs_numbers = set(int(num) for num in re.findall(r'\d+', subs)) #Find numbers in the substituents

    # In order to replace only once:
    for num in subs_numbers:
        subs = subs.replace(str(num), '$'+str(num))

    # loop that updates the cycle numbers in the substituents:
    i = 0
    for num in subs_numbers:
        subs = subs.replace('$'+str(num), str(max_bb_number + i + 1)) # Update cycle numbers in the substituent
        i += 1
        
    return subs

# L-4/4: function to build the combinatory library.
def combinatory_library():
    backbone = import_backbone()
    substitutions = [('(W)', 'W'), ('(X)', 'X'), ('(Y)', 'Y'), ('(Z)', 'Z')]
    substituents_dict = {}

    for insertion_point, prompt in substitutions:
        if insertion_point in backbone:
            filename = input(f"Enter filename for {prompt} substituents: ")
            substituents_dict[insertion_point] = import_substituents(filename)

    comb_lib = [backbone]
    for insertion_point, substituents in substituents_dict.items():
        comb_lib = [
            mol.replace(insertion_point, updateCycles(mol, sub))
            for mol in comb_lib
            for sub in substituents
        ]
    print(comb_lib)

    # Write results to a file
    with open('combinatory_library.csv', 'w',newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["SMILES","molec_form","MW","HBA","HBD","logP","Lip_violation"])
        invalid_smi = 0
        for smi in comb_lib:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                output_row = [smi] + lipinskiRules(smi)
                writer.writerow(output_row)
            else:
                invalid_smi +=1
                continue
            
    
    print(f"Generated {len(comb_lib)-invalid_smi} molecules. {invalid_smi} were invalid SMILES. Results written to combinatory_library.txt")

if __name__ == "__main__":
    combinatory_library()

