import json
import re
from rxnmapper import RXNMapper
rxn_mapper = RXNMapper()


def map_reaction(parent_smiles, children_smiles):
    # Construct the reaction string
    reaction_string = ".".join(children_smiles) + ">>" + parent_smiles
    
    # Assuming rxn_mapper.get_attention_guided_atom_maps is accessible here
    mapped_rxn = rxn_mapper.get_attention_guided_atom_maps([reaction_string])

    results = mapped_rxn[0]['mapped_rxn'].split('>>')
    precursors = results[0]
    precursor_list = precursors.split('.')

    reagents = [molecule for molecule in precursor_list if not re.search(':[0-9]+', molecule)]
    r_p = reaction_string.split('>>')
    rs = r_p[0].split('.')
    reactant_smiles = set(rs) - set(reagents)
    reactant_smiles = [r for r in reactant_smiles if r != '']

    return sorted(reactant_smiles, key=len, reverse=True)

def create_reaction_node(smiles_string):
    if '>>' not in smiles_string:  # Base case: no further reactions
        return {
            "smiles": smiles_string,
            "type": "mol",
            "in_stock": False  # Adjust as needed
        }

    parts = smiles_string.split('>>')
    root_smiles = parts.pop(0)
    root_node = {
        "smiles": root_smiles,
        "type": "mol",
        "in_stock": False,  # Adjust as needed
        "children": []
    }

    for part in parts:
        child_smiles = part.split('.')
        reactant_smiles = map_reaction(root_smiles, child_smiles)
        for smile in reactant_smiles:
            reaction_node = {
                "type": "reaction",
                "children": [create_reaction_node(smile)]
            }
            root_node["children"].append(reaction_node)

    return root_node

# Export to JSON
with open('output_routes.json', 'w') as json_file:
    json.dump(reaction_tree, json_file, indent=4)

print("Data exported to output_routes.json")
