import re
import json
import gzip
from rxnmapper import RXNMapper
rxn_mapper = RXNMapper()

def map_reaction(rxn_smiles):
    mapped_rxn = rxn_mapper.get_attention_guided_atom_maps([rxn_smiles])

    results = mapped_rxn[0]['mapped_rxn'].split('>>')

    precursors = results[0]
    precursor_list = precursors.split('.')

    reagents = [molecule for molecule in precursor_list if not re.search(':[0-9]+', molecule)]

    r_p =  rxn_smiles.split('>>')
    rs = r_p[0].split('.')
    reactant_smiles = set(rs) - set(reagents)
    reactant_smiles = [r for r in reactant_smiles if r != '']

    return sorted(reactant_smiles, key=len, reverse=True)

def filter_reactions(reactions_unfiltered):
    """takes a list of reactions strings 
    and removes the empty list"""
    reactions = [reaction for reaction in reactions_unfiltered if reaction]
    return reactions


def create_rxn_and_child_node(parent_node, reactions):
    """takes in a parent node and a list of reactions
    and recursively itterates down the smiles string 
    until the end of the synthetic sequence is reached"""

    if len(reactions) == 0:  # Base case: no further reactions

        return 
        #     {"smiles": parent_node['smiles'],
        #     "type": "mol",
        #     # You can set 'in_stock' to False by default or add logic to determine its value
        #     "in_stock": False,
        #     "children": []
        # }
        
    else:
        # print(reactions)
        next_reaction = reactions.pop(0)
        # print(next_reaction)
        children = next_reaction.split(".")
        # print(children)
        forward_rxn = next_reaction + ">>" + parent_node["smiles"]
        # print(forward_rxn)
        parent_node_plus_one = ['']
        try:
            parent_node_plus_one = map_reaction(forward_rxn)
            # print(f"this is parent node plus one{parent_node_plus_one}")
        except Exception as e:
            # print(f"Error processing reaction: {e}")
            breakpoint

        reactions_dic = {
                "smiles": '',
                "type": "reaction",
                "children": []
        }
        # print(f"This is children {children}")
        for child in children:
            # print(f'This is child {child}')
            # print(f"This is parent node {parent_node_plus_one}")
            if parent_node_plus_one[0] == child:
                 child_dic = {
                "smiles": child,
                "type": "mol",
                "in_stock": False,
                "children": []}
                #  print("Entering Recursion")
                 result = create_rxn_and_child_node(child_dic, reactions.copy())
                 if result is not None:
                    child_dic["children"].append(result)
                    reactions_dic["children"].append(child_dic)
                 
            else:
                child_dic = {
                    "smiles": child,
                    "type": "mol",
                    "in_stock": False
                }
                reactions_dic["children"].append(child_dic)

        return reactions_dic

trees = []

with open('output_routes.txt', 'r') as file:
    extract_next_line = False  # Initialize a flag to identify when to extract the next line
    for line in file:
        if extract_next_line:
            # This line is the one we want to extract after "Model generation ROOT:"
            extracted_part = line.strip()  # Extract and clean the line
            # print(extracted_part)
            tree['reactions'] = extracted_part  # Store the extracted part in the tree dictionary
            extract_next_line = False  # Reset the flag
            trees.append(tree)
        elif line.startswith("Model generation ROOT: "):
            tree = {}
            # Extract the part after the prefix and assign it to root
            root = line.split("Model generation ROOT: ")[1].strip()
            tree['root'] = root  # Store the root in the tree dictionary
            extract_next_line = True
print(len(trees))
# print(trees[0:10])


transformed_data = []

for tree in trees:
    root_node = {
        "smiles": tree['root'],
        "type": "mol",
        # You can set 'in_stock' to False by default or add logic to determine its value
        "in_stock": False,
        "children": []
    }
    # Split the reactions string into individual reaction steps
    unfiltered_reactions = tree['reactions'].split(">>")
    unfiltered_reactions = [reaction.strip() for reaction in unfiltered_reactions]
    reactions = filter_reactions(unfiltered_reactions)
    # print(f"These are the the reactions that fx is called first time on {reactions}")
    # print(f"This is the root_node that the fx is called first time on {root_node}")
    transformed_tree = create_rxn_and_child_node(root_node, reactions)
    root_node["children"].append(transformed_tree)

    transformed_data.append(root_node)

# Convert the transformed data to JSON format if needed
json_output = json.dumps(transformed_data, indent=4)
# print(json_output)

with open('output_routes.json', 'w') as json_file:
    json.dump(transformed_data, json_file, indent=4)

filename = 'output_file.json.gz'

# Write the JSON data to a compressed file
with gzip.open(filename, 'wt', encoding='utf-8') as f:
    f.write(json_output)
    
print("Data exported to output_routes.json")