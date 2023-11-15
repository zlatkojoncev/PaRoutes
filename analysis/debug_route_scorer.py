from typing import Dict, Any, Set, List, Tuple
from route_distances.utils.type_utils import StrDict
import numpy as np

def route_scorer(routes: List[StrDict]) -> Tuple[List[StrDict], List[float]]:
    """
    Scores and sort a list of routes.
    Returns a tuple of the sorted routes and their costs.

    :param routes: the routes to score
    :return: the sorted routes and their costs
    """
    print(len(routes))
    counter = 0
    for route in routes:
        print(route)
        counter += 1
        if counter > 10:
            break
    scores = np.asarray([route_score(route) for route in routes])
    sorted_idx = np.argsort(scores)
    routes = [routes[idx] for idx in sorted_idx]
    return routes, scores[sorted_idx].tolist()

route1 = {'smiles': 'CC1(C)CCC(C)(C)N1CC(O)COc1cc2ccccc2c(=O)c2ccsc12', 'type': 'mol', 'in_stock': False, 'children': [{'smiles': '', 'type': 'reaction', 'children': [{'smiles': 'O=c1c2ccccc2c(OCC2CO2)c2sccc12', 'type': 'mol', 'in_stock': False}, {'smiles': 'CC1(C)CCC(C)(C)N1', 'type': 'mol', 'in_stock': False}]}]}
print(type(route1))

route_scorer(route1)