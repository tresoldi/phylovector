import itertools
import math
from typing import List, Union

import ete3
import networkx as nx

# TODO: allow different types of branch length (ratio, age, etc.)


def pathgraph2tree(P):
    """
    Converts a path graph to a tree structure using the ETE toolkit.

    This function iteratively removes the edge with the maximum weight in the graph
    and applies the same process to the remaining subgraphs recursively, until all
    nodes are processed and added to the tree. The final tree is then standardized
    to ensure a proper phylogenetic tree structure, where nodes with only one child
    are removed, and multifurcations are resolved.

    Parameters:
    - P (networkx.Graph): A path graph where nodes represent taxa and edges represent
      phylogenetic distances. Each edge must have a 'weight' attribute. Note that the
      graph is modified in-place during the process.

    Returns:
    - ete3.Tree: An ETE tree object representing the phylogenetic tree constructed
      from the path graph.

    Note:
    This function assumes that the input graph is connected and correctly represents
    a phylogenetic tree where the edge with the maximum weight can be considered as
    the most significant bifurcation point at each step of the tree construction process.
    """

    def _max_edge_weight(G: nx.Graph) -> list:
        """
        Find the edge with the largest weight in a given graph and return it.

        Iterates over all edges in the graph, which are expected to have a 'weight' attribute,
        and identifies the edge with the maximum weight. This version returns only the edge,
        not the weight.

        Parameters:
        - G (nx.Graph): A networkx graph object where edges have a 'weight' attribute.

        Returns:
        - list: A list with two elements, representing the edge with the maximum weight.
        """

        return max(
            [(edge, weight) for *edge, weight in G.edges.data("weight")],
            key=lambda x: x[-1],
        )[0]

    # Instantiate the tree to be returned
    tree = ete3.Tree()

    # Iterate over all edges
    if len(P.edges) <= 1:
        # Last case
        for v in P.nodes:
            tree.add_child(name=v)
    else:
        # Remove the edge with the maximum weight
        edge = _max_edge_weight(P)
        P.remove_edge(*edge)

        # Do the same to the remaining subpaths, recursively
        # TODO: investigate if it is possible to do this without making copies
        # of each subgraph
        for comp in nx.connected_components(P):
            component = P.subgraph(comp).copy()
            tree.add_child(pathgraph2tree(component))

    # Process the current tree structure to produce a standardized topology: nodes
    # with only one child are removed and multifurcations are automatically resolved.
    tree.standardize()

    return tree


def vector2tree(
    vector: list[float], leaves: Union[List[str], None] = None
) -> ete3.Tree:
    """
    Converts a vector representing pairwise into a phylogenetic tree.

    Parameters:
    - vector (list[float]): A flat list with the vector representation of the tree,
        where the first `n` elements represent the branch lengths of the tree and the
        remaining `n(n-1)/2` elements represent the pairwise distances between leaves.
    - leaves (list[str]): Names of leaves in the tree. The order of the names must
        correspond to the order of the pairwise distances in the vector.

    Returns:
    - ete3.Tree: An ETE tree object representing the constructed phylogenetic tree.

    Raises:
    - ValueError: If the vector's length does not correspond to a valid complete graph or if there are fewer than 2 leaves.
    """

    # Compute the number of taxa in the tree by solving the equation `n = (sqrt(8*l+1)-1)/2`,
    # first obtained by solving the combinational one for `n` and check if the value
    # is appropriate.
    num_leaves = (math.sqrt(8 * len(vector) + 1) - 1) / 2
    if not num_leaves.is_integer() or num_leaves < 2:
        raise ValueError(
            f"Invalid vector length: {len(vector)}. Vector does not correspond to a valid complete graph with 2 or more leaves."
        )

    num_leaves = int(num_leaves)
    if not leaves:
        leaves = [f"Leaf_{i}" for i in range(num_leaves)]
    elif len(leaves) != num_leaves:
        raise ValueError(
            f"The number of provided leaf names does not match the expected number of leaves ({num_leaves})."
        )

    # Build the graph; note that we don't need to add nodes directly, as
    # `.add_edge()` takes care of it; we also build the ultrametric
    # dictionary in this step (using the offset `n` from above)
    graph = nx.Graph()
    leaf_distances = {}
    for idx, (leave_i, leave_j) in enumerate(itertools.combinations(leaves, 2)):
        # Add the edge to the graph with the weight (i.e., distance) from the vector
        weight = vector[idx + num_leaves]
        graph.add_edge(leave_i, leave_j, weight=weight)

        # Store the pairwise distances for later use; both directions are stored
        # to simplify the process of building the path graph and save the costs of
        # sorting the names
        leaf_distances[leave_i, leave_j] = weight
        leaf_distances[leave_j, leave_i] = weight

    # Build the path graph from the list of nodes
    path_graph = nx.Graph()
    remaining_leaves = leaves[:]  # make a copy
    leaf_i = remaining_leaves.pop()  # start from the last leaf in the set
    while remaining_leaves:
        # Find a leaf node `j` (and its weight) among the set of leaves for which the distance
        # D(i, j) is the minimum, placing (i, j) in the path graph
        leaf_j, weight = min(
            [(leaf, leaf_distances[leaf_i, leaf]) for leaf in remaining_leaves],
            key=(lambda x: x[1]),
        )
        path_graph.add_edge(leaf_i, leaf_j, weight=weight)
        remaining_leaves.remove(leaf_j)

        # Update search node
        leaf_i = leaf_j

    # Build the tree from the path graph; note that this tree is still ultrametric
    tree = pathgraph2tree(path_graph)

    # Iterate over all internal nodes (i.e., everything except leaves and root) and set
    # their lengths (i.e., distance from the ancestor) to the maximum length in the tree
    # (i.e., root to leaves) minus the maximum length between descendants
    root_age = max(leaf_distances.values())
    for leaf in tree.iter_descendants():
        # Get all descendant leaves and their maximum distance, which
        # is the age of the current node
        descendant_leaves = leaf.get_leaf_names()
        distances = [
            leaf_distances[leaf_i, leaf_j]
            for leaf_i, leaf_j in itertools.combinations(descendant_leaves, 2)
        ]

        # If there are no distances, the node is a leaf and we can skip
        if distances:
            node_height = max(distances)
        else:
            node_height = 0

        # Finally set the branch length (i.e., node distance)
        age_anc = tree.get_distance(leaf.up)
        age_curr = root_age - node_height
        leaf.dist = age_curr - age_anc

    # Set all leave branches to lengths
    for leaf_idx, leaf in enumerate(leaves):
        tree_node = tree & leaf
        tree_node.dist -= vector[leaf_idx]

    return tree


def tree2vector(source_tree, components=False, epsilon=1e-6):
    # Make a copy of the tree for manipulation, as we (might) extend the
    # branch leaves to obtain the ultrametric topology. While `ete3`
    # offers copy methods, we take the simpler approach of rebuilding
    # from the Newick representation.
    tree = ete3.Tree(source_tree.write())

    # Obtain list of leaves
    leaves = sorted([leaf.name for leaf in tree.iter_leaves()])

    # Compute the root age (= maximum leaf distance) and extend all branches
    # as necessary, also collecting the `lengths` portion of the vector
    root_age = max([tree.get_distance(tree, leaf) for leaf in leaves])
    lengths = []
    for leaf in leaves:
        # Compute the difference (age/height of node) and append to
        # the lead node distance -- note that this is done even if
        # the node is alive (i.e, difference of zero)
        leaf_node = tree & leaf

        diff = root_age - tree.get_distance(tree, leaf_node)
        leaf_node.dist += diff
        lengths.append(diff)

    # Build ultrametric part of the vector
    ultrametric = []
    for leaf_i, leaf_j in itertools.combinations(leaves, 2):
        comm_anc = tree.get_common_ancestor(leaf_i, leaf_j)
        ultrametric.append(tree.get_distance(comm_anc, leaf_i))

    # For all values in `lengths` and `ultrametric`, if they are very small and close to
    # zero, set them to zero
    lengths = [0 if math.isclose(val, 0, abs_tol=epsilon) else val for val in lengths]
    ultrametric = [
        0 if math.isclose(val, 0, abs_tol=epsilon) else val for val in ultrametric
    ]

    if components:
        return lengths, ultrametric

    return lengths + ultrametric
