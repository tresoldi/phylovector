import itertools
import math

import ete3
import networkx as nx

# TODO: allow different types of branch length (ratio, age, etc.)


def pathgraph2tree(P):
    # Instantiate the tree to be returned
    tree = ete3.Tree()

    def max_edge_weight(G):
        "Return the edge with the largest weight."
        return max(
            [(edge, weight) for *edge, weight in G.edges.data("weight")],
            key=lambda x: x[-1],
        )

    # Iterate over all edges
    if len(P.edges) <= 1:
        # Last case
        for v in P.nodes:
            tree.add_child(name=v)
    else:
        # Remove the edge with the maximum weight
        P = P.copy()
        edge, weight = max_edge_weight(P)
        P.remove_edge(*edge)

        # Do the same to the remaining subpaths, recursively
        for comp in nx.connected_components(P):
            component = P.subgraph(comp).copy()
            tree.add_child(pathgraph2tree(component))

    # Process the current tree structure to produce a standardized topology: nodes
    # with only one child are removed and multifurcations are automatically resolved.
    tree.standardize()

    return tree


# TODO: make it work with nameless leaves
def vector2tree(vector, leaves):

    # Compute the number of taxa in the tree by solving the equation `n = (sqrt(8*l+1)-1)/2`,
    # first obtained by solving the combinational one for `n` and check if the value
    # is appropriate.
    # TODO: check if less than 2 leafs
    n = (math.sqrt(8 * len(vector) + 1) - 1) / 2
    if int(n) != n:
        raise ValueError(f"invalid vector length: {len(vector)}.")
    else:
        n = int(n)

    # Build the graph; note that we don't need to add nodes directly, as
    # `.add_edge()` takes care of it; we also build the ultrametric
    # dictionary in this step (using the offset `n` from above)
    G = nx.Graph()
    ultrametric_d = {}  # TODO: rename later
    for idx, (leave_i, leave_j) in enumerate(itertools.combinations(leaves, 2)):
        G.add_edge(leave_i, leave_j, weight=vector[idx + n])

        ultrametric_d[leave_i, leave_j] = vector[idx + n]
        ultrametric_d[leave_j, leave_i] = vector[idx + n]

    # Build the path graph `P` from the list of nodes
    P = nx.Graph()
    L = set(leaves)  # TODO: sorted list?
    leaf_i = L.pop()  # start for the last leaf in the set
    while L:
        # Find a leaf node `j` (and its weight) among the set of leaves for which the distance
        # D(i, j) is the minimum, placing (i, j) in the path graph
        leaf_j, w = min(
            [(leaf, ultrametric_d[leaf_i, leaf]) for leaf in L], key=(lambda x: x[-1])
        )
        P.add_edge(leaf_i, leaf_j, weight=w)
        L.remove(leaf_j)

        # Update search node
        leaf_i = leaf_j

    # Build the tree from the path graph; note that this tree is still ultrametric
    tree = pathgraph2tree(P)

    # Iterate over all internal nodes (i.e., everything except leaves and root) and set
    # their lengths (i.e., distance from the ancestor) to the maximum length in the tree
    # (i.e., root to leaves) minus the maximum length between descendants
    max_length = max(ultrametric_d.values())
    for leaf in tree.iter_descendants():
        if leaf.is_leaf():
            continue

        # Compute the path distance from the current node to the root; note
        # that ETE3 sets the distance of the root to 1.0, so we need to
        # subtract this
        # TODO: can we just use distance to root?
        path = []
        n = leaf
        while n.up:
            path.append(n.dist)
            n = n.up
        path_length = sum(path) - 1.0

        # Get all descendant leaves and their maximum distance, which
        # is the age of the current node
        node_height = max(
            [
                ultrametric_d[leaf_i, leaf_j]
                for leaf_i, leaf_j in itertools.combinations(leaf.get_leaf_names(), 2)
            ]
        )

        # Finally set the branch length (i.e., node distance)
        leaf.dist = max_length - path_length - node_height

    # Set all leave branches to their minimal distance
    for leaf_idx, leaf in enumerate(leaves):
        tree_node = tree & leaf
        tree_node.dist = vector[leaf_idx]

    return tree


def tree2vector(source_tree):
    # Make a copy of the tree for manipulation, as we (might) extend the
    # branch leaves to obtain the ultrametric topology. While `ete3`
    # offers copy methods, we take the simpler approach of rebuilding
    # from the Newick representation.
    tree = ete3.Tree(source_tree.write())

    # Obtain list of leaves and lengths
    leaves = {leaf.name: leaf.dist for leaf in tree.iter_leaves()}
    lengths = [leaves[key] for key in sorted(leaves)]

    # Compute the root age (= maximum leaf distance) and extend all branches
    # as necessary, keeping track of the original value
    root_age = max([tree.get_distance(tree, leaf) for leaf in leaves])
    for leaf in leaves:
        # TODO: could just add any time, without the checking
        diff = root_age - tree.get_distance(tree, leaf)
        if diff:
            leaf = tree & leaf
            leaf.dist += diff

    # Build ultrametric part of the vector
    vector = lengths + []
    for leaf_i, leaf_j in itertools.combinations(sorted(leaves), 2):
        comm_anc = tree.get_common_ancestor(leaf_i, leaf_j)
        vector.append(tree.get_distance(comm_anc, leaf_i))

    return vector
