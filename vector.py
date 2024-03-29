import itertools
import math
import random

import ete3
import networkx as nx

# TODO: allow different types of branch length (ratio, age, etc.)


def pathgraph2tree(P):
    """
    Converts a path graph to a tree structure.

    This function iteratively removes the edge with the maximum weight in the graph
    and applies the same process to the remaining subgraphs recursively, until all
    nodes are processed and added to the tree. The final tree is then standardized
    to ensure a proper phylogenetic tree structure, where nodes with only one child
    are removed, and multifurcations are resolved.

    Parameters:
    - P (networkx.Graph): A path graph where nodes represent taxa and edges represent
      phylogenetic distances. Each edge must have a 'weight' attribute.

    Returns:
    - ete3.Tree: An ETE tree object representing the phylogenetic tree constructed
      from the path graph.

    Note:
    This function assumes that the input graph is connected and correctly represents
    a phylogenetic tree where the edge with the maximum weight can be considered as
    the most significant bifurcation point at each step of the tree construction process.
    """

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
        edge, _ = max_edge_weight(P)
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
    # TODO: check if less than 2 leaves
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
    # TODO: skip over leaves? (faster)
    # TODO: better solution for [0]? integrate with step below?
    root_age = max(ultrametric_d.values())
    for leaf in tree.iter_descendants():
        # Get all descendant leaves and their maximum distance, which
        # is the age of the current node
        node_height = max(
            [0]
            + [
                ultrametric_d[leaf_i, leaf_j]
                for leaf_i, leaf_j in itertools.combinations(leaf.get_leaf_names(), 2)
            ]
        )

        # Finally set the branch length (i.e., node distance)
        age_anc = tree.get_distance(leaf.up)
        age_cur = root_age - node_height
        leaf.dist = age_cur - age_anc

    # Set all leave branches to lengths
    for leaf_idx, leaf in enumerate(leaves):
        tree_node = tree & leaf
        tree_node.dist -= vector[leaf_idx]

    return tree

# TODO: replace very small numbers with zero, with a tolerance
def tree2vector(source_tree, components=False):
    # Make a copy of the tree for manipulation, as we (might) extend the
    # branch leaves to obtain the ultrametric topology. While `ete3`
    # offers copy methods, we take the simpler approach of rebuilding
    # from the Newick representation.
    tree = ete3.Tree(source_tree.write())

    # Obtain list of leaves and lengths
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

    if components:
        return lengths, ultrametric
    else:
        return lengths + ultrametric

def sorted_newick(newick):
    """
    Build a sorted representation of a Newick string.

    An internal function parses a Newick string by identifying tokens with
    a regular expression, which might fail for complex trees such as
    those carrying information other than branch length and node name.

    :param newick: The Newick tree to be sorted.

    :return: A corresponding but sorted Newick tree.
    """

    tree = ete3.Tree(newick)
    tree.sort_descendants()
    tree.ladderize()

    return tree.write(format=1)

def generate_random_tree(num_leaves=5):
    """
    Generates a simple random tree with a given number of leaves and random branch lengths.
    """
    tree = ete3.Tree()
    tree.populate(num_leaves)
    for node in tree.traverse():
        # Assign random branch lengths, for simplicity we use values between 0.1 and 1.0
        node.dist = round(random.uniform(0.1, 1.0), 2)
    return tree

def compare_trees(original_tree, reconstructed_tree):
    """
    Compares two trees based on their Newick representations after sorting.
    """
    original_newick = sorted_newick(original_tree.write(format=1))
    reconstructed_newick = sorted_newick(reconstructed_tree.write(format=1))
    return original_newick == reconstructed_newick

def main():
    random.seed(13)
    
    # Generate a random tree
    original_tree = generate_random_tree(num_leaves=4)
    original_tree_newick = sorted_newick(original_tree.write(format=1))
    print("Original Tree Newick Representation:", original_tree_newick)

    # Convert the tree to a vector
    vector = tree2vector(original_tree)
    print("\nGenerated Vector:", vector)

    # Reconstruct the tree from the vector
    leaves = [leaf.name for leaf in original_tree.iter_leaves()]
    reconstructed_tree = vector2tree(vector, leaves)
    reconstructed_tree_newick = sorted_newick(reconstructed_tree.write(format=1))
    print("\nReconstructed Tree Newick Representation:", reconstructed_tree_newick)

    # Compare the original and reconstructed trees
    if compare_trees(original_tree, reconstructed_tree):
        print("\nThe reconstructed tree matches the original tree.")
    else:
        print("\nThe reconstructed tree does NOT match the original tree.")

if __name__ == "__main__":
    main()