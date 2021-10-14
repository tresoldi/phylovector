import ngesh
import phylovector


def stress_test1():
    for i in range(10000):
        if (i + 1) % 1000 == 0:
            print("Processing random tree #%i..." % i)

        # Generate a random tree, convert to a vector, and convert
        # back to a tree
        source_tree = ngesh.gen_tree(1.0, 0.5, max_time=2.0, labels="human")
        vector = phylovector.tree2vector(source_tree)
        leaves = sorted([leaf.name for leaf in source_tree.iter_leaves()])
        built_tree = phylovector.vector2tree(vector, leaves)

        # Get newick representations of both trees and sort them
        source_newick = source_tree.write(format=1)
        built_newick = built_tree.write(format=1)
        source_newick = phylovector.sorted_newick(source_newick)
        built_newick = phylovector.sorted_newick(built_newick)

        # Raise an error if they are different
        if source_newick != built_newick:
            print("error: [%s] [%s]" % (source_newick, built_newick))


if __name__ == "__main__":
    stress_test1()
