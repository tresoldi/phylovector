import math
import itertools
import ngesh
import phylovector


def stress_test1():
    INTERVAL = 1000

    for i in range(10000):
        if i % INTERVAL == 0:
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

        if i % INTERVAL == 0:
            print("source", source_newick)
            print("built", built_newick)

def my_minimize():
    def treefun(x):
        ref = ["BAAA", "ABCA", "ABAB", "BABA", "ABAA"]
        for leaf_i, leaf_j in itertools.combinations(range(len(ref)), 2)
            matches = [ref[leaf_i][c] for c in range(
                 
            )]

        n = int((math.sqrt(8 * len(x) + 1) - 1) / 2)

        branch_length = {leaf_idx:x[leaf_idx] for leaf_idx in range(n)}

        ultrametric_d = {}  # TODO: rename later
        root_dist = {leaf_idx:[] for leaf_idx in range(n)} # init to empty
        for idx, (leave_i, leave_j) in enumerate(itertools.combinations(range(n), 2)):
            ultrametric_d[leave_i, leave_j] = vector[idx + n]
            ultrametric_d[leave_j, leave_i] = vector[idx + n]

            root_dist[leave_i].append(vector[idx + n])
            root_dist[leave_j].append(vector[idx + n])

        root_dist = {leaf_idx:max(distances) for leaf_idx, distances in root_dist.items()}

        print(n)
        print(branch_length)
        print(root_dist)
        print(ultrametric_d)

    tree = ngesh.gen_tree(1.0, 0.5, max_time=1.0, labels="human", seed=42)
    vector = phylovector.tree2vector(tree)
    treefun(vector)

    return

    def myfun(x):
        v = [x[0]/x[1], x[0]/x[2], x[0]/x[3], x[0]/x[4]]
        return sum(v)

    from scipy.optimize import minimize, rosen, rosen_der
    x0 =[1.3, 0.7, 0.8, 1.9, 1.2]
    res = minimize(myfun, x0, method="BFGS", #jac=rosen_der,
    options={"gtol":1e-6, "disp":True})
    print(res.x)

    print("..", myfun(res.x))

if __name__ == "__main__":
    #stress_test1()
    my_minimize()