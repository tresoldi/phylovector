import math
import itertools
import ngesh
import phylovector


def stress_test1(iters=1000):
    for i in range(iters):
        if i % (iters/10) == 0:
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

    #    if i % INTERVAL == 0:
    #        print("source", source_newick)
    #        print("built", built_newick)


def my_minimize():
    def treefun(x):
        ref = ["BAAA", "ABCA", "ABAB", "BABA", "ABAA"]
        alm_len = len(ref[0])

        # Collect ultrametric info -- note that it is not correcting by branch length
        n = int((math.sqrt(8 * len(x) + 1) - 1) / 2)
        ultrametric_d = {}
        matches = {}
        for idx, (leaf_i, leaf_j) in enumerate(
            itertools.combinations(range(len(ref)), 2)
        ):
            ultrametric_d[leaf_i, leaf_j] = x[idx + n]
            ultrametric_d[leaf_j, leaf_i] = x[idx + n]

            matches[leaf_i, leaf_j] = len(
                [ci for ci, cj in zip(ref[leaf_i], ref[leaf_j]) if ci == cj]
            )

        # Compute proportion of distance to length
        max_dist = max(ultrametric_d.values())
        max_diff = max(matches.values())
        v = []
        for leaf_i, leaf_j in itertools.combinations(range(len(ref)), 2):
            m = matches[leaf_i, leaf_j]
            d = ultrametric_d[leaf_i, leaf_j]

            m_prop = m / alm_len
            d_prop = 1 - (d / max_dist)

            v.append(abs(m_prop - d_prop))

        return sum(v)

    tree = ngesh.gen_tree(1.0, 0.5, max_time=1.0, labels="human", seed=42)
    tree.show()
    vector = phylovector.tree2vector(tree)
    v = treefun(vector)

    from scipy.optimize import minimize, rosen, rosen_der

    res = minimize(
        treefun,
        vector,
        method="BFGS",  # jac=rosen_der,
        options={"gtol": 1e-6, "disp": True},
    )

    inf_v = res.x
    inf_t = phylovector.vector2tree(inf_v, ["Agu", "E", "Meopo", "Sedu", "Tarami"])
    inf_t.show()

    print(tree.write(format=1))
    print(inf_t.write(format=1))


if __name__ == "__main__":
    stress_test1()
    #my_minimize()
