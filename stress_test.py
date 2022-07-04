import math
import itertools
import ngesh
from numpy.lib.function_base import diff
import phylovector

import ete3


def stress_test1(iters=1000):
    for i in range(iters):
        if i % (iters / 10) == 0:
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


def stress_test2():
    for i in range(10):
        source_tree = ngesh.gen_tree(1.0, 0.5, max_time=2.0, labels="human")
        source_newick = source_tree.write(format=1)
        source_newick = phylovector.sorted_newick(source_newick)

        print(source_newick)


def my_minimize():
    def tree_likelihood(x, alm, age):
        alm_len = len(alm[0])

        # Collect ultrametric info -- note that it is not correcting by branch length
        n = int((math.sqrt(8 * len(x) + 1) - 1) / 2)
        ultrametric_d = {}
        matches = {}
        for idx, (leaf_i, leaf_j) in enumerate(
            itertools.combinations(range(len(alm)), 2)
        ):
            ultrametric_d[leaf_i, leaf_j] = x[idx + n]
            ultrametric_d[leaf_j, leaf_i] = x[idx + n]

            matches[leaf_i, leaf_j] = sum(
                [ci == cj for ci, cj in zip(alm[leaf_i], alm[leaf_j])]
            )

        # Compute proportion of distance to length
        max_dist = max(ultrametric_d.values())
        first_ancestor_dist = {idx: max_dist for idx in range(len(alm))}
        v = []
        for leaf_i, leaf_j in itertools.combinations(range(len(alm)), 2):
            m = matches[leaf_i, leaf_j]
            d = ultrametric_d[leaf_i, leaf_j]

            m_prop = m / alm_len
            d_prop = 1 - (d / max_dist)

            v.append(abs(m_prop - d_prop))

            if d < first_ancestor_dist[leaf_i]:
                first_ancestor_dist[leaf_i] = d
            if d < first_ancestor_dist[leaf_j]:
                first_ancestor_dist[leaf_j] = d

        # high penalty for going back in time (could be a constraint in some methods)
        for idx in range(len(alm)):
            if max_dist - first_ancestor_dist[idx] < x[idx]:
                return 9999

        return sum(v)  # ** abs(age-max_dist)

    ref = ["BAAA", "ABCA", "ABAB", "BABA", "ABAA"]

    tree = ngesh.gen_tree(1.0, 0.8, max_time=1.0, labels="human", seed=42)
    # tree.show()
    vector = phylovector.tree2vector(tree)
    new_tree = phylovector.vector2tree(vector, ["Agu", "E", "Meopo", "Sedu", "Tarami"])
    print("...", vector)
    #tree.show()
    #return
    v = tree_likelihood(vector, ref, 2.5)

    from scipy.optimize import minimize, rosen, rosen_der

    bounds = (
        (0.0, 0.0),  # A
        (0.5, 1.0),  # B
        (0.0, 0.0),  # C
        (0.6, 0.8),  # D
        (0.0, 0.0),  # E
        (0.0, 5.0),  # AB
        (0.0, 5.0),  # AC
        (0.0, 5.0),  # AD
        (0.0, 5.0),  # AE
        (0.0, 5.0),  # BC
        (0.0, 5.0),  # BD
        (0.0, 5.0),  # BE
        (0.0, 5.0),  # CD
        (0.0, 5.0),  # CE
        (0.0, 5.0),  # DE
    )

    res = minimize(
        tree_likelihood,
        vector,
        args=(ref, 2.5),
        bounds=bounds,
        method="Nelder-Mead",
        options={"gtol": 1e-6, "disp": True},
    )

    inf_v = res.x
    inf_t = phylovector.vector2tree(inf_v, ["A", "B", "C", "D", "E"])
    inf_t.show()

    ortree = phylovector.sorted_newick(tree.write(format=1))
    intree = phylovector.sorted_newick(inf_t.write(format=1))

    print(ortree)
    print(intree)
    print(inf_v)


def tiago():
    dm = {
        ("bonobo", "chimpanzee"): 0.09696787751425995,
        ("bonobo", "gorilla"): 0.24977484238967274,
        ("bonobo", "homo"): 0.22095466826778742,
        ("bonobo", "orangutan"): 0.3161212848994296,
        ("chimpanzee", "gorilla"): 0.25788051636145304,
        ("chimpanzee", "homo"): 0.22695887120984692,
        ("chimpanzee", "orangutan"): 0.32182527769438607,
        ("gorilla", "homo"): 0.2668868207745422,
        ("gorilla", "orangutan"): 0.322425697988592,
        ("homo", "orangutan"): 0.31732212548784144,
    }

    bounds = (
        (0.0, 0.1),  # A
        (0.5, 1.0),  # B
        (0.0, 0.1),  # C
        (0.6, 0.8),  # D
        (0.0, 0.1),  # E
        (0.0, 5.0),  # AB
        (0.0, 5.0),  # AC
        (0.0, 5.0),  # AD
        (0.0, 5.0),  # AE
        (0.0, 5.0),  # BC
        (0.0, 5.0),  # BD
        (0.0, 5.0),  # BE
        (0.0, 5.0),  # CD
        (0.0, 5.0),  # CE
        (0.0, 5.0),  # DE
    )

    

    # Define a starting vector
    leaves = ["bonobo", "chimpanzee", "gorilla", "homo", "orangutan"]
    vector = phylovector.start_vector(5, leaves)

    # Infer
    from scipy.optimize import minimize, basinhopping, dual_annealing
    if False:
        res = dual_annealing(phylovector.primate, bounds)
    if True:
        #minimizer_kwargs = {"method": "CG"}
        res = basinhopping(phylovector.primate,
        vector)
        #minimizer_kwargs=minimizer_kwargs)
    if False:
        res = minimize(
            phylovector.method1,
            vector,
            args=(dm,),
            #bounds=bounds,
            method="BFGS",
            options={"gtol": 1e-8, "disp": False},
        )
    inf_v = res.x
    print(inf_v)
    inf_t = phylovector.vector2tree(inf_v, leaves)
    inf_t.show()

def tiago2():
    from skopt import gp_minimize

    # Build bounds
    bounds = [(0.0, 2.0) for n in range(15)]
    #bounds = [(0.0, 0.01) for n in range(10)] + [(0.0, 2.0) for n in range(10)]

    res = gp_minimize(phylovector.primate,                  # the function to minimize
                    bounds,      # the bounds on each dimension of x
                    acq_func="EI",      # the acquisition function
                #    n_calls=15,         # the number of evaluations of f
                    n_random_starts=5,  # the number of random initialization points
                    noise=0.1**2,       # the noise level (optional)
                    random_state=1234)   # the random seed

    leaves = ["bonobo", "chimpanzee", "gorilla", "homo", "orangutan"]
    print(res.x)
    inf_t = phylovector.vector2tree(res.x, leaves)
    inf_t.show()

if __name__ == "__main__":
    # stress_test1()
    # stress_test2()
    # my_minimize()
    tiago()
    #tiago2()
