# Script for supporting development, later to be included as part of
# the paper

import ngesh

import phylovector


def main():
    vector = [3.0, 1.5, 3.0, 4.0, 3.0, 8.0, 8.0, 5.0, 3.0, 3.0, 8.0, 8.0, 8.0, 8.0, 5.0]
    leaves = ["A", "B", "C", "D", "E"]
    tree = phylovector.vector2tree(vector, leaves)

    # print(tree)
    newick = tree.write(format=1)
    print(newick)

    new_vector = phylovector.tree2vector(tree)
    print(new_vector)

    # random
    # t0 = ngesh.gen_tree(1.0, 0.5, max_time=2.0, labels="human")
    t1 = ngesh.gen_tree(1.0, 0.5, max_time=2.0, labels="human", seed=123)
    print(t1)
    v = phylovector.tree2vector(t1)
    print("v", v)

    leaves = sorted([leaf.name for leaf in t1.iter_leaves()])
    t2 = phylovector.vector2tree(v, leaves)

    rf, rf_max, x1, x2, x3, x4, x5 = t1.robinson_foulds(t2)
    print("RF distance is %s over a total of %s" % (rf, rf_max))

    # rf, rf_max, x1, x2, x3, x4, x5 = t1.robinson_foulds(t0)
    # print ("RF distance is %s over a total of %s" %(rf, rf_max))

    print("t1", phylovector.sorted_newick(t1.write(format=1)))
    print("t2", phylovector.sorted_newick(t2.write(format=1)))

    # print ("Partitions in tree2 that were not found in tree1:", edges_t1 - edges_t2)
    # print ("Partitions in tree1 that were not found in tree2:", edges_t2 - edges_t1)


def main2():
    n = "(((Wia:1.13157,(Fisefa:1.03597,(Hebo:0.239366,Aku:0.239366):0.796607):0.0955923):1.30671,((Nenesnu:1.46834,Cebirma:1.46834):0.370696,Nobavta:1.83904):0.599238):0.263869,Bavuli:2.70214);"

    sn = phylovector.sorted_newick(n)
    print(sn)


if __name__ == "__main__":
    main()
#    main2()
