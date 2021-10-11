# Script for supporting development, later to be included as part of
# the paper

import phylovector


def main():
    vector = [3.0, 1.5, 3.0, 4.0, 3.0, 8.0, 8.0, 5.0, 3.0, 3.0, 8.0, 8.0, 8.0, 8.0, 5.0]
    leaves = ["A", "B", "C", "D", "E"]
    tree = phylovector.vector2tree(vector, leaves)

    print(tree)
    print(tree.write(format=1))

if __name__ == "__main__":
    main()
