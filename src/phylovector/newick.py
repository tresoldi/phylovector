import re

# Parsing function following the elegant approach by `tricot` at
# https://stackoverflow.com/a/51375562 (licensed under CC BY-SA 4.0)
def parse(newick):
    tokens = re.finditer(
        r"([^:;,()\s]*)(?:\s*:\s*([\d.]+)\s*)?([,);])|(\S)", newick + ";"
    )

    def recurse(nextid=0, parentid=-1):  # one node
        thisid = nextid
        children = []

        name, length, delim, ch = next(tokens).groups(0)
        if ch == "(":
            while ch in "(,":
                node, ch, nextid = recurse(nextid + 1, thisid)
                children.append(node)
            name, length, delim, ch = next(tokens).groups(0)
            
        return (
            {
                "id": thisid,
                "name": name,
                "length": float(length) if length else None,
                "parentid": parentid,
                "children": children,
            },
            delim,
            nextid,
        )

    return recurse()[0]


def sorted_newick(tree):
    # Parse the tree in newick format
    p = parse(tree)

    # Sort in memory
    # TODO: also sort by number of children
    def mysort(node):
        if node["children"]:
            for child in node["children"]:
                mysort(child)
            node["children"] = sorted(node["children"], key=lambda n: n["name"])

    # Build representation
    def myout(node):
        elms = []
        if not node["children"]:
            return f"{node['name']}:{node['length']}"
        else:
            elms += [myout(n) for n in node["children"]]

        return "(%s)" % (",".join(elms))

    mysort(p)
    return myout(p) + ";"
