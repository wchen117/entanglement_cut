from anytree import Node, RenderTree, LevelOrderGroupIter

def get_all_keys(d):
    for key, value in d.items():
        yield key
        if isinstance(value, dict):
            yield from get_all_keys(value)

# return a list of paths 
def tranverse_tree(tree):

    print(RenderTree(tree))
    return

def build_psuedo_tree():
    tree = Node("root")
    a = Node("a", parent = tree)
    b = Node("b", parent = tree)
    c = Node("c", parent = tree)

    a1 = Node("a1", parent = a)
    a2 = Node("a2", parent = a)
    a3 = Node("a3", parent = a)

    b1 = Node("b1", parent = b)
    b2 = Node("b2", parent = b)
    b3 = Node("b3", parent = b)

    c1 = Node("b1", parent = c)
    c2 = Node("b2", parent = c)
    c3 = Node("b3", parent = c)

    return tree

def auto_leaf(node, leaf_num = 3):

    for idx in range(leaf_num):
        new_node = Node("{idx}".format(idx=idx), parent = node)

    return

def auto_build_tree(leaf_num =3, max_level = 2):

    tree = Node("root")
    
    level = 0;
   
    while (level < max_level):
        if level == 0:
            auto_leaf(tree)
            level = level + 1
        else:
            for leaf in tree.leaves:
                auto_leaf(leaf)
            level = level + 1

    return tree

def buildTree(treeNode):

    # build a tree with three nodes per level
    # and max 5 layers

    
    return treeNode


 
def main():
    #tree = build_psuedo_tree();
    tree = auto_build_tree();
    #RenderTree(tree)
    #print(tree)
    tranverse_tree(tree)
    

    return

if __name__ == '__main__':
    main()
