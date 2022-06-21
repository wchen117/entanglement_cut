from IndexNode import IndexNode
from anytree import RenderTree


def transvers_path(tree):

    """assume the give tree is a root note"""

    big_list = []
    if tree.is_root:
        for leaf in tree.leaves:
            big_list.append(leaf.getIndexPath())


    return big_list


def build_psuedo_tree():
    tree = IndexNode("root")
    a = IndexNode("a", parent = tree)
    b = IndexNode("b", parent = tree)
    c = IndexNode("c", parent = tree)

    a1 = IndexNode("a1", parent = a)
    a2 = IndexNode("a2", parent = a)
    a3 = IndexNode("a3", parent = a)

    b1 = IndexNode("b1", parent = b)
    b2 = IndexNode("b2", parent = b)
    b3 = IndexNode("b3", parent = b)

    c1 = IndexNode("c1", parent = c)
    c2 = IndexNode("c2", parent = c)
    c3 = IndexNode("c3", parent = c)

    return tree

def auto_leaf(node, leaf_num = 3):

    for idx in range(leaf_num):
        new_node = IndexNode("{idx}".format(idx=idx), parent = node)

    return

def auto_build_tree(leaf_num =3, max_level = 3):

    tree = IndexNode("root")
    
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
 
def main():
    #tree = build_psuedo_tree();
    tree = auto_build_tree()
    print(RenderTree(tree))
    #print(tree)
    all_paths = transvers_path(tree)
    print(all_paths)
    

    return

if __name__ == '__main__':
    main()
