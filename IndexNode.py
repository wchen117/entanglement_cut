from anytree import Node, RenderTree

class IndexNode(Node):
    """
       A node implementation of cut indices, based on anytree module 
       (https://anytree.readthedocs.io/en/latest/)
       This implementation of index node contains two additional optional quantities: 
       Index (resid) and Ent (num of entanglement/partial linking number)
    """
    def __init__(self, name, Ent=-1, Index=-1,  parent=None, children=None, **kwargs):
        super().__init__(name, parent, children, **kwargs)
        # set default value to -1 for error handling
        self.Ent = Ent
        self.Index = Index

    def visualize(self):
        """
        print out the tree if the current node is a root node
        """
        if self.is_root:
            print(RenderTree(self))


    def getIndexPath(self):
        """
        for a paticular leaf node, trace to its root and print out 
        the indices (from root to leaf) on its path
        """
       
        index_list = []
        for pa in self.path:
            # omit the root node which by default has an index value of -1
            if (pa.Index == -1 and pa.is_root):
                continue
            index_list.append(pa.Index)
        return index_list



