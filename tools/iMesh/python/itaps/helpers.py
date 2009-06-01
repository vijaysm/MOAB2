class AdjacencyList:
    def __init__(self, adj, offsets):
        self.adj = adj
        self.offsets = offsets

    def __getitem__(self, key):
        if isinstance(key, tuple):
            if self.length(key[0]) <= key[1]:
                raise IndexError
            return self.adj[ self.offsets[key[0]] + key[1] ]
        else:
            return self.adj[ self.offsets[key]:self.offsets[key+1] ]

    def length(self, i = None):
        if i == None:
            return len(self.offsets)-1
        else:
            return self.offsets[i+1] - self.offsets[i]

    def __len__(self):
        return self.length()


class IndexedAdjacencyList:
    def __init__(self, entities, adj, indices, offsets):
        self.entities = entities
        self.adj = adj
        self.indices = indices
        self.offsets = offsets

    def __getitem__(self, key):
        return self.adj[ self.index(key) ]

    def index(self, key):
        if isinstance(key, tuple):
            if self.length(key[0]) <= key[1]:
                raise IndexError
            return self.indices[ self.offsets[key[0]] + key[1] ]
        else:
            return self.indices[ self.offsets[key]:self.offsets[key+1] ]

    def length(self, i = None):
        if i == None:
            return len(self.offsets)-1
        else:
            return self.offsets[i+1] - self.offsets[i]

    def __len__(self):
        return self.length()
