
# class graph
class class_graph:
    number_vertices=0

    # create a molecule graph and init vertexs and edges 
    def __init__(self, nvertices, connect):
        # vertex  - a dict contains vertex laber
        # connect - connectivity information of vertices in upper triangle matrix form 
        #   connect is a dict
        #   for each atom, a list is used to store connectivity
        #   atomindex atom1 bond1 atom2 bond 2
        self.init_graph_method(nvertices,connect)

    # init graph 
    def init_graph_method(self, nvertices, connect):
        # vertex  - a dict contains vertex laber
        # connect - connectivity information of vertices in upper triangle matrix form 
        #   connect is a dict
        #   for each atom, a list is used to store connectivity
        #   atomindex atom1 bond1 atom2 bond 2
        #print connect
        self.atom_number=self.number_vertices
        self.number_vertices = nvertices
        self.edges   ={}
        self.connect ={}

        # edge = 0, not connect
        for i in range(1,self.number_vertices+1):
            self.connect[i]  = connect[i][:]           #原封不动的储存连接情况
            self.edges[i]    = {}
            for j in range(i, self.number_vertices+1):
                self.edges[i][j] = 0 

            leni = len(connect[i])
            if leni == 1:  # only atomindex, connectivity is done    不和其他原子相连
               continue
            else:
               for j in range(1,leni,2):
                   jvertex = int(connect[i][j])
                   self.edges[i][jvertex]=float(connect[i][j+1])    #存第i个和第j个原子之间的键

        # count number of direct connected vertices of each vertex  
        for i in range(1,self.number_vertices+1):
            for j in range(1, i):
                if self.edges[j][i] > 0:
                   self.edges[i][i] = self.edges[i][i] + 1
            for j in range(i+1, self.number_vertices+1):
                if self.edges[i][j] > 0:
                   self.edges[i][i] = self.edges[i][i] + 1         #[i][i]存连接的原子个数

        self.init_graph=self.connected_components_deepthfirstscan()
        self.scan_big_small_mole()

    #初始化 big  samll这两个列表，找到蛋白质和配体的原子序列
    def scan_big_small_mole(self):
        self.full_ord=self.connected_components_deepthfirstscan()
        print(self.full_ord)
        longest_key = max(self.full_ord, key=lambda x: len(self.full_ord[x]))
        self.big_mole_ord=self.full_ord[longest_key]
        self.smal_mole_ord=[]
        for i in self.init_graph:
            if i!=longest_key:
                self.smal_mole_ord.extend(self.init_graph[i])
        
    # cut the edge between ivertex and jvertex
    def cut_edge(self, ivertex, jvertex):
        i  = min(ivertex,jvertex)
        j  = max(ivertex,jvertex)
        ve = self.edges[i][j]
        if ve < 1 : # small value treate as zero
           #print("Vetex %4d and %4d do not connect ..." % (ivertex,jvertex))
            pass     
        else:
           # reset vertex connective information
           self.edges[i][i]=self.edges[i][i]-1
           self.edges[j][j]=self.edges[j][j]-1
           # cut edge
           self.edges[i][j]=0
           # reset connectivity
           #print self.connect[i] 
           for k in range(1,len(self.connect[i]),2):
               #print k
               if self.connect[i][k] == j:
                   # we label these elements to NULL and delete them later
                   self.connect[i][k]="NULL"
                   self.connect[i][k+1]="NULL"
                   break

           #print self.connect[i]
           for k in range(len(self.connect[i])-1,0,-1):
               if self.connect[i][k] == "NULL":
                   del self.connect[i][k]

           self.connect[i][0]=self.connect[i][0]-1 
           #print self.connect[i]

    def add_edge(self, ivertex, jvertex,bond_weight):
        i = min(ivertex, jvertex)
        j = max(ivertex, jvertex)
        ve = self.edges[i][j]
        if ve >= 1:  # if the edge already exists
            pass
        else:
            # set vertex connective information
            self.edges[i][i] = self.edges[i][i] + 1
            self.edges[j][j] = self.edges[j][j] + 1
            self.edges[i][j] = bond_weight
            # set connectivity
            self.connect[i].append(j)
            self.connect[i].append(bond_weight)  # assuming the weight of the edge is 1
            self.connect[i][0] = self.connect[i][0] + 1 


    #只记录大的，没有的话还是显示连接节点数
    def refine_connectivity(self):
        for i in range(1, self.number_vertices+1):
             if self.connect[i][0] == 0: # may connect with nodes whose indice small than current one
                 for j in range(1,i):
                     if self.edges[j][i] > 0:
                          self.connect[i][0] = self.connect[i][0] + 1 

    #只记录比当前大的节点，没有的话连接数为0
    def compress_connectivity(self):
        # for i in range(1, self.number_vertices+1):
        for i in range(1, len(self.connect)+1):
            if len(self.connect[i]) > 1:
                temp=[" "]
                for j in range(1, len(self.connect[i]), 2):
                    if i<self.connect[i][j]:  #双键的问题不知道解决没有
                        # for _ in range(int(self.connect[i][j+1])):
                        #     temp.append(self.connect[i][j])
                        #     temp.append(self.connect[i][j+1])
                        temp.append(self.connect[i][j])
                        temp.append(self.connect[i][j+1])
                temp[0]=int((len(temp)-1)/2)
                self.connect[i] = temp
        for i in range(1,len(self.connect)+1):
            if len(self.connect[i])==2:
                self.connect[i]=[0]

    #完全信息记录
    def extend_connectivity(self):
        # for i in range(1, self.number_vertices+1):
        for i in range(1, len(self.connect)+1):
             if len(self.connect[i]) == 1:
                 self.connect[i][0] = 0
             for j in range(1,i):
                 for k in range(1,len(self.connect[j]),2):
                      if i == self.connect[j][k]:
                          self.connect[i].append(j)
                          self.connect[i].append(self.connect[j][k+1])
                          self.connect[i][0]=self.connect[i][0]+1
                          break

    #找到与给定点直接相连的点    这个函数目前没有使用
    def find_connected_vertices(self, vertex, tags):
        connected_vertices = []
        #print self.connect[vertex],vertex
        if len(self.connect[vertex]) == 1: 
            tags[vertex]=1
            return connected_vertices

        # number of connect nodes
        nc = self.connect[vertex][0] 
        for ii in range(1,2*nc+1,2):
            inode = self.connect[vertex][ii]

            if tags[inode] == 1:
                continue
            tags[inode] = 1

            #print vertex,inode,self.edges[vertex][inode],self.connect[vertex][ii+1]
            connected_vertices.append(inode)
        #print "connected nodes",connected_vertices
        return connected_vertices

    # breadth first scan, from a root of a tree, visiting all son nodes
    def breadth_first_scan(self, vertices, tags):
        sonnodes   = vertices[:]
        leni = len(sonnodes)
        leno = leni - 1
        while (leni != leno):
            lenj = 0 
            for i in range(leno,len(sonnodes)):
                vertex = sonnodes[i]
                #print sonnodes,leni
                #print vertex, self.connect[vertex]
                if len(self.connect[vertex]) == 1:
                    #sonnodes = sonnodes + [inode]
                    continue
                else:
                    for ii in range(1,2*self.connect[vertex][0]+1,2):
                        inode = self.connect[vertex][ii]
                        if tags[inode] != 1 :
                            tags[inode] = 1
                            sonnodes = sonnodes + [inode]
                            lenj = lenj + 1
            leno = leni
            leni = lenj + leni
            #print leni
        #sonnodes = sonnodes + tnode #self.broad_first_scan(tnodes, tags)
        #print sonnodes
        #exit()
        return sonnodes
    #end broad first scan

    # deepth first scan: for a given node, find son nodes until no son 
    def deepth_first_scan(self, vertices, tags):
        sonnodes = vertices[:]
        #print "sonnodes",sonnodes
        for vertex in sonnodes:
            bond = self.edges[vertex][vertex]
            if len(self.connect[vertex]) != 1:
                #print self.connect[vertex]
                tnodes=[]
                for ii in range(1,2*self.connect[vertex][0]+1,2):
                    inode = self.connect[vertex][ii]
                    bond  = bond - 1 
                    if tags[inode] != 1 :
                        tags[inode] = 1
                        sonnodes = sonnodes + self.deepth_first_scan([inode],tags)
            # may connect with some nodes with smaller index
            if bond > 0:
               for j in range(1,vertex):
                   if tags[j] == 1: continue
                   if len(self.connect[j]) == 1: continue
                   for ii in range(1,2*self.connect[j][0]+1,2):
                       jnode = self.connect[j][ii]
                       if jnode == vertex:
                          tags[j] = 1
                          sonnodes = sonnodes + self.deepth_first_scan([j],tags)

        return sonnodes

    # find all connected sub-graphs (in graphic theory - named "connected component")
    # deepth first scan is used in searching 
    def connected_components_deepthfirstscan(self):
        subgraph = {}
        # breadth-first algorithm is used to find a connected component from a given vetice
        node_found = {}
        for inode in range(1, self.number_vertices+1):
            node_found[inode] = 0 

        nsub = 0
        for inode in range(1, self.number_vertices+1):
            if node_found[inode] == 1:
                continue

            #print("inode search",inode)     #先注释了，没发现有什么用

            node_found[inode] = 1
            if self.edges[inode][inode] < 1 :  # no connected vertex
                nsub = nsub+1
                subgraph[nsub]    = [inode]
            else:
                nsub = nsub+1
                subgraph[nsub] = []
                # if number of direct connected node = 1 only connected with a node, skip 
                # if number of direct connected node > 1, we will look for subgraph
                subgraph[nsub] = self.deepth_first_scan([inode], node_found)
            #print sorted(subgraph[nsub])
            #exit()
        
        for i in range(1,len(subgraph)+1):
            if len(subgraph[i])==3:
                j=subgraph[i][0]
                if self.residue_type[j]=="HOH":
                    del subgraph[i]  
        n = len(subgraph)     
        for i in range(1,n+1):
            subgraph[i].sort()
        return subgraph
    #end connected_components_deepthfirstscan

    # find all connected sub-graphs (in graphic theory - named "connected component")
    # broad first scan is used in searching 
    def connected_components_breadthfirstscan(self):
        subgraph = {}
        # breadth-first algorithm is used to find a connected component from a given vetice
        node_found = {}
        for inode in range(1, self.number_vertices+1):
            node_found[inode] = 0 

        nsub = 0
        for inode in range(1, self.number_vertices+1):
        # for inode in range(3,4): #(1, self.number_vertices+1):
            #print "inode search",inode
            if node_found[inode] == 1:
                continue

            node_found[inode] = 1
            if self.edges[inode][inode] < 1 :  # no connected vertex
                nsub = nsub+1
                subgraph[nsub]    = [inode]
            else:
                nsub = nsub+1
                subgraph[nsub] = []
                # if number of direct connected node = 1 only connected with a node, skip 
                # if number of direct connected node > 1, we will look for subgraph
                subgraph[nsub] = self.breadth_first_scan([inode], node_found)
            #print sorted(subgraph[nsub])
            #exit()

        #去除水分子团影响
        # for i in range(1,len(subgraph)+1):
        #     if len(subgraph[i])==3:
        #         j=subgraph[i][0]
        #         if self.AacidName[j]=="HOH":
        #             del subgraph[i]

        return subgraph
    #end connected_components_breadthfirstscan
    
    #创造一个函数，输入我们想寻找的键，可以返回这两个原子的序号
    #固定寻找肽键，返回字典{1:[5,6],2,[7,8]}  或者一个列表，也就是没有键只有值  （我们为什么需要用键值对？）
    #半成品
    def search_peptide_bond(self):
        for i in range(1,self.number_vertices):
            if self.connect[i][0]=="3":
                if "2.0" in self.connect[i]:
                    for j in range(1,len(self.connect[i]),2):
                        aa=int(self.connect[i][j])
                        if self.connect[aa][0]=="3":
                            pass
    
    def is_connected(self,node1,node2):
        nmax = max(node1, node2)
        nmin = min(node1, node2)
        if self.edges[nmin][nmax] >= 1:
            return True