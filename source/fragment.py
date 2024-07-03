import math
import molegraph
import os
# coord={1: ['S', '31.113', '7.238', '-5.638', [2, 41, 2.0, 39, 1.0], 'SG', 'A'], 
#        2: ['S', '30.585', '6.536', '-4.764', [1, 40, 2.0], 'SG', 'B'],}
class class_molecule(molegraph.class_graph):
    natom=0
    def __init__(self,name,coord):
        self.name  = name
        self.natom = len(coord)
        self.label  = {}
        
        self.label2 = {}
        self.chain  = {}
        self.residue_type={}
        
        self.charge = {}
        self.mass   = {}
        self.coord  = {}
        self.moleconnect={}

        self.dangling_bond={}
        self.residue_list={}
        for i in coord:
            self.label[i] =coord[i][0]
            self.coord[i] = {}
            self.coord[i][1]=float(coord[i][1])
            self.coord[i][2]=float(coord[i][2])
            self.coord[i][3]=float(coord[i][3])
            self.moleconnect[i]=coord[i][4]
            self.dangling_bond[i]=0
            if len(coord[i])>5:  #也就是标砖pdb文件
                self.label2[i]=coord[i][5]
                self.chain[i] =coord[i][6]
                self.residue_type[i]=coord[i][7]
                
            
        self.atomdist = self.atomdistance()
        molegraph.class_graph.__init__(self, self.natom, self.moleconnect)

    #计算所有原子的距离
    def atomdistance(self):
        #print self.coord
        dist={}
        for i in range(1,self.natom+1):
           dist[i]={i: 0.0}
           for j in range(1,i):
               v = (self.coord[i][1]-self.coord[j][1])**2 \
                 + (self.coord[i][2]-self.coord[j][2])**2 \
                 + (self.coord[i][3]-self.coord[j][3])**2
               dist[i][j]=math.sqrt(v)

        return dist

    def molecomp_distance(self, molecomp):
        nfrag = len(molecomp)
        dist={}
        r = {}
        for ifrag in range(1, nfrag+1):
            ni = len(molecomp[ifrag])
            for i in range(1,ni+1):
                ia = molecomp[ifrag][i-1]
                r[ia] = [0]
                r[ia].append(float(self.coord[ia][1]))
                r[ia].append(float(self.coord[ia][2]))
                r[ia].append(float(self.coord[ia][3]))

        for ifrag in range(1, nfrag+1):
            ni = len(molecomp[ifrag])
            dist[ifrag] = {}
            for jfrag in range(1,ifrag+1):
                nj = len(molecomp[jfrag])
                
                vmin = 1.0E30
                for i in range(1,ni+1):
                    #print ni,molecomp[ifrag][i-1]
                    ia = molecomp[ifrag][i-1]
                    #print coord[ia]
                    x1 = r[ia][1]
                    y1 = r[ia][2]
                    z1 = r[ia][3]
                    for j in range(1,nj+1):
                        #print nj,j,molecomp[jfrag]
                        ib = molecomp[jfrag][j-1]
                        x2 = r[ib][1]
                        y2 = r[ib][2]
                        z2 = r[ib][3]
                        vd = (x2-x1)**2+(y2-y1)**2+(z2-z1)**2 
                        vd = math.sqrt(vd)
                        if vd < vmin:
                            vmin = vd
                dist[ifrag][jfrag]=vmin
                # put into array

        # print("\nMolecomp distances")
        # print("^"*100)
        # from bdfutil import print_ltmat
        # print_ltmat("Molecomp",nfrag,dist)
        # print("^"*100)
        # exit()
        return dist

    #end of molecomp_distance

    def molecomp_distance_dsmetric(self, molecomp, bdftask, suffix):
        '''
        Similar to molecomp_distance, but use density matrix elements and overlap matrix
        elements as metric, instead of using real-space distances.
        Note that this function overwrites self.atomdist.
        Valid suffixes are '.ioicoulpot' (for the metric that involves both D and S) and
        '.smetric' (for the metric that only involves S)
        '''
        nfrag = len(molecomp)
        dist={}

        # Read metric matrix from file
        from fragutil import portable_float
        with open(bdftask+suffix+'.dsmetric','r') as fn:
            for iatom in range(1, self.natom+1):
                self.atomdist[iatom][iatom] = 0.0
                for jatom in range(iatom+1, self.natom+1):
                    line = fn.readline().strip()
                    self.atomdist[jatom][iatom] = portable_float(line)

        for ifrag in range(1, nfrag+1):
            ni = len(molecomp[ifrag])
            dist[ifrag] = {}
            for jfrag in range(1,ifrag+1):
                nj = len(molecomp[jfrag])
                
                vmin = 1.0E30
                for i in range(1,ni+1):
                    for j in range(1,nj+1):
                        iatom = molecomp[ifrag][i-1]
                        jatom = molecomp[jfrag][j-1]
                        if iatom > jatom:
                            vd = self.atomdist[iatom][jatom]
                        else:
                            vd = self.atomdist[jatom][iatom]
                        if vd < vmin:
                            vmin = vd
                dist[ifrag][jfrag]=vmin

        return dist


    def split_molecule_to_molecomps(self, bondlist):  #删除了linkatom  
        leni=len(self.connect)
        if leni == 0: 
            print ("No connectivity information. Cannot split molecule ...")
            exit()
        print("\nSplitting molecule into molecomps ...\n")
        # Cut bond from input bond list likes
        #   [["C","C",1.0],["C","C",2.0]]
        if bondlist[0] != "None":
            for bond in bondlist:
                # check if the bond is valid
                re = check_bond(bond)
                if not re:
                    print ("Input bond ", bond, " error. Cannot cut bond!")
                    exit()  
            for i in range(1,self.number_vertices+1):
                iatom = self.label[i]
                for j in range(i+1,self.number_vertices+1):
                    jatom = self.label[j]
                    if self.edges[i][j] == 0.0:
                        continue
                    for bond in bondlist:
                        if (bond[0] != iatom or bond[1] != jatom) and \
                           (bond[1] != iatom or bond[0] != jatom):
                            continue
                        if self.edges[i][j] == bond[2]:
                            # print("Cutting bond between %d and %d" % (i,j))
                            self.cut_edge(i,j)
                        break
        self.refine_connectivity()
        molecomp = self.connected_components_deepthfirstscan()

        # wzk20201116: remove molecomps that consist of external linkatoms
        # j = 1
        # n = len(molecomp)
        # for i in range(1,n+1):
        #     if len(molecomp[i])==1 and molecomp[i][0] in linkatom:
        #         continue
        #     molecomp[j] = molecomp[i]
        #     j += 1
        # for i in range(j,n+1):
        #     del molecomp[i]

        n = len(molecomp)
        for i in range(1,n+1):
            molecomp[i].sort()

        return molecomp

#{1: ['molecomp1', '8.67041', '-5.69968', '-1.51405', [1, 2, 1.0]], 2: ['molecomp2', '10.93529', '-3.20107', '-5.64164', [1, 3, 1.0]], 3: ['molecomp3', '13.69965', '1.30594', '-3.68776', [1, 4, 1.0]], 4: ['molecomp4', '7.83971', '1.13818', '-7.65559', [1]]}
#感觉对片段连接处信息的标记有误，上述是一个线性的四肽的例子
def molecomp_centercoord_connectivity(molecule, molecomp):
    nfrag = len(molecomp)
    coord={}
    print(molecomp)
    print("nfrag=",nfrag)
 
    # calculate center
    for ifrag in range(1,nfrag+1):
        name  = "molecomp%d" % ifrag

        #print ifrag,name
        coord[ifrag] = [name]
        #print molecomp[ifrag]
        natom = len(molecomp[ifrag])
        if natom == 1:
            iatom = molecomp[ifrag][0]
            str1  = "%10.5f  %10.5f  %10.5f" % (molecule.coord[iatom][1],molecule.coord[iatom][2],molecule.coord[iatom][3])
            str1  = str1.split()
            coord[ifrag] = coord[ifrag] + str1
        else:
            matom = 0
            x=0.0
            y=0.0
            z=0.0
            for iatom in molecomp[ifrag]:
                # neglect H atom
                #                 str0 = molecule.label[iatom].upper()
                #if str0 == "H":
                #    continue
                matom = matom + 1
                x=x+molecule.coord[iatom][1]
                y=y+molecule.coord[iatom][2]
                z=z+molecule.coord[iatom][3]            
            str1  = "%10.5f  %10.5f  %10.5f" % (x/matom,y/matom,z/matom)
            str1  = str1.split()
            coord[ifrag] = coord[ifrag] + str1

    # generate connectivity information
    edges   = {}
    connect = {}
    for ifrag in range(1,nfrag+1):
        connect[ifrag] = [0]
 
        edges[ifrag]={}
        edges[ifrag][ifrag]=0
        for jfrag in range(ifrag+1,nfrag+1):
            edges[ifrag][jfrag] = 0 
            for iatom in molecomp[ifrag]:
                for jatom in molecomp[jfrag]:
                    if iatom == jatom:
                        print("Error, molecomp %d and %d overlap" % (ifrag,jfrag))
                        exit()
                    i = min(iatom,jatom)
                    j = max(iatom,jatom)
                    #if ifrag == 6 and jfrag == 7:
                    #    print iatom,jatom
                    if molecule.edges[i][j] > 0:
                        edges[ifrag][jfrag] = edges[ifrag][jfrag] + molecule.edges[i][j]
                        edges[ifrag][ifrag] = edges[ifrag][ifrag] + 1 
            if edges[ifrag][jfrag] > 0:
                connect[ifrag][0]=connect[ifrag][0]+1
                connect[ifrag].append(jfrag)
                connect[ifrag].append(edges[ifrag][jfrag])

    # refine connectivity information
    for ifrag in range(1,nfrag+1):
        if connect[ifrag][0] > 0:
            for j in range(1,len(connect[ifrag]),2):
                 jfrag = connect[ifrag][j]
                 if len(connect[jfrag]) == 1:
                     connect[jfrag][0] = connect[jfrag][0] + 1
        coord[ifrag].append(connect[ifrag])
    return coord


#自动猜测初始化分子片段的聚类中心，返回聚类中心列表
#按照   片段包含原子数，邻居距离   排序，再选
def auto_guess_init_clusters_new(molecomp_graph, molecomp, mindist, radcent):
    connected_graph = molecomp_graph.connected_components_deepthfirstscan()
    #complist = []
    #for i in range(1, len(connected_graph)+1):
    #    complist = complist + connected_graph[i]

    clusters = []
    #remaining_molecomp = range(1, len(molecomp)+1)
    remaining_molecomp = list(range(1, len(molecomp)+1))

    # 1. If a component contains more than maxatom atoms, select it as a centroid
    maxatom = 5
    for i in remaining_molecomp:
        if i == 0:
            continue
        if len(molecomp[i]) > maxatom:
            clusters = clusters + [[i]]
            remaining_molecomp[i-1] = 0 # this is more efficient than remove()

    # 2. Sort the molecomps in ascending order according to the number of neighbors. In case
    # of degeneracy, sort in descending order of the distance to the furthest neighbor
    nneighbor = [0]*len(remaining_molecomp)
    maxdist = [0.0]*len(remaining_molecomp)
    for i,ic in enumerate(remaining_molecomp):
        if ic == 0: continue
        for j,jc in enumerate(remaining_molecomp[i+1:]):
            if jc == 0: continue
            if ic > jc:
                dist = mindist[ic][jc]
            else:
                dist = mindist[jc][ic]
            if dist < radcent:
                nneighbor[i] += len(molecomp[j+1])
                nneighbor[j] += len(molecomp[i+1])
                maxdist[i] = max(maxdist[i],dist)
                maxdist[j] = max(maxdist[j],dist)
    tmp = [(nneighbor[i],-maxdist[i],remaining_molecomp[i]) for i in range(len(remaining_molecomp))]
    tmp.sort()
    remaining_molecomp = [c for n,m,c in tmp]

    # 3. If a component is attached to at least nconnect components, select the
    #    former as a centroid and remove the latter components from consideration
    #    Herein nconnect is taken as 4,3,2,1,0 in that order
    for nconnect in range(4,-1,-1):
        for i,ic in enumerate(remaining_molecomp):
            if ic == 0:
                continue
            leni = len(molecomp_graph.connect[ic])
            remaining_connect = []
            for j in range(1,leni,2):
                try:
                    jc = remaining_molecomp.index(molecomp_graph.connect[ic][j])
                    remaining_connect.append(jc)
                except ValueError: # not found
                    pass
            if len(remaining_connect) >= nconnect:
                clusters = clusters + [[ic]]
                remaining_molecomp[i] = 0
                for j in remaining_connect:
                    remaining_molecomp[j] = 0

    # 4. If two centroids are too close to each other, remove the one with more neighbor atoms
    for i,ic in enumerate(clusters):
        if ic == [0]:
            continue
        for j,jc in enumerate(clusters[i+1:]):
            if jc == [0]:
                continue
            if ic[0] > jc[0]:
                dist = mindist[ic[0]][jc[0]]
            else:
                dist = mindist[jc[0]][ic[0]]
            if dist < radcent:
                clusters[j] = [0]
    # remove all zeroes
    #clusters = filter(lambda x: x!=[0], clusters)
    clusters = [x for x in clusters if x!=[0]]

    print("Auto guess cluster centroids:")
    print(clusters)
    return clusters

#对分子进行分割，使用 K-means 聚类算法来创建非重叠的分子片段。
# aggregate 为 True 并且 centroids 有多个，那么会使用 K-means 算法将分子组件聚集到非重叠的中心分子片段中。
def build_none_overlap_molefrag(molecule,molecomp_graph,molecomp,aggregate,nclusters,centroids=None):
     from kmeans import kmeans,repeatedKMeans

     center_molefrag = {}
     if not aggregate:  # do not aggregate, use each molecomp as a center fragment
         number_molecomp = molecomp_graph.natom
         for ifrag in range(1,number_molecomp+1):
             center_molefrag[ifrag]=[]
             center_molefrag[ifrag].append(ifrag)
     elif len(centroids) == 1:
         # wzk20201116: special case - there is only one fragment
         # In this case kmeans will fail, so this case has to be treated separately
         number_molecomp = molecomp_graph.natom
         center_molefrag[1] = range(1,number_molecomp+1)
     else: # aggregate molecule components to none overlap center molecule fragments
         # aggregate cluser by using Kmeans cluster alghrithm
         # init input vectors
         coord={}
         for i in range(1,molecomp_graph.natom+1):
              str1 = "molecomp" + str(i)
              coord[i]=[str1]
              coord[i].append(molecomp_graph.coord[i][1])
              coord[i].append(molecomp_graph.coord[i][2])
              coord[i].append(molecomp_graph.coord[i][3])
         molefrag_old={}
         molefrag_old[1]=["None"]
         maxiter = 10 # maxium number of macro iteration
         niter   = 0
         while(molefrag_old != center_molefrag and niter < maxiter):
             niter = niter+1
             print(" Macro iteration: ", niter)

             molefrag_old=center_molefrag.copy()
             # K-means clustering
             if centroids == None:
                 # user do not define center, perform K-means several times and pick up clusters with smallest withness
                 tmp_molefrag = repeatedKMeans(coord, nclusters, 1) #2*nclusters)
             else:
                 # user has define center
                 tmp_molefrag = kmeans(coord, nclusters, False, centroids)

             #print "not refined",tmp_molefrag["clusters"]
             refined=refine_molefrag_by_connectivity(molecomp_graph, nclusters, tmp_molefrag)
             #print "refined",tmp_molefrag["clusters"]
 
             from kmeans import computeWithinss
             withinss = computeWithinss(coord, tmp_molefrag["clusters"], tmp_molefrag["centroids"])
             #print "refined 0",tmp_molefrag["clusters"]
 
             #withinss = tmp_molefrag["withinss"]       
             print("\n Final cluster withinss %12.4f" % withinss)
             centroids = tmp_molefrag["clusters"]
             for i in range(0, nclusters):
                 center_molefrag[i+1] = tmp_molefrag["clusters"][i]
     #print center_molefrag
     #exit()
     # reture result
     return center_molefrag

#虑分子组件之间的连接性来优化 K-means 聚类的结果
def refine_molefrag_by_connectivity(molecomp_graph, nclusters, tmp_molefrag):
    #print tmp_molefrag["clusters"]
    #print tmp_molefrag["centroids"]
    #print nclusters

    # check connection first
    for icluster in range(0,nclusters):
        molecomplist = tmp_molefrag["clusters"][icluster] # rename, not copy
        disconnect = []
        for icomp in molecomplist:
            connected=False
            for jcomp in molecomplist:
                if icomp == jcomp:
                    continue
                ii = min(icomp, jcomp)
                jj = max(icomp, jcomp)
                if molecomp_graph.edges[ii][jj] > 0.0001:
                    connected = True
                    break
            if not connected:
                disconnect.append(icomp)
        #print disconnect
        #print len(molecomplist)
        if len(disconnect) < 1 or len(molecomplist) < 2: 
            continue

        # move disconnected molecomp in a cluster to nearest connected cluster
        for icomp in disconnect:
            xi = molecomp_graph.coord[icomp][1]
            yi = molecomp_graph.coord[icomp][2]
            zi = molecomp_graph.coord[icomp][3]

            # pick up a nearest cluster 
            leni = len(molecomp_graph.connect[icomp])
            mindis = 1.E30
            iselect= icluster 
            for j in range(1,leni,2):
                jcomp = molecomp_graph.connect[icomp][j]
                for jcluster in range(0, nclusters):
                    if jcluster == icluster: continue
                    if jcomp in tmp_molefrag["clusters"][icluster]:
                        dist = (xi-tmp_molefrag["centroids"][jcluster][1])**2 \
                             + (yi-tmp_molefrag["centroids"][jcluster][2])**2 \
                             + (zi-tmp_molefrag["centroids"][jcluster][3])**2
                        if dist < mindis:
                             mindis  = dist
                             iselect = jcluster
            if iselect == icluster: continue

            #print icomp,iselect 
            #tmp_molefrag["clusters"][icluster].remove(icomp)
            #print "removed",tmp_molefrag["clusters"][icluster]
            #print "old",tmp_molefrag["clusters"][iselect]
            #tmp_molefrag["clusters"][iselect].append(icomp)
            #print "new",tmp_molefrag["clusters"][iselect]
        #print tmp_molefrag["clusters"] 
        # done
    #exit()
    return True


def is_p_block_element(elem):
    from elements import ELEMENTS
    ielem = ELEMENTS[elem].number
    return (5<=ielem<=10) or (13<=ielem<=18) or (31<=ielem<=36) or (49<=ielem<=54) \
           or (81<=ielem<=86) or (113<=ielem<=118)
#end of is_p_block_element

# wzk20200623: added ifpho
def build_buffered_molefrag(molecule,molecomp,molecompdist,cenfrag,molecomp_connect,control,\
                            ifpho):
    natom       = molecule.natom 
    nmolecomp   = len(molecomp)
    ncfrag      = len(cenfrag)

    radii       = control[0]
    onlyconnect = control[1]

    # wzk20201006: there are now two possible formats of radii:
    # (1) give one radius for each fragment
    # (2) use the same radius for all fragments
    # Now we convert format (2) to format (1)
    if len(radii) == 1:
        radii = radii * ncfrag

    # wzk20201012: the closest distance of link/environment atoms from the center region.
    # This is useful for iOI fragment enlarging
    minrenv = [100.] * (ncfrag+1)

    molefrag={}
    for ifrag in range(1,ncfrag+1):
         molefrag[ifrag]={}
        
         iused={}
         for icomp in range(1,nmolecomp+1):
             iused[icomp]=0

         # add center atoms
         cenatom = []
         for icomp in cenfrag[ifrag]:
             iused[icomp]=1
             # add atoms into center fragment
             for iatom in molecomp[icomp]:
                  cenatom.append(iatom)
         # wzk20211214: sort cenatom to prevent problems when nfrag==1
         cenatom.sort()
         #print cenfrag[ifrag],onlyconnect
         #print "cenatom",cenatom

         # wzk20201115: if radii<0, we automatically determine radii assuming that (a) the
         # DSmetric is used and (b) the radii of the last step is -radii
         if radii[ifrag-1] < 0.0:
             vmdist=100.0
             for icomp in range(1,nmolecomp+1): # cycle all molecomp 
                if iused[icomp] == 1:
                    continue
                else:
                    for jcomp in cenfrag[ifrag]:
                        i=max(icomp,jcomp)
                        j=min(icomp,jcomp)
                        dist=molecompdist[i][j]
                        #print i,j,dist
                        if dist < vmdist and dist > -radii[ifrag-1]+1e-10:
                            vmdist = dist
             radii[ifrag-1] = vmdist+1.0

         # add buffer molecomp within a given radii
         bufcomp = []
         for icomp in range(1,nmolecomp+1): # cycle all molecomp 
            if iused[icomp] == 1:
                continue
            else:
                vmdist=1E30 
                for jcomp in cenfrag[ifrag]:
                    i=max(icomp,jcomp)
                    j=min(icomp,jcomp)
                    dist=molecompdist[i][j]
                    #print i,j,dist
                    if dist < vmdist:
                        vmdist = dist
                iused[icomp] = 1 # wzk: redundant?
                if vmdist <= radii[ifrag-1]:
                    bufcomp.append(icomp)
                else:
                    minrenv[ifrag] = min(minrenv[ifrag],vmdist)
         
         bufatom = []
         for icomp in bufcomp:
            linked = True
            if onlyconnect:  # only add direct connected molecomp as buffer
                linked = False
                # if connect with center molecomp
                for ic in cenfrag[ifrag]:
                    i = min(ic,icomp)
                    j = max(ic,icomp)
                    lin = molecomp_connect[i][j]
                    if lin > 0.0001:
                       linked = True
                #print icomp, linked
                #print bufcomp
                # if connect with buffer molecomp
                if not linked:
                    for ic in bufcomp:
                        if ic != icomp:
                            i = min(ic,icomp)
                            j = max(ic,icomp)
                            lin = molecomp_connect[i][j]
                            if lin > 0.0001:
                               linked = True
            #print icomp,linked
            if not linked: 
                continue
            else:    
                # pick up atoms in this molecomp as buffer
                #print icomp,vmdist
                for iatom in molecomp[icomp]:
                    bufatom.append(iatom)
         #print "buffatom", bufatom

         # here should add some self consistent checking for fragment

         # add link atoms, assum all center atoms link with center or buffer atom
         # only add link atom to some buffer atoms
         linkatom=[]
         # wzk20200630: the link atom search must be repeated several times, since
         # environment atoms that are bonded to two or more buffer/center atoms are
         # identified as new buffer atoms, which may have link atoms of its own
         newbuf = True
         while newbuf == True:
             newbuf = False
             #for batom in bufatom:
             # wzk20200630: for very small buffer radii, some link atoms might bond
             # directly to center atoms
             for batom in bufatom + cenatom:
                 leni = len(molecule.connect[batom])
                 #if leni == 1: # wzk20200630
                 #    continue
                 if molecule.edges[batom][batom] > 0:
                     # check atom with index small than present one
                     for catom in range(1,batom):
                         if molecule.edges[catom][batom] > 0 :
                             #print "catom",catom
                             # check if this atom belongs to ceter and buffer atoms
                             if catom in bufatom:
                                 continue
  
                             if catom in cenatom:
                                 continue
  
                             # move link atom to buffer region if it is the link atom
                             # of another atom in the fragment, or if the atom to be
                             # substituted is not a p-block element, or if the bond
                             # to the link atom is not a single bond, or if all environment
                             # atoms bonded to this atom are hydrogens
                             if catom in linkatom:
                                 #touch = True
                                 linkatom.remove(catom)
                                 bufatom.append(catom)
                                 #print "new buf",catom
                                 newbuf = True
                                 continue
                             #if molecule.label[catom] == 'H':
                             if not is_p_block_element(molecule.label[catom]):
                                 bufatom.append(catom)
                                 newbuf = True
                                 continue
                             has_nonH_envatom = False
                             for ia in range(1,len(molecule.connect[catom]),2):
                                 if molecule.connect[catom][ia] == batom:
                                     break
                                 elif molecule.label[molecule.connect[catom][ia]] != 'H':
                                     has_nonH_envatom = True
                             if molecule.connect[catom][ia+1] != 1.0:
                                 bufatom.append(catom)
                                 newbuf = True
                                 continue
                             for ja in range(1,catom):
                                 if molecule.label[ja] != 'H':
                                     for ia in range(1,len(molecule.connect[ja]),2):
                                         if molecule.connect[ja][ia] == catom:
                                             has_nonH_envatom = True
                             if not has_nonH_envatom:
                                 bufatom.append(catom)
                                 newbuf = True
                                 continue
  
                             linkatom.append(catom)
                             #print "catom",catom,batom

                     # check atom with index larger than present one
                     # if this atom only link with one atom, neglect it.
                     for ia in range(1,leni,2):
                         catom = molecule.connect[batom][ia]
                         #print "catom",catom
                         # check if this atom belongs to ceter and buffer atoms
                         if catom in bufatom:
                             continue
 
                         if catom in cenatom:
                             continue
 
                         if catom in linkatom:
                             #touch = True
                             linkatom.remove(catom)
                             bufatom.append(catom)
                             #print "new buf",catom
                             newbuf = True
                             continue
                         if not is_p_block_element(molecule.label[catom]):
                             bufatom.append(catom)
                             newbuf = True
                             continue
                         if molecule.connect[batom][ia+1] != 1.0:
                             bufatom.append(catom)
                             newbuf = True
                             continue
                         has_nonH_envatom = False
                         for ja in range(1,len(molecule.connect[catom]),2):
                             if molecule.label[molecule.connect[catom][ja]] != 'H':
                                 has_nonH_envatom = True
                         for ja in range(1,catom):
                             if molecule.label[ja] != 'H' and ja != batom:
                                 for ka in range(1,len(molecule.connect[ja]),2):
                                     if molecule.connect[ja][ka] == catom:
                                         has_nonH_envatom = True
                         if not has_nonH_envatom:
                             bufatom.append(catom)
                             newbuf = True
                             continue
 
                         linkatom.append(catom)
                         #print "catom",catom,batom

             if newbuf:
                 linkatom = []

         # wzk20200623: for PHO calculations, one should find the environment (MM) atoms
         # that are directly bonded to the boundary atom
         # NOTE: environatom contains the boundary atom itself, as required by PHO
         if ifpho:
             environatom = []
             for latom in linkatom:
                 environatom.append(latom)
                 leni = len(molecule.connect[latom])
                 if leni == 1:
                     continue
                 if molecule.edges[latom][latom] > 0:
                     # check atom with index small than present one
                     for catom in range(1,latom):
                         if molecule.edges[catom][latom] > 0 :
                             if not (catom in cenatom or catom in bufatom or catom in linkatom or catom in environatom):
                                 environatom.append(catom)
                     for ia in range(1,leni,2):
                         catom = molecule.connect[latom][ia]
                         if not (catom in cenatom or catom in bufatom or catom in linkatom or catom in environatom):
                             environatom.append(catom)
                             #print "environ atom",latom,catom

         #print "linkatom ",linkatom
         # check spin multiplicit
         #exit()

         molefrag[ifrag]=[]
         molefrag[ifrag].append(cenatom)
         molefrag[ifrag].append(bufatom)
         molefrag[ifrag].append(linkatom)
         if ifpho: # wzk20200623
             molefrag[ifrag].append(environatom)
         else:
             molefrag[ifrag].append([])
         # wzk20201012: append radbuff and minrenv to molefrag to pass this information out
         molefrag[ifrag].append(radii[ifrag-1])
         molefrag[ifrag].append(minrenv[ifrag])
         #print ifrag, molefrag[ifrag]

    return molefrag 
#end of build_buffered_molefrag

# guess initial cluster by user inputed atom lists
def guess_from_input_atoms(atomlist, molecomp):
    nfrag = len(atomlist)
    print(atomlist)
    centroids = []
    for iatom in atomlist :
        for icomp in range(1,len(molecomp)+1) :
            if int(iatom) in molecomp[icomp]:
                centroids.append([icomp])
    print("\n Initial cluster centroids from input atom list")
    print(" ",centroids)

    return centroids         
#end of guess_from_input_atoms

# auto guess initial clusters
def auto_guess_init_clusters(molecomp_graph, nfragments):
    connected_graph = molecomp_graph.connected_components_deepthfirstscan()
    print("depthfirstscan",connected_graph[1])
    #connected_graph = molecomp_graph.connected_components_breadthfirstscan()
    #print "breadthfirst",connected_graph

    complist = []
    for i in range(1, len(connected_graph)+1):
        complist = complist + connected_graph[i]

    #print complist
    #exit()

    # splits
    leni  = molecomp_graph.natom
    ncomp = leni/nfragments
    
    import random
    #clusters = [complist[i:i+min(ncomp, leni-i),1) for i in range(0, leni, ncomp)]
    clusters = []
    ii = 0
    for i in range(0, nfragments):
        jj = min(ncomp, leni-ii)
        clusters = clusters + [random.sample(complist[ii:ii+jj],1)]   
        ii = ii + jj 

    #if len(clusters) > nfragments:
  
    #connected_graph = molecomp_graph.connected_components_breadthfirstscan()
    #print "breadthfirst",connected_graph
    print(clusters)
    #exit()
    return clusters
#end of auto_guess_init_clusters


def auto_guess_init_clusters_new(molecomp_graph, molecomp, mindist, radcent):
    connected_graph = molecomp_graph.connected_components_deepthfirstscan()
    #complist = []
    #for i in range(1, len(connected_graph)+1):
    #    complist = complist + connected_graph[i]

    clusters = []
    #remaining_molecomp = range(1, len(molecomp)+1)
    remaining_molecomp = list(range(1, len(molecomp)+1))

    # 1. If a component contains more than maxatom atoms, select it as a centroid
    maxatom = 5
    for i in remaining_molecomp:
        if i == 0:
            continue
        if len(molecomp[i]) > maxatom:
            clusters = clusters + [[i]]
            remaining_molecomp[i-1] = 0 # this is more efficient than remove()

    # 2. Sort the molecomps in ascending order according to the number of neighbors. In case
    # of degeneracy, sort in descending order of the distance to the furthest neighbor
    nneighbor = [0]*len(remaining_molecomp)
    maxdist = [0.0]*len(remaining_molecomp)
    for i,ic in enumerate(remaining_molecomp):
        if ic == 0: continue
        for j,jc in enumerate(remaining_molecomp[i+1:]):
            if jc == 0: continue
            if ic > jc:
                dist = mindist[ic][jc]
            else:
                dist = mindist[jc][ic]
            if dist < radcent:
                nneighbor[i] += len(molecomp[j+1])
                nneighbor[j] += len(molecomp[i+1])
                maxdist[i] = max(maxdist[i],dist)
                maxdist[j] = max(maxdist[j],dist)
    tmp = [(nneighbor[i],-maxdist[i],remaining_molecomp[i]) for i in range(len(remaining_molecomp))]
    tmp.sort()
    remaining_molecomp = [c for n,m,c in tmp]

    # 3. If a component is attached to at least nconnect components, select the
    #    former as a centroid and remove the latter components from consideration
    #    Herein nconnect is taken as 4,3,2,1,0 in that order
    for nconnect in range(4,-1,-1):
        for i,ic in enumerate(remaining_molecomp):
            if ic == 0:
                continue
            leni = len(molecomp_graph.connect[ic])
            remaining_connect = []
            for j in range(1,leni,2):
                try:
                    jc = remaining_molecomp.index(molecomp_graph.connect[ic][j])
                    remaining_connect.append(jc)
                except ValueError: # not found
                    pass
            if len(remaining_connect) >= nconnect:
                clusters = clusters + [[ic]]
                remaining_molecomp[i] = 0
                for j in remaining_connect:
                    remaining_molecomp[j] = 0

    # 4. If two centroids are too close to each other, remove the one with more neighbor atoms
    for i,ic in enumerate(clusters):
        if ic == [0]:
            continue
        for j,jc in enumerate(clusters[i+1:]):
            if jc == [0]:
                continue
            if ic[0] > jc[0]:
                dist = mindist[ic[0]][jc[0]]
            else:
                dist = mindist[jc[0]][ic[0]]
            if dist < radcent:
                clusters[j] = [0]
    # remove all zeroes
    #clusters = filter(lambda x: x!=[0], clusters)
    clusters = [x for x in clusters if x!=[0]]

    print("Auto guess cluster centroids:")
    print(clusters)
    return clusters

#end of auto_guess_init_clusters_new

def check_bond(bond):
    if not isinstance(bond[0],str):
        return False
    if not isinstance(bond[1],str):
        return False
    if not isinstance(bond[2],float):
        return False
    return True