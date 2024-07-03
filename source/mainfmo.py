#!/usr/bin/python
import amino_acid
import sys
import fragment

#最初的eMol是将gjf文件识别，在转化成多个bdf输入文件，我是将pdb文件识别，但是这样没有设置电荷和自旋多重度

#
def Read_GaussInput(fn):   

    coord = {}
    # Read fragments information
    lines   =fn.readlines()
    lenlines=len(lines)

    il = -1
    for str1 in lines:
        il = il +1
        str1=str1.strip("\n")
        lines[il]=str1

    iline=-1    
    iset =0
    for str1 in lines:
        iline=iline+1
        # in gaussian input, line starts with % is for system definition
        if str1.startswith("%"):
            continue
        elif str1.startswith("#"):
            # read information from command line
            gausscmd=str1
            break

    bohr = False
    if "bohr" in gausscmd.lower():
        bohr = True
        print(" Unit = Bohr")
    # wzk20200721
    inpconnect = False
    if "connectivity" in gausscmd.lower():
        inpconnect = True

    # decide number of fragment
    str1 = lines[iline+4].strip("\n")
    chargespin = str1
    str2 = str1.split()
    leni = len(str2) + 1
    nfrag=int(leni/2)
    if nfrag == 1:
        print("No fragment information in input ...")
        totcharge = int(str2[0])
        totspin = int(str2[1])
    else:
        nfrag = nfrag-1
        # TODO: constrain individual fragment charges/spins
        totcharge = sum(map(int,str2[0:leni-1:2]))
        totspin = sum(map(int,str2[1:leni-1:2]))
    print("Number of Fragments: ", nfrag)

    molecomp = {}
    # read and process geometry information
    iline  = iline+5
    linkatom = []
    if nfrag == 1:
        # gaussian input does not contain fragment information
        # read geometry only
        icoord = 0
        for il in range(iline, lenlines):
            icoord = icoord+1
            str1 = lines[il]
            print(str1)
            if str1 == "" or len(str1)<3:  # Null line treat as end of geometry input
                break
            coordt = str1.split()
            first_str=coordt[0]
            atomf = first_str[0]
            ct1   = []
            ct1.append(atomf)  # atom label check by Yibo Lei 
            for i in range(2,5):  # cartesian coordinate  ##add a length judgement of string ###########################################################
                if bohr :
                    x1   = float(coordt[i])*0.52917721092
                    str1 = "%20.12f" % x1 
                    ct1.append(str1)
                else:
                    ct1.append(coordt[i])
            coord[icoord]=ct1      
            # wzk20201116: if the global system itself contains link atoms, label them
            if coordt[-1].lower() == 'l':
                linkatom.append(icoord)
        molecomp[1]=[]      
        molecomp[1].append(icoord) 
    else:
        icoord = 0
        for il in range(iline, lenlines):
            icoord = icoord+1
            str1 = lines[il]
            print(str1)
            if str1 == "":  # Null line treat as end of geometry input
                break
            coordt = str1.split()
            atomf = coordt[0].split("(")
            ct1   = []
            ct1.append(atomf[0])
            for i in range(2,5):    ####################1 4->2 5
                if bohr :
                    x1   = float(coordt[i])*0.52917721092
                    str1 = "%20.12f" % x1 
                    ct1.append(str1)
                else:
                    ct1.append(coordt[i])
            coord[icoord]=ct1

            atomf = atomf[1].strip(")")
            str2  = atomf.split("=")
            ifrag = int(str2[1])
            if ifrag in molecomp:
                molecomp[ifrag].append(icoord)
            else:
                molecomp[ifrag]=[]
                molecomp[ifrag].append(icoord)
    natom = icoord-1

    print("^^^^^^^^^^^^^^^^^^^^^^^^")
    print("linkatom=",linkatom)

    iline = iline+natom+1
    # read connectivity information
    icon = 0
    #inpconnect = True
    if inpconnect:
        for il in range(iline,lenlines):
             icon = icon+1
             print(lines[il])
             str1 = lines[il].split()
             if len(str1) == 0:
                 break
             na = int((len(str1)-1)/2)
             ia = int(str1[0])
             if na > 0:
                 cn = []
                 cn.append(na)
                 for i in range(1,2*na+1,2): 
                     cn.append(int(str1[i]))
                     cn.append(float(str1[i+1]))
                 coord[ia].append(cn)
             else:
                 #cn = []
                 #cn.append(1)
                 cn = [0] # wzk20200630
                 coord[ia].append(cn)

    # in gaussian style input, we may input several atoms as inital guess of fragment.
    centeratoms = None
    cutbond     = [] 
    # wzk20200630: the default value is unnecessarily large. Reduce 4.5 to 2.0
    radbuff     = [2.0] # 4.5
    radcent     = 3.0
    nomergelist = []
    # wzk20200701: allow specifying atomic charge and spin
    atomcharge = []
    atomspin = []
    ilinen = iline + icon
    #if iline <= len(lines)-1 :
    for iline in range(ilinen-1,len(lines)):
        str1 = lines[iline].lower()
        if str1.startswith("radbuff"):
            #radbuff = lines[iline+1].split()[0]
            # wzk20201006: now support specifying multiple radbuffs, one for each fragment
            radbuff = lines[iline+1].split()
            for i in range(0,len(radbuff)):
                radbuff[i] = portable_float(radbuff[i])
        elif str1.startswith("radcent"):
            radcent = portable_float(lines[iline+1])
        elif str1.startswith("centeratoms"):
            centeratoms = lines[iline+1].split()
        elif str1.startswith("nomerge"):
            # wzk20201128: do not merge certain subsystems
            nomergelist = map(int,lines[iline+1].split())
            if nomergelist[0] < 0:
                nomergelist = []
        elif str1.startswith("cutbond"):
            print(lines[iline+1])
            list1   = lines[iline+1].split()
            nbonds  = int(list1[0])
            ib = 0
            for ii in range(1,3*nbonds+1,3):
                cutbond.append(list1[ii:ii+3]) 
                cutbond[ib][2]=float(cutbond[ib][2])     
                ib = ib + 1
        elif str1.startswith("charge"):
            atomcharge = map(int,lines[iline+1].split())
        elif str1.startswith("spin"):
            atomspin = map(int,lines[iline+1].split())
    #print cutbond
    #exit()

    if not inpconnect:  # determine connectivity information automaticlly
        print("Determine connectivity information automatically")
        auto_decide_connectivity(coord,atomcharge,atomspin)

    ## read user defined center fragment information
    #ncfrag=0
    #iline = iline+icon+1
    #print iline,lenlines #,lines[iline]
    ##exit()
    #if iline > lenlines:
    #   # user did not set center fragments
    #   for ifrag in range(1,nfrag+1):
    #       ncfrag=ncfrag+1
    #       cenfrag[ncfrag]=[]
    #       cenfrag[ncfrag].append(1)
    #       cenfrag[ncfrag].append(ifrag)
    #else:
    #   # user has set center fragments
    #   for il in (iline,lenlines):
    #       str1=lines[il]
    #       if str1 == "":
    #           break
    #       else:
    #           ncfrag=ncfrag+1     
    #           str1=str1.split()
    #           cenfrag[ifrag]=str1

    #print "Molecule geometry ..."    
    print("Connectivity in gjf style:")
    for il in coord:
        print(' ' + str(il) + ' ' + ' '.join(str(coord[il][4][i]) for i in range(1,len(coord[il][4]))))

    outinfo = {}
    outinfo["coord"]       = coord
    outinfo["centeratoms"] = centeratoms
    outinfo["cutbond"]     = cutbond 
    outinfo["molecomp"]    = molecomp
    outinfo["radbuff"]     = radbuff
    outinfo["radcent"]     = radcent
    outinfo["atomcharge"]  = atomcharge
    outinfo["atomspin"]    = atomspin
    outinfo["totcharge"]   = totcharge
    outinfo["totspin"]     = totspin
    outinfo["linkatom"]    = linkatom
    outinfo["nomergelist"] = nomergelist
    return outinfo 
#end of Read_GaussInput

def read_PDB(fn):    
        coord={}
        lines   =fn.readlines()
        for str1 in lines:
            if str1.startswith("ATOM") or str1.startswith("HETATM"):
                icoord=int(str1.split()[1])    #原子序号
                atomf=str1.split()[-1]     #元素
                ct1   = []
                ct1.append(atomf)  
                for i in range(3,len(str1.split())-1):  
                    if len(str1.split()[i])>4 and str1.split()[i][-4]=="." and str1.split()[i+1][-4]=="." and str1.split()[i+2][-4]==".":  
                        for j in range(i,i+3):
                            ct1.append(str1.split()[j])   #输入坐标
                        break
                coord[icoord]=ct1   
                coord[icoord].append("") 
        for str1 in lines:
            if str1.startswith("ATOM") or str1.startswith("HETATM"):
                icoord=int(str1.split()[1])
                coord[icoord].append(str1.split()[2])
                coord[icoord].append(str1.split()[4]) 
                coord[icoord].append(str1.split()[3])
            if str1.startswith("CONECT"):
                st=str1.split()
                icoord=int(st[1])
                str2=[]
                str2.append("CONNECT")
                str2.append(st[1])
                for ia in st[2:]:
                    if icoord<int(ia):
                        str2.append(ia)
                cn=[]
                cn.append(len(set(str2[2:])))
                i=0
                tmp={}
                for at in set(str2[2:]): 
                    i+=1
                    tmp[i]=at
                for i in range(1,len(set(tmp))+1):
                    n=0
                    for j in range(2,len(str2)):
                        if str2[j]==tmp[i]:
                            n+=1
                    cn.append(int(tmp[i]))
                    cn.append(float(n))    
                coord[icoord][4]=cn

        #在此之前记录的苯环中的碳原子是交替的单双键，此步骤是把单双键全换成1.5键
        # for i in range(1,len(coord)+1):
        #     if coord[i][5]=="CD1" and coord[i+1][5]=="CD2" and coord[i+2][5]=="CE1" and coord[i+3][5]=="CE2" and coord[i+4][5]=="CZ":
        #         ind=coord[i-1][4].index(i)
        #         coord[i-1][4][ind+1]=1.5
        #         ind=coord[i-1][4].index(i+1)
        #         coord[i-1][4][ind+1]=1.5          
        #         ind=coord[i][4].index(i+2)
        #         coord[i][4][ind+1]=1.5
        #         ind=coord[i+1][4].index(i+3)
        #         coord[i+1][4][ind+1]=1.5
        #         ind=coord[i+2][4].index(i+4)
        #         coord[i+2][4][ind+1]=1.5
        #         ind=coord[i+3][4].index(i+4)
        #         coord[i+3][4][ind+1]=1.5        
        # {40: ['C', '31.113', '7.238', '-5.638', '[3, 41, 2.0, 39, 1.0]', 'CB', 'A'，'CYS'], 
        #  41: ['O', '30.585', '6.536', '-4.764', '[1, 40, 2.0]', 'O', 'B','HOH'],}
        
        peptideChain={}
        for str1 in lines:   
            if str1.startswith("SEQRES"):
                parts = str1.split()
                chain_id = parts[2]
                amino_acids = parts[4:]
                if chain_id not in peptideChain:
                    peptideChain[chain_id] = {}
                index = len(peptideChain[chain_id]) + 1
                for aa in amino_acids:
                    peptideChain[chain_id][index] = aa
                    index += 1
        #peptideChain     {'A': {1: 'GLY', 2: 'VAL', 3: 'VAL'}, 
        #                  'B': {1: 'ASN', 2: 'SER', 3: 'LEU'}}
        
        # amino={}
        # chain_n=0
        # label=""
        # chain={}
        # for str1 in lines:   
        #     if str1.startswith("SEQRES"):
        #         if label==str1.split()[2]:
        #             len_lines_split=len(str1.split())
        #             for ires in range(1,len_lines_split-3):
        #                     ordinal+=1
        #                     chain["res"+str(ordinal)]=str1.split()[ires+3]
        #             amino[chain_n]=chain
        #         else: 
        #             chain_n+=1 
        #             chain={}                       #再想第一次空的值，不等于的时候应该怎么补充
        #             label=str1.split()[2]           #要先想一次已经等于A了，相同的时候该怎么做      
        #             chain["chain_label"]=str1.split()[2]
        #             chain["tot_res"]=str1.split()[3]
        #             len_lines_split=len(str1.split())
        #             ordinal=0
        #             for ires in range(1,len_lines_split-3):
        #                 ordinal+=1
        #                 chain["res"+str(ordinal)]=str1.split()[ires+3]
                    # {1: {'chain_label': 'A', 'tot_res': '22', 'res1': 'ASN', 'res2': 'SER'}, 
                    #  2: {'chain_label': 'B', 'tot_res': '24', 'res1': 'ASN', 'res2': 'SER'}}
        #突然感觉这个信息有点鸡肋，用下一个肽链包含的氨基酸信息，就可以得出这个
        
        
        # peptideChain={}
        # temp=[]
        # for str1 in lines:  
        #         if str1.startswith("ATOM") or str1.startswith("HETATM"):
        #             temp.append(str1.split()[4])
        # i=0
        # chain_name={}
        # for a in set(temp):
        #         i+=1
        #         chain_name[i]=a
        #         peptideChain[a]=[]     
        # for str1 in lines:
        #         if str1.startswith("ATOM") or str1.startswith("HETATM"):
        #             for i in range(1,len(set(temp))+1):
        #                 if str1.split()[4]==chain_name[i]:
        #                     peptideChain[chain_name[i]].append(str1.split()[1]) 
        #                 #chain_name[i]  <-> i   字典的键为数字1、2还是字母A、B 
        # #{'A': [2, 3, 4, 5], 'B': [1, 6]}
        radbuff     = [2.0] # 4.5

        outinfo = {}
        outinfo["coord"]       = coord
        outinfo["peptideChain"]=peptideChain
        outinfo["radbuff"]     =radbuff
        return outinfo

def auto_decide_connectivity(coord,charge0,spin0):
    import math
    from elements import ELEMENTS
    # Algorithm:
    # (1) Calculate the Pauling bond order of all atom pairs.
    # (2) If bond order <= 0.5, treat as non-bonded.
    #     NOTE: metal-metal bonds between s, d and f-block metal elements are
    #     always treated as unbonded. This allows for meaningful fragmentation
    #     calculations of metal clusters.
    # (3) Calculate free valence of all B,C,N,O atoms (we assume that only these
    #     atoms can form multiple bonds)
    # (4) If an atom has free valence >0 and, among atoms bonded to it, only one
    #     has free valence, then form a pi bond between them. This should include
    #     most double bonds. Iterating this once again gives most triple bonds
    # (5) Form resonance bonds between all remaining adjacent pairs of atoms that
    #     both have free valence

    natom = len(coord)
    coordnum = [0]*(natom+1)
    for iatom in range(1,natom+1):
        aname = coord[iatom][0]
        irad = ELEMENTS[aname].covrad
        iele = ELEMENTS[aname].number
        cn = [0]
        # Only include bonds I-J where I<J in the connectivity info
        for jatom in range(iatom+1,natom+1):
            aname1 = coord[jatom][0]
            jrad = ELEMENTS[aname1].covrad
            jele = ELEMENTS[aname1].number
            if (3<=iele<=4 or 11<=iele<=12 or 19<=iele<=30 or 37<=iele<=48 or \
                55<=iele<=80 or 87<=iele<=103) and \
                (3<=jele<=4 or 11<=jele<=12 or 19<=jele<=30 or 37<=jele<=48 or \
                55<=jele<=80 or 87<=jele<=103):
                continue
            dist = math.sqrt((float(coord[iatom][1])-float(coord[jatom][1]))**2 \
                + (float(coord[iatom][2])-float(coord[jatom][2]))**2 \
                + (float(coord[iatom][3])-float(coord[jatom][3]))**2)
            rcov = irad+jrad
            bondorder = math.exp((rcov-dist)/(0.3*rcov))
            if dist < 1e-6: # bond order between an atom and itself is 0
                bondorder = 0
            if bondorder <= 0.5:
                continue
            else:
                bondorder = 1.0
            cn[0] += 1
            cn.append(jatom)
            cn.append(bondorder)
            coordnum[iatom] += 1
            coordnum[jatom] += 1
        coord[iatom].append(cn)

    # Calculate max coordination number
    maxcoord = [0]*(natom+1)
    for iatom in range(1,natom+1):
        aname = coord[iatom][0]
        if aname == 'C':
            maxcoord[iatom] = 4
        elif aname == 'O':
            maxcoord[iatom] = 2
        elif aname in ['B','N']:
            maxcoord[iatom] = 3
    for ic in range(0,len(charge0),2):
        iatom = charge0[ic]
        aname = coord[iatom][0]
        if aname == 'C':
            maxcoord[iatom] -= abs(charge0[ic+1])
        elif aname in ['N','O']:
            maxcoord[iatom] += charge0[ic+1]
        elif aname == 'B':
            maxcoord[iatom] -= charge0[ic+1]
    for ic in range(0,len(spin0),2):
        iatom = spin0[ic]
        maxcoord[iatom] -= abs(spin0[ic+1])

    # Calculate free valence
    freeval = [0]*(natom+1)
    for iatom in range(1,natom+1):
        freeval[iatom] = max(maxcoord[iatom] - coordnum[iatom],0)

    # Most double and triple bonds
    updated = True
    while updated: # loop until no new pi bonds are found this way
        updated = False
        neighbor_freeval = [0]*(natom+1) # how many neighbors have free valence
        for iatom in range(1,natom+1):
            leni = len(coord[iatom][4])
            if freeval[iatom] == 0:
                continue
            for ic in range(1,leni,2):
                jatom = coord[iatom][4][ic]
                if freeval[jatom]>0:
                    neighbor_freeval[iatom]+=1
                    neighbor_freeval[jatom]+=1

        for iatom in range(1,natom+1):
            leni = len(coord[iatom][4])
            if freeval[iatom] == 0:
                continue
            for ic in range(1,leni,2):
                jatom = coord[iatom][4][ic]
                if freeval[jatom] == 0:
                    continue
                if (neighbor_freeval[iatom]==1 and neighbor_freeval[jatom]>0) \
                    or (neighbor_freeval[iatom]>0 and neighbor_freeval[jatom]==1):
                    coord[iatom][4][ic+1] += 1 # a new bond
                    freeval[iatom] -= 1
                    freeval[jatom] -= 1
                    neighbor_freeval[iatom] -= 1
                    neighbor_freeval[jatom] -= 1
                    updated = True
                    if freeval[iatom] == 0:
                        break

    # Resonance bonds
    nresonant = [None]
    for iatom in range(1,natom+1):
        nresonant.append([0,0,0])
    for iatom in range(1,natom+1):
        leni = len(coord[iatom][4])
        for ic in range(1,leni,2):
            jatom = coord[iatom][4][ic]
            if freeval[iatom]>0 and freeval[jatom]>0 and \
               neighbor_freeval[iatom]>0 and neighbor_freeval[jatom]>0:
                coord[iatom][4][ic+1] += 0.5
                nresonant[iatom][0] += 1
                nresonant[jatom][0] += 1
                if nresonant[iatom][0] <= 2:
                    nresonant[iatom][nresonant[iatom][0]] = jatom
                if nresonant[jatom][0] <= 2:
                    nresonant[jatom][nresonant[jatom][0]] = iatom

    # If 2n atoms are connected circularly by resonance bonds, and they are not
    # connected to other atoms via resonance bonds (except possibly for only one
    # of the atoms), transform them into Kekule single/double bonds.
    # Of the two Kekule forms possible, the one more consistent
    # with the bond lengths is chosen. This helps keep the orbital localization
    # patterns of adjacent subsystems consistent.
    # This treatment suffices to convert all bonds to single/double/triple bonds
    # for most important biomolecules, even including porphyrins, but still leaves
    # e.g. naphthalene in resonant form
    for iatom in range(1,natom+1):
        if nresonant[iatom][0] != 2: continue
        jatom0 = iatom
        jatom = nresonant[iatom][1]
        ring = [jatom]
        n_external_resonant_bonds = 0
        while jatom != iatom:
            if nresonant[jatom][0] != 2:
                n_external_resonant_bonds = n_external_resonant_bonds+1
                if n_external_resonant_bonds > 1: break
            if nresonant[jatom][1] == jatom0:
                jatom0 = jatom
                jatom = nresonant[jatom][2]
            else:
                jatom0 = jatom
                jatom = nresonant[jatom][1]
            ring.append(jatom)
        if jatom != iatom: continue
        if len(ring)%2 == 1: continue
        # compare the sum of double bond lengths of the two Kekule forms
        oddbond = 0.
        evenbond = 0.
        for i in range(0,len(ring),2):
            jatom = ring[i]
            katom = ring[i+1]
            dist = math.sqrt((float(coord[katom][1])-float(coord[jatom][1]))**2 \
                + (float(coord[katom][2])-float(coord[jatom][2]))**2 \
                + (float(coord[katom][3])-float(coord[jatom][3]))**2)
            oddbond += dist
            jatom = ring[(i+2)%len(ring)]
            dist = math.sqrt((float(coord[katom][1])-float(coord[jatom][1]))**2 \
                + (float(coord[katom][2])-float(coord[jatom][2]))**2 \
                + (float(coord[katom][3])-float(coord[jatom][3]))**2)
            evenbond += dist
        # Zero out nresonant, so that subsequent iatom cycles will not detect this
        # conjugated ring again
        for i in range(len(ring)):
            jatom = ring[i]
            nresonant[jatom][0] = 0
        for i in range(len(ring)):
            jatom = ring[i]
            leni = len(coord[jatom][4])
            for ic in range(1,leni,2):
                 if coord[jatom][4][ic+1] != 1.5: continue
                 katom = coord[jatom][4][ic]
                 if (oddbond < evenbond and i%2==0) or (oddbond >= evenbond and i%2==1):
                     if katom == ring[(i+1)%len(ring)]:
                         coord[jatom][4][ic+1] = 2.0
                     else:
                         coord[jatom][4][ic+1] = 1.0
                 else:
                     if katom == ring[(i+1)%len(ring)]:
                         coord[jatom][4][ic+1] = 1.0
                     else:
                         coord[jatom][4][ic+1] = 2.0


def default_cutbond():
    # Single bonds between the following elements can be cut
    # Note that it is pointless to include halogens because they almost
    # always show up at the end of molecules
    elem = ['Li', 'Be', 'B', 'C', 'N', 'O', \
            'Na', 'Mg', 'Al', 'Si', 'P', 'S', \
            'K', 'Ca', 'Ga', 'Ge', 'As', 'Se', \
            'Rb', 'Sr', 'In', 'Sn', 'Sb', 'Te', \
            'Cs', 'Ba', 'Tl', 'Pb', 'Bi', 'Po', \
            'Fr', 'Ra']
    cutbond = [[str1, str2, 1.0] for str1 in elem for str2 in elem]
    return cutbond

# Assign the charge and multiplicity of each fragment.
def assign_fragment_chargespin(coord, molefrag, charge0, spin0, totcharge, totspin,\
                               connect):
    from elements import ELEMENTS

    nfrag = len(molefrag)
    natom  = len(coord)

    charge = [0]*(nfrag+1)
    spin = [0]*(nfrag+1)
    atomcharge = [0]*(natom+1)
    atomspin = [0]*(natom+1)

    # 0. Count number of bonds bonded to every atom
    #    However, we:
    #    (1) Do not count any bond involving s-block metals;
    #    (2) Zero out the charge of any atom that is connected to a transition metal atom.
    #        This is equivalent to assuming the first coordination sphere of transition metals
    #        is always neutral
    coordnum = [0]*(natom+1)
    is_bonded_to_transition_metal = [False]*(natom+1)
    for iatom in range(1,natom+1):
        leni = len(connect[iatom])
        aname = coord[iatom][0]
        ele = ELEMENTS[aname].number
        is_transition_metal = False
        if 21<=ele<=30 or 39<=ele<=48 or 57<=ele<=80 or 89<=ele<=112:
            is_transition_metal = True
        for ic in range(1,leni,2):
            jatom = connect[iatom][ic]
            if is_transition_metal:
                is_bonded_to_transition_metal[jatom] = True
            bondorder = connect[iatom][ic+1]
            aname2 = coord[iatom][0]
            ele2 = ELEMENTS[aname].number
            if ele2 not in [3,4,11,12,19,20,37,38,55,56,87,88] and \
                ele not in [3,4,11,12,19,20,37,38,55,56,87,88]:
                coordnum[iatom] += bondorder
                coordnum[jatom] += bondorder

    # 1. Deduce atomic charges, temporarily ignoring the total charge constraint
    # Rule:
    # (1) N, P, As, Sb, Bi are cationic when they have 4 bonds attached, and are anionic
    #     when they have 2 bonds attached;
    # (2) O, S, Se, Te, Po are cationic when they have 3 bonds attached, and are anionic
    #     when they have 1 bond attached;
    # (3) F, Cl, Br, I, At are anionic when they have no bonds attached;
    # (4) B, Al, Ga, In, Tl are anionic when they have 4 bonds attached;
    # (5) Li, Na, K, Rb, Cs, Fr always have +1 charge
    # (6) Be, Mg, Ca, Sr, Ba, Ra always have +2 charge
    # (7) The charge of other metals are chosen to neutralize the first coordination sphere.
    # In all of the above, the "number of bonds attached" include free valences
    for ic in range(1,natom+1):
        found_charge0 = False
        for jc in range(0,len(charge0),2):
            if ic == charge0[jc]:
                atomcharge[ic] += charge0[jc+1]
                found_charge0 = True
                break
        if found_charge0: # charge is specified for ic
            continue
        if is_bonded_to_transition_metal[ic]:
            continue
        aname = coord[ic][0]
        ele = ELEMENTS[aname].number
        effcoordnum = coordnum[ic]
        for jc in range(0,len(spin0),2):
            if ic == spin0[jc]:
                #effcoordnum -= abs(spin0[jc+1])
                atomspin[ic] += spin0[jc+1]
                break
        if aname in ['N', 'P', 'As', 'Sb', 'Bi']:
            if effcoordnum == 4:
                atomcharge[ic] += 1
            elif effcoordnum == 2:
                atomcharge[ic] -= 1
        elif aname in ['O', 'S', 'Se', 'Te', 'Po']:
            if effcoordnum == 3:
                atomcharge[ic] += 1
            elif effcoordnum == 1:
                atomcharge[ic] -= 1
        elif aname in ['F', 'Cl', 'Br', 'I', 'At']:
            if effcoordnum == 0:
                atomcharge[ic] -= 1
        elif aname in ['B', 'Al', 'Ga', 'In', 'Tl']:
            if effcoordnum == 4:
                atomcharge[ic] -= 1
        elif aname in ['Li', 'Na', 'K', 'Rb', 'Cs', 'Fr']:
            atomcharge[ic] += 1
        elif aname in ['Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra']:
            atomcharge[ic] += 2

    # 2. Implement total charge and total spin constraints
    diffcharge = sum(atomcharge[1:natom+1])-totcharge
    if diffcharge != 0:
        sign = diffcharge/abs(diffcharge)
        for ic in range(1,natom+1):
            if atomcharge[ic]*sign > 0:
                atomcharge[ic] -= sign
                diffcharge -= sign
            if diffcharge == 0:
                break
    diffspin = sum(atomspin[1:natom+1])-totspin
    if diffspin != 0:
        sign = diffspin/abs(diffspin)
        for ic in range(1,natom+1):
            if atomspin[ic]*sign > 0:
                atomspin[ic] -= sign
                diffspin -= sign
            if diffspin == 0:
                break

    # 3. Generate fragment charges
    for ifrag in range(1,nfrag+1):
        bufcen = molefrag[ifrag][0] + molefrag[ifrag][1]
        for ic in bufcen:
            charge[ifrag] += atomcharge[ic]

    for ifrag in range(1,nfrag+1):
        # 4. Deduce the parity of the sum of nuclear charges in each fragment
        nuclear_charge = 0
        for ic in molefrag[ifrag][0]:
            aname = coord[ic][0]
            ele   = ELEMENTS[aname]
            nuclear_charge = nuclear_charge + ele.number 
        for ic in molefrag[ifrag][1]:
            aname = coord[ic][0]
            ele   = ELEMENTS[aname]
            nuclear_charge = nuclear_charge + ele.number 
        for ic in molefrag[ifrag][2]:
            nuclear_charge = nuclear_charge + 1

        # 5. Correct for atomic charges
        nuclear_charge += charge[ifrag]
 
        # 6. Find fragment multiplicities
        # First, sum up all user-given spins
        spin[ifrag] = 0
        bufcen = molefrag[ifrag][0] + molefrag[ifrag][1]
        for iatom in range(0,len(spin0),2):
            if int(spin0[iatom]) in bufcen:
                spin[ifrag] += int(spin0[iatom+1])
        # Then, round spin or charge to give correct parity
        if (nuclear_charge + spin[ifrag]) % 2 == 1:
            if spin[ifrag] == 0:
                # Keep spin as closed shell, if possible
                if charge[ifrag] > 0:
                    charge[ifrag] -= 1
                else:
                    charge[ifrag] += 1
            elif spin[ifrag] > 0:
                spin[ifrag] -= 1
            else:
                spin[ifrag] += 1
        # Transform number of spins to spin multiplicity
        if spin[ifrag] >= 0:
            spin[ifrag] += 1
        else:
            spin[ifrag] -= 1

    print( "--------------------------------------------------------------------")
    print( "Automatically assigned charges and spin multiplicities of fragments:")
    print( "--------------------------------------------------------------------")
    print( "   Fragment  Total No. of atoms  Charge  SpinMult  SpinAlignment")
    for ifrag in range(1,nfrag+1):
        natomfrag = len(molefrag[ifrag][0]) + len(molefrag[ifrag][1]) + len(molefrag[ifrag][2])
        if spin[ifrag]==1:
           align = 'N.A.'
        elif spin[ifrag]>1:
           align = 'Alpha'
        else:
           align = 'Beta'
        print("   %8d  %18d  %6d  %8d  %13s" % (ifrag, natomfrag, charge[ifrag], abs(spin[ifrag]), align))
    print("--------------------------------------------------------------------")

    return [charge, spin]
#end of assign_fragment_chargespin

def main():
    if len(sys.argv) != 3:
        print("Usage: mfcc.py a.inp b.pdb")
        return

    import time
    start_time=time.time()

    inp_file = sys.argv[1]
    pdb_file = sys.argv[2]

    pdb_info=open(pdb_file,"r")
    outinfo=read_PDB(pdb_info)
    pdb_info.close()
    coord        = outinfo["coord"]
    peptideChain = outinfo["peptideChain"]
    radbuff      = outinfo["radbuff"]   # radius of buffers
    
    inp_open=open(inp_file,"r",encoding='UTF-8')
    inp_info=inp_open.readlines()
    inp_open.close()

    #1，初始化分子图
    molecule = amino_acid.class_aminoAcid("Molecule", coord,peptideChain)
    

    #2，切割分子，返回连通分量列表                  #把水的list给扔了，扔了配体
    # molecomp=molecule.cut_peptides()
    #定义第二种，切割化学键
    cutbond=default_cutbond()
    molecomp = molecule.split_molecule_to_molecomps(cutbond)
    print(molecomp)


    #切割之后要重新初始化一下图，使连接信息恢复
    molecule.init_graph_method(molecule.natom,molecule.moleconnect)


    #3，计算molecomp,片段之间最短的欧氏距离  也就是片段之间原子的最短距离
        #每个片段，记录比自己序号小的片段
    molecomp_mindist = molecule.molecomp_distance(molecomp)

    #4，计算得到每个片段的 中心坐标  连接情况，只用连和没连
    molecompcoord = fragment.molecomp_centercoord_connectivity(molecule, molecomp)

    #5，利用coord信息初始化新的抽象图   其实生成图就一个coord，坐标连接信息就够了
    molecomp_graph = fragment.class_molecule("Molecomp",molecompcoord)
    molecomp_graph.extend_connectivity()
    # aaa=molecomp_graph.connected_components_deepthfirstscan()
    # print("@"*100)
    # print(aaa)

    #6，定义阈值、片段包含的原子数大小，生成聚类中心，并聚类
    centriods = None
    from fragment import auto_guess_init_clusters_new
    radcent=3.0    #聚类时分子片段的中心距离的阈值，小于则聚
    centriods = auto_guess_init_clusters_new(molecomp_graph, molecomp, molecomp_mindist, radcent)   #根据一些规则初猜出聚类中心
    #聚类中心的个数
    nfragments =len(centriods)
    print(nfragments)
    center_molefrag = fragment.build_none_overlap_molefrag(molecule, molecomp_graph, molecomp, True, nfragments, centriods)
    print("NONE buffered fragment")
    print(center_molefrag)

    #7，利用上述中心片段产生有重叠的片段
    connectonly = False  # only add direct connect molecomp as buffer
    controlpar  = []
    controlpar.append(radbuff)
    controlpar.append(connectonly)
    #应该给每个中心片段加buffer应该加的是很多个原子节点，不只是氢
    molefrag = fragment.build_buffered_molefrag(molecule,molecomp,molecomp_mindist,center_molefrag, molecomp_graph.edges,controlpar,False)
    print("buffered fragment")
    print(molefrag)
    #[中心原子][缓冲原子][环境原子][不使用的原子]半径，最小环境距离

    #8，给片段计算电子和自选多多重度
    #eMol用的是gjf文件，里面包含了全电子和全自选多重度，定义了totcharge和totspin
    #同时记录了atomcharge和atomspin
    # fragchargespin = assign_fragment_chargespin(coord, molefrag, outinfo["atomcharge"], \
    #      outinfo["atomspin"], outinfo["totcharge"], outinfo["totspin"], molecule.connect)
    # print(fragchargespin)




    #综合成计算模块

    #end_time=time.time()
    #print("Calculate end,time={}".format(end_time-start_time))

    
if __name__ == "__main__":
    main()