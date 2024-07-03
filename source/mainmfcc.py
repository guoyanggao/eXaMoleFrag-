#!/usr/bin/python
import amino_acid
import sys
import fragment
import copy
#最初的eMol是将gjf文件识别，在转化成多个bdf输入文件，我是将pdb文件识别，但是这样没有设置电荷和自旋多重度
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
                coord[icoord].append(str1.split()[5])
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
        outinfo["peptideChain"]= peptideChain
        outinfo["radbuff"]     = radbuff
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

def main():
    if len(sys.argv) != 3:
        print ("")
        print ("Usage: mfcc.py a.inp b.pdb")
        print ("")
        print ("example:python ./source/main.py .inp  ./test/XXX.pdb")
        print ("")
        print ("Possible command line options:")
        print ("")
        print ("")        
        print ("    -MFCC           Generate the fragments required for the MFCC method")
        print ("")
        print ("    -GMFCC          Generate the fragments required for the GMFCC method")
        print ("")
        print ("    -FLMO           Generate the fragments required for the FLMO method")   
        print ("")
        print ("    -FMO            Generate the fragments required for the FMO method")
        print ("")
        print ("    -x              Generate the fragments required for the MFCC method") 
        print ("                    suitable for BDF software.And The BDF software was ")
        print ("                    called to complete the EE-GMFCC calculation.")
        print ("                    You should enter one more a BDF input file as a template")  
        print ("                    to specify the calculation method and base group and so on. ")   
        print ("                    The usage is as follow:")   
        print ("                    # emf -x  xxx.pdb  template.inp")   
        print ("")   
        exit()

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
    residues_list=molecule.cut_peptide_bond_PDBfile()
    caps_list=molecule.cut_to_make_cap()
    
    molecule.init_graph_method(molecule.natom,molecule.moleconnect)
    #切之后要重新初始化一下图
    residues_coord = fragment.molecomp_centercoord_connectivity(molecule, residues_list)
    residues_graph = fragment.class_molecule("Molecomp",residues_coord)

    #generate caped fragment(dangling bonds not yet hydrogenated)

    # noH_frag_list={}
    # noH_frag_list=copy.deepcopy(residues_list)
    # residues_graph.extend_connectivity()
    # for res_i in residues_graph.connect:
    #     for j in range(1,len(residues_graph.connect[res_i]),2):     #i   and j 是两个frag,要从中找出相连的原子
    #         res_j=residues_graph.connect[res_i][j]
    #         for atom_i in residues_list[res_i]:
    #             for atom_j in residues_list[res_j]:
    #                 if molecule.is_connected(atom_i,atom_j):
    #                     for l in caps_list:
    #                         if atom_j in caps_list[l]:
    #                             noH_frag_list[res_i]=noH_frag_list[res_i]+caps_list[l]
    #                             break
    #                     break
    
    # #generate concaps (dangling bonds not yet hydrogenated)    
    # noH_concaps_list={}
    # for i in range(1,len(noH_frag_list)):
    #     set1=set(noH_frag_list[i])
    #     set2=set(noH_frag_list[i+1])
    #     overlap=list(set1 & set2)
    #     if overlap!=[]:
    #         for j in noH_frag_list[i]:
    #             if j>len(molecule.connect):
    #                 overlap.append(j)
    #         for j in noH_frag_list[i+1]:
    #             if j>len(molecule.connect):
    #                 overlap.append(j)
    #         noH_concaps_list[i]=overlap

    noH_frag_list={}
    for i in range(2,len(residues_list)):
        noH_frag_list[i-1]=residues_list[i-1]+residues_list[i]+residues_list[i+1]

    noH_concaps_list={}
    for i in range(2,len(residues_list)-1):
        noH_concaps_list[i-1]=residues_list[i]+residues_list[i+1]
    
    #generate Gconcap
    disthreshold=4.0            
    print("In GMFCC,if the minimum distance between the two residues is less than a custom threshold,")
    print("these two consitute the generalized concap, which abbreviated as 'Gconcap'.\n")
    print("The threshold you define is {}".format(disthreshold))
    print("Search for Gconcaps...")
    #We chose to store 'Gcaps' as a set, then calculate and store the energy 
    #of each Gcap and extract it when needed, so as not to double count it
    Gcaps=set()
    ngcaps=0
    molecompdis=molecule.molecomp_distance(residues_list)
    for i in residues_list:
        for j in range(i-3,0,-1):
            if molecompdis[i][j]<=disthreshold :
                print(" residues {:<2} and {:<2}  ".format(i,j))
                noH_frag_list["G:"+str(i)+"+"+str(j)]=residues_list[i]+residues_list[j]
                ngcaps+=1
                Gcaps.add(i)
                Gcaps.add(j)
    for i in Gcaps:
        noH_concaps_list["G:"+str(i)]=residues_list[i]
    print("There are a total of {} Gconcaps\n".format(ngcaps))

    #add hydrogen for caped fragment             
    fragment_list = copy.deepcopy(noH_frag_list)
    molecule.extend_connectivity()
    for i in noH_frag_list:
        for j in noH_frag_list[i]:
            if molecule.dangling_bond[j]!=0:
                for ii in range(1,len(molecule.connect[j]),2):
                    l=molecule.connect[j][ii]  
                    if l not in noH_frag_list[i]:
                        fragment_list[i].append(molecule.add_node(j,l))
                        
    #add hydrogen for concaps
    concaps_list=copy.deepcopy(noH_concaps_list) 
    for i in noH_concaps_list:
        for j in noH_concaps_list[i]:
            if molecule.dangling_bond[j]!=0:
                for ii in range(1,len(molecule.connect[j]),2):
                    l=molecule.connect[j][ii]  
                    if l not in noH_concaps_list[i]:
                        concaps_list[i].append(molecule.add_node(j,l))
    molecule.compress_connectivity()

    print("-"*10+"Fragment sequences of GMFCC"+"-"*10)
    for i in fragment_list:
        if not type(i)==str:
            print("  caped fragment {:>3} :{}".format(i,fragment_list[i]))
    for i in concaps_list:
        if not type(i)==str:
            print("      concap     {:>3} :{}".format(i,concaps_list[i]))
    print("All fragment sequences added hydrogen at the dangling bonds")

    have_Gconcaps=False
    for i in fragment_list:
        if type(i)==str:
            print("-"*10+"Gconcaps sequences"+"-"*10)
            have_Gconcaps=True
            break
    if have_Gconcaps:
        for i in fragment_list:
            if type(i)==str:
                print("  Gconcap {:<5} :{}".format(i.split(":")[1],fragment_list[i]))
        for i in concaps_list:
            if type(i)==str:
                print("    Gcap  {:<5} :{}".format(i.split(":")[1],concaps_list[i]))
        print("All Gconcaps sequences added hydrogen at the dangling bonds")
    else:
        print("There are no residues with a minimum distance threshold of less than {}, so there is no Gconcap".format(disthreshold))            
    # Rcharge={1: -0.2234, 2: 0.095, 3: 0.3944, 4: -0.5205, 5: -0.2646, 6: 0.0894, 7: -0.2595, 8: -0.2432, 9: -0.1522, 10: -0.1643, 11: -0.1979, 12: 0.2044, 13: 0.2551, 14: 0.1919, 15: 0.0284, 16: 0.1224, 17: 0.1224, 18: 0.1138, 19: 0.1186, 20: 0.1337, 21: 0.1276, 22: 0.1323, 23: -0.6504, 24: 0.1874, 25: 0.6711, 26: -0.4518, 27: -0.6432, 28: 0.7498, 29: -0.4964, 30: -0.9038, 31: 0.3479, 32: 0.0781, 33: 0.1847, 34: 0.1847, 35: 0.3877, 36: 0.417, 37: -0.7541, 38: 0.2492, 39: 0.5264, 40: -0.4905, 41: -0.3842, 42: -0.16, 43: -0.1479, 44: -0.2717, 45: 0.3596, 46: 0.0798, 47: 0.1652, 48: 0.1652, 49: 0.1099, 50: 0.1099, 51: 0.1218, 52: 0.1218, 53: 0.1218, 54: -0.6027, 55: 0.0, 56: 0.542, 57: -0.4465, 58: -0.0951, 59: -0.4528, 60: 0.8056, 61: -0.5592, 62: -0.9316, 63: 0.3119, 64: 0.1496, 65: 0.0768, 66: 0.0768, 67: 0.1269, 68: 0.1269, 69: 0.4002, 70: 0.4078, 71: -0.0252}
    # Rcharge={1: -0.365, 2: -0.0915, 3: 0.6582, 4: -0.583, 5: 0.2179, 6: -0.6154, 7: 0.1691, 8: 0.261, 9: 0.1754, 10: 0.0977, 11: -0.001, 12: -0.001, 13: 0.4102, 14: -0.892, 15: 0.4102, 16: 0.7159, 17: -0.5847, 18: -0.5068, 19: 0.081, 20: -0.2335, 21: 0.006, 22: -0.3403, 23: 0.179, 24: -0.1745, 25: -0.2958, 26: -0.206, 27: -0.1201, 28: 0.3956, 29: 0.0337, 30: 0.1567, 31: 0.1567, 32: 0.2017, 33: 0.3502, 34: 0.1478, 35: 0.166, 36: 0.1575, 37: 0.1462, 38: -0.7024, 39: -0.0767, 40: 0.628, 41: -0.5232, 42: 0.4555, 43: -0.6514, 44: -0.5766, 45: 0.3297, 46: 0.1097, 47: 0.0404, 48: 0.4358, 49: 0.147, 50: 0.147, 51: 0.147, 52: -0.6471, 53: 0.1841, 54: 0.6878, 55: -0.5352, 56: -0.3429, 57: -0.0523, 58: -0.1614, 59: 0.0154, 60: -0.4119, 61: 0.2562, 62: -0.2271, 63: -0.3714, 64: -0.2056, 65: -0.0989, 66: 0.3342, 67: 0.1131, 68: 0.1363, 69: 0.1363, 70: 0.2004, 71: 0.3651, 72: 0.1582, 73: 0.1762, 74: 0.164, 75: 0.1493, 76: -0.8277, 77: 0.0926, 78: 0.6814, 79: -0.5081, 80: -0.0425, 81: -0.6009, 82: 0.867, 83: -0.6807, 84: -0.5343, 85: 0.3612, 86: 0.1461, 87: 0.0933, 88: 0.0933, 89: 0.1499, 90: 0.1499, 91: -0.6813, 92: 0.1994, 93: 0.5949, 94: -0.5029, 95: -0.5312, 96: 0.8502, 97: -0.5629, 98: -0.9576, 99: 0.3271, 100: 0.0533, 101: 0.1497, 102: 0.1497, 103: 0.4233, 104: 0.4186, 105: -0.5593, 106: -0.1057, 107: 0.5094, 108: -0.4652, 109: 0.3018, 110: 0.1134, 111: 0.1134, 112: 0.5104, 113: -0.0043}
    # Rcharge={1: -0.7498, 2: 0.1904, 3: 0.6094, 4: -0.4959, 5: 0.0727, 6: -0.6308, 7: 0.3682, 8: 0.4257, 9: 0.3271, 10: 0.0576, 11: 0.0431, 12: 0.0431, 13: 0.4383, 14: -0.8458, 15: 0.4757, 16: 0.5583, 17: -0.5648, 18: -0.1894, 19: -0.0718, 20: -0.199, 21: 0.0463, 22: -0.358, 23: 0.195, 24: -0.1437, 25: -0.3318, 26: -0.3043, 27: -0.0937, 28: 0.3776, 29: 0.0332, 30: 0.065, 31: 0.065, 32: 0.1873, 33: 0.3566, 34: 0.0617, 35: 0.1664, 36: 0.191, 37: 0.1375, 38: -0.5185, 39: -0.3201, 40: 0.6465, 41: -0.4478, 42: 0.8021, 43: -0.7286, 44: -0.6935, 45: 0.2269, 46: 0.1089, 47: -0.0944, 48: 0.4562, 49: 0.1648, 50: 0.1648, 51: 0.1648, 52: -0.697, 53: 0.1239, 54: 0.7073, 55: -0.5617, 56: -0.2837, 57: -0.0678, 58: -0.1623, 59: 0.0869, 60: -0.4002, 61: 0.2417, 62: -0.2816, 63: -0.3808, 64: -0.2012, 65: -0.0724, 66: 0.3799, 67: 0.056, 68: 0.1264, 69: 0.1264, 70: 0.2089, 71: 0.3629, 72: 0.2265, 73: 0.1783, 74: 0.1458, 75: 0.1367, 76: -0.7261, 77: 0.0242, 78: 0.7184, 79: -0.472, 80: -0.0834, 81: -0.3635, 82: 0.7408, 83: -0.7315, 84: -0.7229, 85: 0.3721, 86: 0.0979, 87: 0.089, 88: 0.089, 89: 0.0654, 90: 0.0654, 91: -0.6576, 92: 0.1829, 93: 0.5653, 94: -0.4952, 95: -0.5296, 96: 0.8197, 97: -0.557, 98: -0.9401, 99: 0.3275, 100: 0.0408, 101: 0.1566, 102: 0.1566, 103: 0.4192, 104: 0.4177, 105: -0.4716, 106: -0.224, 107: 0.5582, 108: -0.5777, 109: 0.2864, 110: 0.1396, 111: 0.1396, 112: -0.1317, 113: -0.2736, 114: 0.7776, 115: -0.5421, 116: -0.1128, 117: 0.0392, 118: -0.1116, 119: 0.0039, 120: -0.2387, 121: 0.1283, 122: 0.123, 123: 0.0896, 124: 0.0896, 125: 0.0065, 126: 0.0065, 127: 0.0632, 128: 0.0632, 129: 0.067, 130: 0.067, 131: 0.1725, 132: 0.1899, 133: 0.1545, 134: -0.8074, 135: 0.2108, 136: 0.707, 137: -0.5966, 138: -0.3825, 139: -0.0161, 140: -0.2394, 141: 0.094, 142: -0.3036, 143: 0.1144, 144: -0.1665, 145: -0.2863, 146: -0.2303, 147: -0.1379, 148: 0.437, 149: 0.1345, 150: 0.1458, 151: 0.1458, 152: 0.2007, 153: 0.3388, 154: 0.1227, 155: 0.1637, 156: 0.1733, 157: 0.1426, 158: -0.611, 159: -0.0463, 160: 0.5383, 161: -0.3582, 162: 0.382, 163: -0.7086, 164: -0.538, 165: 0.23, 166: 0.1217, 167: 0.0699, 168: 0.4394, 169: 0.1452, 170: 0.1452, 171: 0.1452, 172: -0.5111, 173: 0.068, 174: 0.6741, 175: -0.5624, 176: -0.3096, 177: -0.1181, 178: -0.1982, 179: 0.0602, 180: -0.3704, 181: 0.2411, 182: -0.2826, 183: -0.4493, 184: -0.2629, 185: -0.0994, 186: 0.3361, 187: 0.0476, 188: 0.1468, 189: 0.1468, 190: 0.1886, 191: 0.3405, 192: 0.2138, 193: 0.1666, 194: 0.128, 195: 0.1129, 196: -0.5554, 197: 0.1083, 198: 0.747, 199: -0.6059, 200: -0.1795, 201: -0.0934, 202: -0.3239, 203: 0.0295, 204: -0.249, 205: 0.2273, 206: 0.0189, 207: 0.0659, 208: 0.0659, 209: 0.1146, 210: 0.1146, 211: 0.1386, 212: 0.1386, 213: 0.0771, 214: 0.0771, 215: 0.2347, 216: 0.1692, 217: 0.1658, 218: -0.9055, 219: 0.3868, 220: 0.4111}
    # Rcharge={1: -0.7129, 2: 0.3094, 3: 0.1203, 4: -0.2971, 5: 0.2903, 6: -0.6629, 7: 0.3617, 8: 0.0975, 9: -0.0004, 10: -0.0004, 11: 0.448, 12: -0.2279, 13: -0.2379, 14: 0.1048, 15: -0.28, 16: -0.1813, 17: 0.0136, 18: 0.6615, 19: -0.5062, 20: -0.8736, 21: 0.3064, 22: 0.2478, 23: 0.1154, 24: 0.1154, 25: 0.0071, 26: 0.0071, 27: 0.3805, 28: 0.4306, 29: 0.0331, 30: -0.0776, 31: 0.6705, 32: -0.4577, 33: -0.3865, 34: 0.7775, 35: -0.4565, 36: -0.4717, 37: 0.1531, 38: 0.1096, 39: 0.0973, 40: 0.0973, 41: -0.5176, 42: 0.1668, 43: 0.2428, 44: -0.3634, 45: -0.2602, 46: 0.0546, 47: -0.1511, 48: -0.1454, 49: -0.3378, 50: -0.3587, 51: 0.4482, 52: -0.6153, 53: 0.3047, 54: 0.1013, 55: 0.0927, 56: 0.0927, 57: 0.1886, 58: 0.153, 59: 0.1872, 60: 0.1873, 61: 0.4292, 62: -0.0176, 63: 0.132, 64: 0.4033, 65: -0.4234, 66: -0.542, 67: 0.3553, 68: -0.4533, 69: -0.3787, 70: -0.0565, 71: 0.101, 72: 0.1515, 73: 0.1515, 74: -0.0193, 75: 0.121, 76: 0.121, 77: 0.121, 78: 0.0813, 79: 0.0813, 80: 0.0813, 81: -0.1266, 82: -0.0961, 83: 0.3131, 84: -0.4617, 85: 0.3668, 86: -0.6746, 87: 0.0367, 88: 0.1275, 89: -0.0098, 90: -0.0098, 91: 0.4211, 92: -0.0133, 93: -0.1442, 94: 0.5567, 95: -0.5287, 96: -0.3329, 97: 0.8143, 98: -0.612, 99: -0.6743, 100: -0.0582, 101: 0.1549, 102: 0.1054, 103: 0.1054, 104: -0.2034, 105: -0.3134, 106: 0.4296, 107: -0.3251, 108: -0.0994, 109: 0.1411, 110: -0.3131, 111: 0.1736, 112: 0.0271, 113: -0.3769, 114: 0.2045, 115: 0.1545, 116: 0.0887, 117: 0.0887, 118: 0.3444, 119: -0.1891, 120: 0.144, 121: -0.5073, 122: 0.377, 123: 0.0632, 124: -0.2442, 125: -0.4854, 126: 0.3912, 127: -0.3284, 128: -0.3218, 129: 0.4155, 130: 0.0643, 131: 0.1245, 132: 0.1245, 133: -0.0395, 134: 0.0761, 135: 0.0761, 136: 0.0761, 137: 0.0698, 138: 0.0698, 139: 0.0698, 140: -0.0278, 141: -0.1865, 142: 0.4747, 143: -0.3708, 144: -0.1995, 145: -0.1119, 146: -0.1226, 147: 0.0222, 148: -0.3817, 149: 0.213, 150: -0.1434, 151: -0.3088, 152: -0.2213, 153: -0.0989, 154: 0.1025, 155: 0.14, 156: 0.1255, 157: 0.1255, 158: 0.1895, 159: 0.3599, 160: 0.0976, 161: 0.1533, 162: 0.1548, 163: 0.1348, 164: 0.0695, 165: -0.3491, 166: 0.4562, 167: -0.4338, 168: 0.1317, 169: -0.0731, 170: 0.1998, 171: -0.6745, 172: 0.7714, 173: -0.7448, 174: -0.7786, 175: 0.1823, 176: 0.0868, 177: 0.0078, 178: 0.0078, 179: 0.0509, 180: 0.0509, 181: 0.0072, 182: 0.0072, 183: 0.4038, 184: 0.3934, 185: 0.3594, 186: 0.3785, 187: 0.4159, 188: 0.0104, 189: 0.0616, 190: 0.517, 191: -0.5235, 192: -0.3373, 193: -0.2104, 194: 0.0693, 195: 0.105, 196: 0.105, 197: 0.105, 198: -0.1554, 199: 0.1389, 200: 0.359, 201: -0.492, 202: -0.2231, 203: 0.5131, 204: -0.4292, 205: -0.5326, 206: -0.1716, 207: 0.0715, 208: 0.0564, 209: 0.0564, 210: -0.1528, 211: 0.0913, 212: 0.0913, 213: 0.0913, 214: 0.1148, 215: 0.1148, 216: 0.1148, 217: -0.1225, 218: -0.231, 219: 0.66, 220: -0.5438, 221: -0.2306, 222: 0.6217, 223: -0.5569, 224: -0.6305, 225: 0.2875, 226: 0.1823, 227: 0.0896, 228: 0.0896, 229: 0.2789, 230: 0.3347, 231: -0.2923, 232: 0.0371, 233: 0.681, 234: -0.3198, 235: -0.5746, 236: -0.4212, 237: 0.101, 238: 0.1251, 239: 0.1624, 240: 0.1624, 241: 0.1624}
    # Rcharge={1: -0.4552, 2: 0.0903, 3: 0.432, 4: -0.5104, 5: 0.3185, 6: 0.2909, 7: 0.3277, 8: 0.0426, 9: 0.0426, 10: -0.2268, 11: -0.3573, 12: 0.6743, 13: -0.5351, 14: 0.5356, 15: -0.5116, 16: -0.5288, 17: 0.1178, 18: 0.1272, 19: -0.0485, 20: 0.1234, 21: 0.1234, 22: 0.1234, 23: 0.1053, 24: 0.1053, 25: 0.1053, 26: -0.731, 27: -0.4428, 28: 0.7203, 29: -0.5378, 30: -0.1192, 31: 0.0066, 32: 0.1052, 33: -0.177, 34: -0.2477, 35: 0.4327, 36: 0.2632, 37: 0.0923, 38: 0.0923, 39: -0.0017, 40: -0.0017, 41: 0.0377, 42: 0.0377, 43: 0.1078, 44: 0.1078, 45: 0.2312, 46: 0.2953, 47: 0.2853, 48: -0.4394, 49: -0.1949, 50: 0.553, 51: -0.518, 52: 0.0829, 53: -0.5203, 54: 0.3145, 55: 0.1444, 56: 0.071, 57: 0.071, 58: 0.4047, 59: -0.3625, 60: 0.0968, 61: 0.5031, 62: -0.5046, 63: -0.1161, 64: -0.5557, 65: 0.3896, 66: 0.0528, 67: 0.1179, 68: 0.1179, 69: 0.4398, 70: -0.5266, 71: -0.1648, 72: 0.5096, 73: -0.6524, 74: 0.114, 75: -0.3587, 76: 0.7853, 77: -0.6147, 78: -0.6762, 79: 0.1963, 80: 0.1566, 81: 0.0064, 82: 0.0064, 83: 0.0818, 84: 0.0818, 85: -0.0667, 86: -0.3437, 87: 0.7987, 88: -0.5561, 89: 0.4063, 90: -0.5705, 91: -0.5716, 92: -0.1236, 93: 0.1578, 94: 0.0507, 95: 0.4082, 96: 0.141, 97: 0.141, 98: 0.141, 99: -0.6466, 100: -0.1535, 101: 0.1929, 102: -0.3279, 103: 0.3142, 104: -0.5632, 105: -0.3461, 106: 0.3503, 107: 0.1347, 108: 0.0401, 109: 0.3692, 110: 0.1083, 111: 0.1083, 112: 0.1083, 113: 0.0743, 114: 0.0347, 115: 0.4182, 116: -0.4086, 117: -0.444, 118: 0.2836, 119: -0.3017, 120: -0.4293, 121: 0.2394, 122: 0.0516, 123: 0.1424, 124: 0.1424, 125: 0.0151, 126: 0.0668, 127: 0.0668, 128: 0.0668, 129: 0.1066, 130: 0.1066, 131: 0.1066, 132: -0.1979, 133: -0.252, 134: 0.3691, 135: -0.3987, 136: 0.6411, 137: -0.7336, 138: -0.5431, 139: 0.0578, 140: 0.1822, 141: -0.079, 142: 0.4803, 143: 0.136, 144: 0.136, 145: 0.136, 146: 0.1292, 147: -0.3382, 148: 0.4234, 149: -0.5559, 150: 0.3489, 151: -0.4694, 152: -0.1992, 153: -0.2253, 154: -0.1428, 155: 0.107, 156: 0.0173, 157: 0.0173, 158: 0.2033, 159: 0.2033, 160: 0.1142, 161: 0.1142, 162: 0.1142, 163: -0.1462, 164: 0.2, 165: 0.8928, 166: -0.8322, 167: -0.3952, 168: 0.1508, 169: -0.1029, 170: -0.1798, 171: -0.1507, 172: -0.1578, 173: -0.1555, 174: 0.0117, 175: 0.0232, 176: 0.1289, 177: 0.1289, 178: 0.0159, 179: 0.1332, 180: 0.1435, 181: 0.1395, 182: 0.1467, 183: 0.041, 184: -0.1411, 185: -0.0523, 186: -0.0208, 187: -0.2018, 188: 0.4824, 189: -0.3328, 190: -0.4595, 191: 0.028, 192: 0.1459, 193: 0.0632, 194: 0.0632, 195: -0.0867, 196: 0.0799, 197: 0.0799, 198: 0.0799, 199: 0.1048, 200: 0.1048, 201: 0.1048, 202: 0.0475, 203: -0.1507, 204: 0.0933, 205: -0.2328, 206: 0.0522, 207: -0.1, 208: -0.187, 209: 0.1211, 210: -0.4913, 211: -0.1264, 212: 0.1414, 213: 0.0578, 214: 0.0578, 215: 0.0651, 216: 0.0651, 217: 0.0599, 218: 0.0599, 219: 0.0877, 220: 0.0877, 221: 0.3322, 222: 0.3277, 223: 0.3328, 224: 0.0043, 225: -0.2497, 226: 1.0449, 227: -0.6334, 228: -0.3645, 229: -0.2686, 230: 0.8517, 231: -0.6794, 232: -0.7142, 233: 0.1017, 234: 0.1131, 235: 0.1753, 236: 0.1753, 237: 0.0856, 238: 0.0856, 239: -0.903, 240: 0.3178, 241: 0.4152, 242: -0.4762, 243: -0.1921, 244: -0.2452, 245: -0.1712, 246: -0.2908, 247: 0.4345, 248: 0.0323, 249: 0.084, 250: 0.084, 251: 0.1816, 252: 0.1816, 253: 0.1235, 254: 0.1235, 255: 0.1235, 256: -0.1354, 257: -0.0139, 258: 0.3412, 259: -0.4066, 260: 0.0401, 261: -0.4009, 262: 0.7533, 263: -0.5437, 264: -0.9065, 265: -0.1545, 266: 0.1031, 267: 0.0593, 268: 0.0593, 269: 0.0928, 270: 0.0928, 271: 0.3876, 272: 0.4189, 273: -0.1834, 274: -0.1717, 275: 0.8396, 276: -0.5694, 277: -0.3949, 278: 0.354, 279: -0.4849, 280: -0.3806, 281: 0.0169, 282: 0.2106, 283: 0.1248, 284: 0.1248, 285: -0.0336, 286: 0.1172, 287: 0.1172, 288: 0.1172, 289: 0.0935, 290: 0.0935, 291: 0.0935, 292: -0.7863, 293: 0.0364, 294: 0.4961, 295: -0.5286, 296: -0.1739, 297: 0.1883, 298: -0.1523, 299: -0.0245, 300: -0.3976, 301: 0.5048, 302: 0.1044, 303: 0.0601, 304: 0.0601, 305: -0.032, 306: -0.032, 307: 0.0648, 308: 0.0648, 309: 0.0962, 310: 0.0962, 311: 0.2991, 312: 0.3155, 313: 0.327, 314: -0.1731, 315: -0.2625, 316: 0.6147, 317: -0.5307, 318: 0.1426, 319: 0.1356, 320: 0.1356, 321: -0.4897, 322: 0.1102, 323: 0.5161, 324: -0.5508, 325: -0.2392, 326: 0.4705, 327: -0.402, 328: -0.4647, 329: 0.1613, 330: 0.0922, 331: 0.0643, 332: 0.0643, 333: -0.0908, 334: 0.0851, 335: 0.0851, 336: 0.0851, 337: 0.0838, 338: 0.0838, 339: 0.0838, 340: -0.1998, 341: -0.0433, 342: 0.8019, 343: -0.6392, 344: -0.0257, 345: -0.2283, 346: -0.1598, 347: -0.6635, 348: 0.0482, 349: 0.0275, 350: 0.0275, 351: 0.0807, 352: 0.0807, 353: 0.1492, 354: 0.1492}
    # molecule.EEMFCC_iter_charge(fragment_list,concaps_list,inp_info,molecompdis,4)
    # molecule.EEMFCC_haverespcharge(fragment_list,concaps_list,inp_info,molecompdis,Rcharge)
    # molecule.MFCC_bigcap(fragment_list,concaps_list,inp_info)
    molecule.GMFCC_bigcap(fragment_list,concaps_list,inp_info)
    end_time=time.time()

    
if __name__ == "__main__":
    main()