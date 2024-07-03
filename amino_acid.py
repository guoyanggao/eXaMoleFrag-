import fragment
import os
import shutil
import numpy as np

#peptideChain  {'A': {1: 'GLY', 2: 'VAL', 3: 'VAL'}, 'B': {1: 'ASN', 2: 'SER', 3: 'LEU'}}
class class_aminoAcid(fragment.class_molecule):
    
    def __init__(self,name,coord,peptideChain):
        super().__init__(self,coord)
        # print("constructor called")
        self.peptideChain=peptideChain          

        # coord_list=[]
        # for i in coord_list:
        #     list1=np.array(list(self.coord[i].values()))
        #     coord_list.append(list1)
        # coord_array  = np.array(coord_list)
        # center_of_coord=np.mean(coord_array, axis=0)
        # print(coord_array)
        # print(center_of_coord)
        # return center_of_coord

    #mainly used for adding hydrogen
    def add_node(self,a,b):   #a是被加氢的原子  
        self.atom_number+=1
        i=self.atom_number   #不能是原子数目，因为对于有一些pdb文件，原子序数不连续
        self.label[i]="H"
        x=self.coord[b][1]-self.coord[a][1]
        y=self.coord[b][2]-self.coord[a][2]
        z=self.coord[b][3]-self.coord[a][3]
        length=(x*x+y*y+z*z)**0.5

        self.coord[i]={}
        self.coord[i][0]="H"
        self.coord[i][1]=round(self.coord[a][1]+x/length*1.09,3)
        self.coord[i][2]=round(self.coord[a][2]+y/length*1.09,3)
        self.coord[i][3]=round(self.coord[a][3]+z/length*1.09,3)
        
        # print("add {} H to {} atom ".format(i,a))
        return i

    #excise peptide and disulfide bonds to get isolated amino acid residues
    def cut_peptide_bond_PDBfile(self):   
        self.extend_connectivity()
        for i in range(1,self.natom+1):
            #There will also be disulfide bonds between CYS, cutting disulfide bonds
            if self.label2[i]=="SG" and self.connect[i][0]==2 and self.residue_type[i]=="CYS" and {self.label[self.connect[i][1]],self.label[self.connect[i][3]]}=={"S","C"} :  #找到S原子
                for j1 in range(1,len(self.connect[i]),2):
                    i2=self.connect[i][j1]
                    if self.label2[i2]=="SG" and self.connect[i2][0]==2 and self.residue_type[i]=="CYS" and {self.label[self.connect[i2][1]],self.label[self.connect[i2][3]]}=={"S","C"}:  #找到S原子
                        self.cut_edge(i,i2)
                        self.dangling_bond[i]+=1
                        self.dangling_bond[i2]+=1
                        print("Cutting SG-SG bond between {:>4} and {:>4}".format(i,i2))
                        break
            elif self.label2[i]=="C" and self.connect[i][0]==3 and {self.label[self.connect[i][1]],self.label[self.connect[i][3]],self.label[self.connect[i][5]]}=={"C","N","O"}:  #找到碳原子
                for j1 in range(1,len(self.connect[i]),2):
                    i2=self.connect[i][j1]
                    if self.label2[i2]=="N" and self.residue_type[i2]!="NH2":
                        self.cut_edge(i,i2)
                        self.dangling_bond[i]+=1
                        self.dangling_bond[i2]+=1
                        print("Cutting C-N peptide bonds between {:>4} and  {:>4}".format(i,i2))
                        break
        self.compress_connectivity()
        self.amino_acids_list=self.connected_components_deepthfirstscan()
        print("\n"+"-"*10+"Residues in protein"+"-"*10)
        for i in self.amino_acids_list:
            print("residue-{}:  {}".format(i,self.amino_acids_list[i]))
        return self.amino_acids_list

    #In the amino acid residue，Cut αC and C of the carbonyl group
    #The two parts formed serve as the left and right caps, respectively.
    def cut_to_make_cap(self):
        # self.cut_peptide_bond_PDBfile()
        self.extend_connectivity()
        for i in range(1,self.natom+1):
            if self.label2[i]=="CA" and self.connect[i][0]==4 :
                for j1 in range(1,len(self.connect[i]),2):
                    i2=self.connect[i][j1]
                    if self.label2[i2]=="C":
                        self.cut_edge(i,i2)
                        self.dangling_bond[i]+=1
                        self.dangling_bond[i2]+=1
                        print("Cutting CA-C bond between {:>4} and {:>4} ".format(i,i2))
                        break
            # Actually,in this cap creation method,the ring of PRO need not cut.
            elif self.residue_type[i]=="PRO" and self.label[i]=="N" and self.connect[i][0]==3:
                for j1 in range(1,len(self.connect[i]),2):
                    i2=self.connect[i][j1]
                    if self.label2[i2]=="CD" and self.connect[i2][0]==4:
                        self.cut_edge(i,i2)
                        self.dangling_bond[i]+=1
                        self.dangling_bond[i2]+=1
                        print("Cutting N-CD bond in PRO between {:>4} and {:>4}".format(i,i2))
                        break
        self.compress_connectivity()
        amino_acids_list=self.connected_components_deepthfirstscan()
        print(amino_acids_list)
        return amino_acids_list

    #前段时间改的，有分枝，处理pdb和非pdb文件     <不能正常工作>
    def cut_peptides(self):
        #gaussi保存文件不记录双键，把氧原子变成双键
        # for i in range(1,self.natom+1):
        #     if self.label[i]=="O" and self.connect[i][0]==1 and self.label[self.connect[i][1]]=="C":
        #         self.connect[i][2]=2.0
        print("输入文件为标准pdb格式文件")
        self.extend_connectivity()
        if self.residue_type!={} and self.residue_type!={}:
            for i in range(1,self.natom+1):
                #切割二硫键，注意切割前可能要加一个compress
                if self.label2[i]=="SG" and self.residue_type[i]=="CYS" and self.connect[i][0]==2:  
                    for ii in range(1,len(self.connect[i]),2):
                        j=self.connect[i][ii]
                        if self.label2[j]=="SG":
                            print("二硫键{},{}".format(i,j))
                            self.compress_connectivity()
                            self.cut_edge(i,j)
                            self.extend_connectivity()
                            break
                #切割PRO环
                elif self.residue_type[i]=="PRO" and self.label2[i]=="N":
                    for ii in range(1,len(self.connect[i]),2):
                        j=self.connect[i][ii]
                        if self.residue_type[j]=="PRO" and self.label2[j]!="CA":
                            self.cut_edge(i,j)
                            break
                #切割肽键                    
                elif self.label2[i]=="CA":
                    i+=1
                    for ii in range(1,len(self.connect[i]),2):
                        j=self.connect[i][ii]
                        if self.label[j]=="N": 
                            self.cut_edge(i,j)
                            break
        else:  #还要优化，目前无法正常使用
            for i in range(1,self.natom+1):
                if self.label[i]=="S" and self.connect[i][0]==2:  #找到S原子
                    print(i)
                    for ii in range(1,len(self.connect[i]),2):
                        j=self.connect[i][ii]
                        if self.label[j]=="S":
                            print("二硫键{},{}".format(i,j))
                            self.cut_edge(i,j)
                            print(self.connect[j])
                            break
            for i in range(1,self.natom+1):
                if self.label[i]=="C" and self.connect[i][0]==3 and {self.label[self.connect[i][1]],self.label[self.connect[i][3]],self.label[self.connect[i][5]]}=={"C","N","O"}:  #找到碳原子
                    for j1 in range(1,len(self.connect[i]),2):
                        i2=self.connect[i][j1]
                        if self.label[i2]=="N" and self.connect[i2][0]==3 and ((self.label[self.connect[i2][1]]=="H" and self.label[self.connect[i2][3]]=="C" and self.label[self.connect[i2][5]]=="C") or (self.label[self.connect[i2][1]]=="C" and self.label[self.connect[i2][3]]=="H" and self.label[self.connect[i2][5]]=="C") or (self.label[self.connect[i2][1]]=="C" and self.label[self.connect[i2][3]]=="C" and self.label[self.connect[i2][5]]=="H")): 
                            self.cut_edge(i,i2)
                            break
            for i in range(1,self.natom+1):
                if self.label[i]=="N" and self.connect[i][0]==3 and {self.label[self.connect[i][1]],self.label[self.connect[i][3]],self.label[self.connect[i][5]]}=={"C"} and (self.connect[self.connect[i][1]][0]==3 or self.connect[self.connect[i][3]][0]==3 or self.connect[self.connect[i][5]][0]==3):  #找 N原子
                    for ii in range(1,len(self.connect[i]),2):
                        j=self.connect[i][ii]
                        if self.label[j]=="C" and self.connect[j][0]==4 and {self.label[self.connect[j][1]],self.label[self.connect[j][3]],self.label[self.connect[j][5]],self.label[self.connect[j][7]]}=={"N","C","H"} and (2.0 not in self.connect[self.connect[j][1]] or 2.0 not in self.connect[self.connect[j][3]] or 2.0 not in self.connect[self.connect[j][5] or 2.0 not in self.connect[self.connect[j][7]]]):
                            print("脯氨酸{}和{}".format(i,j))
                            self.cut_edge(i,j)
                            break

        self.compress_connectivity()
        self.full_ord=self.connected_components_deepthfirstscan()
        return self.full_ord

    def genxyzfile(self,atomlist,filename):
        input="    "+str(len(atomlist))+"\n"
        input+=" Generated by eMaMoleFrag \n"
        for i in atomlist:            
            input+="{}      {:>10}   {:>9}   {:>9}\n".format(self.label[i],str(self.coord[i][1]),str(self.coord[i][2]),str(self.coord[i][3]))
        if ".xyz" not in filename:
            filename+=".xyz"
        with open(filename,"w")as writ:
            writ.write(input)

    #all fragment is generated as .xyz format file
    def generate_fragments_file_of_GMFCC(self,fraglist,caplist):
        print("In the 'GMFCC_fragment_xyz_file' folder in the current directory, the xyz format  \
              file for each fragment required by the GMFCC method has been generated.")
        print("NOTICE:\n")
        print(" Generate fragments .xyz format file ....")

        os.system("mkdir GMFCC_fragment_xyz_file")
        os.chdir("GMFCC_fragment_xyz_file")
        for i in fraglist:
            if not type(i)==str:
                self.genxyzfile(fraglist[i],"caped_fragment_"+str(i))
            else:
                self.genxyzfile(fraglist[i],"Gconcap_"+str(i.split(":")[1]))
        for i in caplist:
            if not type(i)==str:
                self.genxyzfile(caplist[i],"concaps_"+str(i))    
            else:
                self.genxyzfile(caplist[i],"Gcap_"+str(i.split(":")[1]))    
        # print("-"*10+"End eXaMoleFrag"+"-"*10)
        print("fragment .xyz format file generation done, total elapsed time 0.02 s")
  
    
    #NOTICE:The following four functions call the BDF program to test the correctness and usefulness of the fragments produced by eXaMoleFrag program.
    #version 3,date 5.21  生成bdf的inp文件  #该函数需要输入inp文件的模板   &  可以再写一个默认的设置参数，不需要inp文件模板
    def genbdfeasyinput(self,atomlist,filename,module_file):
        input=""
        for line in module_file:
            if line!="\n":
                input+=line
                if line.startswith("Geometry"):
                    for i in atomlist:            
                        input+="{}  {:>7}  {:>6}  {:>6}\n".format(self.label[i],str(self.coord[i][1]),str(self.coord[i][2]),str(self.coord[i][3]))
        if ".inp" not in filename:
            filename+=".inp"
        with open(filename,"w")as writ:
                writ.write(input)
    
    def genbdfextchargefile(self,chargelist,charfilena):
        input="Title \n"+str(len(chargelist))+"\n"
        for i in chargelist:
            input+=str(self.label[i])+" "+str(self.respcharge[i])+" "+str(self.coord[i][1])+" "+str(self.coord[i][2])+" "+str(self.coord[i][3])+"\n"
        with open(charfilena,"w")as writ:
            writ.write(input)
    

    #普通的MFCC 
    def MFCC_bigcap(self,fraglist,caplist,inp_info):
        os.system("mkdir EE-GMFCC_calculFile")
        os.chdir("EE-GMFCC_calculFile")
        E_tot=0
        for i in caplist: 
            if not type(i)==str:
                dirname="concap-"+str(i)
                os.system("mkdir "+dirname)
                os.chdir(dirname)
                filename=dirname+".inp"
                self.genbdfeasyinput(caplist[i],filename,inp_info)
                os.system("/home/gao/software/bdf-pkg/build-gnu13/bdf-pkg-full/sbin/run.sh "+filename)
                with open(filename+'.out', 'r')as file:
                    Ea=float(next(l for l in file if 'E_tot =' in l).split()[2])
                E_tot-=Ea   
                os.chdir("..")
        for i in fraglist: 
            if not type(i)==str:
                dirname="frag-"+str(i)
                os.system("mkdir "+dirname)
                os.chdir(dirname)
                filename=dirname+".inp"
                self.genbdfeasyinput(fraglist[i],filename,inp_info)
                os.system("/home/gao/software/bdf-pkg/build-gnu13/bdf-pkg-full/sbin/run.sh "+filename)
                with open(filename+'.out', 'r')as file:
                    Ea=float(next(l for l in file if 'E_tot =' in l).split()[2])
                E_tot+=Ea
                os.chdir("..")
        E_ele=0
        print("No.1 iteration energy={}".format(E_tot-E_ele))
        
        # print("最终能量-   {}".format(E1-E2-E_ele))     

    #GMFCC
    def GMFCC_bigcap(self,fraglist,caplist,inp_info):

        os.system("mkdir EE-GMFCC_calculFile")
        os.chdir("EE-GMFCC_calculFile")
        E_tot=0
        E_gcap={}
        for i in caplist: 
            dirname="concap-"+str(i)
            os.system("mkdir "+dirname)
            os.chdir(dirname)
            filename=dirname+".inp"
            self.genbdfeasyinput(caplist[i],filename,inp_info)
            os.system("/home/gao/software/bdf-pkg/build-gnu13/bdf-pkg-full/sbin/run.sh "+filename)
            with open(filename+'.out', 'r')as file:
                Ea=float(next(l for l in file if 'E_tot =' in l).split()[2])
                if not type(i)==str:
                    E_tot-=Ea   
                else:
                    E_gcap[int(i.split(":")[1])]=Ea
            os.chdir("..")
        for i in fraglist: 
            dirname="frag-"+str(i)
            os.system("mkdir "+dirname)
            os.chdir(dirname)
            chargefilename=dirname+".extcharge"
            filename=dirname+".inp"
            self.genbdfeasyinput(fraglist[i],filename,inp_info)
            os.system("/home/gao/software/bdf-pkg/build-gnu13/bdf-pkg-full/sbin/run.sh "+filename)
            with open(filename+'.out', 'r')as file:
                Ea=float(next(l for l in file if 'E_tot =' in l).split()[2])
                E_tot+=Ea
                if type(i)==str:
                    E_tot-=E_gcap[int(i.split("+")[1])]
                    E_tot-=E_gcap[int(i.split("+")[0].split(":")[1])]
            os.chdir("..")
        
        print("No.1 iteration energy={}".format(E_tot))

    #直接加左右两个残基为cap 
    def EEGMFCC_iter_charge_big(self,fraglist,caplist,inp_info,molecompdis,niter=1):
        self.respcharge={}
        for i in self.init_graph[1]:
            self.respcharge[i]=0.0

        os.system("mkdir EE-GMFCC_calculFile")
        os.chdir("EE-GMFCC_calculFile")
        E_tot=0
        E_gcap={}
        for i in caplist: 
            dirname="concap-"+str(i)
            os.system("mkdir "+dirname)
            os.chdir(dirname)
            chargefilename=dirname+".extcharge"
            filename=dirname+".inp"
            self.genbdfeasyinput(caplist[i],filename,inp_info)
            chargelist=set(self.init_graph[1])-set(caplist[i])
            self.genbdfextchargefile(chargelist,chargefilename)
            os.system("/home/gao/software/bdf-pkg/build-gnu13/bdf-pkg-full/sbin/run.sh "+filename)
            with open(filename+'.out', 'r')as file:
                Ea=float(next(l for l in file if 'E_tot =' in l).split()[2])
                if not type(i)==str:
                    E_tot-=Ea   
                else:
                    E_gcap[int(i.split(":")[1])]=Ea
            os.system("export OMP_NUM_THREADS=1")                    
            os.system("/home/gao/software/bdf-pkg/build-gnu13/bdf-pkg-full/bin/molden2resp.x "+dirname+".scf.molden > "+dirname+".charge")
            with open(dirname+'.charge', 'r')as file:
                next(l for l in file if 'Atomic charges:' in l)
                for j in caplist[i]:
                    self.respcharge[j]=float(next(file).split()[1])
            os.chdir("..")
        for i in fraglist: 
            dirname="frag-"+str(i)
            os.system("mkdir "+dirname)
            os.chdir(dirname)
            chargefilename=dirname+".extcharge"
            filename=dirname+".inp"
            self.genbdfeasyinput(fraglist[i],filename,inp_info)
            chargelist=set(self.init_graph[1])-set(fraglist[i])
            self.genbdfextchargefile(chargelist,chargefilename)
            os.system("/home/gao/software/bdf-pkg/build-gnu13/bdf-pkg-full/sbin/run.sh "+filename)
            with open(filename+'.out', 'r')as file:
                Ea=float(next(l for l in file if 'E_tot =' in l).split()[2])
                E_tot+=Ea
                if type(i)==str:
                    E_tot-=E_gcap[int(i.split("+")[1])]
                    E_tot-=E_gcap[int(i.split("+")[0].split(":")[1])]
            os.system("export OMP_NUM_THREADS=1")        
            os.system("/home/gao/software/bdf-pkg/build-gnu13/bdf-pkg-full/bin/molden2resp.x "+dirname+".scf.molden > "+dirname+".charge")
            with open(dirname+'.charge', 'r')as file:
                next(l for l in file if 'Atomic charges:' in l)
                for j in fraglist[i]:
                    self.respcharge[j]=float(next(file).split()[1])
            os.chdir("..")
        
        E_ele=0
        for i in range(2,len(self.amino_acids_list)-1):
            for j in range(i+2,len(self.amino_acids_list)+1):
                if not molecompdis[j][i-1]<=4.0 :
                    for k in self.amino_acids_list[i-1]:
                        for l in self.amino_acids_list[j]:
                            E_ele+=float(self.respcharge[l])*float(self.respcharge[k])/(self.atomdist[l][k]/0.52918)
        print("Column energy%f" % E_ele)
        print("No.1 iteration energy={}".format(E_tot-E_ele))
        
        if niter>=2:
            for ii in range(2,niter+1):
                E_tot=0
                E_gcap={}
                for i in caplist: 
                    dirname="concap-"+str(i)
                    os.chdir(dirname)
                    filename=dirname+".inp"
                    chargefilename=dirname+".extcharge"
                    chargelist=set(self.init_graph[1])-set(caplist[i])
                    self.genbdfextchargefile(chargelist,chargefilename)
                    os.system("/home/gao/software/bdf-pkg/build-gnu13/bdf-pkg-full/sbin/run.sh "+filename)
                    with open(filename+'.out', 'r')as file:
                        Ea=float(next(l for l in file if 'E_tot =' in l).split()[2])
                        if not type(i)==str:
                            E_tot-=Ea   
                        else:
                            E_gcap[int(i.split(":")[1])]=Ea
                    os.system("export OMP_NUM_THREADS=1")                    
                    os.system("/home/gao/software/bdf-pkg/build-gnu13/bdf-pkg-full/bin/molden2resp.x "+dirname+".scf.molden > "+dirname+".charge")
                    with open(dirname+'.charge', 'r')as file:
                        next(l for l in file if 'Atomic charges:' in l)
                        for j in caplist[i]:
                            self.respcharge[j]=float(next(file).split()[1])
                    os.chdir("..")
                for i in fraglist: 
                    dirname="frag-"+str(i)
                    os.chdir(dirname)
                    filename=dirname+".inp"
                    chargefilename=dirname+".extcharge"
                    chargelist=set(self.init_graph[1])-set(fraglist[i])
                    self.genbdfextchargefile(chargelist,chargefilename)
                    os.system("/home/gao/software/bdf-pkg/build-gnu13/bdf-pkg-full/sbin/run.sh "+filename)
                    with open(filename+'.out', 'r')as file:
                        Ea=float(next(l for l in file if 'E_tot =' in l).split()[2])
                        E_tot+=Ea
                        if type(i)==str:
                            E_tot-=E_gcap[int(i.split("+")[1])]
                            E_tot-=E_gcap[int(i.split("+")[0].split(":")[1])]
                    os.system("export OMP_NUM_THREADS=1")                    
                    os.system("/home/gao/software/bdf-pkg/build-gnu13/bdf-pkg-full/bin/molden2resp.x "+dirname+".scf.molden > "+dirname+".charge")
                    with open(dirname+'.charge', 'r')as file:
                        next(l for l in file if 'Atomic charges:' in l)
                        for j in fraglist[i]:
                            self.respcharge[j]=float(next(file).split()[1])
                    os.chdir("..")

                E_ele=0
                for i in range(2,len(self.amino_acids_list)-1):
                    for j in range(i+2,len(self.amino_acids_list)+1):
                        if not molecompdis[j][i-1]<=4.0 :
                            for k in self.amino_acids_list[i-1]:
                                for l in self.amino_acids_list[j]:
                                    E_ele+=self.respcharge[l]*self.respcharge[k]/(self.atomdist[l][k]/0.52918)

                print("Column energy%f" % E_ele)
                print("No."+str(ii)+" iteration energy={}".format(E_tot-E_ele))

        # print("最终能量-   {}".format(E1-E2-E_ele))     

    
    #以下两个函数均用的我的加cap方法
    
    #already have RESP charge
    def EEMFCC_haverespcharge(self,fraglist,caplist,inp_info,molecompdis,Rcharge):
        E_ele=0
        for i in range(2,len(self.amino_acids_list)-1):
                    for j in range(i+2,len(self.amino_acids_list)+1):
                        if not molecompdis[j][i-1]<=4.0 :
                            for k in self.amino_acids_list[i-1]:
                                for l in self.amino_acids_list[j]:
                                    E_ele+=Rcharge[l]*Rcharge[k]/(self.atomdist[l][k]/0.52918)
        print("Column energy%f" % E_ele)
        
        os.system("mkdir EE-GMFCC_calculFile")
        os.chdir("EE-GMFCC_calculFile")
        E_tot=0
        E_gcap={}
        for i in caplist: 
                dirname="concap-"+str(i)
                os.system("mkdir "+dirname)
                os.chdir(dirname)
                charfilena=dirname+".extcharge"
                filename=dirname+".inp"
                self.genbdfeasyinput(caplist[i],filename,inp_info)
                chargelist=set(self.init_graph[1])-set(caplist[i])
                chargetxt="Title \n"+str(len(chargelist))+"\n"
                for j in chargelist:
                    chargetxt+=str(self.label[j])+" "+str(Rcharge[j])+" "+str(self.coord[j][1])+" "+str(self.coord[j][2])+" "+str(self.coord[j][3])+"\n"
                with open(charfilena,"w")as writ:
                    writ.write(chargetxt)
                os.system("/home/gao/software/bdf-pkg/build-gnu13/bdf-pkg-full/sbin/run.sh "+filename)
                with open(filename+'.out', 'r')as file:
                    Ea=float(next(l for l in file if 'E_tot =' in l).split()[2])
                    if not type(i)==str:
                        E_tot-=Ea   
                    else:
                        E_gcap[int(i.split(":")[1])]=Ea
                os.chdir("..")
        for i in fraglist: 
                dirname="frag-"+str(i)
                os.system("mkdir "+dirname)
                os.chdir(dirname)
                charfilena=dirname+".extcharge"
                filename=dirname+".inp"
                self.genbdfeasyinput(fraglist[i],filename,inp_info)
                chargelist=set(self.init_graph[1])-set(fraglist[i])
                chargetxt="Title \n"+str(len(chargelist))+"\n"
                for j in chargelist:
                    chargetxt+=str(self.label[j])+" "+str(Rcharge[j])+" "+str(self.coord[j][1])+" "+str(self.coord[j][2])+" "+str(self.coord[j][3])+"\n"
                with open(charfilena,"w")as writ:
                    writ.write(chargetxt)
                os.system("/home/gao/software/bdf-pkg/build-gnu13/bdf-pkg-full/sbin/run.sh "+filename)
                with open(filename+'.out', 'r')as file:
                    Ea=float(next(l for l in file if 'E_tot =' in l).split()[2])
                    E_tot+=Ea
                    if type(i)==str:
                        E_tot-=E_gcap[int(i.split("+")[1])]
                        E_tot-=E_gcap[int(i.split("+")[0].split(":")[1])]
                os.chdir("..")
        print("No.1 iteration energy={}".format(E_tot-E_ele))
       
        # if niter>=2:
        #     for ii in range(2,niter+1):
        #         E_tot=0
        #         E_gcap={}
        #         for i in caplist: 
        #                 dirname="concap-"+str(i)
        #                 os.chdir(dirname)
        #                 filename=dirname+".inp"
        #                 charfilena=dirname+".extcharge"
        #                 chargelist=set(self.init_graph[1])-set(caplist[i])
        #                 chargetxt="Title \n"+str(len(chargelist))+"\n"
        #                 for j in chargelist:
        #                     chargetxt+=str(self.label[j])+" "+str(Rcharge[j])+" "+str(self.coord[j][1])+" "+str(self.coord[j][2])+" "+str(self.coord[j][3])+"\n"
        #                 with open(charfilena,"w")as writ:
        #                     writ.write(chargetxt)
        #                 os.system("/home/gao/software/bdf-pkg/build-gnu13/bdf-pkg-full/sbin/run.sh "+filename)
        #                 with open(filename+'.out', 'r')as file:
        #                     Ea=float(next(l for l in file if 'E_tot =' in l).split()[2])
        #                     if not type(i)==str:
        #                         E_tot-=Ea   
        #                     else:
        #                         E_gcap[int(i.split(":")[1])]=Ea
        #                 os.chdir("..")
        #         for i in fraglist: 
        #                 dirname="frag-"+str(i)
        #                 os.chdir(dirname)
        #                 filename=dirname+".inp"
        #                 charfilena=dirname+".extcharge"
        #                 chargelist=set(self.init_graph[1])-set(fraglist[i])
        #                 chargetxt="Title \n"+str(len(chargelist))+"\n"
        #                 for j in chargelist:
        #                     chargetxt+=str(self.label[j])+" "+str(Rcharge[j])+" "+str(self.coord[j][1])+" "+str(self.coord[j][2])+" "+str(self.coord[j][3])+"\n"
        #                 with open(charfilena,"w")as writ:
        #                     writ.write(chargetxt)
        #                 os.system("/home/gao/software/bdf-pkg/build-gnu13/bdf-pkg-full/sbin/run.sh "+filename)
        #                 with open(filename+'.out', 'r')as file:
        #                     Ea=float(next(l for l in file if 'E_tot =' in l).split()[2])
        #                     E_tot+=Ea
        #                     if type(i)==str:
        #                         E_tot-=E_gcap[int(i.split("+")[1])]
        #                         E_tot-=E_gcap[int(i.split("+")[0].split(":")[1])]
        #                 os.chdir("..")
        #         E_ele=0
        #         for i in range(2,len(self.amino_acids_list)-1):
        #             for j in range(i+2,len(self.amino_acids_list)+1):
        #                 if not molecompdis[j][i-1]<=4.0 :
        #                     for k in self.amino_acids_list[i-1]:
        #                         for l in self.amino_acids_list[j]:
        #                             E_ele+=Rcharge[l]*Rcharge[k]/(self.atomdist[l][k]/0.52918)

        #         print("Column energy%f" % E_ele)
        #         print("No."+str(ii)+" iteration energy={}".format(E_tot-E_ele))

    #self-iterating  RESP charge   
    def EEMFCC_iter_charge(self,fraglist,caplist,inp_info,molecompdis,niter=1):
        self.respcharge={}
        for i in self.init_graph[1]:
            self.respcharge[i]=0.0

        os.system("mkdir EE-GMFCC_calculFile")
        os.chdir("EE-GMFCC_calculFile")
        E_tot=0
        E_gcap={}
        for i in caplist: 
            dirname="concap-"+str(i)
            os.system("mkdir "+dirname)
            os.chdir(dirname)
            chargefilename=dirname+".extcharge"
            filename=dirname+".inp"
            self.genbdfeasyinput(caplist[i],filename,inp_info)
            chargelist=set(self.init_graph[1])-set(caplist[i])
            self.genbdfextchargefile(chargelist,chargefilename)
            os.system("/home/gao/software/bdf-pkg/build-gnu13/bdf-pkg-full/sbin/run.sh "+filename)
            with open(filename+'.out', 'r')as file:
                Ea=float(next(l for l in file if 'E_tot =' in l).split()[2])
                if not type(i)==str:
                    E_tot-=Ea   
                else:
                    E_gcap[int(i.split(":")[1])]=Ea
            os.system("export OMP_NUM_THREADS=1")                    
            os.system("/home/gao/software/bdf-pkg/build-gnu13/bdf-pkg-full/bin/molden2resp.x "+dirname+".scf.molden > "+dirname+".charge")
            with open(dirname+'.charge', 'r')as file:
                next(l for l in file if 'Atomic charges:' in l)
                for j in caplist[i]:
                    self.respcharge[j]=float(next(file).split()[1])
            os.chdir("..")
        for i in fraglist: 
            dirname="frag-"+str(i)
            os.system("mkdir "+dirname)
            os.chdir(dirname)
            chargefilename=dirname+".extcharge"
            filename=dirname+".inp"
            self.genbdfeasyinput(fraglist[i],filename,inp_info)
            chargelist=set(self.init_graph[1])-set(fraglist[i])
            self.genbdfextchargefile(chargelist,chargefilename)
            os.system("/home/gao/software/bdf-pkg/build-gnu13/bdf-pkg-full/sbin/run.sh "+filename)
            with open(filename+'.out', 'r')as file:
                Ea=float(next(l for l in file if 'E_tot =' in l).split()[2])
                E_tot+=Ea
                if type(i)==str:
                    E_tot-=E_gcap[int(i.split("+")[1])]
                    E_tot-=E_gcap[int(i.split("+")[0].split(":")[1])]
            os.system("export OMP_NUM_THREADS=1")        
            os.system("/home/gao/software/bdf-pkg/build-gnu13/bdf-pkg-full/bin/molden2resp.x "+dirname+".scf.molden > "+dirname+".charge")
            with open(dirname+'.charge', 'r')as file:
                next(l for l in file if 'Atomic charges:' in l)
                for j in fraglist[i]:
                    self.respcharge[j]=float(next(file).split()[1])
            os.chdir("..")
        
        E_ele=0
        for i in range(2,len(self.amino_acids_list)-1):
            for j in range(i+2,len(self.amino_acids_list)+1):
                if not molecompdis[j][i-1]<=4.0 :
                    for k in self.amino_acids_list[i-1]:
                        for l in self.amino_acids_list[j]:
                            E_ele+=float(self.respcharge[l])*float(self.respcharge[k])/(self.atomdist[l][k]/0.52918)
        print("Column energy%f" % E_ele)
        print("No.1 iteration energy={}".format(E_tot-E_ele))
        
        if niter>=2:
            for ii in range(2,niter+1):
                E_tot=0
                E_gcap={}
                for i in caplist: 
                    dirname="concap-"+str(i)
                    os.chdir(dirname)
                    filename=dirname+".inp"
                    chargefilename=dirname+".extcharge"
                    chargelist=set(self.init_graph[1])-set(caplist[i])
                    self.genbdfextchargefile(chargelist,chargefilename)
                    os.system("/home/gao/software/bdf-pkg/build-gnu13/bdf-pkg-full/sbin/run.sh "+filename)
                    with open(filename+'.out', 'r')as file:
                        Ea=float(next(l for l in file if 'E_tot =' in l).split()[2])
                        if not type(i)==str:
                            E_tot-=Ea   
                        else:
                            E_gcap[int(i.split(":")[1])]=Ea
                    os.system("export OMP_NUM_THREADS=1")                    
                    os.system("/home/gao/software/bdf-pkg/build-gnu13/bdf-pkg-full/bin/molden2resp.x "+dirname+".scf.molden > "+dirname+".charge")
                    with open(dirname+'.charge', 'r')as file:
                        next(l for l in file if 'Atomic charges:' in l)
                        for j in caplist[i]:
                            self.respcharge[j]=float(next(file).split()[1])
                    os.chdir("..")
                for i in fraglist: 
                    dirname="frag-"+str(i)
                    os.chdir(dirname)
                    filename=dirname+".inp"
                    chargefilename=dirname+".extcharge"
                    chargelist=set(self.init_graph[1])-set(fraglist[i])
                    self.genbdfextchargefile(chargelist,chargefilename)
                    os.system("/home/gao/software/bdf-pkg/build-gnu13/bdf-pkg-full/sbin/run.sh "+filename)
                    with open(filename+'.out', 'r')as file:
                        Ea=float(next(l for l in file if 'E_tot =' in l).split()[2])
                        E_tot+=Ea
                        if type(i)==str:
                            E_tot-=E_gcap[int(i.split("+")[1])]
                            E_tot-=E_gcap[int(i.split("+")[0].split(":")[1])]
                    os.system("export OMP_NUM_THREADS=1")                    
                    os.system("/home/gao/software/bdf-pkg/build-gnu13/bdf-pkg-full/bin/molden2resp.x "+dirname+".scf.molden > "+dirname+".charge")
                    with open(dirname+'.charge', 'r')as file:
                        next(l for l in file if 'Atomic charges:' in l)
                        for j in fraglist[i]:
                            self.respcharge[j]=float(next(file).split()[1])
                    os.chdir("..")

                E_ele=0
                for i in range(2,len(self.amino_acids_list)-1):
                    for j in range(i+2,len(self.amino_acids_list)+1):
                        if not molecompdis[j][i-1]<=4.0 :
                            for k in self.amino_acids_list[i-1]:
                                for l in self.amino_acids_list[j]:
                                    E_ele+=self.respcharge[l]*self.respcharge[k]/(self.atomdist[l][k]/0.52918)

                print("Column energy%f" % E_ele)
                print("No."+str(ii)+" iteration energy={}".format(E_tot-E_ele))

        # print("最终能量-   {}".format(E1-E2-E_ele))     

