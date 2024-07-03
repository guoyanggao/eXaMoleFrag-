"""
 BDF tools 
 Author: Bingbing Suo
 Version 1. 2013-5-2
 Calculation correlation energy by using molecular fragments

"""

# Calculate correlation energy of a single molecular fragment or a pair of fragments

def portable_float(str1):
    'Converts string to float, but also respects Fortran-type floats like 1.D-6.'
    return float(str1.lower().replace('d','e'))

# find a string in a file: small file
def findstrinfile(filename, lookup):
    return lookup in open(filename,'rt').read()

# find a string in a file: large file
def findstrinlargefile(filename, lookup):
    handle=open(filename,"r")
    for ln in handle:
        if lookup in ln:
            print (ln)
            handle.close()
            return ln 
    else:
       handle.close()
       return False
    #with open(filename, 'rt') as handle:
    #    for ln in handle:
    #        if lookup in ln:
    #            print ln
    #            return ln 
    #    else:
    #       return False

# fagment nmr, generate information for NMR subsystem calculation 
def fragment_nmr(taskname,nopair=False, cimfrag=False,inplines=None):
    import glob            # list files by names
    import shutil          # copy file
    import sys             # read argument 
    import os              # get and set evironment variable
    import subprocess      # execute command
    import re              # keyword search
    #import datetime            # count time
    import time 

    BDFhome    = os.environ.get("BDFHOME")
    BDFtmpdir  = os.environ.get('BDF_TMPDIR')
    BDFworkdir = os.environ.get("BDF_WORKDIR")

    filedat  = BDFworkdir + "/" + taskname + ".datapunch"
    fdat     = open(filedat,"r+")
    lines    = fdat.readlines()
    #print lines

    # read basis and ri basis name, decide RI mode
    notouch = True
    ribasis ="NULL"
    rimode  = 0
    for str1 in lines:
        #print str1
        if str1.startswith("BASIS"):
           list1 = str1.split()
           basis = list1[1].strip()
           if len(list1) > 2:  # BASIS NAME 5d/6d
              basis = basis + "    " + list1[2].strip()
           notouch = False
        if str1.startswith("RIBASIS"):
           list1 = str1.split()
           ribasis = list1[1].strip()
           ribasis = ribasis + "    " + list1[2].strip()
        # RI approaximation mode
        # 0 none RI 1 RI for all  2 RI for fragment pairs 3 RI for strong pairs
        if "RIMODE" in str1:
           list1  = str1.split()
           rimode = int(list1[1].strip())
    if notouch:
        print ("Could not find basis set! Stop ...")
        exit()
    print ("Ribasis",ribasis," rimode",rimode)

    # read number of fragments 
    notouch = True
    for str1 in lines:
        if "MOLEFRAG" in str1:
           list1     = str1.split()
           nfragment = int(list1[1])
           notouch = False 
           break

    print (" nfragment =",nfragment )
    if notouch:
        print ("Could not find MOLEFRAG! Stop ...")
        exit()
  
    # read fragment list 
    notouch = True
    fraglist=[]
    i = -1
    for str1 in lines:
        i = i+1
        if "FRAGMENTLIST" in str1:
           notouch = False
           list1   = str1.split()
           str2 = " "
           k = i
           for j in range(0,nfragment,20):
              k = k+1
              print (lines[k])
              str2 += lines[k].strip("\n")
           break

    if notouch and not nopair:
        print ("Could not find FRAGMENTLIST! Stop ...")
        exit()
    if not notouch:
        fraglist=str2.split()

    print (basis)
    print (nfragment)
    print (fraglist)
    #print pairlist
    #exit()

    starttime=time.time() #datetime.datetime.now()

    # generate nmr input according to user input lines
    from bdfutil import extract_bdf_module_input
    extract_bdf_module_input(inplines,"nmr")

    nmrfrag=True
    print (fraglist)
    #define a list
    energy_monomer=[]
    for j in range(0, nfragment):
        i=int(fraglist[j])
        print ("\n Subsystem NMR - Calculating Fragment Monomer %-4i" % i)
        # copy punch file into scratch dir
        if cimfrag:
            srcfile = BDFworkdir + '/'+ taskname +'.cpunflmo' + str(i) 
        elif nmrfrag:
            srcfile = BDFworkdir + '/'+ taskname +'.nmrsubmole' + str(i) 
            from os import path
            if not path.exists(srcfile):
                srcfile = BDFworkdir + '/'+ taskname +'.spunflmo' + str(i) 
        else:
            srcfile = BDFworkdir + '/'+ taskname +'.spunflmo' + str(i) 

        disfile = BDFtmpdir  + '/fragpunch'
        shutil.copyfile(srcfile, disfile)

        # Prepare the file storing the global coulpot for fragment calculation. 
        from os import path
        from shutil import copy
        denpot = BDFtmpdir + "/" + taskname + ".denpot"
        if path.exists(denpot):
            denpotsub = BDFtmpdir + "/" + taskname+ "_sub"  + ".denpot"
            copy(denpot,denpotsub)
        gridinfo = BDFtmpdir + "/" + taskname + ".gridinfo"
        if path.exists(gridinfo):
            gridinfosub = BDFtmpdir + "/" + taskname+ "_sub"  + ".gridinfo"
            copy(gridinfo,gridinfosub)
        coulinfo = BDFtmpdir + "/" + taskname + ".coulinfo"
        if path.exists(coulinfo):
            coulinfosub = BDFtmpdir + "/" + taskname+ "_sub"  + ".coulinfo"
            copy(coulinfo,coulinfosub)

        # perform tddft calcualtion
        re = frag_nmr_energy(taskname, "nmr", basis, ribasis, rimode,str(i+1))

        # backup interface file

 
    endtime =time.time() #datetime.datetime.now()

    monotime=endtime-starttime #float((endtime-starttime).seconds)

    #print energy_dimmer
   
    print ("\n")
    print ("\n")

    # close datapunch
    fdat.close()    
    return 0


# fagment tddft, generate information for renormlized excitation calculation
def fragment_tddft(taskname,nopair=False, cimfrag=False):
    import glob            # list files by names
    import shutil          # copy file
    import sys             # read argument 
    import os              # get and set evironment variable
    import subprocess      # execute command
    import re              # keyword search
    #import datetime            # count time
    import time 

    BDFhome    = os.environ.get("BDFHOME")
    BDFtmpdir  = os.environ.get('BDF_TMPDIR')
    BDFworkdir = os.environ.get("BDF_WORKDIR")

    filedat  = BDFworkdir + "/" + taskname + ".datapunch"
    fdat     = open(filedat,"r+")
    lines    = fdat.readlines()
    #print lines

    # read basis and ri basis name, decide RI mode
    notouch = True
    ribasis ="NULL"
    rimode  = 0
    for str1 in lines:
        #print str1
        if str1.startswith("BASIS"):
           list1 = str1.split()
           basis = list1[1].strip()
           if len(list1) > 2:  # BASIS NAME 5d/6d
              basis = basis + "    " + list1[2].strip()
           notouch = False
        if str1.startswith("RIBASIS"):
           list1 = str1.split()
           ribasis = list1[1].strip()
           ribasis = ribasis + "    " + list1[2].strip()
        # RI approaximation mode
        # 0 none RI 1 RI for all  2 RI for fragment pairs 3 RI for strong pairs
        if "RIMODE" in str1:
           list1  = str1.split()
           rimode = int(list1[1].strip())
    if notouch:
        print ("Could not find basis set! Stop ...")
        exit()
    print ("Ribasis",ribasis," rimode",rimode)

    # read number of fragments 
    notouch = True
    for str1 in lines:
        if "MOLEFRAG" in str1:
           list1     = str1.split()
           nfragment = int(list1[1])
           notouch = False 
           break
 
    if notouch:
        print ("Could not find MOLEFRAG! Stop ...")
        exit()
  
    # read number of fragment pairs
    notouch = True
    pairlist=[]
    i = -1
    for str1 in lines:
        i = i+1
        if "FRAGPAIRLIST" in str1:
           notouch = False
           list1      = str1.split()
           nfragpairs = int(list1[1])
           str2 = " "
           k = i
           for j in range(0,2*nfragpairs,20):
              k = k+1
              str2 += lines[k].strip("\n")
           break

    if notouch and not nopair:
        print ("Could not find FRAGPAIRLIST! Stop ...")
        exit()
    if not notouch:
        pairlist=str2.split()

    #print basis
    #print nfragment
    #print nfragpairs
    #print pairlist
    #exit()

    starttime=time.time() #datetime.datetime.now()

    #define a list
    energy_monomer=[]
    
    for i in range(0, nfragment):
        print ("\n Calculate Fragment Monomer %-4i" % i)
        # copy punch file into scratch dir
        if cimfrag:
            srcfile = BDFworkdir + '/'+ taskname +'.cpunflmo' + str(i+1) 
        else:
            srcfile = BDFworkdir + '/'+ taskname +'.spunflmo' + str(i+1) 

        disfile = BDFtmpdir  + '/fragpunch'
        shutil.copyfile(srcfile, disfile)

        # perform tddft calcualtion
        re = frag_tddft_energy(taskname, "tddft", basis, ribasis, rimode,str(i+1))

        # backup interface file

 
    endtime =time.time() #datetime.datetime.now()

    monotime=endtime-starttime #float((endtime-starttime).seconds)


    starttime=time.time() #datetime.datetime.now()
    if not nopair: 
        #Calculate fragment pair interaction energy
        print (pairlist)
        for k in range(0, 2*nfragpairs, 2):
            ifrag=pairlist[k]
            jfrag=pairlist[k+1]

            print ("\n Calculate Fragment Dimmer %-4i  %-4i" % (int(ifrag),int(jfrag)))
 
            srcfile = BDFworkdir + '/'+ taskname +'.ppunflmo_' + ifrag + '_' + jfrag
            disfile = BDFtmpdir  + '/fragpunch' 
            shutil.copyfile(srcfile, disfile)

            # calcualte correlation energy
            re = frag_tddft_energy(taskname, "tddft", basis, ribasis, rimode,ifrag,jfrag)

            # backup interface files

 
    endtime =time.time() #datetime.datetime.now()
    dimmtime=endtime-starttime #float((endtime-starttime).seconds)

    #print energy_dimmer
   
    print ("\n")
    print ("\n")

    # close datapunch
    fdat.close()    
    return 0


def fragment_energy(taskname,nopair=False, cimfrag=False):
    import glob            # list files by names
    import shutil          # copy file
    import sys             # read argument 
    import os              # get and set evironment variable
    import subprocess      # execute command
    import re              # keyword search
    #import datetime            # count time
    import time 

    BDFhome    = os.environ.get("BDFHOME")
    BDFtmpdir  = os.environ.get('BDF_TMPDIR')
    BDFworkdir = os.environ.get("BDF_WORKDIR")

    filedat  = BDFworkdir + "/" + taskname + ".datapunch"
    fdat     = open(filedat,"r+")
    lines    = fdat.readlines()
    #print lines

    # read basis and ri basis name, decide RI mode
    notouch = True
    ribasis ="NULL"
    rimode  = 0
    for str1 in lines:
        #print str1
        if str1.startswith("BASIS"):
           list1 = str1.split()
           basis = list1[1].strip()
           if len(list1) > 2:  # BASIS NAME 5d/6d
              basis = basis + "    " + list1[2].strip()
           notouch = False
        if str1.startswith("RIBASIS"):
           list1 = str1.split()
           ribasis = list1[1].strip()
           ribasis = ribasis + "    " + list1[2].strip()
        # RI approaximation mode
        # 0 none RI 1 RI for all  2 RI for fragment pairs 3 RI for strong pairs
        if "RIMODE" in str1:
           list1  = str1.split()
           rimode = int(list1[1].strip())
    if notouch:
        print ("Could not find basis set! Stop ...")
        exit()
    print ("ribasis",ribasis," rimode",rimode)

    # read number of fragments 
    notouch = True
    for str1 in lines:
        if "MOLEFRAG" in str1:
           list1     = str1.split()
           nfragment = int(list1[1])
           notouch = False 
           break
 
    if notouch:
        print ("Could not find MOLEFRAG! Stop ...")
        exit()
  
    # read number of fragment pairs
    notouch = True
    pairlist=[]
    i = -1
    for str1 in lines:
        i = i+1
        if "FRAGPAIRLIST" in str1:
           notouch = False
           list1      = str1.split()
           nfragpairs = int(list1[1])
           str2 = " "
           k = i
           for j in range(0,2*nfragpairs,20):
              k = k+1
              str2 += lines[k].strip("\n")
           break

    if notouch and not nopair:
        print ("Could not find FRAGPAIRLIST! Stop ...")
        exit()
    if not notouch:
        pairlist=str2.split()

    #print basis
    #print nfragment
    #print nfragpairs
    #print pairlist
    #exit()

    starttime=time.time() #datetime.datetime.now()

    #define a list
    energy_monomer=[]
    
    for i in range(0, nfragment):
        print ("\n Calculate Fragment Monomer %-4i" % i)
        # copy punch file into scratch dir
        if cimfrag:
            srcfile = BDFworkdir + '/'+ taskname +'.cpunflmo' + str(i+1) 
        else:
            srcfile = BDFworkdir + '/'+ taskname +'.spunflmo' + str(i+1) 

        disfile = BDFtmpdir  + '/fragpunch'
        shutil.copyfile(srcfile, disfile)

        enelist=[]
        enelist.append("SFRAG")
        enelist.append(str(i+1))

        # calcualte correlation energy
        re = frag_correlation_energy(taskname, "mp2", basis, ribasis, rimode)

        ## bbs test debug tijab - rename tijab file
        #print "Rename tijab to calculate wave function overlap. Should remove! BBS DEBUG"
        #fname = BDFtmpdir + "/" + taskname + ".tijab"
        #oname = BDFtmpdir + "/" + taskname + ".tijab" + str(i+1) 
        #os.rename(fname,oname) 
   
        # rewind fdat, then find last ecorrelation
        fdat.seek(0,0)
        lines = fdat.readlines()
        j = len(lines)-1
        ecore=0.0
        ecorc=0.0
        ecord=0.0
        for k in range(j,-1,-1):
            if "ECORRELATION" in lines[k]:
                print (lines[k])
                list1 = lines[k].strip("\n")
                list1 = list1.split()
                ecore = list1[1]
                if cimfrag:
                    ecorc = list1[2]
                    ecord = list1[3]
                    ecora = list1[4]
                break
        if cimfrag :
            enelist.append([ecore,ecorc,ecord,ecora])
            energy_monomer.append(enelist)
        else:
            enelist.append(ecore)
            energy_monomer.append(enelist)
 
    endtime =time.time() #datetime.datetime.now()

    monotime=endtime-starttime #float((endtime-starttime).seconds)

    print (energy_monomer)

    starttime=time.time() #datetime.datetime.now()
    energy_dimmer=[] 
    if not nopair: 
        #Calculate fragment pair interaction energy
        print (pairlist)
        for k in range(0, 2*nfragpairs, 2):
            ifrag=pairlist[k]
            jfrag=pairlist[k+1]

            print ("\n Calculate Fragment Dimmer %-4i  %-4i" % (int(ifrag),int(jfrag)))
 
            srcfile = BDFworkdir + '/'+ taskname +'.ppunflmo_' + ifrag + '_' + jfrag
            disfile = BDFtmpdir  + '/fragpunch' 
            shutil.copyfile(srcfile, disfile)

            enelist=[]
            enelist.append("PFRAG")
            enelist.append(ifrag)
            enelist.append(jfrag)

            # calcualte correlation energy
            re = frag_correlation_energy(taskname, "mp2", basis, ribasis, rimode)
 
            #print "Rename tijab to calculate wave function overlap. Should remove! BBS DEBUG"
            #fname = BDFtmpdir + "/" + taskname + ".tijab"
            #oname = BDFtmpdir + "/" + taskname + ".tijab_" + ifrag + "_" + jfrag 
            #os.rename(fname,oname) 
            
            # rewind fdat, then find last ecorrelation
            fdat.seek(0,0)
            lines = fdat.readlines()
            j = len(lines)-1
            ecore = 0.0
            for k in range(j,-1,-1):
                if "ECORRELATION" in lines[k]:
                    print (lines[k])
                    list1 = lines[k].strip("\n")
                    list1 = list1.split()
                    ecore = list1[1]
                    break

            enelist.append(ecore)
            energy_dimmer.append(enelist)

    endtime =time.time() #datetime.datetime.now()
    dimmtime=endtime-starttime #float((endtime-starttime).seconds)

    #print energy_dimmer
   
    print ("\n")
    print ("\n")

    print ("|***********************************************************************|")
    print (" Summary fragment correlation energy \n")
    list1=[]
    #print and sum over correlation energies
    enemono=0.0
    enedim1=0.0
    enedim2=0.0
    totene1=0.0
    # print monomer energy
    for line in energy_monomer:
        if "SFRAG" in line[0]:
            if cimfrag:
                list1=line[2]
                print ("  %s  %-4i    %14.8f %14.8f %14.8f %14.8f" % (line[0],int(line[1]),float(list1[0]),float(list1[1]),float(list1[2]),float(list1[3])))
                totene1 = totene1+float(list1[0])  
                enemono = enemono+float(list1[1])  
                enedim1 = enedim1+float(list1[2])  
                enedim2 = enedim2+float(list1[3])                      
            else:
                print ("  %s  %-4i    %14.8f" % (line[0],int(line[1]),float(line[2])))
                enemono = enemono+float(line[2])
 
    # print dimmer energy
    enedimm=0.0
    if not nopair:
       for line in energy_dimmer:
           if "PFRAG" in line[0]:
               print ("  %s  %-4i %-4i %14.8f" % (line[0],int(line[1]),int(line[2]),float(line[3])))
               enedimm = enedimm+float(line[3])

    # summary energy
    if cimfrag:
        enedimm=(enedim1+enedim2)/2.0
        totene = enemono+enedimm
        if abs(totene - totene1) > 1.E-8:
            print (" CIM energy may error! TotE1 = %14.8f TotE2= %14.8f" % (totene, totene1))
        if totene < 1.E-3: totene=totene1
    else:
        totene = enemono+enedimm
   
    #fn=open("")

    list1=[]

    print ("\n")
    print ("  Intra-fragment energy: %14.8f" % enemono)
    if cimfrag:
        print ("  Inter-fragment energy: %14.8f %14.8f %14.8f \n" % (enedimm,enedim1,enedim2))
    else:
        print ("  Inter-fragment energy: %14.8f" % enedimm)
    print ("  Total          energy: %14.8f" % totene)
    print ("\n")
    print ("  Wall time to calculate single fragments: %14.2f S" % monotime)
    print ("  Wall time to calculate  fragment-pairs : %14.2f S" % dimmtime)
    tottime = monotime + dimmtime
    print ("  Total wall time: %14.2f S" % tottime) #monotime+dimmtime)
    print ("|***********************************************************************|")
    print ("\n")

    fdat.seek(0, 2) # move file pointer to tail for adding information
    fdat.write("\n")
    str1 = "  Intra-fragment energy: %14.8f \n" % enemono
    fdat.write(str1)
    if cimfrag:
        str1 = "  Inter-fragment energy: %14.8f %14.8f %14.8f \n" % (enedimm,enedim1,enedim2)
    else: 
        str1 = "  Inter-fragment energy: %14.8f \n" % enedimm
    fdat.write(str1)
    str1 = "  Total          energy: %14.8f \n" % totene
    fdat.write(str1)
    fdat.write("\n")
    
    # close datapunch
    fdat.close()    
    return 0

# Fragment optimization in correlation calculation
def fragopt(finp, taskname, tmpdir, workdir, bdfhome):
    from os import system,environ

    # read all input lines into a list
    ifdecfrag= False
    ifcimfrag= False
    fn=open(finp,"r")
    inplines=[]
    i=-1
    for str1 in fn.readlines():
        str1=str1.strip("\n")
        str1=str1.strip(" ")
        # line starts with # and * is comment line
        if str1.startswith("#") or str1.startswith("*"):
            # skip
            j=1 
        elif str1.strip()=='':
            j=1
        elif str1.startswith("$end"):
            j=1  # skip $end
        elif str1.startswith("decfrag"):
            inplines.append(str1)   
            ifdecfrag = True     
        elif str1.startswith("cimfrag"):
            inplines.append(str1)   
            ifcimfrag = True     
        else:
            i=i+1
            if not str1.startswith("%"):
                str1=str1.lower()
            inplines.append(str1) 
    fn.close()
    print (inplines)

    # set radii to default value, unit is angstrom
    if ifdecfrag:
        revir = 2.0  # Radius of Unoccupied EOS
        rbocc = 1.0  # Radius Occupied buffers
        rbvir = 0.0  # Radius Unoccupied buffers
    elif ifcimfrag:
        revir = 3.0  # occupied buff 
        rbocc = 1.5  # enviromental buff
        rbvir = 0.0  # none
    else:
        revir = 3.0
        rbocc = 15.0
        rbvir = 0.0
    #rstep = 3*0.52977 # step size 3 au
    rstep = 1.0  # step size 1.0 au

    maxstep = 20   # maxium allowed steps
    threshfrag = 0.0001 # default threshhold for changing of the fragment correlation energy

    # check if user set inital value 
    ifsetrad     = False   
    ifsetithfrag = False
    radii        = [] 
    unitbohr     = False
    noenergy     = False
    nopairenergy = False
    rdistpair    = 0.0
    i=-1
    for str1 in inplines:
       i=i+1
       # user has set default values, use them
       if str1.startswith("threshrad"):
           ifsetrad = True
           radii    = inplines[i+1].split()
           radii[0] = float(radii[0])
           radii[1] = float(radii[1])
           radii[2] = float(radii[2])
           revir    = radii[0]
           rbocc    = radii[1]
           rbvir    = radii[2]
       elif str1.startswith("radstep"):
           rstep = float(inplines[i+1])
       elif str1.startswith("optfrag"):
           if inplines[i+1].startswith("all"):
               j=0
           else:
               fraglist = inplines[i+1].split()
               ifsetithfrag = True
       elif str1.startswith("threshfrag"):
           threshfrag=float(inplines[i+1])
       elif str1.startswith("maxstep"):
           maxstep=int(inplines[i+1])
       elif str1.startswith("unit"):
           if inplines[i+1].startswith("bohr"):
               unitbohr = True
       elif str1.startswith("noenergy"):
           noenergy = True
       elif str1.startswith("nopair"):
           nopairenergy = True
       elif str1.startswith("distpair"):
           rdistpair=float(inplines[i+1])

    ## if user prefer to use bohr as units, change to angstrom
    #if unitbohr:
    #    revir = revir/0.529177  
    #    rbocc = rbocc/0.529177  
    #    rbvir = rbvir/0.529177  
    #    rstep = 3  # default step size 3 a.u.

    # user do not set inital radii
    if not ifsetrad:
        radii=[]
        radii.append(revir)
        radii.append(rbocc)
        radii.append(rbvir)
        inplines.append("threshrad") 
        str1 = str(radii[0]) + ' ' + str(radii[1]) + ' ' + str(radii[2])
        inplines.append(str1)
 
    radii_start=radii[:]
    # end initial radii

    ncompinafrag={}
    if not ifsetithfrag:
        # read number of fragments form database
        nfrag = 0
        ncomp = 0
        datafil = tmpdir +"/"+"database"
        fndat   = open(datafil, "r")
        tlines  = fndat.readlines()
        fndat.close()
        i=-1
        for str1 in tlines:
            i = i+1
            if str1.startswith("molefrag"):
                nfrag = int(tlines[i+1])     
                for j in range(i+2,i+2+nfrag):
                    tmp = tlines[j].split()
                    ncompinafrag[j-i-2+1]=int(tmp[0]) 
                #break
            if str1.startswith("molecomp"):
                ncomp = int(tlines[i+1])

        if nfrag == 0:
            print (" Could not find molefrag defination! Assume each molecomp is a center molefrag")
            nfrag = ncomp
            for j in range(0,nfrag):
                #ncompinafrag.append(1)
                ncompinafrag[j+1]= 1

        fraglist=[]
        for i in range(0, nfrag):
            fraglist.append(str(i+1))
    else:
        # read number of molecomp of ithfrag
        ncomp = 0
        ntfra = 0
        datafil = tmpdir +"/"+"database"
        fndat   = open(datafil, "r")
        tlines  = fndat.readlines()
        fndat.close()
        i=-1
        for str1 in tlines:
            i = i+1
            if str1.startswith("molefrag"):
                ntfra = int(tlines[i+1]) 
                tmpnc = {}   
                for k in range(0,ntfra):
                    j = i + 2 + k
                    tmp = tlines[j].split()
                    ncompinafrag[k+1]=int(tmp[0]) 
                      
                #break
            if str1.startswith("molecomp"):
                ncomp = int(tlines[i+1])
        if ntfra == 0:
            for j in range(0,ncomp):
                #ncompinafrag.append(1)
                ncompinafrag[j+1]=1 


    print ("\n")
    print (" optimizing fragment lists")
    print (fraglist)
    print (ncompinafrag)
    print (" threshfrag:  %f" % threshfrag)
    print (" maxiter:     %i" % maxstep)
    print (radii)
    print ("\n")

    # cycle fragment
    optoutinfo = [] # optimization information
    norb=[]
    for ithfrag in fraglist:
        radii = []
        radii = radii_start[:]  # copy radii_start to radii
        print ("radii_start %6.2f %6.2f %6.2f" % (radii_start[0],radii_start[1],radii_start[2]))
        print ("radii %6.2f %6.2f %6.2f" % (radii[0],radii[1],radii[2]))
 
        #str1 = "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
        #optoutinfo.append(str1)    
        str1 = "Optimizing fragment %d with FOT %f \n" % (int(ithfrag), threshfrag*ncompinafrag[int(ithfrag)])
        optoutinfo.append(str1)  
        if ifdecfrag:
            str1 = "Iter  istep          Ecor             DeltaEcor   Radii_1   Radii_2   Radii_3 \n"
        else:  
            str1 = "Iter          Ecor             DeltaEcor   Radii  \n"
        optoutinfo.append(str1)    
 
        # now we begin optimize fragment
        norb     = "ORBNUM 0 0 0 0"
        environ["ORBNUM"]=norb
        Ecor_old = 0.0
        
        converge = False
        simplified_dec = False #True
 
        if ifdecfrag: 
            for it in range(0,maxstep+1): 
                if converge: break

                if simplified_dec: 
                    maxk = 1
                else:
                    maxk = 2
                for istep in range(0, maxk):  # 0 enlarge unoccpied EOS, 1 enlarge occupied AOS
                    if it == 0:
                        if istep == 1: break
                    else:
                        if simplified_dec:
                            radii[0] = radii[0] + rstep 
                            radii[1] = radii[1] + rstep 
                            radii[2] = radii[2] + 0.0    #rstep
                        else:
                            # step forward - enlarge 
                            if istep == 0:   # enlarge  unoccpied EOS, 
                                radii[0] = radii[0] + rstep 
                                radii[1] = radii[1] 
                                radii[2] = radii[2] + 0.0    #rstep
                            elif istep == 1: # enlarge occupied AOS
                                radii[0] = radii[0]  
                                radii[1] = radii[1] + rstep
                                radii[2] = radii[2] + 0.0    #rstep
                        #    print "radii",radii,it,istep 

                    environ["OPTFRAGITER"]=str(it)
                    environ["OPTFRAGSTEP"]=str(istep)
                    updated=fragopt_kernal(finp, inplines, ithfrag, radii, taskname, tmpdir, workdir, norb)

                    if updated:
                        # read newest correlation energy
                        fdat = workdir+"/"+taskname+".datapunch"
                        fn=open(fdat,"r")
                        lines = fn.readlines()
                        j = len(lines)-1
                        Ecor_new = 0.0
                        for k in range(j,-1,-1):
                            if "ECORRELATION" in lines[k]:
                                print (lines[k])
                                list1 = lines[k].strip("\n")
                                list1 = list1.split()
                                Ecor_new = list1[1]
                                break
                        fn.close()

                        Ecor_new  =float(Ecor_new)
                        Delta_Ecor=Ecor_new-Ecor_old
                        print ("summery energy ",Ecor_old,Ecor_new,Delta_Ecor)

                        if it == 0:
                            str1 = "%4d  %4d  %12.8f   %12.8f %8.3f %8.3f %8.3f \n" % (it, istep, Ecor_new, Ecor_new,  radii[0],radii[1],radii[2])
                        else:
                            str1 = "%4d  %4d  %12.8f   %12.8f %8.3f %8.3f %8.3f \n" % (it, istep, Ecor_new, Delta_Ecor,radii[0],radii[1],radii[2])
                        optoutinfo.append(str1)

                        Ecor_old = Ecor_new
                        if it == 0:
                            break

                        if abs(Delta_Ecor) < threshfrag:
                            print (" FOT converge. Stop fragment optimization!")
                            converge = True
                            break
                    else:
                        # orbital space do not update, write old correlation energy into punch file
                        fdat = workdir+"/"+taskname+".datapunch"
                        fn=open(fdat,"a")
                        str1="ECORRELATION    "+str(Ecor_old)+"\n"
                        fn.write(str1) 
                        fn.close()
                        str1 = "%4d  %4d  %12.8f   %12.8f %8.3f %8.3f %8.3f skipped \n" % (it, istep, Ecor_old, 0.0,radii[0],radii[1],radii[2])
                        optoutinfo.append(str1)


        elif ifcimfrag:
            for it in range(0,maxstep+1): 
                if converge: break

                for istep in range(0, 1):  # 0 enlarge occupied buffer, 1 enlarge occupied enviromental
                    if it == 0:
                        if istep == 1: break
                    else:
                        if istep == 0:
                            radii[0] = radii[0] + rstep 
                            radii[1] = radii[1] 
                            radii[2] = radii[2] + 0.0    #rstep
                        elif istep == 1:
                            radii[0] = radii[0]  
                            radii[1] = radii[1] + rstep
                            radii[2] = radii[2] + 0.0    #rstep
                        #    print "radii",radii,it,istep 

                    environ["OPTFRAGITER"]=str(it)
                    environ["OPTFRAGSTEP"]=str(istep)
                    updated=fragopt_kernal(finp, inplines, ithfrag, radii, taskname, tmpdir, workdir, norb, ifcimfrag)

                    if updated:
                        # read newest correlation energy
                        fdat = workdir+"/"+taskname+".datapunch"
                        fn=open(fdat,"r")
                        lines = fn.readlines()
                        j = len(lines)-1
                        Ecor_new = 0.0
                        for k in range(j,-1,-1):
                            if "ECORRELATION" in lines[k]:
                                print (lines[k])
                                list1 = lines[k].strip("\n")
                                list1 = list1.split()
                                Ecor_new = list1[1]
                                break
                        fn.close()

                        Ecor_new  =float(Ecor_new)
                        Delta_Ecor=Ecor_new-Ecor_old
                        print ("summery energy ",Ecor_old,Ecor_new,Delta_Ecor)

                        if it == 0:
                            str1 = "%4d  %4d  %12.8f   %12.8f %8.3f %8.3f %8.3f \n" % (it, istep, Ecor_new, Ecor_new,  radii[0],radii[1],radii[2])
                        else:
                            str1 = "%4d  %4d  %12.8f   %12.8f %8.3f %8.3f %8.3f \n" % (it, istep, Ecor_new, Delta_Ecor,radii[0],radii[1],radii[2])
                        optoutinfo.append(str1)

                        Ecor_old = Ecor_new
                        if it == 0:
                            break

                        if abs(Delta_Ecor) < threshfrag*ncompinafrag[int(ithfrag)]:
                            print (" FOT converge. Stop fragment optimization!")
                            converge = True
                            break
                    else:
                        # orbital space do not update, write old correlation energy into punch file
                        fdat = workdir+"/"+taskname+".datapunch"
                        fn=open(fdat,"a")
                        str1="ECORRELATION    "+str(Ecor_old)+"\n"
                        fn.write(str1) 
                        fn.close()
                        str1 = "%4d  %4d  %12.8f   %12.8f %8.3f %8.3f %8.3f skipped \n" % (it, istep, Ecor_old, 0.0,radii[0],radii[1],radii[2])
                        optoutinfo.append(str1)


        else:
            for it in range(0,maxstep+1):
                if converge : break
 
                if it != 0:
                    # step forward, one parameter
                    radii[0] = radii[0] + rstep 

                environ["OPTFRAGITER"]=str(it) 
                updated=fragopt_kernal(finp, inplines, ithfrag, radii, taskname, tmpdir, workdir, norb)
 
                if updated:
                    # read newest correlation energy
                    fdat = workdir+"/"+taskname+".datapunch"
                    fn=open(fdat,"r")
                    lines = fn.readlines()
                    j = len(lines)-1
                    Ecor_new = 0.0
                    for k in range(j,-1,-1):
                        if "ECORRELATION" in lines[k]:
                            print (lines[k])
                            list1 = lines[k].strip("\n")
                            list1 = list1.split()
                            Ecor_new = list1[1]
                            break
                    fn.close()

                    Ecor_new=float(Ecor_new)
                    Delta_Ecor=Ecor_new-Ecor_old
                    if it == 0:
                        str1 = "%4d    %12.8f   %12.8f %8.3f \n" % (it, Ecor_new, Ecor_new,  radii[0])
                    else:
                        str1 = "%4d    %12.8f   %12.8f %8.3f \n" % (it, Ecor_new, Delta_Ecor,radii[0])
                    optoutinfo.append(str1)
                    Ecor_old = Ecor_new

                    if abs(Delta_Ecor) < threshfrag*ncompinafrag[int(ithfrag)]:
                        print (" FOT converge. Stop fragment optimization!")
                        converge = True
                        break
                else:
                    fdat = workdir+"/"+taskname+".datapunch"
                    fn=open(fdat,"a")
                    str1="ECORRELATION    "+str(Ecor_old)+"\n"
                    fn.write(str1) 
                    fn.close()
                    str1 = "%4d  %4d  %12.8f   %12.8f %8.3f  skipped \n" % (it, istep, Ecor_old, 0.0,radii[0])
                    optoutinfo.append(str1)


    fdat = workdir+"/"+taskname+".datapunch"
    fn=open(fdat,"a")
 
    print ("\n|--------------------------------------------------------------------------|")
    str1 = "Summary of the fragments' optimization ...\n"
    print (str1.strip("\n"))
    fn.write(str1)   
    for str1 in optoutinfo:
        print (str1.strip("\n"))
        fn.write(str1)  
    print ("|--------------------------------------------------------------------------|\n")
    str1 =  "|--------------------------------------------------------------------------|\n"
    fn.write(str1) 
    fn.write("\n")
    fn.close()

    # fragment have been optimized, generate fragment pairs
    if not ifsetithfrag and not nopairenergy:
       print (" Generate fragment pairs ...")
       re=gen_fragpairs(finp, bdfhome, taskname, tmpdir, rdistpair)
       if re !=0:
           exit()

       if not noenergy: 
           # if fragment pairs have been generated, perform fragment mp2
           from fragutil import fragment_energy
           print (" Envoke fragment correlation calculation automaticlly ...")
           re=fragment_energy(taskname, nopairenergy, False)
 
    return 0
#end of fragopt

#For given optimized fragments, genrate fragment pairs
def gen_fragpairs(finp, bdfhome, taskname, tmpdir, rdistpair):
    from os import environ
    from bdfutil import bdf_exec
 
    finp = tmpdir + "/" + taskname + ".genfrag" + ".inp"
    #execute compass
    fn=open(finp, "w+")
    fn.write("$genfrag\n")
    fn.write("fragpairs\n")
    fn.write("energy\n")
    fn.write("mp2\n")
    if rdistpair > 0.0:
        fn.write("distpair\n")
        str1 = "%10.4f \n" % (rdistpair)
        fn.write(str1)
    fn.write("$end\n")
    fn.close()
 
    shcmd = bdfhome + "/bin/genfrag.x" + "<" + finp
    re=bdf_exec("genfrag", finp, None, True) #system(shcmd)
    if re !=0:
        print ("bdf module genfrag failed! STOP...")
        return re

    return re
#end of gen_fragpairs

# 
def fragopt_kernal(finp, inplines, ithfrag, radii, taskname, tmpdir, workdir, norb, ifcimfrag=False):
     from os import system,environ
     from bdfutil import bdf_exec
 
     re = 0
     # genarate genfrag input and generate fragment 
     orbupdate=genfraginp(finp, inplines, ithfrag, radii, taskname, tmpdir, norb)
     if not orbupdate:
        print ("Orbital do not update! Skip this step ...")
        return False
     else:
        # calculate Ecor for a give punch file
        if ifcimfrag:
            punchfile = taskname + "." + "cpunflmo"+ithfrag
        else:
            punchfile = taskname + "." + "spunflmo"+ithfrag
        shcmd = "cp " + workdir + "/" + punchfile + " " + tmpdir +  "/fragpunch"
        re = system(shcmd)
        if re != 0:
            print ("Copy punch file error! STOP")
            exit()

        # read basis name
        filedat  = workdir + "/" + taskname + ".datapunch"
        fdat     = open(filedat,"r+")
        lines    = fdat.readlines()
        fdat.close()

        # read basis and ri basis name, decide RI mode
        notouch = True
        ribasis ="NULL"
        rimode  = 0
        for str1 in lines:
            #print str1
            if str1.startswith("BASIS"):
               list1 = str1.split()
               basis = list1[1].strip()
               if len(list1) > 2:  # BASIS NAME 5d/6d
                  basis = basis + "    " + list1[2].strip()
               notouch = False
            if str1.startswith("RIBASIS"):
               list1 = str1.split()
               ribasis = list1[1].strip()
               ribasis = ribasis + "    " + list1[2].strip()
            # RI approaximation mode
            # 0 none RI 1 RI for all  2 RI for fragment pairs 3 RI for strong pairs
            if "RIMODE" in str1:
               list1  = str1.split()
               rimode = int(list1[1].strip())
        if notouch:
            print ("Could not find basis set! Stop ...")
            exit()
        print ("ribasis",ribasis," rimode",rimode)

        re = frag_correlation_energy(taskname, "mp2", basis, ribasis, rimode)
        if re != 0:
            print ("Calculate fragment correlation energy error! STOP")
            exit()
        return True
 
#end of fragopt

#Generate input file for module genfrag
def genfraginp(finp, inplines, ithfrag, radii, taskname, tmpdir, norb):
    #finp=tmpdir + "/" + taskname + "." + "genfrag" + ".inp"
    from os import system, environ
    from bdfutil import check_module_result,bdf_exec

    bdfhome = environ.get("BDFHOME")
    workdir = environ.get("BDF_WORKDIR")
    ithiter = environ.get("OPTFRAGITER")
    ithiter = int(ithiter)  # present iter, start from 0

    # get old number of orbitals
    norb=environ.get("ORBNUM")
    print ("norb old",norb)

    fn=open(finp,"w+")
    i=-1
    for str1 in inplines:
        i=i+1
        # write present line into input file
        fn.write(str1)
        fn.write("\n")    
 
        if str1.startswith("optfrag"):
            inplines[i+1]=str(ithfrag)
        elif str1.startswith("threshrad"):
            str2 = str(radii[0]) + ' ' +str(radii[1]) + ' ' + str(radii[2])
            inplines[i+1]=str2
 
    print (inplines)

    # add $end
    str1="$end"
    fn.write(str1)
    fn.write("\n")
    fn.close()

    str1 = "cat " + finp
    system(str1)  
      
    # perform genfrag.x
    bdfcmd = bdfhome + "/bin/" + "genfrag.x" + "<" + finp  #tmpdir + "/" + taskname + "." + "genfrag" + ".inp" 
    re = bdf_exec("genfrag", finp, None, True) #system(bdfcmd)
    if re != 0:
        print ("Generate fragment error! STOP")
        exit()
    re = check_module_result("genfrag")

    list1=norb.split()
    list2=list1
    #print list1
    #print list2

    fdat = workdir+"/"+taskname+".datapunch"
 
    shcmd = "cat " + fdat
    system(shcmd)

    fn=open(fdat,"r+")
    lines = fn.readlines()
    nline = len(lines)
    for i in range(nline,0,-1):
        str1 = lines[i-1]
        if str1.startswith("ORBNUM"):
             norb =str1
             list2=str1.split()
             break

    # move pointer to end of datapunch
    fn.seek(0,2)
    str1="RADII " + str(radii[0]) + " " + str(radii[1]) + " " + str(radii[2])
    fn.write(str1)
    fn.write("\n")
    print ("norb new",norb)

    if ithiter == 0:
        iforbupdated = True
    else:
        iforbupdated = False
        for i in range(1,5):
            print (i)
            if list1[i] != list2[i]:
                iforbupdated= True            
                break
        fn.close()
        print ("\n  compare orbital number")
        print ("  Orb old",list1)
        print ("  Orb new",list2)

    # set to new number of orbital
    environ["ORBNUM"]=norb.strip()

    if not iforbupdated:
        return False # orb not update
    else:
        return True  # orb updated


# For the given radii: revir, rbocc,rbvir, calculate correlation energy
def frag_correlation_energy(taskname, posthf, basis, ribas="NULL", rimode=0):
    from os   import environ, system, remove, rename
    from glob import glob
    from bdfutil import check_module_result,bdf_exec

    tmpdir  = environ.get("BDF_TMPDIR")
    workdir = environ.get("BDF_WORKDIR")
    bdfhome = environ.get("BDFHOME")

    punfile = tmpdir  + '/fragpunch'
 
    # rimode 0, none RI 1 RI all 2 RI For pair, 3 RI for strong pair
    rimp2 = False
    strongpair = False
    if ribas != "NULL" and rimode > 0:
        if rimode == 1: # all RI
           rimp2 = True
        if rimode > 1:
           rimp2 = False
           fn=open(punfile,"r")
           lines=fn.readlines()
           for str1 in lines:
               if str1.startswith("FRAGTYPE"):
                   list1=str1.split()
                   if int(list1[1].strip()) == 1:
                       break
                   else:
                      if rimode == 2:
                          rimp2 = True
                          break
               if str1.startswith("STRONGPAIR"):
                   list1=str1.split()
                   if list1[1].strip() == "TRUE":
                       strongpair = True
                       if rimode == 3:
                           rimp2 = False
                       else:
                           rimp2 = True
                       break
                   else:
                       strongpair = False
                       rimp2 = True
                       break
           fn.close()

    print ("rimp2",rimp2, "  strongpair",strongpair, " ribas",ribas)
    print ("basis",basis)
    #return 0

    # generate input files for fragment correlation calculation
    # set BDFTASK to a temparory name
    tmptask=taskname + "_sub" 

    environ["BDFTASK"]  = tmptask
    environ["XESTTASK"] = tmptask

    finp = tmpdir + "/" + tmptask + ".inp"
    #execute compass
    fn=open(finp, "w+")
    fn.write("$compass\n")
    fn.write("Fragpunch\n")
    fn.write("Basis\n")
    fn.write(basis.strip())
    fn.write("\n")
    print (ribas)
    if rimp2 :
       list1=ribas.split()
       fn.write(list1[0])
       fn.write("\n")
       fn.write(list1[1])
    fn.write("\n")
    fn.write("Skeleton\n")
    fn.write("nosymm\n")
    fn.write("$end\n")
    fn.close()

    print ("\n***** Print input file for checking ******")
    shcmd = "cat " + finp   
    re=system(shcmd)
    print ("***** End of input file ****")

    #shcmd = bdfhome + "/bin/compass.x" + "<" + finp
    re=bdf_exec("compass", finp, None, True) #system(shcmd)
    if re !=0:
        print ("bdf module compass failed! STOP...")
        exit()
    re = check_module_result("compass")

    #execute xuanyuan
    fn=open(finp, "w+")
    fn.write("$Xuanyuan\n")
    fn.write("Direct\n")
    fn.write("Schwarz\n")
    fn.write("Maxmem\n")
    fn.write(" 2GW\n")
    fn.write("$End\n")
    fn.close()
    #shcmd = "cat " + finp
    print ("\n***** Print input file for checking ******")
    shcmd = "cat " + finp   
    re=system(shcmd)
    print ("***** End of input file ****")
    #system(shcmd)
    #shcmd = bdfhome + "/bin/xuanyuan.x" + "<" + finp
    re=bdf_exec("xuanyuan", finp, None, True) #system(shcmd)
    if re !=0:
        print ("bdf module xuanyuan failed! STOP...")
        exit()
    re = check_module_result("xuanyuan")


    #execute mp2/ccsd
    fn=open(finp, "w+")
    fn.write("$Mp2\n")
    fn.write("PCMO\n")
    fn.write("$End\n")
    fn.close()
    print ("\n***** Print input file for checking ******")
    shcmd = "cat " + finp   
    re=system(shcmd)
    print ("***** End of input file ****")
    #shcmd = bdfhome + "/bin/mp2.x" + "<" + finp
    re=bdf_exec("mp2", finp, None, True) #system(shcmd)
    if re !=0:
        print ("bdf module mp2 failed! STOP...")
        exit()
    re = check_module_result("mp2")

    
    # add content of tmptask.datapunch to taskname.datapunch
    srcfile = workdir + "/" + tmptask  + ".datapunch"
    disfile = workdir + "/" + taskname + ".datapunch"

    fns=open(srcfile,"r")
    fnd=open(disfile,"a")
    line1=fns.readlines()
    for str1 in line1:
      fnd.write(str1)

    fnd.close()
    fns.close()

    # delete scratch file
    remove(srcfile)
    srcfile = workdir + "/" + tmptask  + ".chkfil"
    remove(srcfile)

    # bbs test debug tijab - rename tijab file
    print ("Rename tijab to calculate wave function overlap. Should remove! BBS DEBUG")
    fname = tmpdir + "/" + tmptask + ".tijab"
    oname = tmpdir + "/" + taskname + ".tijab" 
    rename(fname,oname) 

    fname=tmpdir+"/"+tmptask + "*"
    filelist=glob(fname)
    #print filelist
    for fil in filelist:
         remove(fil)

    # set back taskname
    environ["BDFTASK"]  = taskname
    environ["XESTTASK"] = taskname

    print (" Finish fragment correlation energy calculation ..." )
    return 0

#Input for renormilized excition, calcualte sub-molecule and sub-molecule pairs via tddft
def frag_tddft_energy(taskname, posthf, basis, ribas="NULL", rimode=0,ifrag="NULL",jfrag="NULL"):
    from os   import environ, system, remove, rename
    from glob import glob
    from bdfutil import check_module_result,bdf_exec

    tmpdir  = environ.get("BDF_TMPDIR")
    workdir = environ.get("BDF_WORKDIR")
    bdfhome = environ.get("BDFHOME")

    punfile = tmpdir  + '/fragpunch'
 
    # rimode 0, none RI 1 RI all 2 RI For pair, 3 RI for strong pair
    rimp2      = False
    strongpair = False
    ifdimmer   = False
    fn         = open(punfile,"r")
    if ribas != "NULL" and rimode > 0:
        if rimode == 1: # all RI
           rimp2 = True
        if rimode > 1:
           rimp2 = False
           lines=fn.readlines()
           for str1 in lines:
               if str1.startswith("FRAGTYPE"):
                   list1=str1.split()
                   if int(list1[1].strip()) == 1:
                       ifdimmer = False
                       break
                   else:
                       ifdimmer = True
                       if rimode == 2:
                           rimp2 = True
                           break
               if str1.startswith("STRONGPAIR"):
                   list1=str1.split()
                   if list1[1].strip() == "TRUE":
                       strongpair = True
                       if rimode == 3:
                           rimp2 = False
                       else:
                           rimp2 = True
                       break
                   else:
                       strongpair = False
                       rimp2 = True
                       break
    else:
           lines=fn.readlines()
           for str1 in lines:
               if str1.startswith("FRAGTYPE"):
                   list1=str1.split()
                   if int(list1[1].strip()) == 1:
                       ifdimmer = False
                   else:
                       ifdimmer = True
    fn.close()

    print ("rimp2",rimp2, "  strongpair",strongpair, " ribas",ribas)
    print ("basis",basis)
    #return 0

    # generate input files for fragment correlation calculation
    # set BDFTASK to a temparory name
    tmptask=taskname + "_sub" 

    environ["BDFTASK"]  = tmptask
    environ["XESTTASK"] = tmptask

    finp = tmpdir + "/" + tmptask + ".inp"
    #execute compass
    fn=open(finp, "w+")
    fn.write("$compass\n")
    fn.write("Fragpunch\n")
    fn.write("Basis\n")
    fn.write(basis.strip())
    fn.write("\n")
    print (ribas)
    if rimp2 :
       list1=ribas.split()
       fn.write(list1[0])
       fn.write("\n")
       fn.write(list1[1])
    fn.write("\n")
    fn.write("Skeleton\n")
    fn.write("nosymm\n")
    fn.write("$end\n")
    fn.close()

    print ("\n***** Print input file for checking ******")
    shcmd = "cat " + finp   
    re=system(shcmd)
    print ("***** End of input file ****")

    #shcmd = bdfhome + "/bin/compass.x" + "<" + finp
    re=bdf_exec("compass", finp, None, True) #system(shcmd)
    if re !=0:
        print ("bdf module compass failed! STOP...")
        exit()
    re = check_module_result("compass")

    #execute xuanyuan
    fn=open(finp, "w+")
    fn.write("$Xuanyuan\n")
    fn.write("Direct\n")
    fn.write("Schwarz\n")
    fn.write("Maxmem\n")
    fn.write(" 2GW\n")
    fn.write("$End\n")
    fn.close()
    #shcmd = "cat " + finp
    print ("\n***** Print input file for checking ******")
    shcmd = "cat " + finp   
    re=system(shcmd)
    print ("***** End of input file ****")
    #system(shcmd)
    #shcmd = bdfhome + "/bin/xuanyuan.x" + "<" + finp
    re=bdf_exec("xuanyuan", finp, None, True) #system(shcmd)
    if re !=0:
        print ("bdf module xuanyuan failed! STOP...")
        exit()
    re = check_module_result("xuanyuan")

    import sys
    from pprint import pprint
    pprint(sys.path)

    import remnju
    reminfo = remnju.class_remnju()

    #execute mp2/ccsd
    fn=open(finp, "w+")
    fn.write("$tddft\n")
    fn.write("imethod\n")
    fn.write(" 1\n")
 
    fn.write("itda\n")
    str = " " + reminfo.input['itda'] + "\n"
    fn.write(str)

    fn.write("idiag\n")
    fn.write(" 1\n")
#
#    fn.write(" 0\n")

    fn.write("iexit\n")
    lists = reminfo.input["nstates"].split()
    nstatmon = lists[0]
    nstatdim = lists[1]
    if ifdimmer :
       str = nstatdim + "\n"
       fn.write(str)
    else:
       str = nstatmon + "\n"
       fn.write(str)
#       fn.write(" 2\n")

    fn.write("istore\n")
    fn.write(" 1\n")
#hjz-jul10.2017
    if reminfo.input["dft"] != "none":
        fn.write("DFT\n")
        str = " " + reminfo.input["dft"] + "\n"
        fn.write(str)
#zhangy: mpec+cosx
    if reminfo.input["mpec+cosx"] == "yes":
        #fn.write("mpec+cosx\n")
        fn.write("coulpot+cosx\n")

#zhangy: print transition density to cube file
    if ifdimmer :
       if reminfo.input["trdcubed"] == "yes":
          fn.write("trdcube\n")
    else:
       if reminfo.input["trdcubem"] == "yes":
          fn.write("trdcube\n")

    if reminfo.input["cubenxyz"] != "none":
        fn.write("cubenxyz\n")
        str = " " + reminfo.input["cubenxyz"] + "\n"
        fn.write(str)

    if reminfo.input["cubeorigin"] != "none":
        fn.write("cubeorigin\n")
        str = " " + reminfo.input["cubeorigin"] + "\n"
        fn.write(str)

    if reminfo.input["cubetxyz"] != "none":
        fn.write("cubetxyz\n")
        str = " " + reminfo.input["cubetxyz"] + "\n"
        fn.write(str)

    fn.write("grid\n")
    str = " " + reminfo.input["grid"] + "\n"
    fn.write(str)

    fn.write("memjkop\n")  # memory control
    str = " " + reminfo.input["memjkop"] + "\n"
    fn.write(str)
 
    fn.write("aokxc\n")
    fn.write("directgrid\n")
    fn.write("iguess\n")
    fn.write("0 \n")
    fn.write("lmotocmo\n")
    #fn.write("fragtd\n") # fragment TDDFT
    #fn.write("remstore\n")
    fn.write("$End\n")
    fn.close()
    print ("\n***** Print input file for checking ******")
    shcmd = "cat " + finp   
    re=system(shcmd)
    print ("***** End of input file ****")
    #shcmd = bdfhome + "/bin/mp2.x" + "<" + finp
    re=bdf_exec("tddft", finp, None, True) #system(shcmd)
    if re !=0:
        print ("bdf module tddft failed! STOP...")
        exit()
    re = check_module_result("tddft")

    if jfrag == "NULL" :
        f1name=tmpdir + "/" + tmptask + ".tdinfo1"
        f2name=tmpdir + "/" + taskname + ".tdinfo" + "_" + ifrag
    else :
        f1name=tmpdir + "/" + tmptask + ".tdinfo1"
        f2name=tmpdir + "/" + taskname + ".tdinfo" + "_" + ifrag + "_" + jfrag

    # rename file
    print ("\n")
    str1 = " Fragment TDDFT : rename file " + f1name + " to" + f2name
    print (str1)
    shcmd = "mv" + " " + f1name + " " + f2name
    re=system(shcmd)
    print ("\n")
 
    # add content of tmptask.datapunch to taskname.datapunch
    #srcfile = workdir + "/" + tmptask  + ".datapunch"
    #disfile = workdir + "/" + taskname + ".datapunch"

    #fns=open(srcfile,"r")
    #fnd=open(disfile,"a")
    #line1=fns.readlines()
    #for str1 in line1:
    #  fnd.write(str1)

    #fnd.close()
    #fns.close()

    # rename data files

    # delete scratch file
    #remove(srcfile)
    srcfile = workdir + "/" + tmptask  + ".chkfil"
    remove(srcfile)

    fname=tmpdir+"/"+tmptask + "*"
    filelist=glob(fname)
    #print filelist
    for fil in filelist:
         remove(fil)

    # set back taskname
    environ["BDFTASK"]  = taskname
    environ["XESTTASK"] = taskname

    print (" Finish fragment tddft energy calculation ..." )
    return 0


#Input for subsystem NMR calculation
def frag_nmr_energy(taskname, posthf, basis, ribas="NULL", rimode=0,ifrag="NULL",jfrag="NULL"):
    from os   import environ, system, remove, rename, path
    from glob import glob
    from bdfutil import check_module_result,bdf_exec

    tmpdir  = environ.get("BDF_TMPDIR")
    workdir = environ.get("BDF_WORKDIR")
    bdfhome = environ.get("BDFHOME")

    punfile = tmpdir  + '/fragpunch'
 
    # rimode 0, none RI 1 RI all 2 RI For pair, 3 RI for strong pair
    rimp2      = False
    strongpair = False
    ifdimmer   = False
    fn         = open(punfile,"r")
    if ribas != "NULL" and rimode > 0:
        if rimode == 1: # all RI
           rimp2 = True
        if rimode > 1:
           rimp2 = False
           lines=fn.readlines()
           for str1 in lines:
               if str1.startswith("FRAGTYPE"):
                   list1=str1.split()
                   if int(list1[1].strip()) == 1:
                       ifdimmer = False
                       break
                   else:
                       ifdimmer = True
                       if rimode == 2:
                           rimp2 = True
                           break
               if str1.startswith("STRONGPAIR"):
                   list1=str1.split()
                   if list1[1].strip() == "TRUE":
                       strongpair = True
                       if rimode == 3:
                           rimp2 = False
                       else:
                           rimp2 = True
                       break
                   else:
                       strongpair = False
                       rimp2 = True
                       break
    else:
           lines=fn.readlines()
           for str1 in lines:
               if str1.startswith("FRAGTYPE"):
                   list1=str1.split()
                   if int(list1[1].strip()) == 1:
                       ifdimmer = False
                   else:
                       ifdimmer = True
    fn.close()

    print ("rimp2",rimp2, "  strongpair",strongpair, " ribas",ribas)
    print ("basis",basis)
    #return 0

    # generate input files for fragment correlation calculation
    # set BDFTASK to a temparory name
    tmptask=taskname + "_sub" 

    environ["BDFTASK"]  = tmptask
    environ["XESTTASK"] = tmptask

    finp = tmpdir + "/" + tmptask + ".inp"
    #execute compass
    fn=open(finp, "w+")
    fn.write("$compass\n")
    fn.write("Fragpunch\n")
    fn.write("Basis\n")
    fn.write(basis.strip())
    fn.write("\n")
    print (ribas)
    if rimp2 :
       list1=ribas.split()
       fn.write(list1[0])
       fn.write("\n")
       fn.write(list1[1])
    fn.write("\n")
    fn.write("Skeleton\n")
    fn.write("nosymm\n")
    fn.write("$end\n")
    fn.close()

    print ("\n***** Print input file for checking ******")
    shcmd = "cat " + finp   
    re=system(shcmd)
    print ("***** End of input file ****")

    #shcmd = bdfhome + "/bin/compass.x" + "<" + finp
    re=bdf_exec("compass", finp, None, True) #system(shcmd)
    if re !=0:
        print ("bdf module compass failed! STOP...")
        exit()
    re = check_module_result("compass")

    #execute xuanyuan
    fn=open(finp, "w+")
    fn.write("$Xuanyuan\n")
    fn.write("Direct\n")
    fn.write("Schwarz\n")
    fn.write("Maxmem\n")
    fn.write(" 2GW\n")
    fn.write("$End\n")
    fn.close()
    #shcmd = "cat " + finp
    print ("\n***** Print input file for checking ******")
    shcmd = "cat " + finp   
    re=system(shcmd)
    print ("***** End of input file ****")
    #system(shcmd)
    #shcmd = bdfhome + "/bin/xuanyuan.x" + "<" + finp
    re=bdf_exec("xuanyuan", finp, None, True) #system(shcmd)
    if re !=0:
        print ("bdf module xuanyuan failed! STOP...")
        exit()
    re = check_module_result("xuanyuan")


    #execute mp2/ccsd
    ftinp = tmpdir + "/" + taskname + ".nmr.inp" 
    if path.exists(finp): 
      shcmd = "cp " + ftinp + " " + finp
      re=system(shcmd)
    else:
        fn=open(finp, "w+")
        fn.write("$nmr\n")
        fn.write("igiao\n")
        fn.write("$End\n")
        fn.close()


    print ("\n***** Print input file for checking ******")
    shcmd = "cat " + finp   
    re=system(shcmd)
    print ("***** End of input file ****")
    #shcmd = bdfhome + "/bin/mp2.x" + "<" + finp
    re=bdf_exec("nmr", finp, None, True) #system(shcmd)
    if re !=0:
        print ("bdf module nmr failed! STOP...")
        exit()
    re = check_module_result("nmr")

    #if jfrag == "NULL" :
    #    f1name=tmpdir + "/" + tmptask + ".tdinfo1"
    #    f2name=tmpdir + "/" + taskname + ".tdinfo" + "_" + ifrag
    #else :
    #    f1name=tmpdir + "/" + tmptask + ".tdinfo1"
    #    f2name=tmpdir + "/" + taskname + ".tdinfo" + "_" + ifrag + "_" + jfrag

    ## rename file
    #print "\n"
    #str1 = " Fragment TDDFT : rename file " + f1name + " to" + f2name
    #print str1
    #shcmd = "mv" + " " + f1name + " " + f2name
    #re=system(shcmd)
    #print "\n"
 
    # add content of tmptask.datapunch to taskname.datapunch
    #srcfile = workdir + "/" + tmptask  + ".datapunch"
    #disfile = workdir + "/" + taskname + ".datapunch"

    #fns=open(srcfile,"r")
    #fnd=open(disfile,"a")
    #line1=fns.readlines()
    #for str1 in line1:
    #  fnd.write(str1)

    #fnd.close()
    #fns.close()

    # rename data files

    # delete scratch file
    #remove(srcfile)
    srcfile = workdir + "/" + tmptask  + ".chkfil"
    remove(srcfile)

    fname=tmpdir+"/"+tmptask + "*"
    filelist=glob(fname)
    #print filelist
    for fil in filelist:
         remove(fil)

    # set back taskname
    environ["BDFTASK"]  = taskname
    environ["XESTTASK"] = taskname

    print (" Finish subsystem nmr calculation of fragment",int(ifrag)-1)
    return 0


