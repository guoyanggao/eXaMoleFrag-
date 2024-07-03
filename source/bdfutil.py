"""
#
#  Utils for bdfdrv
#
"""
# formatted ouput of a low triangle matrix
def print_lowtri_mat(label,nn,ltmat,forminput):
   
    lenl = max(len(label),4)

    # check format,
    formin=forminput
    formin=formin.replace(" ","")
    form  =formin.lower()
       
    if "d" in form:
        str0  = formin.replace("d","")
    elif "f" in form:
        str0  = formin.replace("f","")
    
    str0  = str0.replace("%","")
    list1 = str0.split(".")
    leni  = max(int(list1[0]),4)
    str0 = "%"
    if int(list1[0]) < 4:
        if "d" in form:
           formin = "%d" % leni
           formin = str0 + formin + "d  "
        elif "f" in form:
           formin = "%d" % leni
           formin = formin + "." + list1[1]
           formin = str0 + formin + "f  "
    else:
        formin = formin + "  " 

    if leni < 1:
        print("print_ltmat: Print format error ...")
        exit()
    leni   = leni+2
    colum  = int(100/leni) 

    # lable colum    
    strt = ""
    for i in range(0,lenl+2):
         strt = strt + " "
    form0 = "%" + str(lenl)   + "d  "
    form1 = "%" + str(leni-2) + "d  "

    #print colum,formin,form0
    for i in range(1,nn+1,colum):
        if i == 1:
            str1 = label + "  "
        else:
            str1 = strt
        n = min(colum, nn+1-i)
        for j in range(0,n):
            str0 = form1 % int(j+i)
            str1 = str1 + str0
        print(str1)

    for i in range(1, nn+1):
        for j in range(1,i+1, colum):
            if j == 1:
               str1 = form0 % i
            else:
               str1 = strt

            n = min(colum,i+1-j)
            for k in range(0,n):
              str0 = formin % ltmat[i][k+j]
              str1 = str1 + str0
            print(str1)
#end of print_lowtri_mat

# formatted ouput of a low triangle matrix
def print_uptri_mat(label,nn,upmat,formin):
    # for upper tri matrix, we tranform it to lower triangle and pint
    ltmat={}
    for i in range(1,nn+1):
        ltmat[i]={}
        for j in range(1,i+1):
            ltmat[i][j]=upmat[j][i]
    print_lowtri_mat(label,nn,ltmat,formin)
#end of print_upmat



# formatted ouput of a low triangle matrix
def print_ltmat(label,nn,ltmat):

    colum=15

    lent = len(label)
    lenm = lent-4
    strt = ""
    for i in range(0,lenm):
         strt = strt+" "

    str1 = label + strt
    for i in range(1,nn+1):
        str0 = "%4d" % i 
        str1 = str1 + str0
    print(str1)

    for i in range(1, nn+1):
        str1 = "%4d" % i
        str1 = str1+strt
        for j in range(1,i+1, colum):
            if j != 1:
               str1 = "    "
            n = min(colum,i+1-j)
            for k in range(0,n):
              str0 = "%8.3f" % ltmat[i][k+j]
              str1 = str1 + str0
            print(str1)
#end of print_ltmat

