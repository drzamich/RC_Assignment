import Elements
reload(Elements)
from Elements import *


#this creates a data from the input file

def DataInput(fileName):
    MList = []
    LoadList = []
    ElemList = []
    BoundList = []
    NodeList = []
    ff = open(fileName,'r')
    z1 = ff.readline()  #reading line by line the whole line
    z2 = z1.strip()
    z3 = z2.split(',')   #1st element of the splitted list
    while z1<>"":           #
        if z3[0]=="*node":                              # key NODE
            IType = "node"  #define  out item type to the string "node"
        elif z3[0]=="*element":                         # key ELEMENT
            IType = "element"
            for i in z3:
                if i.find("eltype")>-1: eltype=i.split("=")[1]# element type
        elif z3[0]=="*material":                         # key ELEMENT
            IType = "mat"
        elif z3[0]=="*loading":                         # key ELEMENT
            IType = "load"
        elif z3[0]=="*boundary":                         # key ELEMENT
            IType = "bound"
        elif z3[0]=="" or z3[0].find("**")>-1:          # blank line / comment line
            pass
        elif IType<>"":
            if IType=="node":                           # data NODE
                NodeList += [Node( int(z3[0]), float(z3[1]), float(z3[2]), float(z3[3]))]
            elif IType=="element":                      # data ELEMENT
                if eltype=='T1D2':
                    ElemList += [T1D2( int(z3[0]), [int(z3[1]),int(z3[2])], float(z3[3]), NodeList)]
                if eltype=='T2D2':
                    ElemList += [T2D2( int(z3[0]), [int(z3[1]),int(z3[2])], float(z3[3]), NodeList)]
                if eltype=='CPS4':
                    ElemList += [CPS4( int(z3[0]), [int(z3[1]),int(z3[2]),int(z3[3]),int(z3[4])], 1.0, True, NodeList)]
                if eltype =='CPE4':
                    ElemList += [CPS4(int(z3[0]), [int(z3[1]), int(z3[2]), int(z3[3]), int(z3[4])], float(z3[5]), False, NodeList)]
            elif IType=="mat":
                MList += [float(z3[0]), float(z3[1])]  #youngs modulus and poissons ratio
            elif IType=="load":
                i0 = FindIndexByLabel(NodeList, int(z3[0]))
                if i0>-1: LoadList += [[ i0, int(z3[1]), float(z3[2])]]
                else: raise NameError ("load node not found")
            elif IType=="bound":
                i0 = FindIndexByLabel(NodeList, int(z3[0]))
                if i0>-1: BoundList += [[ i0, int(z3[1]), float(z3[2]) ]]
                else: raise NameError ("boundary node not found")
        z1 = ff.readline()  #we go one line forward
        z2 = z1.strip()
        z3 = z2.split(',')
    ff.close()
    return MList,BoundList,LoadList, NodeList, ElemList

def DataOut(fileName, uu ):
    ff = open(fileName,'w')
    for i in xrange(len(uu)):   #it's written in the script how to get stresse
         ff.write('%3i%15.10f\n'%(i, uu[i]))  #its not necessary to create a graph
    ff.close()

def DataOutStresses(fileName, uu ):
    ff = open(fileName,'w')
    for i in xrange(len(uu)/3):   #it's written in the script how to get stresse
        # ff.write('%3i%15.10f\n' % (i, uu[i]))
        ff.write('Element %3i sigma_x: %e\n' % (i, uu[i*3]))
        ff.write('Element %3i sigma_y: %e\n' % (i, uu[i*3+1]))
        ff.write('Element %3i tau_xy:  %e\n' % (i, uu[i*3+2]))
    ff.close()
