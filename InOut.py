import Elements
reload(Elements)
from Elements import *

def DataInput(fileName):
    MList = []
    LoadList = []
    ElemList = []
    BoundList = []
    NodeList = []
    ff = open(fileName,'r')
    z1 = ff.readline()
    z2 = z1.strip()gh
    z3 = z2.split(',')
    while z1<>"":           #
        if z3[0]=="*node":                              # key NODE
            IType = "node"
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
            elif IType=="mat":
                MList += [float(z3[0]), float(z3[1])]
            elif IType=="load":
                i0 = FindIndexByLabel(NodeList, int(z3[0]))
                if i0>-1: LoadList += [[ i0, int(z3[1]), float(z3[2])]]
                else: raise NameError ("load node not found")
            elif IType=="bound":
                i0 = FindIndexByLabel(NodeList, int(z3[0]))
                if i0>-1: BoundList += [[ i0, int(z3[1]), float(z3[2]) ]]
                else: raise NameError ("boundary node not found")
        z1 = ff.readline()
        z2 = z1.strip()
        z3 = z2.split(',')
    ff.close()
    return MList,BoundList,LoadList, NodeList, ElemList

def DataOut(fileName, uu ):
    ff = open(fileName,'w')
    for i in xrange(len(uu)):
        ff.write('%7i%9.4f\n'%(i, uu[i]))
    ff.close()
