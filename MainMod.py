from time import time
from numpy import *
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse.linalg.dsolve import linsolve
from scipy.sparse.linalg import aslinearoperator

import Elements

reload(Elements)  #a little bit tricky, just accept it, as the professor said such
from Elements import *
import InOut

reload(InOut)
from InOut import *

SamplePoints = {}  #this is needed for higher degree elements, in the script there should be some table with the values
SampleWeight = {}
# Sample points and weighting factors for Gauss Quadrature 1D   - two nodes on one "line" - linear shape function
SamplePoints[0, 0, 0] = [0., 0., 0.]
SampleWeight[0, 0, 0] = 2.
SamplePoints[0, 1, 0] = [-0.577350269189626, 0., 0.]  #he mentioned this one
SamplePoints[0, 1, 1] = [0.577350269189626, 0., 0.]
SampleWeight[0, 1, 0] = 1.
SampleWeight[0, 1, 1] = 1.
SamplePoints[0, 2, 0] = [-0.774596669241483, 0., 0.]
SamplePoints[0, 2, 1] = [0., 0., 0.]
SamplePoints[0, 2, 2] = [0.774596669241483, 0., 0.]
SampleWeight[0, 2, 0] = 0.555555555555556
SampleWeight[0, 2, 1] = 0.888888888888888
SampleWeight[0, 2, 2] = 0.555555555555556

# Sample points and weighting factors for Gauss Quadrature 2D
SamplePoints[1, 0, 0] = [0., 0., 0.]
SampleWeight[1, 0, 0] = 4.
SamplePoints[1, 1, 0] = [-0.577350269189626, -0.577350269189626, 0.]
SamplePoints[1, 1, 1] = [-0.577350269189626, 0.577350269189626, 0.]
SamplePoints[1, 1, 2] = [0.577350269189626, -0.577350269189626, 0.]
SamplePoints[1, 1, 3] = [0.577350269189626, 0.577350269189626, 0.]
SampleWeight[1, 1, 0] = 1.
SampleWeight[1, 1, 1] = 1.
SampleWeight[1, 1, 2] = 1.
SampleWeight[1, 1, 3] = 1.


def MatC(PlSt, Emod, nu, dim):  #out of this we have to make it osothropic - this we have to edit for task 2
    if dim == 1:  #for truss or a beam
        MatM = array([[Emod, 0], [0, 0]])
    elif dim == 2:  #same as we did in our excercies last time
        if PlSt:
            MatM = Emod / (1 - nu ** 2) * array([[1, nu, 0], \
                                                 [nu, 1, 0], \
                                                 [0, 0, (1 - nu) / 2]])  # plane stress
        else:
            MatM = Emod * (1 - nu) / ((1 + nu) * (1 - 2 * nu)) * array([[1, nu / (1 - nu), 0],       #this is what we've done in the
                                                                        [nu / (1 - nu), 1, 0],
                                                                        [0, 0, (1 - 2 * nu) / (
                                                                        2 * (1 - nu))]])  # plane strain
    return MatM  #returns the matrix for the material


def AssignGlobalDof(NodeList, ElList):  # assign dof indices to nodes and elements
    Index = 0
    for i in NodeList:  # loop over nodes
        i.GlobDofStart = Index  # assign current global index
        Index += len(i.DofT)  # Node.DofT has been filled during data input / element initialization
    for i in ElList:  # loop over all elements to assign global dof index to element table if element shares this dof (which is not mandatory)
        for j in xrange(i.nNod):  # loop over nodes of element
            set1 = i.DofT[j]  # set of dof types of element node
            Node = NodeList[i.Inzi[j]]  # global node of element
            set2 = Node.DofT  # set of dof types of global node - not necessarily the same as the element node
            iterator1 = iter(set1)  #
            for k in xrange(len(set1)):  # loop over number of dofs of element node
                k0 = iterator1.next()  # dof type of element node
                iterator2 = iter(set2)  # must be here to make next loop start again with 0
                for l in xrange(len(set2)):  # loop over number of dofs of global node
                    l0 = iterator2.next()  # dof type of global node
                    if k0 == l0:  # element node dof found as global node dof
                        i.DofI[j, k] = Node.GlobDofStart + l  # assign global dof index to element table
                        break
    return Index  # total number of dofs in system


def FindGlobalDof(node, dofT):
    iterator = iter(node.DofT)  # iterator for set of dof types of node
    for j in xrange(len(node.DofT)):  # loop over all dofs of node
        j0 = iterator.next()  # dof type
        if j0 == dofT:  # dofT has dof type to be loaded
            break  # dof type equals prescribed dof type
    return node.GlobDofStart + j  # loaded dof global index


if __name__ == "__main__":   #here the program starts
    Name = "input3.txt"      #name of the input file
    # ~ Name = "../Data/inT2D2.txt"
    # ~ Name = "../Data/SimplePlate.txt"
    MList, BoundList, LoadList, NodeList, ElemList = DataInput(Name)  #lists are extracted from this function
    stime = time()
    NE = len(ElemList)
    NN = len(NodeList)
    N = AssignGlobalDof(NodeList, ElemList)  # N total number of dofs in system
    print NE, NN, N
    EE, nu = MList[0], MList[1]
    Sparse = False  #"cheap" python programming
    if Sparse:
        KK = sparse.lil_matrix((N, N))  # sparse stiffness matrix initialization
    else:
        KK = zeros((N, N), dtype=float)
    pp = zeros((N), dtype=float)

    # build stiffness matrix
    for i in ElemList:   #i is the element in the element list
        MatM = MatC(i.PlSt, EE, nu, i.dim)  #what we got from material matrix function is now named MatM
        KL = zeros((i.DofE, i.DofE), dtype=float)
        #
        for j in xrange(i.nIntL):  #more complicated and generalized that in our excercise example
            r = SamplePoints[i.IntT, i.nInt, j][0]
            s = SamplePoints[i.IntT, i.nInt, j][1]
            t = SamplePoints[i.IntT, i.nInt, j][2]
            BB, JJ = i.FormB(r, s, t)
            f = JJ * i.Geom * SampleWeight[i.IntT, i.nInt, j]
            KL = KL + f * dot(transpose(BB), dot(MatM, BB))
        #
        ndof0 = 0                   #local stiffness matrix to global stiffness matrix
        for j0 in xrange(i.nNod):  # assemble rows with loop over node rows -- nNod: number of nodes per element
            ndof1 = 0
            for j1 in xrange(i.nNod):  # assemble columns with loop over column rows
                for k in xrange(i.DofN[j0]):  # loop over row dofs -- DofN: number of degrees of freedom of element node
                    kk = i.DofI[j0, k]  # global row index of dof
                    for l in xrange(i.DofN[j1]):  # loop over column dofs
                        ll = i.DofI[j1, l]  # global column index of dof
                        KK[kk, ll] = KK[kk, ll] + KL[
                            ndof0 + k, ndof1 + l]  # add components of element stiffness matrix to global stiffness matrix
                ndof1 = ndof1 + i.DofN[j1]  # update entry for dof column index
            ndof0 = ndof0 + i.DofN[j0]  # update entry for dof row index

    # fill load vector
    for i in LoadList:
        k = FindGlobalDof(NodeList[i[0]], i[1])
        pp[k] = pp[k] + i[2]  # global load vector, i[2] has load value

    # assign displacement boundary conditions
    for i in BoundList:
        k = FindGlobalDof(NodeList[i[0]], i[1])
        for j in xrange(N):
            pp[j] = pp[j] - KK[j, k] * i[2]
            KK[j, k] = 0.
            KK[k, j] = 0.
        KK[k, k] = 1.
        pp[k] = i[2]

    # solve system
    if Sparse:
        K_LU = linsolve.splu(KK.tocsc(), permc_spec=3)  # triangulization of stiffness matrix
        uu = K_LU.solve(pp)  # solution of K*u=R -> displacement increment
    else:
        uu = linalg.solve(KK, pp)
    print time() - stime
    print uu  #displacement vector that we've calculated and shall be writted in the output file
    DataOut("results.txt", uu)
    if not Sparse:
        plt.matshow(KK)
        plt.grid()
        plt.show()
    print 'finish'