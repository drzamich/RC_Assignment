from numpy import zeros, array, double, sqrt


def FindIndexByLabel(NodeList, Key):  # find index from label
    if len(NodeList) == 0: raise NameError("no nodes defined")
    for i in range(len(NodeList)):  # loop over all nodes
        if NodeList[i].Label == Key: return i
    return -1  # no node found


class Node(object):
    def __init__(self, Label, XCo, YCo, ZCo):
        self.Label = Label
        self.XCo = XCo
        self.YCo = YCo
        self.ZCo = ZCo
        self.DofT = set([])  # extendible table (set) for indicators for type of dofs of this node
        self.GlobDofStart = 0  # will be later be filled with first index of global degrees of freedom for this node


class Element(object): #daughter of class object
    """ holds element parameters.
    """

    def __init__(self, TypeVal, nNodVal, DofEVal, IntTVal, nIntVal, nIntLVal, DofTVal, DofNVal, dimVal, Label, InzList,
                 NoList):
        self.Type = TypeVal  # element type
        self.nNod = nNodVal  # number of nodes
        self.dim = dimVal  # material type dimension
        self.DofT = DofTVal  # tuple, type of dof for every node of this element. 1 -> u_x, 2->u_y, 3->u_z, 4->phi around x, 5->phi around y, 6->phi around z
        self.DofN = DofNVal  # tuple, number of dofs ( -> len(DofT[i]) ) for every node of this element -- should be the same for all nodes of a particular element
        self.DofE = DofEVal  # number of dofs for whole element
        self.IntT = IntTVal  # Integration type 0: 1D Gaussian, 1: 2D Gaussian, 2: 3D Gaussian
        self.nInt = nIntVal - 1  # integration order (--> index for points and weights = nInt-1)
        self.nIntL = nIntLVal  # total number of integration points
        self.Label = Label  # element number in input
        self.Inzi = []
        for i in InzList:
            self.Inzi += [FindIndexByLabel(NoList, i)]  # indices of nodes belonging to this element
        i_ = 0
        for i in self.Inzi:
            NoList[i].DofT = NoList[i].DofT.union(
                self.DofT[i_])  # carry dof types from element to nodes/NodeList belonging to this element
            i_ += 1
        self.DofI = zeros((self.nNod, self.DofN[0]), dtype=int)  # table for global indices of dofs of every element
        self.PlSt = None  # later relevant for CPS4 only


class T1D2(Element):  #daughter of class element
    def __init__(self, Label, InzList, AA, NoList):
        Element.__init__(self, "T1D2", 2, 2, 0, 1, 1, (set([1]), set([1])), (1, 1), 1, Label, InzList, NoList)  #here we call the init function of the mother
        self.LL = sqrt((NoList[self.Inzi[1]].XCo - NoList[self.Inzi[0]].XCo) ** 2)  # element length
        self.Geom = AA

    def FormN(self, r, s, t):
        N = array([0.5 * (1 - r), 0.5 * (1 + r)])
        return N

    def FormB(self, r, s, t):
        L = self.LL
        B = array([[-1. / L, 1. / L], [0, 0]])
        return B, L / 2.


class T2D2(Element): #its the daughter class of class 'element'
    def __init__(self, Label, InzList, AA, NoList):
        Element.__init__(self, "T2D2", 2, 4, 0, 1, 1, (set([1, 2]), set([1, 2])), (2, 2), 1, Label, InzList, NoList)
        L = sqrt((NoList[self.Inzi[1]].XCo - NoList[self.Inzi[0]].XCo) ** 2 + (
        NoList[self.Inzi[1]].YCo - NoList[self.Inzi[0]].YCo) ** 2)  # element length
        self.LL = L
        self.Geom = AA
        self.sinA = (NoList[self.Inzi[1]].YCo - NoList[self.Inzi[0]].YCo) / L
        self.cosA = (NoList[self.Inzi[1]].XCo - NoList[self.Inzi[0]].XCo) / L

    def FormN(self, r, s, t):
        N = array([[0.5 * (1 - r), 0., 0.5 * (1 + r), 0.],
                   [0., 0.5 * (1 - r), 0., 0.5 * (1 + r)]])
        return N

    def FormB(self, r, s, t):
        L = self.LL
        B = array([[-self.cosA / L, -self.sinA / L, self.cosA / L, self.sinA / L], [0, 0, 0, 0]])
        return B, L / 2.


class CPS4(Element):
    def __init__(self, Label, InzList, thickness, PlSt, NoList):
        Element.__init__(self, "CPS4", 4, 8, 1, 2, 4, (set([1, 2]), set([1, 2]), set([1, 2]), set([1, 2])),
                         (2, 2, 2, 2), 2, Label, InzList, NoList)
        self.PlSt = PlSt  # flag for plane stress (True->plane stress, False->plane strain)
        self.Geom = thickness
        self.X0 = NoList[self.Inzi[0]].XCo
        self.Y0 = NoList[self.Inzi[0]].YCo
        self.X1 = NoList[self.Inzi[1]].XCo
        self.Y1 = NoList[self.Inzi[1]].YCo
        self.X2 = NoList[self.Inzi[2]].XCo
        self.Y2 = NoList[self.Inzi[2]].YCo
        self.X3 = NoList[self.Inzi[3]].XCo
        self.Y3 = NoList[self.Inzi[3]].YCo
        self.AA = 0.5 * (
        (self.X1 - self.X3) * (self.Y2 - self.Y0) + (self.X2 - self.X0) * (self.Y3 - self.Y1))  # element area
        if self.AA <= 0.: raise NameError("Something is wrong with this CPE4-element")

    def FormN(self, r, s, t):
        N = array([[(1 + r) * (1 + s) * 0.25, 0, (1 - r) * (1 + s) * 0.25, 0, (1 - r) * (1 - s) * 0.25, 0,
                    (1 + r) * (1 - s) * 0.25, 0],
                   [0, (1 + r) * (1 + s) * 0.25, 0, (1 - r) * (1 + s) * 0.25, 0, (1 - r) * (1 - s) * 0.25, 0,
                    (1 + r) * (1 - s) * 0.25]])
        return N

    def FormB(self, r, s, t):
        h00 = (1 + s) * 0.25
        h01 = (1 + r) * 0.25
        h10 = -(1 + s) * 0.25
        h11 = (1 - r) * 0.25
        h20 = (-1 + s) * 0.25
        h21 = (-1 + r) * 0.25
        h30 = (1 - s) * 0.25
        h31 = -(1 + r) * 0.25
        JAC00 = h00 * self.X0 + h10 * self.X1 + h20 * self.X2 + h30 * self.X3
        JAC01 = h00 * self.Y0 + h10 * self.Y1 + h20 * self.Y2 + h30 * self.Y3
        JAC10 = h01 * self.X0 + h11 * self.X1 + h21 * self.X2 + h31 * self.X3
        JAC11 = h01 * self.Y0 + h11 * self.Y1 + h21 * self.Y2 + h31 * self.Y3
        det = JAC00 * JAC11 - JAC01 * JAC10
        deti = 1. / det
        a1 = self.Y0 * h01 + self.Y1 * h11 + self.Y2 * h21 + self.Y3 * h31
        a2 = self.Y0 * h00 + self.Y1 * h10 + self.Y2 * h20 + self.Y3 * h30
        b1 = self.X0 * h01 + self.X1 * h11 + self.X2 * h21 + self.X3 * h31
        b2 = self.X0 * h00 + self.X1 * h10 + self.X2 * h20 + self.X3 * h30
        B00 = deti * (h00 * a1 - h01 * a2)
        B10 = deti * (h10 * a1 - h11 * a2)
        B20 = deti * (h20 * a1 - h21 * a2)
        B30 = deti * (h30 * a1 - h31 * a2)
        B01 = deti * (-h00 * b1 + h01 * b2)
        B11 = deti * (-h10 * b1 + h11 * b2)
        B21 = deti * (-h20 * b1 + h21 * b2)
        B31 = deti * (-h30 * b1 + h31 * b2)
        B = array([[B00, 0, B10, 0, B20, 0, B30, 0],
                   [0, B01, 0, B11, 0, B21, 0, B31],
                   [B01, B00, B11, B10, B21, B20, B31, B30]])
        return B, det  #B matrix and Jacobian

