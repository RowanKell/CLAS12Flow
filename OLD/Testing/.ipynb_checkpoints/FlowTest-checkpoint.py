import hipopy.hipopy as hip
import numpy as np

class TVector2:
    def __init__(self, X, Y):
        self.X = X
        self.Y = Y
    def X(self):
        return self.X
    def Y(self):
        return self.Y
class TVector3:
    def __init__(self, X, Y, Z):
        self.X = X
        self.Y = Y
        self.Z = Z
    def X(self):
        return self.X
    def Y(self):
        return self.Y
    def Z(self):
        return self.Z
class TLorentzVector:
    def __init__(self):
        self.px = 0
        self.py = 0
        self.pz = 0
        self.E = 0
        self.M = 0
    def M(self):
        return E * E - Px * Px - Py * Px - Pz * Pz
    def SetPxPyPzE(self, Px, Py, Pz, E):
        self.px = Px
        self.py = Py
        self.pz = Pz
        self.E = E
        self.M = E * E - Px * Px - Py * Py - Pz * Pz
    def Px(self):
        return self.px
    def Py(self):
        return self.py
    def Pz(self):
        return self.pz
    def E(self):
        return E


def Pfunc(Px,  Py, Pz):
    return np.sqrt(Px*Px + Py*Py + Pz*Pz)

def Efunc(M, P):
    return np.sqrt(M * M + P * P)

def Ptfunc(Px, Py):
    return np.sqrt(Px*Px + Py*Py)

def PtVectfunc(lv):
    Px = lv.Px()
    Py = lv.Py()
    Pt = TVector2(Px, Py)
    return Pt

class MCParticle:
    def __init__(self):
#         super(MCParticle, self).__init__()
        #Lund bank variables
        self.pid = 0;
        self.myid = 0;
        self.px = 0;
        self.py = 0;
        self.pz = 0;
        self.daughter = 0;
        self.parent = 0;
        self.mass = 0;
        self.P = 0;
        self.E = 0;
        self.vz = 0;
        self.mytype = 0;
        #TLorentzVector
        self.lv = TLorentzVector();

        #Calculations
        self.Pt = 0;
        self.PtVect = TVector2(0,0);
        
    #fillparticles
    #setvectors
    #update
    def fillParticle(self, _id, _pid, _px, _py, _pz, _daughter, _parent, _mass, _vz, _type):
        self.myid = _id;
        self.pid = _pid;
        self.px = _px;
        self.py = _py;
        self.pz = _pz;
        self.daughter = _daughter;
        self.parent = _parent;
        self.mass = _mass;
        self.vz = _vz;
        self.mytype = _type;

        self.P = Pfunc(self.px, self.py, self.pz);
        self.E = Efunc(self.mass, self.P);

        self.lv.SetPxPyPzE(self.px,self.py,self.pz,self.E);

        self.Pt = Ptfunc(self.px, self.py);
        self.PtVect = PtVectfunc(self.lv);
    def test(self, px):
        self.px = px
        
    def setVectors(self):
        self.lv.SetPxPyPzE(self.px,self.py,self.pz,self.E);

        self.P = Pfunc(self.px, self.py, self.pz);
        self.E = Efunc(self.mass, self.P);

        self.Pt = Ptfunc(self.px, self.py);
        self.PtVect = PtVectfunc(self.lv);
        
class MultiParticle(MCParticle):
    def __init__(self):
        MCParticle.__init__(self)
        
        self.v_id = []
        self.v_pid = []
        self.v_px = []
        self.v_py = []
        self.v_pz = []
        self.v_daughter = []
        self.v_parent = []
        self.v_mass = []
        self.v_vz = []
        self.v_name = []
        self.v_type = []
    
    
    def update(self, _id,  _pid, _px, _py, _pz,  _daughter,  _parent, _mass, _vz,  _type):
        self.v_id.append(_id)
        self.v_pid.append(_pid)
        self.v_px.append(_px)
        self.v_py.append(_py)
        self.v_pz.append(_pz)
        self.v_daughter.append(_daughter)
        self.v_parent.append(_parent)
        self.v_mass.append(_mass)
        self.v_vz.append(_vz)
        self.v_type.append(_type)
    def testin(self, px, py, pz, E):
        self.px = px
        self.py = py
        self.pz = pz
        self.E = E

a = MultiParticle()
a.testin(3, 4, 5, 6)
a.setVectors()
print(a.lv.Px())