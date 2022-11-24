#Main file for HipopyFlow
#
# Written by Rowan Kelleher
#
# Uses Hipopy package by Matthew McEneaney
# 
#

# MCParticle Class for organizing all particle data from MC::Lund
import hipopy.hipopy as hip
import numpy as np
import shutil


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
        return self.E * self.E - self.px * self.px - self.py * self.px - self.pz * self.pz
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
        return self.E


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

class Pidi(MultiParticle):
    def __init__(self):
        MultiParticle.__init__(self)
        self.select_id1 = -999
        self.select_id2 = -999
        self.count = 0
        self.exist = False
class Quark(MultiParticle):
    def __init__(self):
        MultiParticle.__init__(self)
        self.pair_id = -999
        self.final_id = -999
        
class Diquark(MultiParticle):
    def __init__(self):
        MultiParticle.__init__(self)
        self.select_id  = -999
    
    def diquarkReset():
        self.v_id.clear()
        self.v_pid.clear()
        self.v_px.clear()
        self.v_py.clear()
        self.v_pz.clear()
        self.v_daughter.clear()
        self.v_parent.clear()
        self.v_mass.clear()
        self.v_vz.clear()
        
def main():
    jupyterfile = False;
    userpid1 = 211
    userpid2 = -211
    userpidN = 1 #1 is 2 in 0-index
    
    do_dihadronflow = True;
    jupyterfile = True;
    
    pidselect_at_least = True # If you want at least the hadrons specified
    pidselect_exact = False #If you want exact hadron endstate
    
    
    filename = "/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV/45nA_job_3051_0.hipo"
    file = hip.open(filename,mode='r')
    PID_map = {
        # Quarks
        -6 : "Anti-Top",
        -5 : "Anti-Bottom",
        -4 : "Anti-Charm",
        -3 : "Anti-Strange",
        -2 : "Anti-Up",
        -1 : "Anti-Down",
        1 : "Down",
        2 : "Up",
        3 : "Strange",
        4 : "Charm",
        5 : "Bottom",
        6 : "Top",
        # Leptons
        11 : "Electron",
        13 : "Muon",
        15 : "Tau",
        22 : "Photon",
        # Bosons
        -24 : "W-",
        24 : "W+",
        23 : "Z",
        # Baryons
        -3122 : "Anti-Lambda",
        -211 : "Pi-",
        111 : "Pi0",
        211 : "Pi+",
        1114 : "Delta-",
        2114 : "Delta0",
        2212 : "Proton",
        2214 : "Delta+",
        2224 : "Delta++",
        3112 : "Sigma-",
        3114 : "Sigma*-",
        3122 : "Lambda",
        3214 : "Sigma*0",
        3222 : "Sigma+",
        3224 : "Sigma*+",
        3312 : "Xi-",
        3324 : "Xi*0",
        # Mesons
        -323 : "K*-",
        -313 : "K*0",
        -213 : "Rho-",
        113 : "Rho0",
        213 : "Rho+",
        221 : "Eta",
        223 : "Omega",
        310 : "Ks0",
        313 : "K*0",
        323 : "K*+",
        331 : "Eta",
        333 : "Phi",
        # Event Gen
        91 : "Gen=91",
        92 : "Gen=92",
    }
    # Quarks


    vdiquarklist = [1103, 2101, 2103, 2203, 3101, 3103, 3201, 3203, 3303, 4101, 4103, 4201, 4203, 4301, 4303, \
                    4403, 5101, 5103, 5201, 5203, 5301, 5303, 5401, 5403, 5503];
    vquarklist = [-6,-5,-4,-3,-2,-1,1,2,3,4,5,6];
    vhadronlist = [-3122, -211, 111, 211, 1114, 2114, 2212, 2214, 2224, 3112, 3114, 3122, 3214, 3222, 3224, 3312, \
                   3324, -323, -313, -213, 113, 213, 221, 223, 310, 313, 323, 331, 333];
    for num in vdiquarklist:
        PID_map[num] = "Diquark"
        
    #copy the input file so we can run the code many times

    original = r'FlowDev/pythondiagram.py'
    target = r'FlowDev/copypythondiagram.py'

    shutil.copyfile(original, target)
    #Set base input file
    pyfile = open("FlowDev/copypythondiagram.py", "a")
    
    # Flowchart variables
    equalshadron = " = Compute("
    equalsmhadron = " = Aurora("
    equalsfphoton = " = Neptune("
    endparan = ")\n"
    finalst_cluster = "with Cluster(\"Final State\"):"
    final_hadron_cluster = "with Cluster(\"Hadrons\"):"
    mhadron_cluster = "with Cluster(\"Hadrons\"):"
    notmeasured_cluster = "with Cluster(\"Not Measured\"):"
    final_photon_cluster = "with Cluster(\"Decay Photon Pairs\"):"
    post_collision_cluster = "with Cluster(\"Post Collision\"):"
    connect_right = " >> "
    connect_left = " << "
    connect = " - "
    endline = "\n"
    tab = "    "
    quote = "\""
        
    event_count = 0
    for event in file:
        event_count += 1
        v_id = []
        v_pid = []
        v_px = []
        v_py = []
        v_pz = []
        v_parent = []
        v_daughter = []
        v_vz = []
        v_type = []
        v_id = file.getBytes('MC::Lund','index')
        v_pid = file.getInts('MC::Lund','pid')
        v_px = file.getFloats('MC::Lund','px')
        v_py = file.getFloats('MC::Lund','py')
        v_pz = file.getFloats('MC::Lund','pz')
        v_parent = file.getBytes('MC::Lund','parent')
        v_daughter = file.getBytes('MC::Lund','daughter')
        v_vz = file.getFloats('MC::Lund','vz')
        v_type = file.getBytes('MC::Lund','type')
        v_mass = file.getFloats('MC::Lund','mass')

        electron = MCParticle()
        proton = MCParticle()
        Lund = MCParticle()

        photon = Pidi()

        quark = Quark()

        diquark = Pidi()
        
        MidHadron = Pidi()
        EndHadron = Pidi()
        
        for i in range(len(v_pid)):
            id_ = v_id[i]
            pid = v_pid[i]
            px = v_px[i]
            py = v_py[i]
            pz = v_pz[i]
            parent = v_parent[i]
            daughter = v_daughter[i]
            vz = v_vz[i]
            mass = v_mass[i]
            type_ = v_type[i]
            
            if(pid==11 and parent==1):
                electron.fillParticle(id_, pid, px, py, pz, daughter, parent, mass, vz, type_);
                electron.setVectors();
            #all quarks
            elif(pid in vquarklist):
                quark.fillParticle(id_, pid, px, py, pz, daughter, parent, mass, vz, type_);
                quark.update(id_, pid, px, py, pz, daughter, parent, 
                              mass, vz, type_);
                quark.v_name.append(" ");

            #MCParticle
            elif(pid==92 or pid == 91):
                Lund.fillParticle(id_, pid, px, py, pz, daughter, parent, mass, vz, type_);
                Lund.setVectors();
            
            #Diquark
            elif(pid in vdiquarklist):
                diquark.fillParticle(id_, pid, px, py, pz, daughter, parent, mass, vz, type_);
                diquark.update(id_, pid, px, py, pz, daughter, parent, 
                              mass, vz, type_);
                diquark.exist = True;
            
            #Photon
            if(pid == 22 and type_ == 1):
                photon.fillParticle(id_, pid, px, py, pz, daughter, parent, mass, vz, type_);
                photon.update(id_, pid, px, py, pz, daughter, parent, 
                              mass, vz, type_);
            
            #Proton target
            elif(id == 2):
                proton.fillParticle(id_, pid, px, py, pz, daughter, parent, mass, vz, type_);
                proton.setVectors();
            
            # Mid Hadrons
            elif((pid in vhadronlist) and (parent == 2)):
                MidHadron.fillParticle(id_, pid, px, py, pz, 
                                       daughter, parent, mass, vz, type_);
                MidHadron.update(id_, pid, px, py, pz, daughter, parent, 
                              mass, vz, type_);
                MidHadron.v_name.append(" ");
            
            #End Hadrons
            elif((pid in vhadronlist) and (parent != 2)):
                EndHadron.fillParticle(id_, pid, px, py, pz, 
                                       daughter, parent, mass, vz, type_);
                EndHadron.update(id_, pid, px, py, pz, daughter, parent, 
                              mass, vz, type_);
                EndHadron.v_name.append(" ");
        
        #-------------------------------------
        #Back in event loop, not particle loop
        #-------------------------------------
        #Selecting pions that come from Lund particle
        
        #piit is pion iterator
        piit = 0
        while piit < len(EndHadron.v_id):
            if(EndHadron.v_parent[piit] == Lund.myid and EndHadron.v_pid[piit] == userpid1):
                EndHadron.select_id1 = piit
                EndHadron.count += 1;

            if(EndHadron.v_parent[piit] == Lund.myid and EndHadron.v_pid[piit] == userpid2):
                EndHadron.select_id2 = piit
                EndHadron.count += 1;
            piit += 1
        
        #Selecting final quark and pair quark
        
        #qit is quark iterator
        qit = 0
        while qit < len(quark.v_id):
            if(quark.v_parent[qit] == 0):
                quark.final_id = qit
            elif(quark.v_parent[qit] == 2):
                quark.pair_id = qit
            qit += 1

        
        #Selecting diquark
        
        #dit is diquark iterator
        dit = 0
        while dit < len(diquark.v_id):
            if(diquark.v_parent[dit] == 2):
                diquark.select_id1 = dit
            dit += 1

        #Selecting midhadron
                
        #mit is midhadron iterator
        mit = 0
        if(len(MidHadron.v_id) > 0):
            while mit < len(MidHadron.v_id):
                if(MidHadron.v_parent[mit] == 2):
                    MidHadron.select_id1 = mit
                mit += 1

        
        #Checking for having 2 pions

        if(pidselect_at_least == True):
            if(EndHadron.select_id1 != -999 and EndHadron.select_id2 != -999):
                pid_at_least = True
            else:
                pid_at_least = False
                continue
        elif(pidselect_exact == True):
            # checking to see if every hadron is specified by user, and if amount is correct
            if(EndHadron.select_id1 != -999 and EndHadron.select_id2 != -999 and EndHadron.count == userpidN):
                pid_exact = True
            else:
                pid_exact = False 
                continue

        # this is for exact inputs
        if(pidselect_exact == True and \
            pid_exact == True \
        ):
            good_dihadronflow = True
        #this is for at least inputs
        elif(pidselect_at_least == True and \
                pid_at_least == True \
               ):
            good_dihadronflow = True
        else:
            good_dihadronflow = False
        
        if(do_dihadronflow == True and good_dihadronflow == True):
            #assign names from pid for end state hadrons
            ehadit = 0
            while(ehadit < len(EndHadron.v_id)):
                EndHadron.v_name[ehadit] = PID_map[EndHadron.v_pid[ehadit]]
                ehadit += 1
            #assign names to midhadrons
            mhadit = 0
            while(mhadit < len(MidHadron.v_id)):
                MidHadron.v_name[mhadit] = PID_map[MidHadron.v_pid[mhadit]]
                mhadit += 1
            #assign quark names
            qnameit = 0
            while(qnameit < len(quark.v_id)):
                quark.v_name[qnameit] = PID_map[quark.v_pid[qnameit]]
                qnameit += 1
                
        
            # Code to create the right number of diagram variable names
            
            # Diagram variable names for endstate hadrons
            vendsthadron_init_names = []
            endhadnameit = 0
            while endhadnameit < len(EndHadron.v_id):
                vendsthadron_init_names.append("hadron" + str(endhadnameit))
                endhadnameit += 1
                
            #diagram variable names for final state photons
            vfphoton_init_names = []
            vfphotonit = 0
            while(vfphotonit < len(photon.v_id)):
                vfphoton_init_names.append("fphoton" + str(vfphotonit))
                vfphotonit += 1
            
            vmhadron_init_names = []
            #diagram names for midhadrons
            if(len(MidHadron.v_name) > 0):
                mhadnameit = 0
                while(mhadnameit < len(MidHadron.v_name)):
                    vmhadron_init_names.append("mhadron" + str(mhadnameit))
                    
            #------------------------------------------------------------
            # CPP CODE
            #------------------------------------------------------------
            #Module for writing when there is a diquark, quark and meson
            if(diquark.exist and len(MidHadron.v_id) > 0):
                pyfile.write(tab + tab + "quark = AMI(\"" + quark.v_name[quark.final_id] + quote + endparan + "\n")
                pyfile.write(tab + tab + "double_quark = AMI(\"")
                if(quark.v_pid[quark.final_id] == 1):
                    pyfile.write("Up, Up\")\n")
                elif(quark.v_pid[quark.final_id] == 2):
                    pyfile.write("Up, Down\")\n")
                pyfile.write(tab + "Lund = Compute(\"Lund\")" + "\n")
#                pyfile << tab << "quark = AMI(\"" << quark.v_name[quark.final_id] << quote << endparan << "\n";
                pyfile.write(tab + post_collision_cluster + "\n")
                pyfile.write(tab + tab + "struck_quark = AMI(\"Struck " + quark.v_name[quark.final_id] + quote + endparan + "\n")
                pyfile.write(tab + tab + "diquark = Compute(\"diquark " + str(diquark.v_pid[diquark.select_id1]) + quote + endparan + "\n")
                mwriteit = 0
                while(mwriteit < len(MidHadron.v_id)):
                    pyfile.write(tab + tab + vmhadron_init_names[mwriteit] + equalsmhadron + quote + MidHadron.v_name[mwriteit] + quote + endparan + "\n")
            #Module for writing when quark and diquark no meson
            elif(diquark.exist and len(MidHadron.v_id) == 0 ):
                pyfile.write(tab + tab + "diquark = Compute(\"diquark " + str(diquark.v_pid[diquark.select_id1]) + quote + endparan + "\n")
                pyfile.write(tab + tab + "quark = AMI(\"" + quark.v_name[quark.final_id] + quote + endparan + "\n")
                pyfile.write(tab + "Lund = Compute(\"Lund\")" + "\n")
                pyfile.write(tab + "struck_quark = AMI(\"Struck " + quark.v_name[quark.final_id] + quote + endparan + "\n")
            
            # Module for writing when quark anti-quark pair and baryon
            elif(quark.v_id.size() == 4 and len(MidHadron.v_id) != 0):
                pyfile.write(tab + tab + "quark = AMI(\"" + quark.v_name[quark.final_id] + quote + endparan + "\n")
                pyfile.write(tab + tab + "pairquark = AMI(\"" + quark.v_name[quark.pair_id] + quote + endparan + "\n")
                pyfile.write(tab + tab + "baryon = Compute(\"" + MidHadron.v_name[MidHadron.select_id1] + quote + endparan + "\n")
                pyfile.write(tab + "struck_quark = AMI(\"Struck " + quark.v_name[quark.final_id] + quote + endparan + "\n")
                pyfile.write(tab + "Lund = Compute(\"Lund\")" + "\n")
            #Final State writing
            pyfile.write(tab + finalst_cluster + "\n")
            pyfile.write(tab + tab + final_hadron_cluster + "\n")
            ewriteit = 0
            while(ewriteit < len(EndHadron.v_id)):
                pyfile.write(tab + tab + tab + vendsthadron_init_names[ewriteit] + equalshadron + quote + EndHadron.v_name[ewriteit] + quote + endparan + "\n")
                ewriteit += 1
                

            #photon endstates:
            if(len(photon.v_id) > 0):
                pyfile.write(tab + tab + final_photon_cluster + "\n")
                photonwriteit = 0
                while(photonwriteit < len(photon.v_id)):
                    pyfile.write(tab + tab + tab + vfphoton_init_names[i] + equalsfphoton + quote + PID_map[22] + quote + endparan + "\n")
                    photonwriteit += 1

            pyfile.write(tab + tab + "scattered_electron = Chime(\"Scattered Electron\")\n")
            #connecting decay particles to their photons
            if(len(EndHadron.v_id) > 0 and len(photon.v_id) > 0):
                endhadpit = 0
                while(photonhadit < len(EndHadron.v_id)):
                    photonhit = 0
                    while(hadphotonit < len(photon.v_id)):
                        if(EndHadron.v_id[endhadpit] == photon.v_parent[photonhit]):
                            pyfile.write(tab + vendsthadron_init_names[endhadpit] + connect_right + vfphoton_init_names[photonhit + 1] + "\n")
                            pyfile.write(tab + vfphoton_init_names[photonhit + 1] + connect + vfphoton_init_names[photonhit] + "\n")
                        endhadpit += 1
                        photonhit += 2

            #connecting decay particles to their decays
            decaycheckit = 0
            while(decaycheckit < len(EndHadron.v_id)):
                if(EndHadron.v_type[decaycheckit] != 1):
                    decayit = 0
                    while(decayit < len(EndHadron.v_id)):
                        decayit2 = 0
                        while(decayit2 < len(EndHadron.v_id)):
                            if(EndHadron.v_id[decayit] == EndHadron.v_parent[decayit2]):
                                pyfile.write(tab + vendsthadron_init_names[decayit] + connect_right + vendsthadron_init_names[decayit2] + "\n")
                            decayit2 += 1
                        decayit += 1
                    break
                decaycheckit += 1
            
            if(diquark.exist and len(MidHadron.v_id) == 0):
                pyfile.write(tab + "proton" + connect_right + "diquark" + connect_right + "Lund" + "\n")
                pyfile.write(tab + "proton" + connect_right + "quark" + connect_right + "collision" + connect_right  + "struck_quark" + connect_right + "Lund \n")
            if(diquark.exist and len(MidHadron.v_id) > 0):
                pyfile.write(tab + "proton" + connect_right + "double_quark" + "\n")
                pyfile.write(tab + "proton" + connect_right + "quark" + connect_right + "collision" + connect_right  + "struck_quark" + connect_right + "Lund \n")
                mhadit = 0
                while(mhadit < len(MidHadron.v_id)):
                    pyfile.write(tab + "double_quark" + connect_right + vmhadron_init_names[mhadit] + connect_right + "Lund" + "\n")
                    mhadit += 1
                pyfile.write(tab + "double_quark" + connect_right + "diquark" + connect_right + "Lund" + "\n")

            if(len(quark.v_id) == 4 and len(MidHadron.v_id) != 0):
                pyfile.write(tab + "proton" + connect_right + "quark" + connect_right + "collision" + connect_right + "struck_quark" + connect_right + "Lund \n")
                pyfile.write(tab + "proton" + connect_right + "pairquark" + connect_right + "Lund" + "\n")
                pyfile.write(tab + "proton" + connect_right + "baryon" + connect_right + "Lund" + "\n")
            #connect mid state hadrons to lund
            #Connect all endstate hadrons to lund
            hadlundit = 0
            while(hadlundit < len(EndHadron.v_id)):
                if(EndHadron.v_parent[hadlundit] == Lund.myid):
                    pyfile.write(tab + "Lund" + connect_right + vendsthadron_init_names[hadlundit] + "\n")
                hadlundit += 1
            #--------------------------------------------------------------------------------
            #Broken currently: need to make this connect midstate hadrons to endstate hadrons
            #--------------------------------------------------------------------------------
            midlundit = 0
            while(midlundit < len(MidHadron.v_id)):
                    endlundit = 0
                    while(endlundit < len(EndHadron.v_id)):
                        if(EndHadron.v_id[lundit1] == MidHadron.v_parent[lundit2]):
                            pyfile.write(tab + vmidsthadron_init_names[midlundit] + connect_right + vendsthadron_init_names[endlundit] + "\n")
                        endlundit += 1
                    midlundit += 1
            pyfile.write(tab + "electron" + connect_right + "vphoton" + connect_right + "collision" + "\n")
            pyfile.write(tab + "electron" + connect_right + "scattered_electron" + "\n")
            pyfile.close()
            break
    
main()
