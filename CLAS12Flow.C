#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <string>
#include <sstream>
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

using namespace std;

double Pfunc(double Px, double Py, double Pz)
{
    return sqrt(Px*Px + Py*Py + Pz*Pz);
}

double Efunc(double M, double P)
{
    return sqrt(M * M + P * P);
}

double Ptfunc(double Px, double Py)
{
    return sqrt(Px*Px + Py*Py);
}

TVector2 PtVectfunc(TLorentzVector lv)
{
    TVector2 Pt;
    double Px = lv.Px();
    double Py = lv.Py();
    Pt.SetX(Px);
    Pt.SetY(Py);
    return Pt;
}
class MCParticle
{
    public:
    
    //Lund bank variables
    int pid = 0;
    int id = 0;
    double px = 0;
    double py = 0;
    double pz = 0;
    int daughter = 0;
    int parent = 0;
    double mass = 0;
    double P = 0;
    double E = 0;
    double vz = 0;
    int type = 0;
    //TLorentzVector
    TLorentzVector lv;
    
    //Calculations
    double Pt = 0;
    TVector2 PtVect;
    
    void inputPxPyPzM(double _px, double _py, double _pz, double _m);
    
    void SetParentDaughter(double _parent,double _daughter);
    
    void fillParticle(int _id, int _pid, double _px, double _py, double _pz, int _daughter, int _parent, double _mass, double _vz, int _type);
    
    void Calculate();
    
    void setVectors();
};

void MCParticle::inputPxPyPzM(double _px, double _py, double _pz, double _m)
{
    px = _px;
    py = _py;
    pz = _pz;
    mass = _m;
    
    P = Pfunc(px, py, pz);
    E = Efunc(mass, P);
    
    Pt = Ptfunc(px, py);
    PtVect = PtVectfunc(lv);
    
}
void MCParticle::SetParentDaughter(double _parent,double _daughter)
{
    parent = _parent;
    daughter = _daughter;
}
void MCParticle::Calculate()
{
    P = Pfunc(px, py, pz);
    E = Efunc(mass, P);
    
    lv.SetPxPyPzE(px,py,pz,E);
    
    Pt = Ptfunc(px, py);
    PtVect = PtVectfunc(lv);
}
void MCParticle::setVectors()
{
    lv.SetPxPyPzE(px,py,pz,E);
    
    P = Pfunc(px, py, pz);
    E = Efunc(mass, P);
    
    Pt = Ptfunc(px, py);
    PtVect = PtVectfunc(lv);
}
void MCParticle::fillParticle(int _id, int _pid, double _px, double _py, double _pz, int _daughter, int _parent, double _mass, double _vz, int _type)
{
    id = _id;
    pid = _pid;
    px = _px;
    py = _py;
    pz = _pz;
    daughter = _daughter;
    parent = _parent;
    mass = _mass;
    vz = _vz;
    type = _type;
    
    P = Pfunc(px, py, pz);
    E = Efunc(mass, P);
    
    lv.SetPxPyPzE(px,py,pz,E);
    
    Pt = Ptfunc(px, py);
    PtVect = PtVectfunc(lv);
}

class MultiParticle : public MCParticle
{
    public:
    
    vector<int> v_id;
    vector<int> v_pid;
    vector<double> v_px;
    vector<double> v_py;
    vector<double> v_pz;
    vector<int> v_daughter;
    vector<int> v_parent;
    vector<double> v_mass;
    vector<double> v_vz;
    vector<string> v_name;
    vector<int> v_type;
    
    
    void update(int _id, int _pid, double _px, double _py, double _pz, int _daughter, int _parent, double _mass, double _vz, int _type)
    {
        v_id.push_back(_id);
        v_pid.push_back(_pid);
        v_px.push_back(_px);
        v_py.push_back(_py);
        v_pz.push_back(_pz);
        v_daughter.push_back(_daughter);
        v_parent.push_back(_parent);
        v_mass.push_back(_mass);
        v_vz.push_back(_vz);
        v_type.push_back(_type);
    }
    
};

class Pidi : public MultiParticle
{
    public:
    
    int select_id1 = -999;
    int select_id2 = -999;
    int count = 0;
    bool exist = false;
    
};
class Quark : public MultiParticle
{
    public:
    
    int initial_id = -999;
    int final_id = -999;
};
class Diquark : public MultiParticle
{
    public:
    
    int select_id  = -999;
    
    void diquarkReset()
    {
        v_id.clear();
        v_pid.clear();
        v_px.clear();
        v_py.clear();
        v_pz.clear();
        v_daughter.clear();
        v_parent.clear();
        v_mass.clear();
        v_vz.clear();
    }
};


int CLAS12Flow() {
    gROOT->ProcessLine("#include <vector>");
    
    auto hipofile = "/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV/45nA_job_3301_3.hipo";
    HipoChain chain;
    
    chain.Add(hipofile);
    auto config_c12 = chain.GetC12Reader();
    
    //
    // For dev version, Pi+ Pi- Dihadron endstates are used
    // User can adjust for different endstates using PID codes
    // See 'https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf' for PIDs
    // syntax: config_c12->addExactPid(<PID>,<# of particles>)
    //
    int userpid1; // These are user set values for hadrons they want
    int userpid2;
//    int userpid3;
//    int userpid4;
//    int userpid5;
//    int userpid6;
    
    userpid1 = 221;
    userpid2 = 221;
//    userpid3 = ;
//    userpid4 = ;
//    userpid5 = ;
//    userpid6 = ;
    
    bool do_dihadronflow = true;
    bool jupyterfile = true;
    
    bool pidselect_at_least = true; // If you want at least the hadrons specified
    bool pidselect_exact = false; //If you want exact hadron endstate
    
    std::vector<int> vuserpid;
    
    vuserpid = {userpid1, 
                          userpid2,
/*                          userpid3,
                          userpid4,
                          userpid5,
                          userpid6
*/                          };

    int userpidN = vuserpid.size();
//    config_c12->addExactPid(11,1);    //exactly 1 electron
//    config_c12->addExactPid(userpid1,1);    //exactly 1 pi+
//    config_c12->addExactPid(userpid2,1);    //exactly 1 pi-
//    config_c12->addExactPid(2212,1);    //exactly 1 proton
    
    // Initializations for map
    std::vector<int> vdiquarklist;
    std::vector<int> vquarklist;
    std::vector<int> vhadronlist;
    int vdiquarksize;
    // Defining particle PIDs - look into how to use libraries/dictionaries
    // Vectors will suffice for now
    
    map<int, string> PID_map;
    // Quarks
    PID_map.insert(pair<int,string>(-6, "Anti-Top"));
    PID_map.insert(pair<int,string>(-5, "Anti-Bottom"));
    PID_map.insert(pair<int,string>(-4, "Anti-Charm"));
    PID_map.insert(pair<int,string>(-3, "Anti-Strange"));
    PID_map.insert(pair<int,string>(-2, "Anti-Up"));
    PID_map.insert(pair<int,string>(-1, "Anti-Down"));
    PID_map.insert(pair<int,string>(1, "Down"));
    PID_map.insert(pair<int,string>(2, "Up"));
    PID_map.insert(pair<int,string>(3, "Strange"));
    PID_map.insert(pair<int,string>(4, "Charm"));
    PID_map.insert(pair<int,string>(5, "Bottom"));
    PID_map.insert(pair<int,string>(6, "Top"));
    // Leptons
    PID_map.insert(pair<int,string>(11, "Electron"));
    PID_map.insert(pair<int,string>(13, "Muon"));
    PID_map.insert(pair<int,string>(15, "Tau"));
    PID_map.insert(pair<int,string>(22, "Photon"));
    // Bosons
    PID_map.insert(pair<int,string>(-24, "W-"));
    PID_map.insert(pair<int,string>(24, "W+"));
    PID_map.insert(pair<int,string>(23, "Z"));
    // Baryons
    PID_map.insert(pair<int,string>(-3122, "Anti-Lambda"));
    PID_map.insert(pair<int,string>(-211, "Pi-"));
    PID_map.insert(pair<int,string>(111, "Pi0"));
    PID_map.insert(pair<int,string>(211, "Pi+"));
    PID_map.insert(pair<int,string>(1114, "Delta-"));
    PID_map.insert(pair<int,string>(2114, "Delta0"));
    PID_map.insert(pair<int,string>(2212, "Proton"));
    PID_map.insert(pair<int,string>(2214, "Delta+"));
    PID_map.insert(pair<int,string>(2224, "Delta++"));
    PID_map.insert(pair<int,string>(3112, "Sigma-"));
    PID_map.insert(pair<int,string>(3114, "Sigma*-"));
    PID_map.insert(pair<int,string>(3122, "Lambda"));
    PID_map.insert(pair<int,string>(3214, "Sigma*0"));
    PID_map.insert(pair<int,string>(3222, "Sigma+"));
    PID_map.insert(pair<int,string>(3224, "Sigma*+"));
    PID_map.insert(pair<int,string>(3312, "Xi-"));
    PID_map.insert(pair<int,string>(3324, "Xi*0"));
    // Mesons
    PID_map.insert(pair<int,string>(-323, "K*-"));
    PID_map.insert(pair<int,string>(-313, "K*0"));
    PID_map.insert(pair<int,string>(-213, "Rho-"));
    PID_map.insert(pair<int,string>(113, "Rho0"));
    PID_map.insert(pair<int,string>(213, "Rho+"));
    PID_map.insert(pair<int,string>(221, "Eta"));
    PID_map.insert(pair<int,string>(223, "Omega"));
    PID_map.insert(pair<int,string>(310, "Ks0"));
    PID_map.insert(pair<int,string>(313, "K*0"));
    PID_map.insert(pair<int,string>(323, "K*+"));
    PID_map.insert(pair<int,string>(331, "Eta"));
    PID_map.insert(pair<int,string>(333, "Phi"));
    // Event Gen
    PID_map.insert(pair<int,string>(91, "Gen=91")); //Lund Cluster
    PID_map.insert(pair<int,string>(92, "Gen=92")); //Lund String
    map<int, string>::iterator it;
    
    // Particle list vectors for categorization
    vdiquarklist = {1103, 2101, 2103, 2203, 3101, 3103, 3201, 3203, 3303, 4101, 4103, 4201, 4203, 4301, 4303, 4403, 5101, 5103, 5201, 5203, 5301, 5303, 5401, 5403, 5503};
    vquarklist = {-6,-5,-4,-3,-2,-1,1,2,3,4,5,6};
    vhadronlist = {-3122, -211, 111, 211, 1114, 2114, 2212, 2214, 2224, 3112, 3114, 3122, 3214, 3222, 3224, 3312, 3324, -323, -313, -213, 113, 213, 221, 223, 310, 313, 323, 331, 333};
    
    // Diquarks - we don't usually care what they are, but we can label them all Diquark this way
    vdiquarksize = vdiquarklist.size();

    for (int i = 0; i < vdiquarksize; i++) {
        PID_map.insert(pair<int,string>(vdiquarklist[i], "Diquark"));
    }
    
    // Initializing variables needed from MC::Lund
    int pid;
    int id;
    int daughter;
    int parent;
    int type;
    
    double px; 
    double py; 
    double pz; 
    double mass; 
    double vz; 
    double P;
    double E;
    
    int event_count = -1;
    
    // Add MC::Lund bank for taking Lund data
    auto idx_MCLund= config_c12->addBank("MC::Lund");
    
    // Add needed bank items
    auto iPid=config_c12->getBankOrder(idx_MCLund,"pid");           // I THINK THIS IS ALL UNUSED
    auto idaughter=config_c12->getBankOrder(idx_MCLund,"daughter");
    auto iparent=config_c12->getBankOrder(idx_MCLund,"parent");
    
    
    // Name variables for storing hadron names from map
    std::vector<string> vhadronname;
    std::vector<string> vmhadronname;
    

    
    // Selection
    bool good_dihadronflow = false;
    bool pid_exact = false;
    bool pid_at_least = true;
    
    // Flowchart file creation
    string inputfile;
    string outputfile = "flowdiagram.py";
    if(jupyterfile == true) {
        inputfile = "FlowDev/jupyterdiagram.py";
    } else {inputfile = "FlowDev/pythondiagram.py";}
    std::ifstream  src(inputfile, std::ios::binary);
    std::ofstream  dst(outputfile,   std::ios::binary);
    dst << src.rdbuf();
    dst.close();
    
    // Flowchart variables
    char equalshadron [] = " = Compute(";
    char equalsmhadron [] = " = Aurora(";
    char equalsfphoton [] = " = Neptune(";
    char endparan [] = ")\n";
    char finalst_cluster [] = "with Cluster(\"Final State\"):";
    char final_hadron_cluster [] = "with Cluster(\"Hadrons\"):";
    char mhadron_cluster [] = "with Cluster(\"Hadrons\"):";
    char notmeasured_cluster [] = "with Cluster(\"Not Measured\"):";
    char final_photon_cluster [] = "with Cluster(\"Decay Photon Pairs\"):";
    char post_collision_cluster [] = "with Cluster(\"Post Collision\"):";
    char connect_right [] = " >> ";
    char connect_left [] = " << ";
    char connect [] = " - ";
    char endline [] = "\n";
    char tab [] = "    ";
    char quote [] = "\"";

    
    // This comes from AnalysisWithExtraBanks.C in CLAS12root documentation
    auto& c12 = chain.C12ref();
    
    // Loop over all events in Hipo file
    while (chain.Next()==true && good_dihadronflow == false) {
        event_count += 1;
        if (c12->getDetParticles().empty())
            continue;
        MCParticle electron;
        MCParticle proton;
        MCParticle Lund;

        Pidi photon;

        Quark quark;

        Pidi diquark;
        
        Pidi MidHadron;
        Pidi EndHadron;
        
        // Hadron name vector
        vhadronname.clear();
        // Loop over MC::Lund entries in this event using index -> ID = idx_MCLund
        for (auto imc = 0; imc < c12->getBank(idx_MCLund)->getRows(); imc++) {
            auto mcparticles = c12->mcparts();
            id = mcparticles->getIndex(imc);
            pid = mcparticles->getPid(imc);
            px = mcparticles->getPx(imc);
            py = mcparticles->getPy(imc);
            pz = mcparticles->getPz(imc);
            daughter = mcparticles->getDaughter(imc);
            parent = mcparticles->getParent(imc);
            mass = mcparticles->getMass(imc);
            P = Pfunc(px,py,pz);
            E = Efunc(mass,P);
            vz = mcparticles->getVz(imc);
            type = mcparticles->getType(imc);
            
            if(pid==11 && parent==1){
                electron.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz, type);
                electron.setVectors();
            }
            //all quarks
            else if(std::count(vquarklist.begin(), vquarklist.end(), pid)){
                quark.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz, type);
                quark.update(id, pid, px, py, pz, daughter, parent, 
                              mass, vz, type);
                quark.v_name.push_back(" ");
            }
            //MCParticle
            else if(pid==92 || pid == 91){
                Lund.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz, type);
                Lund.setVectors();
            }
            
            //Diquark
            else if(std::count(vdiquarklist.begin(), vdiquarklist.end(), pid)){
                diquark.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz, type);
                diquark.update(id, pid, px, py, pz, daughter, parent, 
                              mass, vz, type);
                diquark.exist = true;
            }
            //Photon
            if(pid == 22 && type == 1){
                photon.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz, type);
                photon.update(id, pid, px, py, pz, daughter, parent, 
                              mass, vz, type);
            }
            //Proton target
            else if(id == 2){
                proton.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz, type);
                proton.setVectors();
            }
            // Mid Hadrons
            else if((std::count(vhadronlist.begin(), vhadronlist.end(), pid))
                    && (parent == 2)) {
                MidHadron.fillParticle(id, pid, px, py, pz, 
                                       daughter, parent, mass, vz, type);
                MidHadron.update(id, pid, px, py, pz, daughter, parent, 
                              mass, vz, type);
                MidHadron.v_name.push_back(" ");
            }
            //End Hadrons
            else if((std::count(vhadronlist.begin(), vhadronlist.end(), pid))
                    && (parent != 2)) {
                EndHadron.fillParticle(id, pid, px, py, pz, 
                                       daughter, parent, mass, vz, type);
                EndHadron.update(id, pid, px, py, pz, daughter, parent, 
                              mass, vz, type);
                EndHadron.v_name.push_back(" ");
            }
        }
        
        //Selecting pions that come from Lund particle
        for(int i = 0; i < EndHadron.v_id.size(); i++) {
            if(EndHadron.v_parent[i] == Lund.id && EndHadron.v_pid[i] == userpid1) {
                EndHadron.select_id1 = i;
                EndHadron.count += 1;
            }
            if(EndHadron.v_parent[i] == Lund.id && EndHadron.v_pid[i] == userpid2) {
                EndHadron.select_id2 = i;
                EndHadron.count += 1;
            }
        }
        
        //Selecting final quark
        for(int i = 0; i < quark.v_id.size(); i++) {
            if(quark.v_parent[i] == 0) {
                quark.final_id = i;
            }
        }
        
        //Selecting diquark
        for(int i = 0; i < diquark.v_id.size(); i++) {
            if(diquark.v_parent[i] == 2) {
                diquark.select_id1 = i;
            }
        }

        //NEED TO SELECT OTHER HADRONS?
        if(pidselect_at_least == true) {
            for(int i = 0; i < userpidN; i++) {
                if(EndHadron.select_id1 != -999 && EndHadron.select_id2 != -999) {
                    pid_at_least = true;
                }
                else {
                    pid_at_least = false;
                    break;
                }
            }
        }
        else if(pidselect_exact == true) {
            for(int i = 0; i < EndHadron.v_pid.size(); i++) {
                if(EndHadron.select_id1 != -999 && EndHadron.select_id2 != -999 && EndHadron.count == userpidN) { //checking to see if every hadron is specified by user, and if amount is correct
                    pid_exact = true;
                }
                else {
                    pid_exact = false ;
                    break;
                }
            }
        }
        
        if(pidselect_exact == true && // this is for exact inputs
            pid_exact == true
        ) {good_dihadronflow = true;}
        else if(pidselect_at_least == true && //this is for at least inputs
                pid_at_least == true
               ) {good_dihadronflow = true;}
        else{good_dihadronflow = false;}
        
        if(do_dihadronflow == true && good_dihadronflow == true) {
            
           for(it=PID_map.begin(); it!=PID_map.end(); ++it){

               for(int i = 0; i < EndHadron.v_id.size(); i++) {
                   if(EndHadron.v_pid[i] == it->first) {
                      EndHadron.v_name[i] = it->second;
                    }
                }
                for(int i = 0; i < MidHadron.v_id.size(); i++) {
                   if(MidHadron.v_pid[i] == it->first) {
                      MidHadron.v_name[i] = it->second;
                        }
                    }
                for(int i = 0; i < quark.v_id.size(); i++) {
                    if(quark.v_pid[i] == it->first) {
                        quark.v_name[i] = it->second;
                    }
                }
           }
//            endsthadron_init_names
            
            // Code to create the right number of diagram variable names
            
            // Diagram variable names for endstate hadrons
            std::vector<string> vendsthadron_init_names(EndHadron.v_id.size());
            stringstream vhadronss;
            for(int i = 0; i < EndHadron.v_id.size(); i++) {
                vhadronss << "hadron";
                vhadronss << i;
                vhadronss >> vendsthadron_init_names[i];
                vhadronss.str("");
                vhadronss.clear();
            }
            //Diagram variable names for photons
            std::vector<string> vfphoton_init_names(photon.v_id.size());
            stringstream vfphotonss;
            if(photon.v_id.size() > 0) {
                for(int i = 0; i < photon.v_id.size(); i++) {
                    vfphotonss << "fphoton";
                    vfphotonss << i;
                    vfphotonss >> vfphoton_init_names[i];
                    vfphotonss.str("");
                    vfphotonss.clear();
                }
            }
            std::vector<string> vmhadron_init_names(MidHadron.v_id.size());
            if(MidHadron.v_id.size() > 0) {
                    //Diagram variable names for mhadrons
                
                stringstream vmhadronss;
                for(int i = 0; i < MidHadron.v_id.size(); i++) {
                    vmhadronss << "mhadron";
                    vmhadronss << i;
                    vmhadronss >> vmhadron_init_names[i];
                    vmhadronss.str("");
                    vmhadronss.clear();
                }
            }
            
            

            
            //
            // Code to insert python code into file
            //
            ofstream file(outputfile, ios::app);
            
            //Module for writing when there is a diquark, quark and meson
            if(diquark.exist && MidHadron.v_id.size() > 0) {
                file << tab << tab << "quark = AMI(\"" << quark.v_name[quark.final_id] << quote << endparan << "\n";
                file << tab << tab << "double_quark = AMI(\"";
                if(quark.v_pid[quark.final_id] == 1) {
                    file << "Up, Up\")\n";
                }
                else if(quark.v_pid[quark.final_id] == 2) {
                    file << "Up, Down\")\n";
                }
                file << tab << "Lund = Compute(\"Lund\")" << "\n";
//                file << tab << "quark = AMI(\"" << quark.v_name[quark.final_id] << quote << endparan << "\n";
                file << tab << post_collision_cluster << "\n";
                file << tab << tab << "struck_quark = AMI(\"Struck " << quark.v_name[quark.final_id] << quote << endparan << "\n";
                file << tab << tab << "diquark = Compute(\"diquark " << diquark.v_pid[diquark.select_id1] << quote << endparan << "\n";
                for(int i = 0; i < MidHadron.v_id.size(); i++) {
                    file << tab << tab << vmhadron_init_names[i] << equalsmhadron << quote << MidHadron.v_name[i] << quote << endparan << "\n";
                }
                
                
            }
            // Module for writing when quark and diquark no meson
            else if(diquark.exist && MidHadron.v_id.size() == 0 ) {
                file << tab << tab << "diquark = Compute(\"diquark " << diquark.pid << quote << endparan << "\n";
                file << tab << tab << "quark = AMI(\"" << quark.v_name[quark.final_id] << quote << endparan << "\n";
                file << tab << "Lund = Compute(\"Lund\")" << "\n";
                file << tab << "struck_quark = AMI(\"Struck " << quark.v_name[quark.final_id] << quote << endparan << "\n";
            }
            //MidHadron writing
/*            if(MidHadron.v_id.size() > 0) {
                file << tab << mhadron_cluster << "\n";
                for(int i = 0; i < MidHadron.v_id.size(); i++) {
                    file << tab << tab << vmhadron_init_names[i] << equalsmhadron << quote << MidHadron.v_name[i] << quote << endparan << "\n";
                }
            }*/
            //Final State writing
            file << tab << finalst_cluster << "\n";
            file << tab << tab << final_hadron_cluster << "\n";
            for(int i = 0; i < EndHadron.v_id.size(); i++) {
                file << tab << tab << tab << vendsthadron_init_names[i] << equalshadron << quote << EndHadron.v_name[i] << quote << endparan << "\n";
            }
            
            //photon endstates:
            if(photon.v_id.size() > 0) {
                file << tab << tab << final_photon_cluster << "\n";
                for(int i = 0; i < photon.v_id.size(); i++) {
                    file << tab << tab << tab << vfphoton_init_names[i] << equalsfphoton << quote << PID_map[22] << quote << endparan << "\n";
                }
            }
            file << tab << tab << "scattered_electron = Chime(\"Scattered Electron\")\n";
            //connecting decay particles to their photons
            if(EndHadron.v_id.size() > 0 && photon.v_id.size() > 0) {
                for(int i = 0; i < EndHadron.v_id.size(); i++) {
                    for(int j = 0; j < photon.v_id.size(); j += 2) {
                        if(EndHadron.v_id[i] == photon.v_parent[j]) {
                            file << tab << vendsthadron_init_names[i] << connect_right << vfphoton_init_names[j + 1] << "\n";
                            file << tab << vfphoton_init_names[j + 1] << connect << vfphoton_init_names[j] << "\n";
                        }
                    }
                }
            }
            //connecting decay particles to their decays
            for(int k = 0; k < EndHadron.v_id.size(); k++) {
                
            if (EndHadron.v_type[k] != 1) {
                for(int i = 0; i < EndHadron.v_id.size(); i++) {
                    for(int j = 0; j < EndHadron.v_id.size(); j++) {
                        if(EndHadron.v_id[i] == EndHadron.v_parent[j]) {
                            file << tab << vendsthadron_init_names[i] << connect_right << vendsthadron_init_names[j] << "\n";
                        }
                    }
                }
            break;}
            }
            
            if(diquark.exist && MidHadron.v_id.size() == 0) {
                file << tab << "proton" << connect_right << "diquark" << connect_right << "Lund" << "\n";
                file << tab << "proton" << connect_right << "quark" << connect_right << "collision" << connect_right  << "struck_quark" << connect_right << "Lund" "\n";
            }
            if(diquark.exist && MidHadron.v_id.size() > 0) {
                file << tab << "proton" << connect_right << "double_quark" << "\n";
                file << tab << "proton" << connect_right << "quark" << connect_right << "collision" << connect_right  << "struck_quark" << connect_right << "Lund" "\n";
                for(int i = 0; i < MidHadron.v_id.size(); i++) {
                file << tab << "double_quark" << connect_right << vmhadron_init_names[i] << connect_right << "Lund" << "\n";
                }
                file << tab << "double_quark" << connect_right << "diquark" << connect_right << "Lund" << "\n";
            }
            /*for(int i = 0; i < vendsthadron_init_names.size(); i++) {
                if(EndHadron.v_type[i] != 1){
                    file << tab << "Lund" << connect_right << vendsthadron_init_names[i] << "\n";
                }
            }*/
            for(int k = 0; k < EndHadron.v_id.size(); k++) {
                
            if (EndHadron.v_type[k] != 1) {
                for(int i = 0; i < EndHadron.v_id.size(); i++) {
                    for(int j = 0; j < EndHadron.v_id.size(); j++) {
                        if(EndHadron.v_id[i] == EndHadron.v_parent[j]) {
                            file << tab << "Lund" << connect_right << vendsthadron_init_names[i] << "\n";
                        break;}
                    }
                }
            break;}
            }
            file << tab << "electron" << connect_right << "vphoton" << connect_right << "collision" << "\n";
            file << tab << "electron" << connect_right << "scattered_electron" << "\n";
            if(jupyterfile == true) {
                file << "web";
            }
            file.close();
            break;
        }
    }
    return 0;
}