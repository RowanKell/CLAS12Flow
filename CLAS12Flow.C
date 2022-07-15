#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <string>
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

using namespace std;

int CLAS12Flow()
{
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
    
    config_c12->addExactPid(11,1);    //exactly 1 electron
    config_c12->addExactPid(211,1);    //exactly 1 pi+
    config_c12->addExactPid(-211,1);    //exactly 1 pi-
    config_c12->addExactPid(2212,1);    //exactly 1 proton
    
    // Initializations for map
    std::vector<int> vdiquarklist;
    std::vector<int> vquarklist;
    std::vector<int> vhadronlist;
    int vdiquarksize;
    // Defining particle PIDs - look into how to use libraries/dictionaries
    // Vectors will suffice for now
    
    map<int, string> PID_map;
    //Quarks
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
    //Leptons
    PID_map.insert(pair<int,string>(11, "Electron"));
    PID_map.insert(pair<int,string>(13, "Muon"));
    PID_map.insert(pair<int,string>(15, "Tau"));
    PID_map.insert(pair<int,string>(22, "Photon"));
    //Bosons
    PID_map.insert(pair<int,string>(-24, "W-"));
    PID_map.insert(pair<int,string>(24, "W+"));
    PID_map.insert(pair<int,string>(23, "Z"));
    //Baryons
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
    //Mesons
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
    //Event Gen
    PID_map.insert(pair<int,string>(91, "Gen=91")); //Lund Cluster
    PID_map.insert(pair<int,string>(92, "Gen=92")); //Lund String

    //Particle list vectors for categorization
    vdiquarklist = {1103, 2101, 2103, 2203, 3101, 3103, 3201, 3203, 3303, 4101, 4103, 4201, 4203, 4301, 4303, 4403, 5101, 5103, 5201, 5203, 5301, 5303, 5401, 5403, 5503};
    vquarklist = {-6,-5,-4,-3,-2,-1,1,2,3,4,5,6};
    vhadronlist = {-3122, -211, 111, 211, 1114, 2114, 2212, 2214, 2224, 3112, 3114, 3122, 3214, 3222, 3224, 3312, 3324, -323, -313, -213, 113, 213, 221, 223, 310, 313, 323, 331, 333};
    
    //Diquarks - we don't usually care what they are, but we can label them all Diquark this way
    vdiquarksize = vdiquarklist.size();
    for(int i = 0;i < vdiquarksize;i++){
    PID_map.insert(pair<int,string>(vdiquarklist[i], "Diquark"));
    }
    
    // Initializing variables needed from MC::Lund
    int pid;
    int id;
    int daughter;
    int parent;
    int type;
    
    //Initializing counting variables
    int qcount;
    int MCGenindex;
    
    //Add MC::Lund bank for taking Lund data
    auto idx_MCLund= config_c12->addBank("MC::Lund");
    
    //Add needed bank items
    auto iPid=config_c12->getBankOrder(idx_MCLund,"pid");           // I THINK THIS IS ALL UNUSED
    auto idaughter=config_c12->getBankOrder(idx_MCLund,"daughter");
    auto iparent=config_c12->getBankOrder(idx_MCLund,"parent");
    
    //Vectors for storing particle info for each event
    std::vector<float> vpid;
    std::vector<float> vdaughter;
    std::vector<float> vparent;
    
    std::vector<float> vquarkindex;
    std::vector<float> vquarkpid;    
    std::vector<float> vquarkparent;
    std::vector<float> vquarkdaughter;
    
    std::vector<float> vhadronindex;
    std::vector<float> vhadronpid;
    std::vector<float> vhadronparent;
    std::vector<float> vhadrondaughter;
    
    std::vector<float> vmhadronindex;
    std::vector<float> vmhadronpid;
    std::vector<float> vmhadronparent;
    std::vector<float> vmhadrondaughter;
    
    // MC Gen data
    int MCGenparent;
    int MCGendaughter;
    int MCGenpid;
    int MCGenindex;
    
    // Quark data
    int quarkindex;
    int quarkpid;
    int quarkparent;
    int quarkdaughter;
    
    // Diquark data
    int diquarkindex;
    int diquarkpid;
    int diquarkparent;
    int diquarkdaughter;
    
    // Virtual poton data;
    int vphotonindex;
    int vphotonpid;
    int vphotonparent;
    int vphotondaughter;
    
    // Beam Electron;
    int belectronindex;
    int belectronpid;
    int belectronparent;
    int belectrondaughter;
        
    // Scattered Electron
    int selectronindex;
    int selectronpid;
    int selectronparent;
    int selectrondaughter;
    
    // Target
    int targetindex;
    int targetpid;
    int targetparent;
    int targetdaughter;
    
    //This comes from AnalysisWithExtraBanks.C in CLAS12root documentation
    auto& c12 = chain.C12ref();
    
    //Loop over all events in Hipo file
    while(chain.Next()==true){
        if(c12->getDetParticles().empty())
            continue;
        
        qcount = 0;
        MCGenindex = 0;
        
        // General Vectors
        vpid.clear();
        vdaughter.clear();
        vparent.clear();
        
        // Quark vectors
        vquarkindex.clear();
        vquarkpid.clear();
        vquarkparent.clear();
        vquarkdaughter.clear();
        
        // Final State Hadron vectors
        vhadronindex.clear();
        vhadronparent.clear();
        vhadrondaughter.clear();
        vhadronpid.clear();
        
        // Mid-State Hadron Vectors
        vmhadronindex.clear();
        vmhadronparent.clear();
        vmhadronpid.clear();
        vmhadrondaughter.clear();
        
        //Loop over MC::Lund entries in this event using index -> ID = idx_MCLund
        for(auto imc = 0; imc < c12->getBank(idx_MCLund)->getRows(); imc++){
            auto mcparticles = c12->mcparts();
            id = mcparticles->getIndex(imc);
            pid = mcparticles->getPid(imc);
            daughter = mcparticles->getDaughter(imc);
            parent = mcparticles->getParent(imc);
            type = mcparticles->getType(imc);
            
            // Identifying Lund String
            if(pid==92 or pid==91){ 
                MCGenindex += 1;
                MCGenparent = parent;
                MCGendaughter = daughter;
                MCGenindex = id;
                MCGenpid = pid;
            }
            // Identifying pre-fragment quark candidates
            else if(std::count(vquarklist.begin(), vquarklist.end(), pid)){
                qcount += 1;
                vquarkindex.push_back(id);
                vquarkparent.push_back(parent);
                vquarkdaughter.push_back(daughter);
                vquarkpid.push_back(pid);
            }
            // Identifying diquark
            else if(std::count(vdiquarklist.begin(), vdiquarklist.end(), pid) && parent == 2){
                diquarkindex = id;
                diquarkpid = pid;
                diquarkparent = parent;
                diquarkdaughter = daughter;
            }
            //Identifying endstate hadrons
            else if(type == 1 && std::count(vhadronlist.begin(), vqhadronlist.end(), pid)){
                vhadronpid.push_back(pid);
                vhadronparent.push_back(parent);
                vhadronindex.push_back(id);
                vhadrondaughter.push_back(daughter);
            }
            //Identifying mid-state hadrons
            else if(type != 1 && std::count(vhadronlist.begin(), vqhadronlist.end(), pid)){
                vmhadronpid.push_back(pid);
                vmhadronparent.push_back(parent);
                vmhadronindex.push_back(id);
                vmhadrondaughter.push_back(daughter);
            }
            //Identifying virtual photon
            else if(pid == 22 && parent == 1) {
                vphotonindex = id;
                vphotonpid = pid;
                vphotonparent = parent;
                vphotondaughter = daughter;
            }
            //Identifying Beam Electron
            else if(id == 1) {
                belectronindex = id;
                belectronpid = pid;
                belectronparent = parent;
                belectrondaughter = daughter;
            }
            //Identifying Scattered Electron
            else if(pid == 11 && parent == 1) {
                selectronindex = id;
                selectronpid = pid;
                selectronparent = parent;
                selectrondaughter = daughter;
            }
            //Identifying Target Proton
            else if(id == 2) {
                targetindex = id;
                targetpid = pid;
                targetparent = parent;
                targetdaughter = daughter;
            }
        }
        //Identifying pre-fragment quark
        for(int i = 0; i<qcount; i++) 
        {
            if(vquarkindex[i] == MCGenparent){
                quarkindex; =  vquarkindex[i];
                quarkpid; vquarkpid[i];
                quarkparent; vquarkparent[i];
                quarkdaughter; vquarkdaughter[i];
          }
        }
    }
    return 0;
}