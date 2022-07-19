#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <string>
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

using namespace std;

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
    
    userpid1 = 211;
    userpid2 = -211;
//    userpid3 = ;
//    userpid4 = ;
//    userpid5 = ;
//    userpid6 = ;
    
    bool do_dihadronflow = true;
    
    bool pidselect_at_least = false; // If you want at least the hadrons specified
    bool pidselect_exact = true; //If you want exact hadron endstate
    
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
    
    // Initializing counting variables
    int qcount;
    int qusecount;
    int diquarkcount;
    int MCGenindex;
    int elec_count;
    int event_count;
    
    // Add MC::Lund bank for taking Lund data
    auto idx_MCLund= config_c12->addBank("MC::Lund");
    
    // Add needed bank items
    auto iPid=config_c12->getBankOrder(idx_MCLund,"pid");           // I THINK THIS IS ALL UNUSED
    auto idaughter=config_c12->getBankOrder(idx_MCLund,"daughter");
    auto iparent=config_c12->getBankOrder(idx_MCLund,"parent");
    
    // Vectors for storing particle info for each event
    std::vector<int> vpid;
    std::vector<int> vdaughter;
    std::vector<int> vparent;
    
    std::vector<int> vquarkindex;
    std::vector<int> vquarkpid;    
    std::vector<int> vquarkparent;
    std::vector<int> vquarkdaughter;
    
    std::vector<int> vhadronindex;
    std::vector<int> vhadronpid;
    std::vector<int> vhadronparent;
    std::vector<int> vhadrondaughter;
    
    std::vector<int> vmhadronindex;
    std::vector<int> vmhadronpid;
    std::vector<int> vmhadronparent;
    std::vector<int> vmhadrondaughter;
    
    std::vector<int> vdiquarkindex;
    std::vector<int> vdiquarkpid;
    std::vector<int> vdiquarkdaughter;
    std::vector<int> vdiquarkparent;
    
    std::vector<int> vusequarkindex;
    std::vector<int> vusequarkpid;
    std::vector<int> vusequarkdaughter;
    std::vector<int> vusequarkparent;
    
    //making these have 6 0s
    vusequarkindex.push_back(0);
    vusequarkpid.push_back(0);
    vusequarkdaughter.push_back(0);
    vusequarkparent.push_back(0);
    
    vusequarkindex.push_back(0);
    vusequarkpid.push_back(0);
    vusequarkdaughter.push_back(0);
    vusequarkparent.push_back(0);
    
    vusequarkindex.push_back(0);
    vusequarkpid.push_back(0);
    vusequarkdaughter.push_back(0);
    vusequarkparent.push_back(0);
    
    vusequarkindex.push_back(0);
    vusequarkpid.push_back(0);
    vusequarkdaughter.push_back(0);
    vusequarkparent.push_back(0);
    
    vusequarkindex.push_back(0);
    vusequarkpid.push_back(0);
    vusequarkdaughter.push_back(0);
    vusequarkparent.push_back(0);
    
    vusequarkindex.push_back(0);
    vusequarkpid.push_back(0);
    vusequarkdaughter.push_back(0);
    vusequarkparent.push_back(0);
    
    std::vector<string> vhadronname;
    
    std::vector<string> vendsthadron_init_names;
    vendsthadron_init_names.push_back("hadron0");
    vendsthadron_init_names.push_back("hadron1");
    vendsthadron_init_names.push_back("hadron2");
    vendsthadron_init_names.push_back("hadron3");
    vendsthadron_init_names.push_back("hadron4");
    vendsthadron_init_names.push_back("hadron5");
    vendsthadron_init_names.push_back("hadron6");
    vendsthadron_init_names.push_back("hadron7");
    vendsthadron_init_names.push_back("hadron8");
    
    // MC Gen data
    int MCGenparent;
    int MCGendaughter;
    int MCGenpid;
    
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
    
    // Selection
    bool good_dihadronflow = false;
    bool pid_exact = false;
    bool pid_at_least = false;
    
    // Flowchart file creation
    std::ifstream  src("FlowDev/test.py", std::ios::binary);
    std::ofstream  dst("flowtest.py",   std::ios::binary);
    dst << src.rdbuf();
    dst.close();
    
    // Flowchart variables
    char equalscompute [] = " = compute(";
    char endparan [] = ")\n";
    char finalst_cluster [] = "with Cluster(\"Final State\"):";
    char connect_right [] = " >> ";
    char connect_left [] = " << ";
    char endline [] = "\n";
    char tab [] = "    ";
    char quote [] = "\"";

    
    // This comes from AnalysisWithExtraBanks.C in CLAS12root documentation
    auto& c12 = chain.C12ref();
    
    // Loop over all events in Hipo file
    while (chain.Next()==true && good_dihadronflow == false) {
        if (c12->getDetParticles().empty())
            continue;
        
        qcount = 0;
        qusecount = 0;
        diquarkcount = 0;
        MCGenindex = 0;
        elec_count = 0;
        event_count += 1;
        
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
        
        // Diquark vector
        vdiquarkindex.clear();
        vdiquarkpid.clear();
        vdiquarkdaughter.clear();
        vdiquarkparent.clear();
        
        // Hadron name vector
        vhadronname.clear();
        
        // Loop over MC::Lund entries in this event using index -> ID = idx_MCLund
        for (auto imc = 0; imc < c12->getBank(idx_MCLund)->getRows(); imc++) {
            auto mcparticles = c12->mcparts();
            id = mcparticles->getIndex(imc);
            pid = mcparticles->getPid(imc);
            daughter = mcparticles->getDaughter(imc);
            parent = mcparticles->getParent(imc);
            type = mcparticles->getType(imc);
            
            // Identifying Lund String
            if (pid==92 or pid==91) { 
                MCGenindex += 1;
                MCGenparent = parent;
                MCGendaughter = daughter;
                MCGenindex = id;
                MCGenpid = pid;
            }
            // Identifying pre-fragment quark candidates
            else if (std::count(vquarklist.begin(), vquarklist.end(), pid)) {
                qcount += 1;
                vquarkindex.push_back(id);
                vquarkparent.push_back(parent);
                vquarkdaughter.push_back(daughter);
                vquarkpid.push_back(pid);
            }
            // Identifying diquark
            else if (std::count(vdiquarklist.begin(), vdiquarklist.end(), pid) && parent == 2) {
                diquarkcount += 1;
                vdiquarkpid.push_back(pid);
                vdiquarkparent.push_back(parent);
                vdiquarkindex.push_back(id);
                vdiquarkdaughter.push_back(daughter);
            }
            // Identifying endstate hadrons
            else if (type == 1 && parent != 2 && std::count(vhadronlist.begin(), vhadronlist.end(), pid)) {
                vhadronpid.push_back(pid);
                vhadronparent.push_back(parent);
                vhadronindex.push_back(id);
                vhadrondaughter.push_back(daughter);
                vhadronname.push_back(" ");
            }
            // Identifying mid-state hadrons
            else if (type != 1 && std::count(vhadronlist.begin(), vhadronlist.end(), pid)) {
                vmhadronpid.push_back(pid);
                vmhadronparent.push_back(parent);
                vmhadronindex.push_back(id);
                vmhadrondaughter.push_back(daughter);
            }
            // Identifying virtual photon
            else if (pid == 22 && parent == 1) {
                vphotonindex = id;
                vphotonpid = pid;
                vphotonparent = parent;
                vphotondaughter = daughter;
            }
            // Identifying Beam Electron
            else if (id == 1) {
                elec_count += 1;
                belectronindex = id;
                belectronpid = pid;
                belectronparent = parent;
                belectrondaughter = daughter;
            }
            // Identifying Scattered Electron
            else if (pid == 11 && parent == 1) {
                elec_count += 1;
                selectronindex = id;
                selectronpid = pid;
                selectronparent = parent;
                selectrondaughter = daughter;
            }
            // Identifying Target Proton
            else if (id == 2) {
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
                vusequarkindex[qusecount] =  vquarkindex[i];
                vusequarkpid[qusecount] = vquarkpid[i];
                vusequarkdaughter[qusecount] = vquarkparent[i];
                vusequarkparent[qusecount] = vquarkdaughter[i];

                qusecount += 1; //counting number of useful quarks that contribute to MC string
          }

        }
        if(pidselect_at_least == true) {
            for(int i = 0; i < userpidN; i++) {
                if(std::count(vhadronpid.begin(), vhadronpid.end(), vuserpid[i])) {
                    pid_at_least = true;
                }
                else {
                    pid_at_least = false;
                    break;
                }
            }
        }
        else if(pidselect_exact == true) {
            for(int i = 0; i < vhadronpid.size(); i++) {
                if(std::count(vuserpid.begin(),vuserpid.end(), vhadronpid[i]) && vhadronpid.size() == userpidN) { //checking to see if every hadron is specified by user, and if amount is correct
                    pid_exact = true;
                }
                else {
                    pid_exact = false ;
                    break;
                }
            }
        }
        
        if(pidselect_exact == true && // this is for exact inputs
            pid_exact == true &&
//            pid_at_least == true &&
            qcount ==2 &&
            MCGenpid == 92 &&
            elec_count == 2
        ) {good_dihadronflow = true;}
        else if(pidselect_at_least == true && //this is for at least inputs
                pid_at_least == true &&
                qcount == 2 && 
                MCGenpid == 92 &&
                elec_count == 2
               ) {good_dihadronflow = true;}
        else{good_dihadronflow = false;}
        
        if(do_dihadronflow == true && good_dihadronflow == true) {
               for(it=PID_map.begin(); it!=PID_map.end(); ++it){
                  if(vhadronpid[0] == it->first) {
                      vhadronname[0] = it->second;
                  }
                  if(vhadronpid[1] == it->first) {
                      vhadronname[1] = it->second;
                  }
               }
//            endsthadron_init_names
            ofstream file("flowtest.py", ios::app);
            file << tab << finalst_cluster << "\n";
            for(int i = 0; i < vhadronpid.size(); i++) {
                file << tab << tab << vendsthadron_init_names[i] << equalscompute << quote << vhadronname[i] << quote << endparan << "\n";
                file << "#hadron id: " << vhadronindex[i] << "\n";
                file << "#hadron parent: " << vhadronparent[i] << "\n";
                file << "#hadron daughter: " << vhadrondaughter[i] << "\n";
            }
            file << "\n" << "#MCGen index: " << MCGenindex << "\n";
            file << "\n" << "#MCGen parent: " << MCGenparent << "\n";
            file << "\n" << "#MCGen pid: " << MCGenpid << "\n";
            file << "\n" << "#MCGen daughter: " << MCGendaughter << "\n";
            file << "#event count: " << event_count << "\n";
            file << "#hadron count: " << vhadronpid.size() << "\n";
            file.close();
//            break;
        }
    }
    return 0;
}