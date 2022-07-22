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
    bool jupyterfile = true;
    
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
    
    std::vector<int> vfphotonindex;
    std::vector<int> vfphotonpid;
    std::vector<int> vfphotonparent;
    std::vector<int> vfphotondaughter;
    
    std::vector<int> vmhadronindex;
    std::vector<int> vmhadronpid;
    std::vector<int> vmhadronparent;
    std::vector<int> vmhadrondaughter;
    
    std::vector<int> vdiquarkindex;
    std::vector<int> vdiquarkpid;
    std::vector<int> vdiquarkdaughter;
    std::vector<int> vdiquarkparent;
    
    std::vector<int> vusequarkindex(6);
    std::vector<int> vusequarkpid(6);
    std::vector<int> vusequarkdaughter(6);
    std::vector<int> vusequarkparent(6);
    
    // Name variables for storing hadron names from map
    std::vector<string> vhadronname;
    std::vector<string> vmhadronname;
    
    // Vectors for connecting pions to photons
    std::vector<std::array<int,2>> vphoton_pion;

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
    char mhadron_cluster [] = "with Cluster(\"Hadrons\"):";
    char notmeasured_cluster [] = "with Cluster(\"Not Measured\"):";
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
        
        // Final state photons
        vfphotonpid.clear();
        vfphotonparent.clear();
        vfphotonindex.clear();
        vfphotondaughter.clear();
        
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
            // Identifying final state photons
            else if (type == 1 && pid == 22) {
                vfphotonpid.push_back(pid);
                vfphotonparent.push_back(parent);
                vfphotonindex.push_back(id);
                vfphotondaughter.push_back(daughter);
            }
            // Identifying mid-state hadrons
            else if (type != 1 && id != 2 && std::count(vhadronlist.begin(), vhadronlist.end(), pid)) {
                vmhadronpid.push_back(pid);
                vmhadronparent.push_back(parent);
                vmhadronindex.push_back(id);
                vmhadrondaughter.push_back(daughter);
                vmhadronname.push_back(" ");
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
                   for(int i = 0; i < vhadronpid.size(); i++) {
                       if(vhadronpid[i] == it->first) {
                          vhadronname[i] = it->second;
                      }
                    for(int i = 0; i < vmhadronpid.size(); i++) {
                       if(vmhadronpid[i] == it->first) {
                          vmhadronname[i] = it->second;
                       }
                      }
                    }
               }
//            endsthadron_init_names
            
            // Code to create the right number of diagram variable names
            
            // Diagram variable names for endstate hadrons
            std::vector<string> vendsthadron_init_names(vhadronname.size());
            stringstream vhadronss;
            for(int i = 0; i < vhadronname.size(); i++) {
                vhadronss << "hadron";
                vhadronss << i;
                vhadronss >> vendsthadron_init_names[i];
                vhadronss.str("");
                vhadronss.clear();
            }
            
            //Diagram variable names for photons
            std::vector<string> vfphoton_init_names(vfphotonpid.size());
            stringstream vfphotonss;
            for(int i = 0; i < vfphotonpid.size(); i++) {
                vfphotonss << "fphoton";
                vfphotonss << i;
                vfphotonss >> vfphoton_init_names[i];
                vfphotonss.str("");
                vfphotonss.clear();
            }
            
            //Diagram variable names for mhadrons
            std::vector<string> vmhadron_init_names(vmhadronname.size());
            stringstream vmhadronss;
            for(int i = 0; i < vmhadronname.size(); i++) {
                vmhadronss << "mhadron";
                vmhadronss << i;
                vmhadronss >> vmhadron_init_names[i];
                vmhadronss.str("");
                vmhadronss.clear();
            }
            

            
            //
            // Code to insert python code into file
            //
            ofstream file(outputfile, ios::app);
            file << tab << mhadron_cluster << "\n";
            for(int i = 0; i < vmhadronpid.size(); i++) {
                file << tab << tab << vmhadron_init_names[i] << equalsmhadron << quote << vmhadronname[i] << quote << endparan << "\n";
//               file << "# id: " << vmhadronindex[i] << "\n";
//               file << "# parent: " << vmhadronparent[i] << "\n";
//               file << "# daughter: " << vmhadrondaughter[i] << "\n";
            }
            file << tab << finalst_cluster << "\n";
            for(int i = 0; i < vhadronpid.size(); i++) {
                file << tab << tab << vendsthadron_init_names[i] << equalshadron << quote << vhadronname[i] << quote << endparan << "\n";
//                file << "# id: " << vhadronindex[i] << "\n";
//                file << "# parent: " << vhadronparent[i] << "\n";
//                file << "# daughter: " << vhadrondaughter[i] << "\n";
            }
            for(int i = 0; i < vfphotonpid.size(); i++) {
                file << tab << tab << vfphoton_init_names[i] << equalsfphoton << quote << PID_map[22] << quote << endparan << "\n";
//                file << "# id: " << vfphotonindex[i] << "\n";
//                file << "# parent: " << vfphotonparent[i] << "\n";
//                file << "# daughter: " << vfphotondaughter[i] << "\n";
            }
            file << tab << notmeasured_cluster << "\n";
//            for(int i = 0;) // NEED TO WORK ON THIS
//            file << "\n" << "#MCGen index: " << MCGenindex << "\n";
//            file << "\n" << "#MCGen parent: " << MCGenparent << "\n";
//            file << "\n" << "#MCGen pid: " << MCGenpid << "\n";
//            file << "\n" << "#MCGen daughter: " << MCGendaughter << "\n";
//            file << "#event count: " << event_count << "\n";
//            file << "#hadron count: " << vhadronpid.size() << "\n";
//            file << "\n" << "testing:" << "\n";
//            file << "photon count" << vfphoton_init_names.size() << "\n";
            
            for(int i = 0; i < vmhadronindex.size(); i++) {
                for(int j = 0; j < vfphotonindex.size(); j++) {
                    if(vmhadronindex[i] == vfphotonparent[j]) {
                        file << tab << tab << vmhadron_init_names[i] << connect_right << vfphoton_init_names[j] << "\n";
                    }
                }
            }
            if(jupyterfile == true) {
                file << "web";
            }
            file.close();
        }
    }
    return 0;
}