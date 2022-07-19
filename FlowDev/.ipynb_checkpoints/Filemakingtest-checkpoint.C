#include <fstream>
#include <iostream>
using namespace std;
int Filemakingtest()
{
    std::ifstream  src("test.py", std::ios::binary);
    std::ofstream  dst("copied.py",   std::ios::binary);

    dst << src.rdbuf();
    dst.close();
    ofstream file("copied.py", ios::app);
//    file << "\nimport os\nos.environ[\"PATH\"] += os.pathsep + 'C:/Program Files/Graphviz/bin/'\n";
//    char my_with [] = "\nwith Diagram(\"Particle Collision\", show=False) as web:\n";
//    char my_imports [] = "from diagrams import Diagram\nfrom diagrams import Cluster\nfrom diagrams.aws.compute import *\nfrom diagrams.aws.database import RDS\nfrom diagrams.aws.network import ELB, Route53\nfrom diagrams.aws.analytics import *\nfrom diagrams.aws.business import *\n";
//    file << my_imports;
//    file << my_with;
    
    char firstcluster [] = "\n    with Cluster(\"Initial Particles\"):\n        proton = A4B(\"proton\") \n        electron = Chime(\"electron\")\n    collision = ComputeOptimizer(\"collision\")\n    LundString = Compute(\"Lund String\")\n    vphoton = Workmail(\"virtual photon\")";
    file << firstcluster;
    file.close();
    return 0;
}