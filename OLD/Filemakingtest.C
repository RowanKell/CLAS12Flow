#include <fstream>
#include <iostream>
using namespace std;
int Filemakingtest()
{
    equalscompute = "= compute(";
    endparan = ")";
    finalst_cluster = "with Cluster(\"Final State\"):";
    connect_right = ">>";
    connect_left = "<<";
    
    
    std::ifstream  src("test.py", std::ios::binary);
    std::ofstream  dst("copied.py",   std::ios::binary);

    dst << src.rdbuf();
    dst.close();
    ofstream file("copied.py", ios::app);
    
    char firstcluster [] = 
    file << firstcluster;
    file.close();
    return 0;
}