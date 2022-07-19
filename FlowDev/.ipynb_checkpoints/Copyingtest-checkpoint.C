#include <fstream>
#include <iostream>

int Copyingtest()
{
    std::ifstream  src("test.py", std::ios::binary);
    std::ofstream  dst("copied.py",   std::ios::binary);

    dst << src.rdbuf();
    return 0;
}