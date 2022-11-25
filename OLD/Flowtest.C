//TODO - ADD CONDITION OPTIONS, DIRECTIONS, COLONS

#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <string>
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

using namespace std;

int Flowtest() {
    
    //
    // Defining General flowchart.js tools
    //
    
    // Types
    char op;
    char cond;
    char e;
    char st;
    // Syntax
    char def_arrow;
    char connect_arrow;
    // Names
    char operation;
    char condition;
    char end;
    char start;
    
    
    //
    // Defining variables for particles
    //
    
    // Initializations
    char init_elec;
    char target;
    // Initial Electron
    init_elec = "st1=>start: Electron";
    // Target Proton
    target = "op1=>operation:";
    // Types
    op = "op";
    cond = "cond";
    e = "e";
    st = "st";
    // Syntax
    def_arrow = "=>";
    connect_arrow = "->";
    //Names
    operation = "operation";
    condition = "condition";
    start = "start";
    end = "end";
    return 0;
}