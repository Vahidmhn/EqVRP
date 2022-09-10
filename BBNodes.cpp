#include "BBNodes.h"
#include "Declarations.h"



void BBNode1::Copy(BBNode1& NewNode)
{
    NewNode.DB = DB;
    NewNode.PB = PB;
    NewNode.NotBranchedDem = NotBranchedDem;
    NewNode.Required = Required;
    NewNode.Tabu = Tabu;
    std::copy(RequiredCap,RequiredCap+nV,NewNode.RequiredCap);
    
}


BBNode1::BBNode1(){
    PB = 0;
    DB = 0;    
    RequiredCap = new double[nV];
    for(int v=0; v<nV; v++){
        RequiredCap[v] = 0;
    }
    NotBranchedDem = std::vector<int>(nD,0);
    
}
BBNode1::BBNode1(double pb, double db){
    PB = db;
    DB = pb;
    RequiredCap = new double[nV];
    for(int v=0; v<nV; v++){
        RequiredCap[v] = 0;
    }
    NotBranchedDem =  std::vector<int>(nD,0);
}

BBNode1::~BBNode1(){
    delete[] RequiredCap ;
    NotBranchedDem.clear();
    Required.clear();
    Tabu.clear();
}