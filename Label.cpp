#include "Label.h"

Label::Label(int n, int node, double cp)
{
    Cost = 0;
    P = 0;
    Node = node;
    CP = cp;
    Delta = 0;
    C = 0;
    Cap = 0;
    nVisited = 0;
    Visits = new bool[n]; 
    for (int i=0; i<n;i++){
        Visits[i] = 0;
    }
}

Label::Label(int n, int node, double cp, double cost)
{
    Cost = cost;
    P = 0;
    Node = node;
    CP = cp;
    Delta = 0;
    C = 0;
    Cap = 0;
    nVisited = 1;
    Visits = new bool[n];
    for (int i=0; i<n;i++){
        Visits[i] = 0;
    }
    
    Visits[0] = 1;
    
}

Label::Label(int n, int node, double cp,double c, double cost,
        std::vector <int> &Tabu, int VehEnd)
{
    Cost = cost;
    P = 0;
    Node = node;
    CP = cp;
    Delta = 0;
    C = c;
    Cap = 0;
    nVisited = 1;
    Visits = new bool[n];
    for (int i=0; i<n;i++){
        Visits[i] = 0;
    }
    
    Visits[0] = 1;
    for (int i=0; i<nD;i++){
        if(!IsInSet(i+1, Tabu)){
            Valid.push_back(i+1);
        }
    }
    Valid.push_back(VehEnd);
    
}


//1
Label::Label(int n, int node, double p, double delta, double c, double cp, double cost,
        std::vector <int> &Tabu, int VehEnd)
{
    Cost = cost;
    P = p;
    Node = node;
    CP = cp;
    Delta = delta;
    C = c;
    Cap = 0;
    nVisited = 1;
    Visits = new bool[n];
    for (int i=0; i<n;i++){
        Visits[i] = 0;
    }
    
    Visits[0] = 1;
    for (int i=0; i<nD;i++){
        if(!IsInSet(i+1, Tabu)){
            Valid.push_back(i+1);
        }
    }
    Valid.push_back(VehEnd);
    
}

Label::Label(int n, int node, double p, double delta, double c, double cp, double cost, double cap,
     Label* from){
    Cost = cost;
    P = p;
    Node = node;
    CP = cp;
    Delta = delta;
    C = c;
    Cap = cap;
    nVisited = from->nVisited;
    From = from;
    
    Visits = new bool[n];
    std::copy(from->Visits,from->Visits+nD+1, Visits);
    
    Valid = from->Valid;
}


Label::~Label()
{
    delete[] Visits;
    Valid.clear();
}

