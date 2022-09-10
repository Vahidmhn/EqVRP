#ifndef Label_H
#define Label_H
#include "Declarations.h"
#include <vector>
#include <bitset>
#include "MyLib.h"
class Label
{
public:
    double Cost;
    double Cap ;
    double CP;
    double C;
    double P;
    double Delta;
    int Node;
    int nVisited;
    bool* Visits;
    Label* From;
    std::vector <int> Valid;
    
    Label(int n, int node, double cp);
    Label(int n, int node, double cp, double cost );
    Label(int n, int node, double cp,double c, double cost,
        std::vector <int>& Tabu, int VehEnd);

    
    
    Label(int n, int node, double p, double delta, double c, double cp, double cost,
        std::vector <int> &Tabu, int VehEnd);
    Label(int n, int node, double p, double delta, double c, double cp, double cost, double cap,
     Label* from);
    
    
    virtual ~Label();
    
};

#endif

