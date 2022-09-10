

#ifndef BBNODES_H
#define BBNODES_H

#include "Declarations.h"
#include "MyLib.h"
#include <vector>
#include <map>

class BBNode1{
public:
    double PB;
    double DB;
    double* RequiredCap;
    std::map<int, std::vector <int>> Required;
    std::map<int, std::vector <int>> Tabu;
    std::vector<int> NotBranchedDem;
    
    
    BBNode1();
    BBNode1(double, double);
    virtual ~BBNode1();
    
    
    void Copy(BBNode1& NewNode);
};

#endif /* BBNODES_H */

