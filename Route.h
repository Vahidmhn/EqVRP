
#ifndef ROUTE_H
#define ROUTE_H

#include "CommonLibs.h"
#include <vector>


class Route{
public:
    double Cost;
    double P;
    std::vector<int> Seq;
    
    Route();
    virtual ~Route();
};



#endif /* ROUTE_H */

