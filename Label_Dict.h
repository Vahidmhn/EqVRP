#ifndef Label_Dict_H
#define Label_Dict_H


//#include"CommonLibs.h"
#include "Label.h"
#include "MyLib.h"
#include <map>
#include <vector>
#include <iostream>
class Label_Dict
{
public:
    std::map<int, std::vector<Label*>> U;
    std::map<int, std::vector<Label*>> Q;
    std::map<int, std::vector<Label*>> Qp;
//    Label_Dict();
    virtual ~Label_Dict();
    
    long get_size();
    long get_size(int Key);
    long get_size_wo(int Key);
    
    Label* PopLabel(int Node, int nNode);
    Label* PopLabel(int Node, int label, int nNode);
    
    int SelectNode();
    int SelectNode(int End);
    
    // For Zstar Model
    void AddLabel(Label* Sel_Lab, double NewCost, double NewCap, double NewC, double NewCP, double NewP, int NewNode,
                    int NodeIndexInValid, int Veh, int Min_Or_Max);
    
    // For NSW Model
    void AddLabel(Label* Sel_Lab, double NewP, double NewDelta, double NewC, double NewCP,
        double NewCost, double NewCap, int NewNode, int NodeIndexInValid, int Veh, int Min_Or_Max);
    
    
    
    //----------------------------------
    void AddLabel_NSW(Label* Sel_Lab, double NewP, double NewDelta, double* Demand_Dual, double NewC, double NewCP,
        double NewCost, double NewCap, int NewNode, int NodeIndexInValid, int Min_Or_Max,
    double Veh_Dual, double CompanyCost_Dual, double Hypotenus_Dual_Coeff, int Veh);
    
    int Dominates_NSW(Label* Lab1, Label* Lab2,int n, int Dir, double , double, double, int,double*,int);
    
    int Find_iStar(Label*, int, int );
    
    bool HasNeg(double*, Label*,int, int);
    //----------------------------------
    
    int Dominates(Label* Lab1, Label* Lab2, int n, int Veh, int Dir);
    
    
  
    
};

#endif

