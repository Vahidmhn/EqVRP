#ifndef COLUMN_GEN_H
#define COLUMN_GEN_H
#include "CommonLibs.h"
#include "MyLib.h"
#include "InitColumns.h"
#include "Solve_ESPPRC.h"
#include "BBNodes.h"
#include <map>
#include <vector>
#include "Route.h"

ILOSTLBEGIN;

typedef IloArray<IloNumVarArray> D2NumVarMat;
class Column_Gen1{
public:
    
    std::vector<Route> BestSol;
    double glb, gub, iter;
    
    
    // Constructors
    Column_Gen1();
    
    virtual ~Column_Gen1();
    
    // Methods
    void Solve_Zstar2();
    double PrintBestSol();
    void GetRv(std::map<int,std::vector<Route*>>& NewRv, std::vector<int>& Vstar);
    void KillRv();
private:
    
    IloEnv Env;
    IloModel RMP;
    IloCplex Solver;
    IloObjective Obj;
    D2NumVarMat Xvr;
    IloNumArray Demand_RhS;
    IloNumArray Veh_RhS;
    IloRangeArray Demand_Cons;
    IloRangeArray Veh_Cons;
    
    
    // Algorithm variables ------------------------------------------
    long* nRoute_v;
    double* Demand_Dual;
    double* Veh_Dual;
    double* Reduced_Cost;
    /// List of the routes or columns
    std::map<int,std::vector<Route*> > Rv;
    std::vector<BBNode1*> Tree;
    bool CutsChanged=false;


    // Methods
    int Solve_Master(BBNode1& node);
    void Update_RedCost1(int v);
    double Solve_Int_RMP(std::vector<Route>& BestSol, double UB);

    void AddColumns2Rv( std::vector <Route*>& NewColVector, int v);
    void AddColumns2RMP(std::vector <Route*>& NewColVector, int v);
    
    bool FindBranchedDem_Veh(BBNode1& node, int& Dem, int& Veh);
    void Branch1(BBNode1 &node, int Ind_Dem, int Veh );
    int Add2Tree(BBNode1 *node);
    double UpdateGDB();
    
    bool IsSolInteger();
    void Update_Duals();
    void SaveBestSol();
    
    void FreeTree();
    
    BBNode1* PopNode(int i);

};


class Column_Gen2{
public:
    
    std::vector<Route> BestSol;
    double glb, gub , iter;
    
    // Constructors
    Column_Gen2(double Zbar, std::vector<int>& Vstar);
    Column_Gen2(double Zbar, std::vector<int>& Vstar, std::map<int,
                std::vector<Route*> > &Rv, std::vector<Route> InitSol);
    Column_Gen2(double Zbar, std::vector<int>& Vstar, std::map<int,
                std::vector<Route*> > &Rv, std::vector<Route> InitSol, double GUB);
    
    virtual ~Column_Gen2();
    
    // Methods
    void Solve_NSW();
    double PrintBestSol();
    void GetRv(std::map<int,std::vector<Route*>>& NewRv, std::vector<int>& Vstar);
    void KillRv();
private:
    
    IloEnv Env;
    IloModel RMP;
    IloCplex Solver;
    IloObjective Obj;
    D2NumVarMat Xvr;
    IloNumArray Demand_RhS;
    IloNumArray Veh_RhS;

    IloRangeArray Demand_Cons;
    IloRangeArray Veh_Cons;
    IloRange CompanyCost;
    IloRangeArray HypotenusCut;
    
    
    // Algorithm variables ------------------------------------------
    long* nRoute_v;
    long HypoCutCounter;
    double* Demand_Dual;
    double* Veh_Dual;
    double CompanyCost_Dual;
    std::vector<double> Hypotenus_Dual;
    /// List of the routes or columns
    std::map<int,std::vector<Route*> > Rv;
    std::vector<BBNode1*> Tree;
    std::vector<std::vector<Route>> Xp;
    bool CutsChanged=false;

    // Problem Parameters
    std::vector<int> Vset ;
    std::vector<int> IntSol;

    // Methods
    int Solve_Master(BBNode1& node);
    double Solve_Int_RMP(std::vector<Route>& BestSol, double LB);
    double Solve_MaxMin_RMP(double LB);
//
    void AddColumns2Rv(std::vector <Route*>& NewColVector, int v);
    void AddColumns2RMP(std::vector <Route*>& NewColVector, int v);
    void AddHypoCut();
//    
    bool FindBranchedDem_Veh(BBNode1& node, int& Dem, int& Veh);
    void Branch(BBNode1 &node, int Ind_Dem, int Veh );
    int Add2Tree(BBNode1* node);
    double UpdateGDB();
//    
    bool IsSolInteger();
    void Update_Duals();
    void SaveBestSol();
    
    
    
    void PrintFracSol();
    void FreeTree();
    
    BBNode1* PopNode(int i);

};




#endif /* COLUMN_GEN_H */

