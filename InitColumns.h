#ifndef INITCOLUMNS_H
#define INITCOLUMNS_H

#include <numeric>
#include <vector>
#include <ilcplex/ilocplex.h>
#include "ortools/linear_solver/linear_expr.h"
#include "ortools/linear_solver/linear_solver.h"

#include "ortools/constraint_solver/routing.h"
#include "ortools/constraint_solver/routing_enums.pb.h"
#include "ortools/constraint_solver/routing_index_manager.h"
#include "ortools/constraint_solver/routing_parameters.h"


#include "Route.h"
#include "MyLib.h"


void Add_Col_2_List(std::map<int,std::vector<Route*> >& list, Route* col, int v);
void Add_OD_Columns(std::map<int,std::vector<Route*> >& Rv);
Route* Find_Shortest_Path(int Veh, std::vector<int>& Dems);

void Gen_Rand_Cols(std::map<int,std::vector<Route*> >& Rv, int nCols);
double Add_Heuristic_Col(std::map<int,std::vector<Route*> >& Rv, std::vector<Route>& BestSol, double GUB);

double Add_Google_VRP_Cols(std::map<int,std::vector<Route*> >& Rv, std::vector<Route>& BestSol);

double Gen_Rand_VRPs(std::map<int,std::vector<Route*> >& Rv, std::vector<Route>& BestSol, double GUB,  int Iteration);
double Gen_Rand_VRPs2(std::map<int,std::vector<Route*> >& Rv,
                    std::vector<Route>& BestSol, double GPB,
                    std::vector<int>& VehList,
                    int Iteration);


bool CheckPositiveDemand(std::vector <int> &Seq);        

#endif /* INITCOLUMNS_H */

