#include "InitColumns.h"
#include "CommonLibs.h"
#include "Declarations.h"
#include "ortools/constraint_solver/routing.h"
#include "ortools/constraint_solver/routing_enums.pb.h"
#include "ortools/constraint_solver/routing_index_manager.h"
#include "ortools/constraint_solver/routing_parameters.h"

      
void Add_Col_2_List(std::map<int,std::vector<Route*> >& list, Route* col, int v)
{
    if (list[v].size() == 0)
    {
        list[v].push_back(col);
        return;
    }
    for(int i=0; i<list[v].size(); i++)
    {
        if (col->Cost <= list[v].at(i)->Cost)
        {
            if (abs(col->Cost - list[v].at(i)->Cost)<=EqualityEPS && IsEqual(col->Seq, list[v].at(i)->Seq)){
                col->~Route();
                return;
            }else{
                list[v].insert(list[v].begin()+i, col);
                return;
            }
        }
    }
    list[v].push_back(col);
}



void Add_OD_Columns(std::map<int,std::vector<Route*> >& Rv)
{
    for(int v=0; v<nV;v++)
    {
        Route* NewCol = new Route();
        NewCol->Seq.push_back(Sv[v]);
        NewCol->Seq.push_back(Ev[v]);
        NewCol->Cost = 0;
        Add_Col_2_List(Rv,NewCol,v);
    }
}

void Gen_Rand_Cols(std::map<int,std::vector<Route*> >& Rv, int nCols)
{

    std::vector<int> DemList(nD);
    std::vector<int> TempDemList;
    std::vector<int> Sel_DemList;
    int Sel_Veh = 0;
    int Sel_Dem = 0;
    
    for(int i=0; i<nD;i++)
    {
        DemList.at(i) = i+1;
    }
    
    
    for(int i=0; i<nCols;i++)
    {
        
        Sel_Veh = Randi(0,nV-1);
        double Cap = Qv[Sel_Veh];
        TempDemList = DemList;
        for(;;)
        {
            Sel_Dem = Randi(0,TempDemList.size()-1);
            
            
            if (Cap - Qj[TempDemList.at(Sel_Dem)-1]<0 )
            {
                if(Sel_DemList.size()>0){
                    break;
                }else{
                    TempDemList.erase(TempDemList.begin()+Sel_Dem);
                }
            }else{
                Cap -= Qj[TempDemList.at(Sel_Dem)-1];
                Sel_DemList.push_back(TempDemList.at(Sel_Dem));
                TempDemList.erase(TempDemList.begin()+Sel_Dem);
            }
        }
        Route* New_Col = Find_Shortest_Path(Sel_Veh, Sel_DemList);    
        double P = 0;
        bool Key = true;
        for(int i = 0; i<New_Col->Seq.size()-1; i++){
            int node1 = New_Col->Seq.at(i);
            int node2 = New_Col->Seq.at(i+1);
            if(node2>=1 && node2<=nD){
                P += Pj[node2-1];
                if (Pj[node2-1] +Cij[node1*nN+node2]-CPij[node1*nN+node2]<0){
                    Key=false;
                    break;
                }
            }
            
        }
        
        if (Key && New_Col->Cost/(1-BETA) <= LIMv[Sel_Veh] && P+New_Col->Cost*(1-1/(1-BETA))-Cij[Sv[Sel_Veh]*nN+Ev[Sel_Veh]] > 0){
            New_Col->Cost = New_Col->Cost - Cij[Sv[Sel_Veh]*nN+Ev[Sel_Veh]];
            Add_Col_2_List(Rv, New_Col,Sel_Veh);
        }
        Sel_DemList.clear();
    }

}



double Add_Heuristic_Col(std::map<int,std::vector<Route*> >& Rv, std::vector<Route>& BestSol, double GUB)
{
    std::vector<std::vector <int> > Assign(nV);
    std::vector<double> Load_v(nV,0);
    for(int i = 0; i<nD; i++){
        // Find the closest node
        double Min = INF;
        int Selected_v = -1;
        for(int v=0; v<nV; v++){
            
            std::vector<int> TempAssign(Assign.at(v).begin(), Assign.at(v).end());
            TempAssign.push_back(i+1);
            Route* New_Col = Find_Shortest_Path(v, TempAssign);
            
            if(New_Col->Cost/(1-BETA) <= LIMv[v] && Load_v[v]+Qj[i] < Qv[v] && CheckPositiveDemand(New_Col->Seq)){
                New_Col->Cost -= Cij[Sv[v]*nN+Ev[v]];
                
                
                if(New_Col->Cost < Min){
                    Min = New_Col->Cost;
                    Selected_v = v;
                }                
            }
            New_Col->~Route();
        }
        
        if (Selected_v == -1){
            std::cout << "No initial column is found." << std::endl;
            return -1;
        }
        
        Assign.at(Selected_v).push_back(i+1);
        Load_v[Selected_v] += Qj[i]; 
    }
    
    
    
    double UB = 0;
    std::vector<Route> Sol(nV);
    for (int v = 0; v < nV; ++v) {
        if(Assign.at(v).size()>0){
            Route* New_Col = Find_Shortest_Path(v, Assign.at(v));
            New_Col->Cost -= Cij[Sv[v]*nN+Ev[v]];
            
            for(int j = 0; j < New_Col->Seq.size()-1; j++)
            {
                New_Col->P += Cij[New_Col->Seq.at(j)*nN+New_Col->Seq.at(j+1)] 
                            - CPij[New_Col->Seq.at(j)*nN+New_Col->Seq.at(j+1)];
                if(New_Col->Seq.at(j)<=nD && New_Col->Seq.at(j)>= 1){
                    New_Col->P += Pj[New_Col->Seq.at(j)-1]; 
                }
            }
            New_Col->P -= Cij[Sv[v]*nN+Ev[v]];
            
            
            Add_Col_2_List(Rv, New_Col,v);
            UB += New_Col->Cost;
            Sol.at(v).Cost = New_Col->Cost;
            Sol.at(v).Seq = New_Col->Seq;
        }
    }
    
    if(UB<GUB){
        BestSol = Sol;
        GUB = UB;
    }
    
    return GUB;
}

bool CheckPositiveDemand(std::vector <int> &Seq){
    double C = 0;
    double CP = 0;
    double P = 0;
    for(int k = 0; k < Seq.size()-1; k++){
        int i = Seq.at(k);
        int j = Seq.at(k+1);
        C += Cij[i*nN+j];
        CP += CPij[i*nN+j];
        if (j>=1 && j<=nD){
            P += Pj[j-1];
            if (Pj[j-1] - CPij[i*nN+j] + Cij[i*nN+j] < 0){
                return false;
            }
        }
    }
    if (P-CP+C-Cij[Seq.at(0)*nN+Seq.at(Seq.size()-1)] <= 0) {return false;}
    
    
    return true;
}

double Gen_Rand_VRPs(std::map<int,std::vector<Route*> >& Rv, std::vector<Route>& BestSol, double GPB, int Iteration)
{
    std::vector<int> DemList(nD);
    std::vector<int> VehList(nV);
    std::vector<int> Sel_DemList;
    int Sel_Veh = 0;
    int Sel_Dem = 0;
    
    for(int i=0; i<nD;i++)
    {
        DemList.at(i) = i+1;
    }
    for(int i=0; i<nV;i++)
    {
        VehList.at(i) = i;
    }
    
    
    for(int i=0; i<Iteration;i++)
    {
        std::vector<int> TempDemList(DemList.begin(), DemList.end());
        std::vector<int> TempVehList(VehList.begin(), VehList.end());
        std::vector<Route> TempBestSol(nV);
        
        double PB = 0;
        bool Key = true;
        while(TempDemList.size()>0 && TempVehList.size()>0)
        {     
            int Sel_Veh_Ind= Randi(0,TempVehList.size()-1);
            Sel_Veh = TempVehList.at(Sel_Veh_Ind);
            TempVehList.erase(TempVehList.begin()+Sel_Veh_Ind);
            
            double Cap = Qv[Sel_Veh];
            for(;;)
            {
                
                Sel_Dem = Randi(0,TempDemList.size()-1);
                
                if (Cap - Qj[TempDemList.at(Sel_Dem)-1]<0)
                {
                    break;
                }else{
                    Cap -= Qj[TempDemList.at(Sel_Dem)-1];
                    Sel_DemList.push_back(TempDemList.at(Sel_Dem));
                    TempDemList.erase(TempDemList.begin()+Sel_Dem);
                    if (TempDemList.size()==0){break;}
                }
            }
            // Add the Column to the list
            Route* New_Col = Find_Shortest_Path(Sel_Veh, Sel_DemList);
            if (New_Col->Cost/(1-BETA) <= LIMv[Sel_Veh] && CheckPositiveDemand(New_Col->Seq)){
                New_Col->Cost -= Cij[Sv[Sel_Veh]*nN+Ev[Sel_Veh]];
                
                
                for(int j = 0; j < New_Col->Seq.size()-1; j++)
                {
                    New_Col->P += Cij[New_Col->Seq.at(j)*nN+New_Col->Seq.at(j+1)] 
                                - CPij[New_Col->Seq.at(j)*nN+New_Col->Seq.at(j+1)];
                    if(New_Col->Seq.at(j)<=nD && New_Col->Seq.at(j)>= 1){
                        New_Col->P += Pj[New_Col->Seq.at(j)-1]; 
                    }
                }
                New_Col->P -= Cij[Sv[Sel_Veh]*nN+Ev[Sel_Veh]];
                
                
                
                
                TempBestSol.at(Sel_Veh).Seq = New_Col->Seq;

                TempBestSol.at(Sel_Veh).Cost = New_Col->Cost;
                //Update the primal bound
                PB += New_Col->Cost;
            
                Add_Col_2_List(Rv, New_Col,Sel_Veh);
            }else{
                Key = false;
                New_Col->~Route();
            }
            Sel_DemList.clear();
        }
        if(TempDemList.size() == 0 && Key == true && PB < GPB)
        {
            GPB = PB;
            BestSol = TempBestSol;
        }
    }

    return GPB;
}


double Gen_Rand_VRPs2(std::map<int,std::vector<Route*> >& Rv, std::vector<Route>& BestSol, double GPB,
                    std::vector<int>& VehList, int Iteration)
{

    std::map<int, int> V2ind;
    for(int i = 0; i<VehList.size(); i++)
    {
        V2ind[VehList.at(i)] = i;
    }
    std::vector<int> DemList(nD);
    std::vector<int> Sel_DemList;
    int Sel_Veh = 0;
    int Sel_Dem = 0;
    
    for(int i=0; i<nD;i++)
    {
        DemList.at(i) = i+1;
    }
    
    
    for(int i=0; i<Iteration;i++)
    {
        std::vector<int> TempDemList(DemList.begin(), DemList.end());
        std::vector<int> TempVehList(VehList.begin(), VehList.end());
        std::vector<Route> TempBestSol(VehList.size());
        
        double PB = 0;
        bool Key = true;
        while(TempDemList.size()>0 && TempVehList.size()>0)
        {     
            int Sel_Veh_Ind = Randi(0,TempVehList.size()-1);
            Sel_Veh = TempVehList.at(Sel_Veh_Ind);
            TempVehList.erase(TempVehList.begin()+Sel_Veh_Ind);
            
            double Cap = Qv[Sel_Veh];
            for(;;)
            {
                Sel_Dem = Randi(0,TempDemList.size()-1);
                
                if (Cap - Qj[TempDemList.at(Sel_Dem)-1]<0)
                {
                    break;
                }else{
                    Cap -= Qj[TempDemList.at(Sel_Dem)-1];
                    Sel_DemList.push_back(TempDemList.at(Sel_Dem));
                    TempDemList.erase(TempDemList.begin()+Sel_Dem);
                    if (TempDemList.size()==0){break;}
                }
            }
            // Add the Column to the list
            Route* New_Col = Find_Shortest_Path(Sel_Veh, Sel_DemList);
            
            if (New_Col->Cost/(1-BETA) <= LIMv[Sel_Veh]){
                for(int j = 0; j < New_Col->Seq.size()-1; j++)
                {
                    New_Col->P += Cij[New_Col->Seq.at(j)*nN+New_Col->Seq.at(j+1)] 
                                - CPij[New_Col->Seq.at(j)*nN+New_Col->Seq.at(j+1)];
                    if(New_Col->Seq.at(j)<=nD && New_Col->Seq.at(j)>= 1){
                        New_Col->P += Pj[New_Col->Seq.at(j)-1]; 
                    }
                }
                New_Col->P -= Cij[Sv[Sel_Veh]*nN+Ev[Sel_Veh]];
                if(New_Col->P>0){  
                    TempBestSol.at(V2ind[Sel_Veh]).Seq = New_Col->Seq;
                    TempBestSol.at(V2ind[Sel_Veh]).Cost = New_Col->Cost;
                    TempBestSol.at(V2ind[Sel_Veh]).P = New_Col->P;
                    // Update the Primal bound
                    PB += log(New_Col->P);
            
            
                    Add_Col_2_List(Rv, New_Col,Sel_Veh);
                }else{
                    New_Col->~Route();
                    Key = false;
                }
            }else{
                New_Col->~Route();
                Key = false;
            }  
           
            Sel_DemList.clear();
        }
        if(TempDemList.size() == 0 && Key == true && PB > GPB)
        {
            GPB = PB;
            BestSol = TempBestSol;
        }
    }

    return GPB;
}

Route* Find_Shortest_Path(int Veh,  std::vector<int>& Dems){
 
    Dems.insert(Dems.begin(),0);
    Dems.push_back(Ev[Veh]);
    int  n = Dems.size();
//    Print(Dems);
    using namespace operations_research;

    const std::vector<RoutingIndexManager::NodeIndex> start{
        RoutingIndexManager::NodeIndex{0}
    };
    const std::vector<RoutingIndexManager::NodeIndex> end{
        RoutingIndexManager::NodeIndex{n-1}
    };

    // Create Routing Index Manager
    RoutingIndexManager IndManager(n, 1, start, end);

    // Create Routing Model.
    RoutingModel routing(IndManager);

    // Create and register a transit callback.
    const int transit_callback_index = routing.RegisterTransitCallback(
        [&Dems,&IndManager](int64 from_index, int64 to_index) -> int64 {
          // Convert from routing variable Index to distance matrix NodeIndex.
          auto from_node = IndManager.IndexToNode(from_index).value();
          auto to_node = IndManager.IndexToNode(to_index).value();
          if (to_node<=nD && to_node>=1 && 
              Pj[to_node-1]+Cij[Dems.at(from_node)*nN+Dems.at(to_node)]-
              CPij[Dems.at(from_node)*nN+Dems.at(to_node)]<0 ){
              return 99999;
          }
          return round(1000*Cij[Dems.at(from_node)*nN+Dems.at(to_node)]);
        });

    // Define cost of each arc.
    routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index);

  // Add Distance constraint.
//  routing.AddDimension(transit_callback_index, 0, 2000,
//                       /*fix_start_cumul_to_zero=*/true, "Distance");
//  routing.GetMutableDimension("Distance")->SetGlobalSpanCostCoefficient(100);

  // Setting first solution heuristic.
    RoutingSearchParameters searchParameters = DefaultRoutingSearchParameters();
    searchParameters.set_first_solution_strategy(
        FirstSolutionStrategy::PATH_CHEAPEST_ARC);

    // Solve the problem.
    const Assignment* solution = routing.SolveWithParameters(searchParameters);

    Route* Sol = new Route();
    Sol->Seq.push_back(Sv[Veh]);

    int64 index = routing.Start(0);
    while (routing.IsEnd(index) == false) {
      Sol->Seq.push_back(Dems.at(IndManager.IndexToNode(index).value()));
      int i = Sol->Seq.at(Sol->Seq.size()-2);
      int j = Sol->Seq.at(Sol->Seq.size()-1);
      Sol->Cost += Cij[i*nN + j ];      
      index = solution->Value(routing.NextVar(index));
    }
    
    
    Sol->Seq.push_back(Dems.at(IndManager.IndexToNode(index).value()));
    int i = Sol->Seq.at(Sol->Seq.size()-2);
    int j = Sol->Seq.at(Sol->Seq.size()-1);
    Sol->Cost += Cij[i*nN + j ];      
    
    return Sol;
}


double Add_Google_VRP_Cols(std::map<int,std::vector<Route*> >& Rv, std::vector<Route>& BestSol){
    using namespace operations_research;

    std::vector<int> Points;
    Points.push_back(0);
    for(int i=1;i<=nD; i++)
    {
        Points.push_back(i);
    }
    for (int v=0; v<nV;v++){
        Points.push_back(Ev[v]);
    };
    int n = 1+nD+nV; // depot + demands + End points 
       
    std::vector<RoutingIndexManager::NodeIndex> starts;
    for (int v=0; v<nV;v++){
        starts.push_back(RoutingIndexManager::NodeIndex{0});
    };
    std::vector<RoutingIndexManager::NodeIndex> ends;
    for (int v=0; v<nV;v++){
        ends.push_back(RoutingIndexManager::NodeIndex{n-nV+v});
    };

    // Create Routing Index Manager
    RoutingIndexManager IndManager(n, nV, starts, ends);

    // Create Routing Model.
    RoutingModel routing(IndManager);

    
    for (int v=0; v<nV; v++)
    {
        const int transit_callback_index = routing.RegisterTransitCallback(
        [v,&Points,&IndManager](int64 from_index, int64 to_index) -> int64 {
          // Convert from routing variable Index to distance matrix NodeIndex.
          auto from_node = IndManager.IndexToNode(from_index).value();
          auto to_node = IndManager.IndexToNode(to_index).value();
          if (from_node==0 && to_node>nD){
              return 0;
          }else if(from_node==0){
              return round(1000*Cij[Sv[v]*nV+0]+1000*Cij[Points.at(from_node)*nN+Points.at(to_node)]);
          }else{
            return round(1000*Cij[Points.at(from_node)*nN+Points.at(to_node)]);
          }
        });
        routing.SetArcCostEvaluatorOfVehicle(transit_callback_index, v);
    }
    
    
    
    
    
    // Add Time Limit constraint.   
    for (int v=0; v<nV; v++)
    {
        const int transit_callback_index = routing.RegisterTransitCallback(
        [v,&Points,&IndManager](int64 from_index, int64 to_index) -> int64 {
          // Convert from routing variable Index to distance matrix NodeIndex.
          auto from_node = IndManager.IndexToNode(from_index).value();
          auto to_node = IndManager.IndexToNode(to_index).value();
          if (from_node==0 && to_node>nD){
              return 0;
          }else if(from_node==0){
              return round(CPij[Sv[v]*nV+0]+CPij[Points.at(from_node)*nN+Points.at(to_node)]);
          }else{
            return round(CPij[Points.at(from_node)*nN+Points.at(to_node)]);
          }
        });
        
        routing.AddDimension(transit_callback_index, 0, (int64)round(LIMv[v]-2),
                             /*fix_start_cumul_to_zero=*/true, "Time_limit");
    }
    
    const int Demand_Callback_Index = routing.RegisterTransitCallback(
        [&Points,&IndManager](int64 from_index, int64 to_index) -> int64 {
          // Convert from routing variable Index to distance matrix NodeIndex.
          auto from_node = IndManager.IndexToNode(from_index).value();
          auto to_node = IndManager.IndexToNode(to_index).value();
          if (to_node>=1 && to_node<<nD){
              return Qj[to_node-1];
          }else{
            return 0;
          }
    });

  // Add Capacity constraint.
    std::vector<int64> vehicle_capacities;
    for (int v=0; v<nV; v++){
        vehicle_capacities.push_back((int64)Qv[v]);
    }
    routing.AddDimensionWithVehicleCapacity(Demand_Callback_Index, 0, vehicle_capacities,
                         /*fix_start_cumul_to_zero=*/true, "Capacity");
    
    

  // Setting first solution heuristic.
    RoutingSearchParameters searchParameters = DefaultRoutingSearchParameters();
    searchParameters.set_first_solution_strategy(
      FirstSolutionStrategy::PATH_CHEAPEST_ARC);
    searchParameters.set_solution_limit(1);
    
    
//    searchParameters.mutable_time_limit()->set_seconds(30);
//    searchParameters.set_solution_limit(50);

    // Solve the problem.
    const Assignment* solution = routing.SolveWithParameters(searchParameters);

    double TotalCost = 0;

    
    for(int v = 0; v<nV; v++)
    {
        Route* Sol = new Route();
        Sol->Seq.push_back(Sv[v]);
        Sol->Cost = 0;

        int64 index = routing.Start(v);
        while (routing.IsEnd(index) == false) {
          Sol->Seq.push_back(Points.at(IndManager.IndexToNode(index).value()));
          int i = Sol->Seq.at(Sol->Seq.size()-2);
          int j = Sol->Seq.at(Sol->Seq.size()-1);
          
          Sol->Cost += Cij[i*nN + j ];  
          
          
          index = solution->Value(routing.NextVar(index));
        }
        Sol->Seq.push_back(Points.at(IndManager.IndexToNode(index).value()));
        int i = Sol->Seq.at(Sol->Seq.size()-2);
        int j = Sol->Seq.at(Sol->Seq.size()-1);
        if (i==0 && j == Ev[v]){
            Sol->~Route();
            continue;            
        }
        Sol->Cost += Cij[i*nN + j ] - Cij[Sv[v]*nN+Ev[v]];
        Add_Col_2_List(Rv, Sol,v);
        TotalCost += Sol->Cost;
    }
    delete solution;
//    std:: cout<< TotalCost << std::endl;
    return TotalCost;
}
