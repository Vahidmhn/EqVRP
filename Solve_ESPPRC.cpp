#include "Solve_ESPPRC.h"
#include "Declarations.h"
#include "MyLib.h"
#include "BBNodes.h"

// Used in Zstar Model
std::vector<Route*> Solve_ESPPRC(int Veh,double Veh_Dual, double* Demand_Dual,
       std::vector<int>& Tabu , long nRoute, int Min_Or_Max){
        
        
    Label_Dict Dict;
    
    Dict.U[0].push_back(new Label(nD+1, 0, CPij[Sv[Veh]*nN+0], Cij[Sv[Veh]*nN+0], Cij[Sv[Veh]*nN+0]-Veh_Dual,Tabu, Ev[Veh]));        
    
    // Main loop
    while(Dict.get_size_wo(Ev[Veh]) > 0)
    {
        
        // Time commands ******************************************
        if (TimeLim-(get_wall_time()-Time)<0){break;}
        // ********************************************************
        
        int Sel_Node = 0;
        Sel_Node = Dict.SelectNode();
//        Label* Sel_Lab = Dict.PopLabel(Sel_Node,nN);
        Label* Sel_Lab = Dict.U[Sel_Node].at(0);
        

         
        // Extend
        /// Extention never happens to the depot. So potential node are >0
        for(int k=0; k<Sel_Lab->Valid.size(); k++)
        {            
            int i = Sel_Lab->Valid.at(k);
            // new node can only be a demand node or Ev[Veh] when Sel_Node!=0
            if ( Sel_Node == 0 && i == Ev[Veh] ){continue;}
                
            // Rules of potential nodes to extend

            // Then extend
            //Calculate the cost
            double NewCost = 0;
            NewCost = Sel_Lab->Cost + Cij[Sel_Node*nN+i]-(i<=nD ? Demand_Dual[i-1] : Cij[Sv[Veh]*nN+Ev[Veh]]);
            
            // Is the cost feasible?
            if(i == Ev[Veh] && NewCost*Min_Or_Max <= EPS){continue;}//

            // if i is a demand point
            // add the demand to cap
            double NewCap = 0;
            NewCap = Sel_Lab->Cap;
            double NewCP = Sel_Lab->CP;
            NewCP += CPij[Sel_Node*nN+i];
            double NewC = Sel_Lab->C;
            NewC += Cij[Sel_Node*nN+i];
            double NewP = 0;
            NewP = Sel_Lab->P;
                    
            if (i<=nD){ 
                NewCap += Qj[i-1];// it shows the node number not demand id
                NewP += Pj[i-1];
                if (Pj[i-1] - (CPij[Sel_Node*nN+i] - Cij[Sel_Node*nN+i]) < 0){continue;}
            }else{
                if (NewP + (NewC-NewCP) - Cij[Sv[Veh]*nN+Ev[Veh]] < 0){continue;}
            }
            
            if (NewCap<=Qv[Veh] && NewCP <= LIMv[Veh] ){
                // Add the new label to the Dict                
                Dict.AddLabel(Sel_Lab, NewCost, NewCap, NewC, NewCP, NewP, i, k, Veh, Min_Or_Max);
            }
                
             
        }
        
//        Sel_Lab->~Label();
        Dict.Q[Sel_Node].push_back(Sel_Lab);
        Dict.U[Sel_Node].erase(Dict.U[Sel_Node].begin());

    }
    
    
    std::vector<Route*> Returned_Route_Vector;
    long nReturnedRoute = nRoute;
    if (nReturnedRoute<0){
        nReturnedRoute = Dict.U[Ev[Veh]].size();
    }else{
        nReturnedRoute = std::min(nReturnedRoute, Dict.get_size(Ev[Veh]));
    }
     
    for (int i=0; i<nReturnedRoute; i++){
        Route *NewRoute = new Route();
        NewRoute->Seq.push_back(Ev[Veh]);
        Label* From = Dict.U[Ev[Veh]].at(0)->From;
        while(true){            
            if (From->Node == 0)
            {
                NewRoute->Seq.insert(NewRoute->Seq.begin(),0);
                NewRoute->Seq.insert(NewRoute->Seq.begin(),Sv[Veh]);
                break;
            }
            NewRoute->Seq.insert(NewRoute->Seq.begin(),From->Node);
            From = From->From;
        }
        
        for(int j = 0; j < NewRoute->Seq.size()-1 ;j++)
        {
            NewRoute->Cost += Cij[NewRoute->Seq.at(j)*nN + NewRoute->Seq.at(j+1)];
        }
        NewRoute->Cost -= Cij[Sv[Veh]*nN + Ev[Veh]];
        NewRoute->P = Dict.U[Ev[Veh]].at(0)->P + Dict.U[Ev[Veh]].at(0)->C - Dict.U[Ev[Veh]].at(0)->CP
                - Cij[Sv[Veh]*nN + Ev[Veh]];
        
        Dict.U[Ev[Veh]].at(0)->~Label();
        Dict.U[Ev[Veh]].erase(Dict.U[Ev[Veh]].begin());     
        
        Returned_Route_Vector.push_back(NewRoute);
    }
    
    
    
    // Free the memory
    for(std::map<int, std::vector<Label*> >::iterator item=Dict.U.begin(); item!=Dict.U.end();item++){
        for(int i = 0; i<item->second.size(); i++){
            Label* temp = item->second.at(i);
            temp->~Label();
            temp = NULL;
        }
        item->second.clear();
    }
    Dict.U.clear();
    for(std::map<int, std::vector<Label*> >::iterator item=Dict.Q.begin(); item!=Dict.Q.end();item++){
        for(int i = 0; i<item->second.size(); i++){
            Label* temp = item->second.at(i);
            temp->~Label();
            temp = NULL;
        }
        item->second.clear();
    }
    Dict.Q.clear();
    return Returned_Route_Vector;
}


// Used in NSW Model
std::vector<Route*> Solve_ESPPRC(int Veh, double* Demand_Dual, double Veh_Dual,
        double CompanyCost_Dual, double Hypotenus_Dual_Coeff, 
        std::vector<int>& Tabu , long nRoute, int Min_Or_Max){
    
    
    Label_Dict Dict;
    
    
    // Parameters from the origin to the depot
    double NewPvr = 0;
    double NewDelta = 0;
    double NewCP = CPij[Sv[Veh]*nN+0];
    double NewC = Cij[Sv[Veh]*nN+0];
    NewPvr = 0 + (NewC-NewCP);
    if (NewPvr > 0){
        Dict.U[0].push_back(new Label(nD+1, 0, 0, NewDelta,  NewC, NewCP, log(NewPvr) - Veh_Dual 
                - NewC * CompanyCost_Dual - NewPvr*Hypotenus_Dual_Coeff, Tabu, Ev[Veh]));
    }else{
        Dict.U[0].push_back(new Label(nD+1, 0, 0, NewDelta,  NewC, NewCP, 
            NewPvr, Tabu, Ev[Veh]));
    }
    

    // Main loop
    while(Dict.get_size_wo(Ev[Veh]) > 0)
    {
        // Time commands ******************************************
        if (TimeLim-(get_wall_time()-Time)<0){break;}
        // ********************************************************
        
        int Sel_Node = 0;
        Sel_Node = Dict.SelectNode();
        Label* Sel_Lab = Dict.U[Sel_Node].at(0);
        
        // Extend
        /// Extention never happens to the depot. So potential node are >0
        for(int k=0; k<Sel_Lab->Valid.size(); k++)
        {            
            int i = Sel_Lab->Valid.at(k);
            // new node can only be a demand node or Ev[Veh] when Sel_Node!=0
            if ( Sel_Node==0 && i == Ev[Veh]){continue;}

            // if i is a demand point
            // add the demand to cap
            double NewCost = 0;
            double NewCap = Sel_Lab->Cap;
            double NewP = Sel_Lab->P;
            double ODcomp = 0;
            
            
            NewDelta = Sel_Lab->Delta;
            NewCP = Sel_Lab->CP + CPij[Sel_Node*nN+i];
            NewC = Sel_Lab->C + Cij[Sel_Node*nN+i]; 
            
            if (i<=nD){ 
                NewCap += Qj[i-1]; // i shows the node number not demand id
                NewP += Pj[i-1];
                NewDelta += Demand_Dual[i-1];
                if ( Pj[i-1] - (CPij[Sel_Node*nN+i] - Cij[Sel_Node*nN+i]) < 0){continue;}
            }else{
                ODcomp = Cij[Sv[Veh]*nN+Ev[Veh]];
            }
            // Check the constraints first
            
            if (NewCap <= Qv[Veh] && NewCP <= LIMv[Veh])
            {
                             
                // Then extend                
                NewPvr = NewP + (NewC-NewCP) - ODcomp;
                if(NewPvr > 0){
                    NewCost = log(NewPvr) - NewDelta - Veh_Dual - (NewC-ODcomp)*CompanyCost_Dual - NewPvr*Hypotenus_Dual_Coeff;
                    
                    // Is the cost feasible?
                    if(i == Ev[Veh] && NewCost * Min_Or_Max <= EPS){continue;}//
                }else{
                    if(i == Ev[Veh]){continue;}
                    NewCost = NewPvr;
                }
                // Add the new label to the Dict
                Dict.AddLabel_NSW(Sel_Lab, NewP, NewDelta, Demand_Dual, NewC, NewCP, NewCost, NewCap, i, k, Min_Or_Max
                ,Veh_Dual,CompanyCost_Dual,Hypotenus_Dual_Coeff, Veh );
            }
        }
        
//        Sel_Lab->~Label();
        if (Sel_Lab->P + (Sel_Lab->C-Sel_Lab->CP) > 0 ){
            Dict.Q[Sel_Node].push_back(Sel_Lab);
        }else{
            Dict.Qp[Sel_Node].push_back(Sel_Lab);
        }
        Dict.U[Sel_Node].erase(Dict.U[Sel_Node].begin());
    }
    std::vector<Route*> Returned_Route_Vector;
    long nReturnedRoute = nRoute;
    if (nReturnedRoute<0){
        nReturnedRoute = Dict.U[Ev[Veh]].size();
    }else{
        nReturnedRoute = std::min(nReturnedRoute, Dict.get_size(Ev[Veh]));
    }
     
    for(int i=0; i<nReturnedRoute; i++){
        Route *NewRoute = new Route();
        NewRoute->Seq.push_back(Ev[Veh]);
        Label* From = Dict.U[Ev[Veh]].at(0)->From;
        while(true){            
            if (From->Node == 0)
            {
                NewRoute->Seq.insert(NewRoute->Seq.begin(),0);
                NewRoute->Seq.insert(NewRoute->Seq.begin(),Sv[Veh]);
                break;
            }
            NewRoute->Seq.insert(NewRoute->Seq.begin(),From->Node);
            From = From->From;
        }    
        NewRoute->P = Dict.U[Ev[Veh]].at(0)->P + (Dict.U[Ev[Veh]].at(0)->C - Dict.U[Ev[Veh]].at(0)-> CP)
                - Cij[Sv[Veh]*nN+Ev[Veh]];
        
        for(int j = 0; j < NewRoute->Seq.size()-1 ;j++)
        {
            NewRoute->Cost += Cij[NewRoute->Seq.at(j)*nN + NewRoute->Seq.at(j+1)];
        }
        NewRoute->Cost -= Cij[Sv[Veh]*nN+Ev[Veh]];
        
        Dict.U[Ev[Veh]].at(0)->~Label();
        Dict.U[Ev[Veh]].erase(Dict.U[Ev[Veh]].begin());
    

        Returned_Route_Vector.push_back(NewRoute);
    }
    
    
//    std:: cout<< "Iteration-------->" << Iteration << std::endl;
    for(std::map<int, std::vector<Label*> >::iterator item=Dict.U.begin(); item!=Dict.U.end();item++){
        for(int i = 0; i<item->second.size(); i++){
            Label* temp = item->second.at(i);
            temp->~Label();
            temp = NULL;
        }
        item->second.clear();
    }
    Dict.U.clear();
    for(std::map<int, std::vector<Label*> >::iterator item=Dict.Q.begin(); item!=Dict.Q.end();item++){
        for(int i = 0; i<item->second.size(); i++){
            Label* temp = item->second.at(i);
            temp->~Label();
            temp = NULL;
        }
        item->second.clear();
    }
    Dict.Q.clear();
    for(std::map<int, std::vector<Label*> >::iterator item=Dict.Qp.begin(); item!=Dict.Qp.end();item++){
        for(int i = 0; i<item->second.size(); i++){
            Label* temp = item->second.at(i);
            temp->~Label();
            temp = NULL;
        }
        item->second.clear();
    }
    Dict.Qp.clear();
    
    
    return Returned_Route_Vector;
}






