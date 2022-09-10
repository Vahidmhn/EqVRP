#include "Label_Dict.h"
#include "Declarations.h"
#include <fstream>

Label_Dict:: ~Label_Dict(){
    std::map<int,std::vector<Label*>>::iterator it;
    for (it=U.begin(); it!= U.end(); ++it)
    {
        for(int i = 0; i < it->second.size(); i++)
        {
            it->second.at(i)->~Label();
        }
    }
    for (it = Q.begin(); it != Q.end(); ++it)
    {
        for(int i = 0; i < it->second.size(); i++)
        {
            it->second.at(i)->~Label();
        }
    }
    for (it = Qp.begin(); it != Qp.end(); ++it)
    {
        for(int i = 0; i < it->second.size(); i++)
        {
            it->second.at(i)->~Label();
        }
    }
    
    Q.clear();
    Qp.clear();
    U.clear();
}

long Label_Dict:: get_size(){
    long size = 0;
    std::map<int,std::vector<Label*>>::iterator it;
    for (it=U.begin(); it!= U.end(); ++it)
    {
        size += it->second.size();
    }
    return size;
}

long Label_Dict::get_size(int Key){
    long size = 0;
    size = U[Key].size();
    return size;
}


long Label_Dict::get_size_wo(int Key){
    long size = 0;
    std::map<int,std::vector<Label*>>::iterator it;
    for (it=U.begin(); it!= U.end(); ++it)
    {
        if (it->first == Key){continue;}
        size += it->second.size();
    }
    return size;
    
}

int Label_Dict::SelectNode()
{
    std::map<int,std::vector<Label*>>::iterator it;
    int node = 0;
    for(it = U.begin();it != U.end();++it)
    {
        if (it->second.size() >=1)
        {
            node = it->first;
            return node;
        }
    }
    return 0;
}

int Label_Dict::SelectNode(int End)
{
    std::map<int,std::vector<Label*>>::iterator it;
    int node = 0;
    for(it = U.begin();it != U.end();++it)
    {
        if (it->second.size() >=1 && it->first!=End && it->second.at(0)->Cost < 0)
        {
            node = it->first;
            return node;
        }
    }
    
    int Max = U.begin()->second.size();
    node = U.begin()->first;
    for(it = U.begin();it != U.end();++it)
    {
        if(it->first == End){continue;}
        if (it->second.size() > Max)
        {
            Max = it->second.size();
            node = it->first;
        }
    }
    return node;
}


Label* Label_Dict::PopLabel(int Node, int nNode)
{
    Label* Lab = new Label(nNode, U[Node].at(0)->Node, U[Node].at(0)->CP);
    Lab->Cost = U[Node].at(0)->Cost;
    Lab->C = U[Node].at(0)->C;
    Lab->CP = U[Node].at(0)->CP;
    Lab->Delta = U[Node].at(0)->Delta;
    Lab->Cap = U[Node].at(0)->Cap;
    Lab->nVisited = U[Node].at(0)->nVisited;
    Lab->From = U[Node].at(0)->From;
    Lab->Valid = U[Node].at(0)->Valid;
    std::copy(U[Node].at(0)->Visits, U[Node].at(0)->Visits+nNode, Lab->Visits);
    U[Node].at(0)->~Label();
    U[Node].erase(U[Node].begin());
    
    return Lab;
}

Label* Label_Dict::PopLabel(int Node, int label, int nNode)
{
    Label* Lab = new Label(nNode,U[Node].at(0)->Node, U[Node].at(0)->CP);
    Lab->Cost = U[Node].at(label)->Cost;
    Lab->C = U[Node].at(0)->C;
    Lab->CP = U[Node].at(0)->CP;
    Lab->Delta = U[Node].at(0)->Delta;
    Lab->Cap = U[Node].at(label)->Cap;
    Lab->nVisited = U[Node].at(label)->nVisited;
    Lab->From = U[Node].at(0)->From;
    Lab->Valid = U[Node].at(label)->Valid;
    std::copy( U[Node].at(label)->Visits,U[Node].at(label)->Visits+nNode, Lab->Visits);
    U[Node].at(label)->~Label();
    U[Node].erase(U[Node].begin()+label);
    
    return Lab;
}

void Label_Dict::AddLabel(Label* Sel_Lab, double NewCost, double NewCap, double C, double CP, double NewP,
                        int NewNode, int NodeIndexInValid, int Veh, int Min_Or_Max){
    Label* NewLabel = new Label(nD+1, NewNode, CP);
    
    
    NewLabel->Cost = NewCost;
    NewLabel->Cap = NewCap;
    NewLabel->P = NewP;
    NewLabel->C = C;
    NewLabel->From = Sel_Lab;
    NewLabel->Valid = Sel_Lab->Valid;
    // Remove NewNode from Valid
    NewLabel->Valid.erase(NewLabel->Valid.begin()+NodeIndexInValid);        
    std::copy(Sel_Lab->Visits,Sel_Lab->Visits+nD+1, NewLabel->Visits);
    if(NewNode<=nD){
        NewLabel->Visits[NewNode] = 1;
        NewLabel->nVisited = Sel_Lab->nVisited + 1;
    }
    
    
    
    
    
    if(U[NewNode].size() > 0){
        // Add if it is not dominated and remove dominated labels
        int i = 0;
        do 
        {
            int Dom_Res = Dominates(U[NewNode].at(i), NewLabel, nD+1, Veh, Min_Or_Max);
            if (Dom_Res == 1){
                NewLabel->~Label(); 
                return;
            }else if(Dom_Res == -1){
                U[NewNode].at(i)->~Label();
                U[NewNode].erase(U[NewNode].begin()+i);
                i--;
            }
            i++;
        }while(i< U[NewNode].size()); 
    }
    if(Q[NewNode].size() > 0){
        // Add if it is not dominated and remove dominated labels
        int i = 0;
        do 
        {
            int Dom_Res = Dominates(Q[NewNode].at(i), NewLabel, nD+1, Veh, Min_Or_Max);
            if (Dom_Res == 1){
                NewLabel->~Label(); 
                return;
            }else if(Dom_Res == -1){
                Q[NewNode].at(i)->~Label();
                Q[NewNode].erase(Q[NewNode].begin()+i);
                i--;
            }
            i++;
        }while(i< Q[NewNode].size()); 
    }
    
    if(NewCost>0){
        U[NewNode].push_back(NewLabel);
    }else{
        U[NewNode].insert(U[NewNode].begin(),NewLabel);
    }
    
}

void Label_Dict::AddLabel(Label* Sel_Lab, double NewP, double NewDelta, double NewC, double NewCP,
        double NewCost, double NewCap, int NewNode, int NodeIndexInValid, int Veh, int Min_Or_Max){
    
    Label* NewLabel = new Label(nD+1, NewNode, NewP, NewDelta, NewC, NewCP, NewCost, NewCap, 
             Sel_Lab);

    // Remove NewNode from Valid
    NewLabel->Valid.erase(NewLabel->Valid.begin()+NodeIndexInValid);
    
    
    if(NewNode<=nD){
        NewLabel->Visits[NewNode] = 1;
        NewLabel->nVisited += 1;
    }
    
    if(NewP+(NewC-NewCP) > 0 && U[NewNode].size() > 0 ){
        // Add if it is not dominated and remove dominated labels
        bool Key = false;
        int i = 0;
        do 
        {
            double Pvr_U = U[NewNode].at(i)->P+(U[NewNode].at(i)->C-U[NewNode].at(i)->CP);
            if(Pvr_U>0){
                int Dom_Res = Dominates(U[NewNode].at(i),NewLabel,nD+1, Veh, Min_Or_Max);
                if (Dom_Res == 1){
                    NewLabel->~Label(); 
                    return;
                }else if(Dom_Res == -1){
                    U[NewNode].at(i)->~Label();
                    U[NewNode].erase(U[NewNode].begin()+i);
                    i--;
                }
            }
            i++;
        }while(i< U[NewNode].size());        
    }
    
    if(Q[NewNode].size() > 0){
        // Add if it is not dominated and remove dominated labels
        int i = 0;
        do 
        {
            int Dom_Res = Dominates(Q[NewNode].at(i), NewLabel, nD+1, Veh, Min_Or_Max);
            if (Dom_Res == 1){
                NewLabel->~Label(); 
                return;
            }else if(Dom_Res == -1){
                Q[NewNode].at(i)->~Label();
                Q[NewNode].erase(Q[NewNode].begin()+i);
                i--;
            }
            i++;
        }while(i < Q[NewNode].size()); 
    }
    if(NewCost<0){
        U[NewNode].push_back(NewLabel);
    }else{
        U[NewNode].insert(U[NewNode].begin(),NewLabel);
    }
    
}


//---------------------------------------------------------------------------------

void Label_Dict::AddLabel_NSW(Label* Sel_Lab, double NewP, double NewDelta, double *Demand_Dual , double NewC, double NewCP,
        double NewCost, double NewCap, int NewNode, int NodeIndexInValid, int Min_Or_Max,
        double Veh_Dual, double CompanyCost_Dual, double Hypotenus_Dual_Coeff, int Veh){
    
    
    
    Label* NewLabel = new Label(nD+1, NewNode, NewP, NewDelta, NewC, NewCP, NewCost, NewCap, 
             Sel_Lab);

    // Remove NewNode from Valid
    NewLabel->Valid.erase(NewLabel->Valid.begin()+NodeIndexInValid);
    
    
    if(NewNode<=nD){
        NewLabel->Visits[NewNode] = 1;
        NewLabel->nVisited += 1;
    }
    
    if(NewP+(NewC-NewCP) > 0 && U[NewNode].size() > 0 ){
        // Add if it is not dominated and remove dominated labels
        bool Key = false;
        int i = 0;
        do 
        {
            double Pvr_U = U[NewNode].at(i)->P+(U[NewNode].at(i)->C-U[NewNode].at(i)->CP);
            
            if(Pvr_U>0){
                int Dom_Res = Dominates_NSW(U[NewNode].at(i),NewLabel,nD+1, Min_Or_Max, Veh_Dual, CompanyCost_Dual,Hypotenus_Dual_Coeff, Veh, Demand_Dual,NewNode );
                if (Dom_Res == 1){
                    NewLabel->~Label(); 
                    return;
                }else if(Dom_Res == -1){
                    U[NewNode].at(i)->~Label();
                    U[NewNode].erase(U[NewNode].begin()+i);
                    i--;
                }
            }
            i++;
        }while(i< U[NewNode].size());        
    }
    
    if(NewP+(NewC-NewCP) > 0 && Q[NewNode].size() > 0){
        // Add if it is not dominated and remove dominated labels
        int i = 0;
        do 
        {
            int Dom_Res = Dominates_NSW(Q[NewNode].at(i), NewLabel, nD+1, Min_Or_Max, Veh_Dual, CompanyCost_Dual,Hypotenus_Dual_Coeff, Veh, Demand_Dual, NewNode);
            if (Dom_Res == 1){
                NewLabel->~Label(); 
                return;
            }else if(Dom_Res == -1){
                Q[NewNode].at(i)->~Label();
                Q[NewNode].erase(Q[NewNode].begin()+i);
                i--;
            }
            i++;
        }while(i < Q[NewNode].size()); 
    }
    if(NewCost<0){
        U[NewNode].push_back(NewLabel);
    }else{
        U[NewNode].insert(U[NewNode].begin(),NewLabel);
    }
    
}

int Label_Dict::Dominates_NSW(Label* Lab1, Label* Lab2,int n, int Dir,
        double Veh_Dual, double CompanyCost_Dual, double Hypotenus_Dual_Coeff, int Veh, double *Demand_Dual, int K){

   if(Lab1->Cost*Dir > Lab2->Cost*Dir && Lab1->nVisited <= Lab2->nVisited 
            && VectorProd(Lab1->Visits,Lab2->Visits,n) == Lab1->nVisited)
   {
       bool HasNeg2 = HasNeg(Demand_Dual, Lab2, K, Veh);
       if(!HasNeg2)
       {
           
            double g1 = - Lab1->Delta - Veh_Dual - Lab1->C*CompanyCost_Dual - (Lab1->P+(Lab1->C-Lab1->CP))*Hypotenus_Dual_Coeff;
            double g2 = - Lab2->Delta - Veh_Dual - Lab2->C*CompanyCost_Dual - (Lab2->P+(Lab2->C-Lab2->CP))*Hypotenus_Dual_Coeff;
            double Ex2EndCost = CPij[Lab1->Node*nN+Ev[Veh]]-Cij[Lab1->Node*nN+Ev[Veh]] ;
   
   
            double Ex2End1 = log(Lab1->P+Lab1->C-Lab1->CP - Ex2EndCost - Cij[Sv[Veh]*nN+Ev[Veh]])
                             - Lab1->Delta - Veh_Dual - (Lab1->C+Cij[Lab1->Node*nN+Ev[Veh]]- Cij[Sv[Veh]*nN+Ev[Veh]])*CompanyCost_Dual
                             - (Lab1->P+(Lab1->C-Lab1->CP) - Ex2EndCost)*Hypotenus_Dual_Coeff;
            
            double Ex2End2 = log(Lab2->P+Lab2->C-Lab2->CP - Ex2EndCost - Cij[Sv[Veh]*nN+Ev[Veh]])
                             - Lab2->Delta - Veh_Dual - (Lab2->C+Cij[Lab1->Node*nN+Ev[Veh]]- Cij[Sv[Veh]*nN+Ev[Veh]])*CompanyCost_Dual 
                             - (Lab2->P+(Lab2->C-Lab2->CP) - Ex2EndCost)*Hypotenus_Dual_Coeff;
           
            if(g2-g1<=0 && Ex2End1 >= Ex2End2){
               return 1;
            }
       }else{
           
            int iStar = Find_iStar(Lab2, K, Veh);//ISTAR[Veh*nD+Lab1->Node-1]; //
            double P_Ext_iStar = Pj[iStar-1]-CPij[Lab1->Node*nN+iStar]+Cij[Lab1->Node*nN+iStar] - 
                                    CPij[iStar*nN+Ev[Veh]] + Cij[iStar*nN+Ev[Veh]] - Cij[Sv[Veh]*nN+Ev[Veh]];
            
            
            
            //-------------------------------
            double Ex2EndCost = CPij[Lab1->Node*nN+Ev[Veh]]-Cij[Lab1->Node*nN+Ev[Veh]] ;
   
   
            double Ex2End1 = log(Lab1->P+Lab1->C-Lab1->CP - Ex2EndCost - Cij[Sv[Veh]*nN+Ev[Veh]])
                             - Lab1->Delta - Veh_Dual - (Lab1->C+Cij[Lab1->Node*nN+Ev[Veh]]- Cij[Sv[Veh]*nN+Ev[Veh]])*CompanyCost_Dual
                             - (Lab1->P+(Lab1->C-Lab1->CP) - Ex2EndCost)*Hypotenus_Dual_Coeff;
            
            double Ex2End2 = log(Lab2->P+Lab2->C-Lab2->CP - Ex2EndCost - Cij[Sv[Veh]*nN+Ev[Veh]])
                             - Lab2->Delta - Veh_Dual - (Lab2->C+Cij[Lab1->Node*nN+Ev[Veh]]- Cij[Sv[Veh]*nN+Ev[Veh]])*CompanyCost_Dual 
                             - (Lab2->P+(Lab2->C-Lab2->CP) - Ex2EndCost)*Hypotenus_Dual_Coeff;                    
            //-------------------------------  
            
            
                
            double P_r1 = Lab1->P+Lab1->C-Lab1->CP;
            double P_r2 = Lab2->P+Lab2->C-Lab2->CP;
           
            double L1 = - Lab1->Delta - Demand_Dual[iStar-1] - Veh_Dual - 
            (Lab1->C + Cij[Lab1->Node*nN+iStar] + Cij[iStar*nN+Ev[Veh]] - Cij[Sv[Veh]*nN+Ev[Veh]])*CompanyCost_Dual - (P_r1 + P_Ext_iStar)*Hypotenus_Dual_Coeff;
           
            double L2 = - Lab2->Delta - Demand_Dual[iStar-1] - Veh_Dual - 
            (Lab2->C + Cij[Lab1->Node*nN+iStar] + Cij[iStar*nN+Ev[Veh]] - Cij[Sv[Veh]*nN+Ev[Veh]])*CompanyCost_Dual - (P_r2 + P_Ext_iStar)*Hypotenus_Dual_Coeff;
           
            if(P_r2 >= P_r1 && P_r1 + P_Ext_iStar >= 0 && log(P_r1 + P_Ext_iStar) + L1>= log(P_r2 + P_Ext_iStar) + L2 && Ex2End1 >= Ex2End2){//
                return 1;
            }
       }
       
   }
   if(Lab1->Cost*Dir < Lab2->Cost*Dir && Lab1->nVisited >= Lab2->nVisited
            && VectorProd(Lab1->Visits,Lab2->Visits,n) == Lab2->nVisited){
       bool HasNeg1 = HasNeg(Demand_Dual, Lab1, K, Veh);
       if(!HasNeg1)
       {   
            double g1 = - Lab1->Delta - Veh_Dual - Lab1->C*CompanyCost_Dual - (Lab1->P+(Lab1->C-Lab1->CP))*Hypotenus_Dual_Coeff;
            double g2 = - Lab2->Delta - Veh_Dual - Lab2->C*CompanyCost_Dual - (Lab2->P+(Lab2->C-Lab2->CP))*Hypotenus_Dual_Coeff;
            double Ex2EndCost = CPij[Lab1->Node*nN+Ev[Veh]] - Cij[Lab1->Node*nN+Ev[Veh]] ;
   
   
            double Ex2End1 = log(Lab1->P+Lab1->C-Lab1->CP - Ex2EndCost - Cij[Sv[Veh]*nN+Ev[Veh]])
                             - Lab1->Delta - Veh_Dual - (Lab1->C+Cij[Lab1->Node*nN+Ev[Veh]]- Cij[Sv[Veh]*nN+Ev[Veh]])*CompanyCost_Dual
                             - (Lab1->P+(Lab1->C-Lab1->CP) - Ex2EndCost)*Hypotenus_Dual_Coeff;
            double Ex2End2 = log(Lab2->P+Lab2->C-Lab2->CP - Ex2EndCost - Cij[Sv[Veh]*nN+Ev[Veh]])
                             - Lab2->Delta - Veh_Dual - (Lab2->C+Cij[Lab1->Node*nN+Ev[Veh]]- Cij[Sv[Veh]*nN+Ev[Veh]])*CompanyCost_Dual 
                             - (Lab2->P+(Lab2->C-Lab2->CP) - Ex2EndCost)*Hypotenus_Dual_Coeff;
           
            if(g1-g2<=0 && Ex2End2 >= Ex2End1){
               return -1;
            }
       }else{
           int iStar = Find_iStar(Lab1, K, Veh);//ISTAR[Veh*nD+Lab1->Node-1];//
           double P_r1 = Lab1->P+Lab1->C-Lab1->CP;
           double P_r2 = Lab2->P+Lab2->C-Lab2->CP;
           double P_Ext_iStar = Pj[iStar-1]-CPij[Lab1->Node*nN+iStar]+Cij[Lab1->Node*nN+iStar] - 
                                    CPij[iStar*nN+Ev[Veh]] + Cij[iStar*nN+Ev[Veh]] - Cij[Sv[Veh]*nN+Ev[Veh]];
           
           
           //-------------------------------
            double Ex2EndCost = CPij[Lab1->Node*nN+Ev[Veh]] - Cij[Lab1->Node*nN+Ev[Veh]] ;
   
   
            double Ex2End1 = log(Lab1->P+Lab1->C-Lab1->CP - Ex2EndCost - Cij[Sv[Veh]*nN+Ev[Veh]])
                             - Lab1->Delta - Veh_Dual - (Lab1->C+Cij[Lab1->Node*nN+Ev[Veh]]- Cij[Sv[Veh]*nN+Ev[Veh]])*CompanyCost_Dual
                             - (Lab1->P+(Lab1->C-Lab1->CP) - Ex2EndCost)*Hypotenus_Dual_Coeff;
            double Ex2End2 = log(Lab2->P+Lab2->C-Lab2->CP - Ex2EndCost - Cij[Sv[Veh]*nN+Ev[Veh]])
                             - Lab2->Delta - Veh_Dual - (Lab2->C+Cij[Lab1->Node*nN+Ev[Veh]]- Cij[Sv[Veh]*nN+Ev[Veh]])*CompanyCost_Dual 
                             - (Lab2->P+(Lab2->C-Lab2->CP) - Ex2EndCost)*Hypotenus_Dual_Coeff;                    
            //------------------------------- 

           double L1 = - Lab1->Delta - Demand_Dual[iStar-1] - Veh_Dual - 
           (Lab1->C + Cij[Lab1->Node*nN+iStar] + Cij[iStar*nN+Ev[Veh]] - Cij[Sv[Veh]*nN+Ev[Veh]])*CompanyCost_Dual - (P_r1 + P_Ext_iStar)*Hypotenus_Dual_Coeff;
           
           double L2 = - Lab2->Delta - Demand_Dual[iStar-1] - Veh_Dual - 
           (Lab2->C + Cij[Lab1->Node*nN+iStar] + Cij[iStar*nN+Ev[Veh]] - Cij[Sv[Veh]*nN+Ev[Veh]])*CompanyCost_Dual - (P_r2 + P_Ext_iStar)*Hypotenus_Dual_Coeff;
           
           if(P_r2 <= P_r1 && P_r2 + P_Ext_iStar >= 0 && log(P_r1 + P_Ext_iStar) + L1 <= log(P_r2 + P_Ext_iStar) + L2 && Ex2End1 <= Ex2End2){//
               return -1;
           }  
       }
       
    } 
   
   return 0;
    
}

int Label_Dict::Find_iStar(Label* Lab,int K,int v)
{
    int iStar = -1;
    double MinCost = INF;
    for(int i = 0; i < Lab->Valid.size();i++){
        if (Lab->Valid.at(i) == K || Lab->Visits[Lab->Valid.at(i)] == 1)
        {
            continue;
        }else{
            int j = Lab->Valid.at(i);
            if(MinCost > Pj[j]-(CPij[Lab->Node*nN+j] - Cij[Lab->Node*nN+j]) - (CPij[j*nN+Ev[v]] - Cij[j*nN+Ev[v]]))
            {
                MinCost = Pj[j]-(CPij[Lab->Node*nN+j] - Cij[Lab->Node*nN+j]) - (CPij[j*nN+Ev[v]] - Cij[j*nN+Ev[v]]);
                iStar = j;      
            }
        }
    }
    return iStar;
}

bool Label_Dict::HasNeg(double* Demand_Dual, Label* Lab,int K, int Veh)
{
    if(K != Ev[Veh]){
    
        for(int i = 0; i < Lab->Valid.size();i++){
            if (Lab->Valid.at(i) == K || Lab->Visits[Lab->Valid.at(i)] == 1)
            {
                continue;
            }else{
                if(Demand_Dual[Lab->Valid.at(i)-1] < 0)
                {
                    return true;
                }
            }
        }
    }
    return false;
}


//-----------------------------------------------------------------------

int Label_Dict::Dominates(Label* Lab1, Label* Lab2,int n, int Veh, int Dir){
    
    double P_Ext_iStar = 0;
    if (Lab1->Node <= nD){
        int iStar = ISTAR[Veh*nD+Lab1->Node-1];
        P_Ext_iStar = std::min(Pj[iStar-1]-CPij[Lab1->Node*nN+iStar]+Cij[Lab1->Node*nN+iStar] 
               - CPij[iStar*nN+Ev[Veh]] + Cij[iStar*nN+Ev[Veh]] - Cij[Sv[Veh]*nN+Ev[Veh]],
               - Cij[Sv[Veh]*nN+Ev[Veh]] );
    }
    
   if(Lab1->Cost*Dir >= Lab2->Cost*Dir && Lab1->nVisited <= Lab2->nVisited 
            && VectorProd(Lab1->Visits,Lab2->Visits,n) == Lab1->nVisited 
            && (Lab1->P + (Lab1->C-Lab1->CP) > Lab2->P + (Lab2->C-Lab2->CP) ||
                Lab1->P + (Lab1->C-Lab1->CP) + P_Ext_iStar > 0)){
     
        return 1; 
    }else if(Lab1->Cost*Dir <= Lab2->Cost*Dir && Lab1->nVisited >= Lab2->nVisited
            && VectorProd(Lab1->Visits,Lab2->Visits,n) == Lab2->nVisited
            && ( Lab1->P + (Lab1->C-Lab1->CP) < Lab2->P + (Lab2->C-Lab2->CP) ||
            Lab2->P + (Lab2->C-Lab2->CP) + P_Ext_iStar > 0)){ 
        return -1;
    }else{
        return 0;
    }
    
}