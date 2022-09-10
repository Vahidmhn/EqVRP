#include "Column_Gen.h"
#include "Declarations.h"

void Column_Gen1::Solve_Zstar2(){ 
    std::string Name;
    /// B&B parameters
    double GLB = 0;
    double GUB = INF;
    
    
    //----------------------------------------------------------------------
    // Initial Columns
    /// Add Obvious O-D columns
    Add_OD_Columns(Rv);
    
    /// Add one set of feasible columns
//    GUB = Add_Feasible_Col(Rv, BestSol, GUB);
    GUB = Add_Heuristic_Col(Rv, BestSol, GUB);
    
    if (GUB < 0){
        std::cout << "The problem is infeasible!!!"<< std::endl;
        return;
    }
    
//    GUB = Gen_Rand_VRPs(Rv,BestSol, GUB, nInit_Col_Rand_VRPs);
    
    
//    GUB = Add_Google_VRP_Cols(Rv, BestSol);

    
    /// Add the initial columns to the model
    for(int v=0; v<nV;v++)
    {
        for(int i=0; i<Rv[v].size(); i++)
        {              
            IloNumColumn NewColExpr = Obj(Rv[v].at(i)->Cost);            
            for(int j=2; j<Rv[v].at(i)->Seq.size(); j++)
            {
                if( Rv[v].at(i)->Seq.at(j) <= nD){
                    NewColExpr += Demand_Cons[Rv[v].at(i)->Seq.at(j)-1](1);
                } 
            }
            NewColExpr += Veh_Cons[v](1);                
            Name = "Xvr(" + Str(v) + ")(" + Str(nRoute_v[v]) + ")";
            nRoute_v[v]++;
            Xvr[v].add(IloNumVar(NewColExpr, 0, IloInfinity, ILOFLOAT, Name.c_str()));
            NewColExpr.end();        
        }
    }
    
//    This is useless here since the GUB returned from Gen_Rand_VRPs is already integer
//    GUB = Solve_Int_RMP(BestSol,GUB);
    
    /// Initial Node
    Tree.push_back(new BBNode1(GLB,GUB));
    Tree.at(0)->NotBranchedDem = Vector(1,nD);
    //----------------------------------------------------------------------
    // Main loop
    bool ExistsNewCol = false;
    bool Non_Optimal_Route = false;
    long Iteration = 0;
    int Stat = 0;
    
    while(!(Tree.size()==0 || GUB-GLB<OptEPS || abs(GUB-GLB)/GUB<OptEPS))
    { 
        Iteration++;
        
        std::cout<< GLB << " -- " << GUB << " -- " << Iteration << std::endl;

        /// Pick a node to Explore
        BBNode1* Sel_Node = PopNode(0);
        
        CutsChanged = true;
        bool IsPruned = false;
        int v=0;
        int LoopCounter=0;
        
        for (;;){        
            
            // Time commands ******************************************
            if (TimeLim-(get_wall_time()-Time)<0){
                glb = GLB;
                gub = GUB;
                iter = Iteration;
                FreeTree();
                return;
            }
            // ********************************************************
            Solver.setParam(IloCplex::TiLim,TimeLim-(get_wall_time()-Time));
            Stat = Solve_Master(*Sel_Node);  
            
            // Time commands ******************************************
            if (TimeLim-(get_wall_time()-Time)<0){
                glb = GLB;
                gub = GUB;
                iter = Iteration;
                FreeTree();
                return;
            }
            // ********************************************************
            
            if (Stat == 0 )
            {
                std:: cout<< "Master Problem is infeasible" << std::endl;
                IsPruned = true;
                break;
            }
            
                 
            bool IsInteger = IsSolInteger();
            if (IsInteger){
                /// A primal bound is found
                try{
                    Sel_Node->PB = Solver.getObjValue();
                }catch(IloException& e){
                    std::cout<<e<< std::endl;
                }
                //Update Global Primal bound
                if (Sel_Node->PB<GUB)
                {
                    GUB = Sel_Node->PB;
                    SaveBestSol();
                }
            }
        
            Update_Duals();
            
            
            //// Solve Subproblem
            std::vector <Route*> NewColumns;
            
            NewColumns =  Solve_ESPPRC(v,Veh_Dual[v], Demand_Dual,
                    Sel_Node->Tabu[v], Max_ESSPRC_nROUTE, -1);
            
           
            // Time commands ******************************************
            if (TimeLim-(get_wall_time()-Time)<0){
                glb = GLB;
                gub = GUB;
                iter = Iteration;
                FreeTree();
                return;
            }
            // ********************************************************


            if(NewColumns.size() !=0 )
            {
                LoopCounter = 0;
                AddColumns2Rv(NewColumns, v);
                /// Update RMP
                AddColumns2RMP(NewColumns, v);
                // Free the memory of NewColumns
//                for (Route* item: NewColumns){
//                    item->~Route();
//                }

            }else{
                LoopCounter++;
            }
            NewColumns.clear();
                
            if(LoopCounter == nV)
            {
                Sel_Node->DB = Solver.getObjValue();
                if(IsInteger){IsPruned = true;}
                break;
            }
            v++;
            if(v == nV){v=0;}
        }// Node is Explored --------------------------------
        
        if(IsPruned)
        {
            // Update global dual bound without this node
            if(Tree.size()>0){
                GLB = UpdateGDB();
            }else{
                GLB = std::max(GLB,Sel_Node->DB);
            }
            // Kill the node
            delete Sel_Node;
            continue;
        }
        // Branch
        int Dem_Ind_NotBranched = 0;
        int BranchedVeh = 0;
        // Find the Branching variable before adding new columns
        bool IsBranchable = FindBranchedDem_Veh(*Sel_Node,Dem_Ind_NotBranched,BranchedVeh);

//        std::cout <<  Sel_Node->NotBranchedDem.at(Dem_Ind_NotBranched) <<
//                "->" << BranchedVeh << std::endl;
        if (!IsBranchable){
            std::cout<< "This must not happen!" << std::endl;
            continue;
        }
        
        // Branch
        Branch1(*Sel_Node, Dem_Ind_NotBranched, BranchedVeh);
        
        // Kill the node
        Sel_Node->~BBNode1();
        // Update the global dual bound with new nodes
        if(Tree.size()>0){
            GLB = UpdateGDB();
        }
    }
    glb = GLB;
    gub = GUB;
    iter = Iteration;
    FreeTree();
}

int Column_Gen1::Solve_Master(BBNode1& node)
{
    if(CutsChanged){
        // Remove all existing branching cuts 
        for (int v=0; v<nV; v++)
        {
            for (int i=0; i<Rv[v].size(); i++)
            {
                Xvr[v][i].setUB(IloInfinity);
            }
        }


        // Add the cuts   
        for (int v=0; v<nV; v++)
        {
             /// Add the Required cuts
//            for(int d=0;d<node.Required[v].size(); d++)
//            {
//                for (int i=0; i<Rv[v].size(); i++)
//                {
//                    if (!IsInSet(node.Required[v].at(d), Rv[v][i]->Seq) ){
//                        Xvr[v][i].setUB(0);
//                    }
//                }
//            }
            /// Add the Tabu cuts
            for(int d=0;d<node.Tabu[v].size(); d++)
            {
                for (int i=0; i<Rv[v].size(); i++)
                {
                    if (IsInSet(node.Tabu[v].at(d), Rv[v][i]->Seq) ){
                        Xvr[v][i].setUB(0);
                    }
                }
            }
        }
        CutsChanged = false;
    }
    Solver.extract(RMP);
//    Solver.exportModel("Test.lp");
    int Stat = 0 ;
    try{
        Stat = Solver.solve();
    }catch(IloException &e){
        std::cout << e << std::endl;
    }
    
    return Stat;
    
}

void Column_Gen1::Update_RedCost1(int v)
{    
    int i,j;
    std::copy(Cij,Cij+(nN*nN) , Reduced_Cost);

    Reduced_Cost[Sv[v]*nN+0] -= Veh_Dual[v];
    for(int j = 1; j<= nD;j++)
    {
        for(int i=0;i<nN; i++){
            Reduced_Cost[i*nN + j] -= Demand_Dual[j-1];
        }
    }
    
}

double Column_Gen1::Solve_Int_RMP(std::vector<Route>& BestSol, double UB){
    std::vector< std::vector<IloConversion>> VarTypeChange(nV);
    
    for (int v=0; v<nV; v++)
    {
        for(int i = 0; i<Rv[v].size(); i++){
            VarTypeChange[v].push_back(IloConversion(Env, Xvr[v][i], ILOBOOL));
            RMP.add(VarTypeChange[v][i]);
        }
        
    } 
    
    Solver.extract(RMP);
    int Stat = Solver.solve();
    
    std::cout <<Solver.getObjValue()<<std::endl;
    
    if(Stat == 1 && Solver.getObjValue() < UB){
        UB = Solver.getObjValue();
        for (int v=0; v<nV; v++)
        {
            for(int i = 0; i<Rv[v].size(); i++)
            {
                try
                {
                    double value = Solver.getValue(Xvr[v][i]);
                    if(value >= 1-EPS) {
                        BestSol.at(v).Seq = Rv[v].at(i)->Seq;
                        BestSol.at(v).Cost = Rv[v].at(i)->Cost;
                    }
                }catch(...){};
            }
        }
    }
    
    for (int v=0; v<nV; v++)
    {
        for(int i = 0; i<Rv[v].size(); i++){
            RMP.remove(VarTypeChange[v][i]);
            VarTypeChange[v][i].end();
        }
        
    } 
    
    return UB;
}

void Column_Gen1::AddColumns2Rv( std::vector <Route*>& NewColVector, int v)
{
    for (int i=0; i<NewColVector.size(); i++)
    {
//        Route* NewRoute = new Route();
//        NewRoute->Cost = NewColVector.at(i)->Cost;
//        NewRoute->Seq = NewColVector.at(i)->Seq;
        Rv[v].push_back(NewColVector.at(i));
    }
        
}

void Column_Gen1::AddColumns2RMP(std::vector <Route*>& NewColVector, int v)
{
    std::string Name;
    for(int col=0; col<NewColVector.size(); col++){
        IloNumColumn NewColExpr = Obj(NewColVector.at(col)->Cost);

        for(int i = 2; i <NewColVector.at(col)->Seq.size() ;i++)
        {
            if( NewColVector.at(col)->Seq.at(i) <= nD){
                NewColExpr += Demand_Cons[NewColVector.at(col)->Seq.at(i)-1](1);
            }
        }  
        NewColExpr += Veh_Cons[v](1);                
        Name = "Xvr(" + Str(v) + ")(" + Str(nRoute_v[v]) + ")";
        nRoute_v[v]++;
        Xvr[v].add(IloNumVar(NewColExpr, 0, IloInfinity, ILOFLOAT, Name.c_str()));
        NewColExpr.end();
                
    }
}

void Column_Gen1::SaveBestSol(){
    double cost=0;
    for (int v=0; v<nV; v++)
    {
        for(int i = 0; i<Rv[v].size(); i++)
        {
            double value = Solver.getValue(Xvr[v][i]);
            if(value>= 1-EPS){
                BestSol.at(v).Seq = Rv[v].at(i)->Seq;
                BestSol.at(v).Cost = Rv[v].at(i)->Cost;  
                BestSol.at(v).P = Rv[v].at(i)->P ;
                cost += BestSol.at(v).Cost;
            }
        }
    }
}

double Column_Gen1::PrintBestSol(){
    double Tcost = 0;
    for (int v=0; v<nV; v++)
    {
        Print(BestSol.at(v).Seq);
        std::cout<<BestSol.at(v).Cost << "--" << BestSol.at(v).P<< std::endl;
        Tcost += BestSol.at(v).Cost;
    }
    return Tcost;
}

bool Column_Gen1::IsSolInteger()
{
    for (int v=0; v<nV; v++)
    {
        for(int i = 0; i<Rv[v].size(); i++)
        {
            try
            {
                double value = Solver.getValue(Xvr[v][i]);
                if(value>= EPS && value<= 1-EPS) {
                    return false;
                }
            }catch(...){};
        }
    }
    return true;
}

bool Column_Gen1::FindBranchedDem_Veh(BBNode1& node, int& Dem, int& Veh){
    if (node.NotBranchedDem.size()==0){return false;};
    for(int v=0; v<nV;v++)
    {
        for(int i=0; i<Rv[v].size(); i++){
            if (Solver.getValue(Xvr[v][i])>EPS && Solver.getValue(Xvr[v][i])<1-EPS)
            {
                for(int d = 0; d < node.NotBranchedDem.size(); d++)
                {
                    if(IsInSet(node.NotBranchedDem.at(d),Rv[v].at(i)->Seq))
                    {
                        Dem = d;
                        Veh = v;
                        return true;                        
                    }
                }
            }
        }
    }
    return false;
}

void Column_Gen1::Branch1(BBNode1 &node, int Ind_Dem, int Veh )
{
    int Dem = node.NotBranchedDem.at(Ind_Dem);
    
    // Add Right node
    BBNode1* Right_Node = new BBNode1(node.DB, node.PB);
    Right_Node->Required = node.Required;
    Right_Node->Tabu = node.Tabu;    
    Right_Node->Tabu[Veh].push_back(Dem);
    
    Right_Node->NotBranchedDem = node.NotBranchedDem;

    std::copy(node.RequiredCap,node.RequiredCap+nV, Right_Node->RequiredCap);
    int ind = Add2Tree(Right_Node);
    
    // Add Left node
    if(node.RequiredCap[Veh]+Qj[Dem-1]<Qv[Veh])
    {
        BBNode1* Left_Node = new BBNode1(node.DB, node.PB);
        Left_Node->Required = node.Required;
        Left_Node->Tabu = node.Tabu;
        Left_Node->Required[Veh].push_back(Dem);
        for(int v=0; v<nV; v++){
           if(v==Veh || IsInSet(Dem, Left_Node->Tabu[v])){continue;}
           Left_Node->Tabu[v].push_back(Dem);
        }
        Left_Node->NotBranchedDem = node.NotBranchedDem;
        Left_Node->NotBranchedDem.erase(Left_Node->NotBranchedDem.begin()+Ind_Dem);
        
        std::copy(node.RequiredCap,node.RequiredCap+nV, Left_Node->RequiredCap);
        Left_Node->RequiredCap[Veh] += Qj[Dem-1];
        Tree.insert(Tree.begin()+ind, Left_Node);
    }
    
    
     
}

int Column_Gen1::Add2Tree(BBNode1* node)
{
    if(Tree.size()==0){
        Tree.push_back(node);
        return 0;
    }
    for(int i=0; i<Tree.size(); i++){
        if (node->DB < Tree.at(i)->DB){
            Tree.insert(Tree.begin()+i,node);
            return i;
        }
    }
    Tree.push_back(node);
    return Tree.size()-1;
}

double Column_Gen1::UpdateGDB()
{
    
    double MIN = Tree.at(0)->DB;
    for(int i=1; i<Tree.size(); i++)
    {
        if(Tree.at(i)->DB < MIN)
        {
            MIN = Tree.at(i)->DB;
        }
    }
    return MIN;
}


BBNode1* Column_Gen1::PopNode(int i){
    BBNode1* Node = new BBNode1(Tree.at(i)->DB,Tree.at(i)->PB);
    std::copy(Tree.at(i)->RequiredCap,Tree.at(i)->RequiredCap+nV,Node->RequiredCap);
    Node->Required = Tree.at(i)->Required;
    Node->Tabu = Tree.at(i)->Tabu;
    Node->NotBranchedDem = Tree.at(i)->NotBranchedDem;
    Tree.at(i)->~BBNode1();
    Tree.erase(Tree.begin()+i);
    return Node;
}

void Column_Gen1::Update_Duals(){      
    for(int i=0; i<nD;i++)
    {
        Demand_Dual[i] = Solver.getDual(Demand_Cons[i]);
    }
    for(int v=0; v<nV; v++)
    {
        Veh_Dual[v] = Solver.getDual(Veh_Cons[v]);;        
    }
}


void Column_Gen1::FreeTree(){
    for(int i=0; i<Tree.size(); i++){
        Tree.at(i)->~BBNode1();
    }
    Tree.clear();
}

Column_Gen1::Column_Gen1()
{
    RMP = IloModel(Env);
    Solver = IloCplex(RMP);
    Obj = IloObjective(Env);
    Xvr = D2NumVarMat(Env,nV);
    for(int v=0; v<nV; v++)
    {
        Xvr[v] = IloNumVarArray(Env,0);        
    } 
    Demand_RhS = IloNumArray(Env, nD);
    for(int i=0; i<nD; i++)
    {
        Demand_RhS[i] = 1;        
    }
    Veh_RhS = IloNumArray(Env, nV);
    for(int v=0; v<nV; v++)
    {
        Veh_RhS[v] = 1;        
    }
    Demand_Cons = IloRangeArray(Env, Demand_RhS, Demand_RhS);//IloInfinity
    Veh_Cons = IloRangeArray(Env, Veh_RhS, Veh_RhS);
    // Building the model ------------------------------------------
    Obj = IloMinimize(Env);
    RMP.add(Obj);
    RMP.add(Demand_Cons);
    RMP.add(Veh_Cons);
    
    // Solver Settings ---------------------------------------------
    Solver.setOut(Env.getNullStream());
    Solver.setParam(IloCplex::Param::Threads, nThreads);
    Solver.setParam(IloCplex::Param::Simplex::Tolerances::Optimality,1e-9);
//    Solver.setParam(IloCplex::Param::MIP::Tolerances::MIPGap,0);
//    Solver.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap,0);
//    Solver.setParam(IloCplex::Param::Preprocessing::Reduce,0);  
    
    
    nRoute_v = new long[nV];
    Veh_Dual = new double[nV];
    for(int v = 0; v<nV; v++)
    {
        nRoute_v[v] = 0;
        Veh_Dual[v] = 0;
        BestSol.push_back(Route());
    }
    Demand_Dual = new double[nD];
    for (int i = 0; i<nD; i++)
    {
        Demand_Dual[i] = 0;
    }
    
    Reduced_Cost = new double[nN*nN];    
    std:: copy(Cij,Cij+(nN*nN),Reduced_Cost);
    
    
}

Column_Gen1::~Column_Gen1()
{
    std::vector<BBNode1*> Tree;
    std::vector<Route> BestSol;
    
    Demand_RhS.end();
    Veh_RhS.end();
    Demand_Cons.end();
    Veh_Cons.end();
    Obj.end();
    Xvr.end();
    Solver.end();
    RMP.end();
    Env.end();
    
    for(int v= 0; v<nV; v++){
        for(int i =0; i<Rv[v].size();i++){
            Rv[v].at(i)->~Route();            
        }
    }
    Rv.clear();
    
    for(int i =0; i<Tree.size();i++){
        Tree.at(i)->~BBNode1();            
    }
    Tree.clear();
    
    delete [] nRoute_v;
    delete [] Veh_Dual;
    delete [] Demand_Dual;
    delete [] Reduced_Cost;
}


void Column_Gen1::GetRv(std::map<int,std::vector<Route*>>& NewRv, std::vector<int>& Vstar){
    
    for(int vv =0; vv<nV; vv++)
    {
        int v = Vstar.at(vv);
        for(int i=0; i<Rv[v].size();i++){
            if(Rv[v].at(i)->Seq.size()>2){
                Route* NewRoute = new Route();
                NewRoute->Cost = Rv[v].at(i)->Cost;
                NewRoute->Seq = Rv[v].at(i)->Seq;
                for(int j = 0; j<Rv[v].at(i)->Seq.size()-1; j++)
                {
                    if (Rv[v].at(i)->Seq.at(j+1) <= nD && Rv[v].at(i)->Seq.at(j+1) >=1)
                    {
                        NewRoute->P += Pj[Rv[v].at(i)->Seq.at(j+1)-1]
                                +Cij[Rv[v].at(i)->Seq.at(j)*nN+Rv[v].at(i)->Seq.at(j+1)]
                                -CPij[Rv[v].at(i)->Seq.at(j)*nN+Rv[v].at(i)->Seq.at(j+1)];
                    }else{
                        NewRoute->P += Cij[Rv[v].at(i)->Seq.at(j)*nN+Rv[v].at(i)->Seq.at(j+1)]
                                -CPij[Rv[v].at(i)->Seq.at(j)*nN+Rv[v].at(i)->Seq.at(j+1)];
                    }
                }
                NewRoute->P -= Cij[Sv[v]*nN+Ev[v]];
                
                if (NewRoute->P > 0){
                    NewRv[v].push_back(NewRoute);
                }else{
                    NewRoute->~Route();
                }
            }
            Rv[v].at(i)->~Route();
        }
        Rv[v].clear();
    }
    
}

void Column_Gen1::KillRv(){
    
    for(std::map<int,std::vector<Route*> >::iterator vv=Rv.begin(); vv != Rv.end(); vv++)
    {
        for(int i=0; i<vv->second.size();i++){
            vv->second.at(i)->~Route();
        }
        vv->second.clear();
    }
    
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void Column_Gen2::Solve_NSW(){
    std::string Name;
    /// B&B parameters
    double GLB = glb;
    double GUB = gub;
    
    
    //----------------------------------------------------------------------
    // Initial Columns
    for(int i = 0; i<BestSol.size(); i++)
    {
        GLB += log(BestSol.at(i).P);
    }
    
    GLB = Solve_MaxMin_RMP(GLB);
    
    
    
    /// Initial Node
    
    Tree.push_back(new BBNode1(GUB,GLB));
    Tree.at(0)->NotBranchedDem = Vector(1,nD);
    //----------------------------------------------------------------------
    // Main loop
    bool ExistsNewCol = false;
    bool Non_Optimal_Route = false;
    long Iteration = 0;
    int Stat = 0;

    while(!(Tree.size()==0 || GUB-GLB < OptEPS || abs(GUB-GLB)/GUB < OptEPS))
    { 
        Iteration++;
        std::cout<< GLB << " -- " << GUB << " -- " << Iteration << std::endl;
        
        /// Pick a node to Explore
        BBNode1* Sel_Node = PopNode(0);
        
        CutsChanged = true;
        bool IsPruned = false;
        int vv=0;
        int v = 0;
        int LoopCounter=0;
        bool Switch = false;
        for (;;){
            v = Vset.at(vv);
            
            //Solve Master Problem
            // Time commands ******************************************
            if (TimeLim-(get_wall_time()-Time)-.01<0){
                glb = GLB;
                gub = GUB;
                iter = Iteration;
                FreeTree();
                return;
            }
            // ********************************************************
            Solver.setParam(IloCplex::TiLim,TimeLim-(get_wall_time()-Time));
            Stat = Solve_Master(*Sel_Node);  
            
            // Time commands ******************************************
            if (TimeLim-(get_wall_time()-Time)<0){
                glb = GLB;
                gub = GUB;
                iter = Iteration;
                FreeTree();
                return;
            }
            // ******************************************************** 
            
            if (Stat == 0)
            {
//                std::cout<< "Master Problem is infeasible" << std::endl;
                IsPruned = true;
                break;
            }
            
            bool IsInteger = IsSolInteger();
            if (IsInteger)
            {
                /// A primal bound is found
                Sel_Node->PB = Solver.getObjValue();
                
                // Update Global Primal bound
                if(Sel_Node->PB > GLB)
                {
                    GLB = Sel_Node->PB;
                    SaveBestSol();
                }
            }
        
            Update_Duals();     
              
            std::vector <Route*> NewColumns;
            
            double Hypotenus_Dual_Coeff=0;
            for(int i = 0;i < Hypotenus_Dual.size(); i++)
            {
                Hypotenus_Dual_Coeff += Hypotenus_Dual.at(i)/(Xp.at(i).at(vv).P);
            }
            
            
            // Solve Subproblem -------------------------------------------
            NewColumns = Solve_ESPPRC(v, Demand_Dual, Veh_Dual[vv], CompanyCost_Dual,
                Hypotenus_Dual_Coeff, Sel_Node->Tabu[v], Max_ESSPRC_nROUTE, 1);            
      
            // Time commands ******************************************
            if (TimeLim-(get_wall_time()-Time)<0){
                glb = GLB;
                gub = GUB;
                iter = Iteration;
                FreeTree();
                return;
            }
            // ********************************************************
            
            if(NewColumns.size() !=0 )
            {
                LoopCounter = 0;
                AddColumns2Rv(NewColumns, v);
                /// Update RMP
                AddColumns2RMP(NewColumns, vv);
            }else{
                LoopCounter++;
            }
            NewColumns.clear();
            
            
            if(LoopCounter == nV)
            {
                Sel_Node->DB = Solver.getObjValue();
                
                if(IsInteger){
                    IsPruned = true;
                    // Add a Hypotenuse cut
                    if(HypoTenuseCutKey)AddHypoCut();
                }
                break;
                
            }
            if(HypoTenuseCutKey && IsInteger){
                // Add a Hypotenuse cut
                AddHypoCut();
            }
            vv++;
            if(vv == nV){vv=0;}
        }// Node is Explored --------------------------------
        
        if(IsPruned)
        {
            // Update global dual bound without this node
            if(Tree.size()>0){
                GUB = UpdateGDB();
            }else{
                GUB = std::min(GUB, Sel_Node->DB);
            }
            // Kill the node
            Sel_Node->~BBNode1();
            continue;
        }
        // Branch
        int Dem_Ind_NotBranched = 0;
        int BranchedVehID = 0;
        // Find the Branching variable before adding new columns
        bool IsBranchable = FindBranchedDem_Veh(*Sel_Node,Dem_Ind_NotBranched,BranchedVehID);
      

        // Branch
        Branch(*Sel_Node, Dem_Ind_NotBranched, BranchedVehID);
        
        // Kill the node
        Sel_Node->~BBNode1();
        
        // Update the global dual bound with new nodes
        if(Tree.size()>0){
            GUB = UpdateGDB();
        }
    }
    glb = GLB;
    gub = GUB;
    iter = Iteration;
    FreeTree();
}

Column_Gen2::Column_Gen2(double Zbar, std::vector<int>& Vstar)
{
    glb = 0;
    gub = INF;
    iter = 0;
            
            
    Vset = Vstar;   

    RMP = IloModel(Env);
    Solver = IloCplex(RMP);
    Obj = IloObjective(Env);
    Xvr = D2NumVarMat(Env,nV);
    for(int v=0; v<nV; v++)
    {
        Xvr[v] = IloNumVarArray(Env,0);        
    } 
    Demand_RhS = IloNumArray(Env, nD);
    for(int i=0; i<nD; i++)
    {
        Demand_RhS[i] = 1;        
    }
    Veh_RhS = IloNumArray(Env, nV);
    for(int v=0; v<nV; v++)
    {
        Veh_RhS[v] = 1;        
    }
    Demand_Cons = IloRangeArray(Env, Demand_RhS, Demand_RhS);
    Veh_Cons = IloRangeArray(Env, Veh_RhS, Veh_RhS);
    CompanyCost = IloRange(Env, 0, Zbar);
    // We can initialize some hypotenus cuts too
    
    // Building the model ------------------------------------------
    Obj = IloMaximize(Env);
    RMP.add(Obj);
    RMP.add(Demand_Cons);
    RMP.add(Veh_Cons);
    RMP.add(CompanyCost);
    RMP.add(HypotenusCut);
    
    // Solver Settings ------------------------------------------

    Solver.setOut(Env.getNullStream());
    Solver.setParam(IloCplex::Param::Threads,nThreads);
    Solver.setParam(IloCplex::Param::Simplex::Tolerances::Optimality,1e-9);
//    Solver.setParam(IloCplex::Param::MIP::Tolerances::MIPGap,0);
//    Solver.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap,0);
//    Solver.setParam(IloCplex::Param::Preprocessing::Reduce,0);  
    
    HypoCutCounter = 0;
    nRoute_v = new long[nV];
    Veh_Dual = new double[nV];
    for (int v = 0; v<nV; v++)
    {
        nRoute_v[v] = 0;
        Veh_Dual[v] = 0;
        BestSol.push_back(Route());
        IntSol.push_back(0);
    }
    Demand_Dual = new double[nD];
    for (int i = 0; i<nD; i++)
    {
        Demand_Dual[i] = 0;
    }
    CompanyCost_Dual = 0;
    
     
}

Column_Gen2::Column_Gen2(double Zbar, std::vector<int>& Vstar,
                 std::map<int,std::vector<Route*> >& InitRv, std::vector<Route> InitSol)
{
    glb = 0;
    gub = INF;
    iter = 0;
    
    Vset = Vstar;   
    Rv = InitRv;
    BestSol = InitSol;
    Xp.push_back(InitSol);
    
    RMP = IloModel(Env);
    Solver = IloCplex(RMP);
    Obj = IloMaximize(Env);
    Xvr = D2NumVarMat(Env,nV);
    for(int vv=0; vv<nV; vv++)
    {
        int v = Vset.at(vv);
        Xvr[vv] = IloNumVarArray(Env,Rv[v].size());        
    } 
    Demand_RhS = IloNumArray(Env, nD);
    for(int i=0; i<nD; i++)
    {
        Demand_RhS[i] = 1;        
    }
    Veh_RhS = IloNumArray(Env, nV);
    for(int v=0; v<nV; v++)
    {
        Veh_RhS[v] = 1;        
    }
    Demand_Cons = IloRangeArray(Env,Demand_RhS, Demand_RhS);
    Veh_Cons = IloRangeArray(Env, Veh_RhS, Veh_RhS);
    CompanyCost = IloRange(Env, 0, Zbar);
    HypoCutCounter = 0;
    HypotenusCut = IloRangeArray(Env, 0);
    HypotenusCut.add(IloRange(Env, nV, IloInfinity));
    HypoCutCounter++;
    Hypotenus_Dual.push_back(0);
    
    // Hypotenus cuts for BestSol
    
    std::string Name;
    nRoute_v = new long[nV];
    for (int vv = 0; vv < nV; vv++)
    {
        int v = Vset.at(vv);
        nRoute_v[vv]= 0;
        for(int i = 0; i<Rv[v].size(); i++)
        {
            IloNumColumn NewColExpr = Obj(log(Rv[v].at(i)->P)); 
            for(int j=2; j < Rv[v].at(i)->Seq.size()-1; j++)
            {
                NewColExpr += Demand_Cons[Rv[v].at(i)->Seq.at(j)-1](1);
            }
            NewColExpr += Veh_Cons[vv](1);  
            NewColExpr += CompanyCost(Rv[v].at(i)->Cost);
            
            
            NewColExpr += HypotenusCut[0]( Rv[v].at(i)->P/
                                       Xp.at(0).at(vv).P );
            
            
            Name = "Xvr(" + Str(v) + ")(" + Str(nRoute_v[vv]) + ")";
            nRoute_v[vv]++;
            Xvr[vv][i]=IloNumVar(NewColExpr, 0, IloInfinity, ILOFLOAT, Name.c_str());
        }
    }
    
        
    // Building the model ------------------------------------------
    
    RMP.add(Obj);
    RMP.add(Demand_Cons);
    RMP.add(Veh_Cons);
    RMP.add(CompanyCost);
    if(HypoTenuseCutKey){
        RMP.add(HypotenusCut);
    }
       
    // Solver Settings ------------------------------------------

    Solver.setOut(Env.getNullStream());
    Solver.setParam(IloCplex::Param::Threads,nThreads);
    Solver.setParam(IloCplex::Param::Simplex::Tolerances::Optimality,1e-9);
//    Solver.setParam(IloCplex::Param::MIP::Tolerances::MIPGap,0);
//    Solver.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap,0);
//    Solver.setParam(IloCplex::Param::Preprocessing::Reduce,0);  
    
    
    Veh_Dual = new double[nV];
    for (int v = 0; v<nV; v++)
    {
        Veh_Dual[v] = 0;
        IntSol.push_back(0);
    }
    Demand_Dual = new double[nD];
    for (int i = 0; i<nD; i++)
    {
        Demand_Dual[i] = 0;
    }
    CompanyCost_Dual = 0;

}

Column_Gen2::Column_Gen2(double Zbar, std::vector<int>& Vstar, std::map<int,
                std::vector<Route*> > &InitRv, std::vector<Route> InitSol, double GUB)
{
    glb = 0;
    gub = GUB;
    iter = 0;
    
    Vset = Vstar;   
    Rv = InitRv;
    BestSol = InitSol;
    Xp.push_back(InitSol);
    
    RMP = IloModel(Env);
    Solver = IloCplex(RMP);
    Obj = IloMaximize(Env);
    Xvr = D2NumVarMat(Env,nV);
    for(int vv=0; vv<nV; vv++)
    {
        int v = Vset.at(vv);
        Xvr[vv] = IloNumVarArray(Env,Rv[v].size());        
    } 
    Demand_RhS = IloNumArray(Env, nD);
    for(int i=0; i<nD; i++)
    {
        Demand_RhS[i] = 1;        
    }
    Veh_RhS = IloNumArray(Env, nV);
    for(int v=0; v<nV; v++)
    {
        Veh_RhS[v] = 1;        
    }
    Demand_Cons = IloRangeArray(Env, Demand_RhS, Demand_RhS);
    Veh_Cons = IloRangeArray(Env, Veh_RhS, Veh_RhS);
    CompanyCost = IloRange(Env, 0, Zbar);
    HypoCutCounter = 0;
    HypotenusCut = IloRangeArray(Env, 0);
    HypotenusCut.add(IloRange(Env, nV, IloInfinity));
    
    Hypotenus_Dual.push_back(0);
    
    // Hypotenus cuts for BestSol
    
    std::string Name;
    nRoute_v = new long[nV];
    for (int vv = 0; vv < nV; vv++)
    {
        int v = Vset.at(vv);
        nRoute_v[vv]= 0;
        for(int i = 0; i<Rv[v].size(); i++)
        {
            IloNumColumn NewColExpr = Obj(log(Rv[v].at(i)->P)); 
            for(int j=2; j < Rv[v].at(i)->Seq.size()-1; j++)
            {
                NewColExpr += Demand_Cons[Rv[v].at(i)->Seq.at(j)-1](1);
            }
            NewColExpr += Veh_Cons[vv](1);  
            NewColExpr += CompanyCost(Rv[v].at(i)->Cost);
            
            NewColExpr += HypotenusCut[0]( Rv[v].at(i)->P/
                                       Xp.at(0).at(vv).P );
            
            
            Name = "Xvr(" + Str(v) + ")(" + Str(nRoute_v[vv]) + ")";
            nRoute_v[vv]++;
            Xvr[vv][i]=IloNumVar(NewColExpr, 0, IloInfinity, ILOFLOAT, Name.c_str());
        }
    }
    
        
    // Building the model ------------------------------------------
    
    RMP.add(Obj);
    RMP.add(Demand_Cons);
    RMP.add(Veh_Cons);
    RMP.add(CompanyCost);
    if(HypoTenuseCutKey){
        HypoCutCounter++;
        RMP.add(HypotenusCut);
    }
       
    // Solver Settings ------------------------------------------

    Solver.setOut(Env.getNullStream());
    Solver.setParam(IloCplex::Param::Threads,nThreads);
    Solver.setParam(IloCplex::Param::Simplex::Tolerances::Optimality,1e-9);
//    Solver.setParam(IloCplex::Param::MIP::Tolerances::MIPGap,0);
//    Solver.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap,0);
//    Solver.setParam(IloCplex::Param::Preprocessing::Reduce,0);  
    
    
    Veh_Dual = new double[nV];
    for (int v = 0; v<nV; v++)
    {
        Veh_Dual[v] = 0;
        IntSol.push_back(0);
    }
    Demand_Dual = new double[nD];
    for (int i = 0; i<nD; i++)
    {
        Demand_Dual[i] = 0;
    }
    CompanyCost_Dual = 0;
    
}

BBNode1* Column_Gen2::PopNode(int i){
    BBNode1* Node = new BBNode1(Tree.at(i)->DB,Tree.at(i)->PB);
    std::copy(Tree.at(i)->RequiredCap,Tree.at(i)->RequiredCap+nV,Node->RequiredCap);
    Node->Required = Tree.at(i)->Required;
    Node->Tabu = Tree.at(i)->Tabu;
    Node->NotBranchedDem = Tree.at(i)->NotBranchedDem;
    Tree.at(i)->~BBNode1();
    Tree.erase(Tree.begin()+i);
    return Node;
}

int Column_Gen2::Solve_Master(BBNode1& node)
{
    if(CutsChanged){
        // Remove all existing branching cuts 
        for (int vv=0; vv<nV; vv++)
        {
            int v = Vset.at(vv);
            for (int i=0; i<Rv[v].size(); i++)
            {
                Xvr[vv][i].setUB(IloInfinity);
            }
        }

        // Add the cuts   
        for (int vv=0; vv<nV; vv++)
        {
            int v = Vset.at(vv);
             /// Add the Required cuts
//            for(int d=0;d<node.Required[v].size(); d++)
//            {
//                for (int i=0; i<Rv[v].size(); i++)
//                {
//                    if (!IsInSet(node.Required[v].at(d), Rv[v][i]->Seq) ){
//                        Xvr[vv][i].setUB(0);
//                    }
//                }
//            }
            /// Add the Tabu cuts
            for(int d=0;d<node.Tabu[v].size(); d++)
            {
                for (int i=0; i<Rv[v].size(); i++)
                {
                    if (IsInSet(node.Tabu[v].at(d), Rv[v][i]->Seq) ){
                        Xvr[vv][i].setUB(0);
                    }
                }
            }
        }
        CutsChanged = false;
    }
    Solver.extract(RMP);
//    Solver.exportModel("Test.lp");
    
    int Stat = 0 ;
    try{
        Stat = Solver.solve();
    }catch(IloException &e){
        std::cout << e << std::endl;
    }
    return Stat;
    
}

double Column_Gen2::Solve_Int_RMP(std::vector<Route>& BestSol, double LB){
    std::vector< std::vector<IloConversion>> VarTypeChange(nV);
    
    for (int vv=0; vv<nV; vv++)
    {
        int v = Vset.at(vv);
        for(int i = 0; i<Rv[v].size(); i++){
            VarTypeChange[vv].push_back(IloConversion(Env, Xvr[vv][i], ILOBOOL));
            RMP.add(VarTypeChange[vv][i]);
        }
    } 
    
    
    Solver.extract(RMP);
    
    int Stat = 0 ;
    try{
        Stat = Solver.solve();
    }catch(IloException &e){
        std::cout << e << std::endl;
    }
    
    if(Stat == 1 && Solver.getObjValue() > LB){
        LB = Solver.getObjValue();
        for (int vv=0; vv<nV; vv++)
        {
            int v = Vset.at(vv);
            for(int i = 0; i<Rv[v].size(); i++)
            {
                try
                {
                    double value = Solver.getValue(Xvr[vv][i]);
                    if(value >= 1-EPS) {
                        BestSol.at(vv).Seq = Rv[v].at(i)->Seq;
                        BestSol.at(vv).Cost = Rv[v].at(i)->Cost;
                    }
                }catch(...){};
            }
        }
    }
    
    for (int v=0; v<nV; v++)
    {
        for(int i = 0; i<Rv[v].size(); i++){
            RMP.remove(VarTypeChange[v][i]);
        }
        
    } 
    
    return LB;
}

double Column_Gen2::Solve_MaxMin_RMP(double LB){
    
    IloModel MaxMin = IloModel(Env);
    IloNumVar ZZ = IloNumVar(Env, 0, IloInfinity, ILOFLOAT);
    IloNumVarArray Z = IloNumVarArray(Env, nV, 0, IloInfinity, ILOFLOAT);
    IloCplex TempSolver(Env);
    IloObjective TempObj = IloMaximize(Env, ZZ);
    
    TempSolver.setOut(Env.getNullStream());
    
    MaxMin.add(TempObj);
    for (int vv=0; vv<nV; vv++)
    {
        int v = Vset.at(vv);
        IloExpr ObjV(Env);
        for(int i = 0; i<Rv[v].size(); i++){
            ObjV += Rv[v].at(i)->P * Xvr[vv][i];            
        }
        MaxMin.add(Z[vv]==ObjV);
        MaxMin.add(Z[vv]>=ZZ);
        ObjV.end();
    } 
    
    MaxMin.add(Demand_Cons);
    MaxMin.add(Veh_Cons);
    MaxMin.add(CompanyCost);
    MaxMin.add(HypotenusCut);
    
    for (int vv=0; vv<nV; vv++)
    {
        int v = Vset.at(vv);
        for(int i = 0; i<Rv[v].size(); i++){
            MaxMin.add(IloConversion(Env, Xvr[vv][i], ILOBOOL));
        }
    } 
    
    TempSolver.extract(MaxMin);
    int Stat = 0 ;
    try{
        Stat = TempSolver.solve();
    }catch(IloException &e){
        std::cout << e << std::endl;
    }
    
    if(Stat ==1){
        double ObjVal = 0;
        for (int vv=0; vv<nV; vv++)
        {
            ObjVal += log(TempSolver.getValue(Z[vv]));
        }
        if(ObjVal >= LB){
            LB = ObjVal;
            for (int vv=0; vv<nV; vv++)
            {
                int v = Vset.at(vv);
                for(int i = 0; i<Rv[v].size(); i++)
                {
                    try
                    {
                        double value = TempSolver.getValue(Xvr[vv][i]);
                        if(value >= 1-EPS) {
                            BestSol.at(vv).Seq = Rv[v].at(i)->Seq;
                            BestSol.at(vv).Cost = Rv[v].at(i)->Cost;
                            BestSol.at(vv).P = Rv[v].at(i)->P;
                        }
                    }catch(...){};
                }
            }
        }
    }
    TempObj.end();
    Z.end();
    ZZ.end();
    MaxMin.end();
    TempSolver.end();
    
    return LB;
}

bool Column_Gen2::IsSolInteger()
{
    for (int vv=0; vv<nV; vv++)
    {
        int v = Vset.at(vv);
        for(int i = 0; i<Rv[v].size(); i++)
        {
            try
            {
                double value = Solver.getValue(Xvr[vv][i]);
                if(value>= EPS && value<= 1-EPS) {
                    return false;
                }
                if(value> 1-EPS){
                    IntSol.at(vv) = i;
                }
            }catch(...){};
        }
    }
    return true;
}

void Column_Gen2::PrintFracSol()
{
    for (int vv=0; vv<nV; vv++)
    {
        int v = Vset.at(vv);
        for(int i = 0; i<Rv[v].size(); i++)
        {
            try
            {
                double value = Solver.getValue(Xvr[vv][i]);
                if(value>= EPS ) {
                    Print(Rv[v].at(i)->Seq);
                    std::cout<< v << "---" << i << " "<<Xvr[vv][i].getName()<<"= "<< value <<std::endl;
                    
                }
            }catch(...){};
        }
    }
}

void Column_Gen2::SaveBestSol(){
//    double cost=0;
    for (int vv=0; vv<nV; vv++)
    {
        int v = Vset.at(vv);
        int i = IntSol.at(vv);
        BestSol.at(vv).Seq = Rv[v].at(i)->Seq;
        BestSol.at(vv).Cost = Rv[v].at(i)->Cost;  
        BestSol.at(vv).P = Rv[v].at(i)->P;
//                cost += BestSol.at(v).Cost;
    }
}

void Column_Gen2::Update_Duals(){      
    for(int i=0; i<nD;i++)
    {
        Demand_Dual[i] = Solver.getDual(Demand_Cons[i]);
    }
    for(int v=0; v<nV; v++)
    {
        Veh_Dual[v] = Solver.getDual(Veh_Cons[v]);      
    }
    CompanyCost_Dual = Solver.getDual(CompanyCost);
    for(int i = 0; i<HypoCutCounter; i++)
    {
        Hypotenus_Dual.at(i) = Solver.getDual(HypotenusCut[i]);
    }
    
}

void Column_Gen2::AddColumns2Rv( std::vector <Route*>& NewColVector, int v)
{
    for (int i=0; i<NewColVector.size(); i++)
    {
//        Route* NewRoute = new Route();
//        NewRoute->Cost = NewColVector.at(i)->Cost;
//        NewRoute->Seq = NewColVector.at(i)->Seq;
//        NewRoute->P = NewColVector.at(i)->P;
        Rv[v].push_back(NewColVector.at(i));
    }
        
}

void Column_Gen2::AddColumns2RMP(std::vector <Route*>& NewColVector, int v)
{
    std::string Name;
    for(int col=0; col<NewColVector.size(); col++){
        IloNumColumn NewColExpr = Obj(log(NewColVector.at(col)->P));

        for(int i = 2; i < NewColVector.at(col)->Seq.size()-1 ;i++)
        {
            NewColExpr += Demand_Cons[NewColVector.at(col)->Seq.at(i)-1](1);
        }  
        NewColExpr += Veh_Cons[v](1);  
          
        NewColExpr += CompanyCost(NewColVector.at(col)->Cost);

        for(int xp=0; xp<HypoCutCounter; xp++){
            NewColExpr += HypotenusCut[xp](NewColVector.at(col)->P/
                                           Xp.at(xp).at(v).P );
        }        
        Name = "Xvr(" + Str(Vset.at(v)) + ")(" + Str(nRoute_v[v]) + ")";
        nRoute_v[v]++;
        Xvr[v].add(IloNumVar(NewColExpr, 0, IloInfinity, ILOFLOAT, Name.c_str()));
        NewColExpr.end();
    }
}

void Column_Gen2::AddHypoCut(){
    std::vector<Route> NewXp;
    
    // Add the solution to Xp
    for(int v=0; v<nV; v++){
        NewXp.push_back(Route());
        NewXp.at(v).P = Rv[Vset.at(v)].at(IntSol.at(v))->P;
    }
    Xp.push_back(NewXp);
    // Add the new coeff to all 
//    HupotenusCut.add(IloRange(Env, nV, IloInfinity));
    HypoCutCounter++;
    Hypotenus_Dual.push_back(0);
    IloExpr NewHypoExpr(Env);
    IloRange NewHypoCut(Env, nV, IloInfinity);
    for (int vv = 0; vv < nV; vv++)
    {
        int v = Vset.at(vv);
        for(int i = 0; i<Rv[v].size(); i++)
        {
            NewHypoExpr +=  Rv[v].at(i)->P/Xp.at(HypoCutCounter-1).at(vv).P * Xvr[vv][i];
        }
    }
    NewHypoCut.setExpr(NewHypoExpr);
    NewHypoExpr.end();
    RMP.remove(HypotenusCut);
    HypotenusCut.add(NewHypoCut);
    RMP.add(HypotenusCut);
    
}

double Column_Gen2::UpdateGDB()
{
    
    double MAX = Tree.at(0)->DB;
    for(int i=1; i<Tree.size(); i++)
    {
        if(Tree.at(i)->DB > MAX)
        {
            MAX = Tree.at(i)->DB;
        }
    }
    return MAX;
}

bool Column_Gen2::FindBranchedDem_Veh(BBNode1& node, int& Dem, int& Veh){
    if (node.NotBranchedDem.size()==0){return false;}
    for(int vv=0; vv<nV;vv++)
    {
        int v = Vset.at(vv);
        for(int i=0; i<Rv[v].size(); i++){
            if (Solver.getValue(Xvr[vv][i])>EPS && Solver.getValue(Xvr[vv][i])<1-EPS)
            {
                for(int d = 0; d < node.NotBranchedDem.size(); d++)
                {
                    if(IsInSet(node.NotBranchedDem.at(d),Rv[v].at(i)->Seq))
                    {
                        Dem = d;
                        Veh = vv;
                        return true;                        
                    }
                }
            }
        }
    }
    return false;
}

void Column_Gen2::Branch(BBNode1 &node, int Ind_Dem, int VehID )
{
    int Dem = node.NotBranchedDem.at(Ind_Dem);
    int Veh = Vset.at(VehID);
    
    // Add the right node
    BBNode1* Right_Node = new BBNode1(node.DB, node.PB);
    Right_Node->Required = node.Required;
    Right_Node->Tabu = node.Tabu;    
    Right_Node->Tabu[Veh].push_back(Dem);
    
    Right_Node->NotBranchedDem = node.NotBranchedDem;

    std::copy(node.RequiredCap,node.RequiredCap+nV, Right_Node->RequiredCap);
    int ind = Add2Tree(Right_Node);
    
    
    // Add the left node
    if(node.RequiredCap[VehID]+Qj[Dem-1]<Qv[VehID])
    {
        BBNode1* Left_Node = new BBNode1(node.DB, node.PB);
        Left_Node->Required = node.Required;
        Left_Node->Tabu = node.Tabu;
        Left_Node->Required[Veh].push_back(Dem);
        for(int v=0; v<nV; v++){
           if(v==VehID || IsInSet(Dem,Left_Node->Tabu[Vset.at(v)])){continue;}
           Left_Node->Tabu[Vset.at(v)].push_back(Dem);
        }
        Left_Node->NotBranchedDem = node.NotBranchedDem;
        Left_Node->NotBranchedDem.erase(Left_Node->NotBranchedDem.begin()+Ind_Dem);
        
        std::copy(node.RequiredCap,node.RequiredCap+nV, Left_Node->RequiredCap);
        Left_Node->RequiredCap[VehID] += Qj[Dem-1];
        Tree.insert(Tree.begin()+ind , Left_Node );
    }
    
    
     
}

int Column_Gen2::Add2Tree(BBNode1* node)
{
    if(Tree.size()==0){
        Tree.push_back(node);
        return 0;
    }
    for(int i=0; i<Tree.size(); i++){
        if (node->DB > Tree.at(i)->DB){
            Tree.insert(Tree.begin()+i,node);
            return i;
        }
    }
    Tree.push_back(node);
    return Tree.size()-1;
}

double Column_Gen2::PrintBestSol(){
    double Tcost = 0;
    for (int v=0; v<nV; v++)
    {
        std::cout<<"P_"<<Vset.at(v)<<" = "<< BestSol.at(v).P << std::endl;
        Print(BestSol.at(v).Seq);
        Tcost += BestSol.at(v).Cost;
    }
    std::cout<< "Total Cost = " << Tcost << std::endl;
    return Tcost;
}


void Column_Gen2::FreeTree(){    
    
    for(int i=0; i<Tree.size(); i++){
        Tree.at(i)->~BBNode1();
    }
    Tree.clear();
    
}

void Column_Gen2::GetRv(std::map<int,std::vector<Route*>>& NewRv, std::vector<int>& Vstar){
    
    for(int vv =0; vv<nV; vv++)
    {
        int v = Vstar.at(vv);
        for(int i=0; i<Rv[v].size();i++){
            if(Rv[v].at(i)->Seq.size()>2){
                Route* NewRoute = new Route();
                NewRoute->Cost = Rv[v].at(i)->Cost;
                NewRoute->Seq = Rv[v].at(i)->Seq;
                NewRoute->P = Rv[v].at(i)->P;
                if (NewRoute->P > 0){
                    NewRv[v].push_back(NewRoute);
                }
                Rv[v].at(i)->~Route();
            }
        }
        Rv[v].clear();
    }
    
}

void Column_Gen2::KillRv(){
    
    for(std::map<int,std::vector<Route*> >::iterator vv = Rv.begin(); vv != Rv.end(); vv++)
    {
        for(int i=0; i<vv->second.size();i++){
            vv->second.at(i)->~Route();
        }
        vv->second.clear();
    }
    
}


Column_Gen2::~Column_Gen2()
{
    std::vector<BBNode1*> Tree;
    std::vector<Route> BestSol;
    
    Demand_RhS.end();
    Veh_RhS.end();
    Demand_Cons.end();
    Veh_Cons.end();
    CompanyCost.end();
    HypotenusCut.end();
    Obj.end();
    Xvr.end();
    Solver.end();
    RMP.end();
    Env.end();
    
    for(int v= 0; v<nV; v++){
        for(int i =0; i<Rv[v].size();i++){
            Rv[v].at(i)->~Route();            
        }
    }
    Rv.clear();
    
    for(int i =0; i<Tree.size();i++){
        Tree.at(i)->~BBNode1();            
    }
    Tree.clear();
    
    delete [] nRoute_v;
    delete [] Veh_Dual;
    delete [] Demand_Dual;
    Hypotenus_Dual.clear();
}
