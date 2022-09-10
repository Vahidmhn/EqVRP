#include"CommonLibs.h"

using namespace std;



// Declarations ------------------------------------------------------
/// Classes
#include "Declarations.h"
#include "Label.h"
#include "Label_Dict.h"
#include "Route.h"
#include "InitColumns.h"
#include <vector>


#include "Solve_ESPPRC.h"
#include "Column_Gen.h"
#include "MyLib.h"

#include <fstream>

fstream Inputs, Outputs;

long nD,nV,nN;

int *Sv, *Ev;
double *Pj, *Qj, *Qv, *LIMv;
int *ISTAR;


double *Cij, *CPij;
int nInit_Col_Rand_VRPs = 500;
double ALPHA, BETA, GAMMA, W_LOAD, W_DIS;
double TimeLim = 3600 * 4;
double Time = 0;
long Max_ESSPRC_nROUTE = INF;
bool HypoTenuseCutKey = false;

// Subroutines -------------------------------------------------------
bool ReadData();
void test();
bool Check_Instance();
double CostCal(vector<int> &Seq);

int main(int argc, char** argv) {
    
    char* FileName = "test.txt";
    ALPHA = .05;    
    BETA = .01;
    GAMMA = 100;
    W_LOAD = 1;
    W_DIS= 0;
    
    if (argc > 1){
        FileName = argv[1];
        ALPHA = atof(argv[2]);   
        BETA = atof(argv[3]);
        GAMMA = atof(argv[4]);
        W_LOAD = atof(argv[5]);
        W_DIS = atof(argv[6]);
    }
    string Instance = FileName;

    bool ReadStat = ReadData();
 
    
    if (!ReadStat )
    {
        cout << "Reading the inputs encountered a problem!" << endl;
        return -1;
    }  
    
    
    
//    bool feas = Check_Instance();
//    Outputs.open("InsFeas.txt", ios_base::app | ios_base::out);
//    if (feas){
//        Outputs << Instance.substr(0,Instance.size()-4)<< ", "<< "Yes" << endl;
//    }else{
//        Outputs << Instance.substr(0,Instance.size()-4)<< ", "<< "No" << endl;
//    }
//    Outputs.close();
//    
//    
//    
//    test();
//    return 0;
//    
//    
//    
//    
//    
    Time = get_wall_time();
    
    Column_Gen1 ZstarModel;
    ZstarModel.Solve_Zstar2();
    
    Time = get_wall_time()-Time;
    double  ZstarVal = ZstarModel.PrintBestSol();
    
    
    
    // Write the routes----------------------------------------
    //---------------------------------------------------------
    int nVstar = 0;
    Outputs.open("Z_Routes.txt", ios_base::app | ios_base::out);
    if(IsFileEmpty(Outputs)){
        Outputs << "Instance, "<<  "[Costs], " << "[Benefits]" <<endl;
    }
    Outputs << Instance.substr(0,Instance.size()-4)<< ", ";
    for(int v=0; v<nV;v++){
        Outputs << ZstarModel.BestSol.at(v).Cost << ", ";
    }
    
    for(int v=0; v<nV;v++){
        Outputs <<  "[ ";
        double P = 0;
        for(int i=0; i < ZstarModel.BestSol.at(v).Seq.size();i++){
            Outputs << ZstarModel.BestSol.at(v).Seq.at(i)<< " ";
            if (ZstarModel.BestSol.at(v).Seq.at(i)>=1 && ZstarModel.BestSol.at(v).Seq.at(i)<=nD){
                P += Pj[ZstarModel.BestSol.at(v).Seq.at(i)-1];
            }
        }
         Outputs << "],  ";
         double tempCost = ZstarModel.BestSol.at(v).Cost + Cij[Sv[v]*nN+Ev[v]];
         ZstarModel.BestSol.at(v).P = P+(1-1/(1-BETA))*tempCost - Cij[Sv[v]*nN+Ev[v]];
         if(P>0) nVstar++;
    }
    for(int v=0; v<nV;v++){
        Outputs << ZstarModel.BestSol.at(v).P << ", ";
    }
    Outputs << endl; 
    Outputs.close();
    
    // Write the Algorithm's stat------------------------------
    //---------------------------------------------------------
     
    Outputs.open("Z_AlgStat.txt", ios_base::app | ios_base::out);
    if(IsFileEmpty(Outputs)){
        Outputs << "Instance, "<<  "ObjVal, " << "Time, " <<
                 "LB, " << "UB, " << "nVstar, "<< "Iteration " << endl;
    }
    
    
    Outputs << Instance.substr(0,Instance.size()-4)<< ", " << ZstarVal << ", " << Time << ",  "
            <<ZstarModel.glb<< ", " << ZstarModel.gub<< ", "<< nVstar<<", " << ZstarModel.iter<< ", " <<endl;

    Outputs.close();
    
    
    
    
    //--------------------------------------------------------------------------
//    return 0;
    
    //--------------------------------------------------------------------------
    // Prepare the parameters for NSW Model
    
    vector <int> Vstar;
    vector <Route> InitSol;
    for (int v = 0; v<nV;v++){
        if(ZstarModel.BestSol.at(v).Seq.size()>2){
            Vstar.push_back(v);
            InitSol.push_back(Route());
            InitSol.back().Cost = ZstarModel.BestSol.at(v).Cost;
            InitSol.back().Seq = ZstarModel.BestSol.at(v).Seq;
            
            for(int i = 0; i<InitSol.back().Seq.size()-1; i++)
            {
                InitSol.back().P += Cij[InitSol.back().Seq.at(i)*nN+InitSol.back().Seq.at(i+1)]
                        -CPij[InitSol.back().Seq.at(i)*nN+InitSol.back().Seq.at(i+1)];
                if(InitSol.back().Seq.at(i+1)>=1 && InitSol.back().Seq.at(i+1)<=nD)
                {
                    InitSol.back().P += Pj[InitSol.back().Seq.at(i+1)-1];
                }
            }
            InitSol.back().P -= Cij[Sv[v]*nN+Ev[v]];
            
        }
    }
    nV = Vstar.size();
    double Zbar =0; 
    double NSW_DB = INF;
    map<int, std::vector<Route*>> Rv;
    ZstarModel.GetRv(Rv,Vstar);
    //------------------------------------------------------------------ 
    
    for(ALPHA = .1; ALPHA >= 0.09; ALPHA -= .024999){
        Zbar = ZstarVal*(1+ALPHA);
        
        
        cout << "-------------------------------->" << ALPHA<< endl;
        Time = get_wall_time();
        
        Column_Gen2 NSWModel(Zbar, Vstar, Rv, InitSol, NSW_DB);
        NSWModel.Solve_NSW();
        double NSWCost = NSWModel.PrintBestSol();
        Time = get_wall_time()-Time;

        
        Rv.clear();
        NSWModel.GetRv(Rv,Vstar);
        NSW_DB = NSWModel.gub ;
        // Write the routes----------------------------------------
        //---------------------------------------------------------
        Outputs.open("NSW_Routes.txt", ios_base::app | ios_base::out);
        if(IsFileEmpty(Outputs)){
            Outputs << "Instance, "<< "ALPHA, "<<  "[Costs], " << "[Benefits]" << endl;
        }
        Outputs << Instance.substr(0,Instance.size()-4)<< ", " << ALPHA<< ", ";
        for(int v=0; v<nV;v++){
            Outputs << NSWModel.BestSol.at(v).Cost << ", ";
        }    
        for(int v=0; v<nV;v++){
            Outputs << NSWModel.BestSol.at(v).P << ", ";
        }
        for(int v=0; v<nV;v++){
            Outputs<< "[ ";
            for(int i=0; i < NSWModel.BestSol.at(v).Seq.size();i++){
                 Outputs << NSWModel.BestSol.at(v).Seq.at(i) << " ";
            }
             Outputs << "], ";
        }
        Outputs << endl;
        Outputs.close();


        // Write the Algorithm's stat------------------------------
        //---------------------------------------------------------

        Outputs.open("NSW_AlgStat.txt", ios_base::app | ios_base::out);
        if(IsFileEmpty(Outputs)){
            Outputs << "Instance, "<< "ALPHA, "<<  "ObjVal, " << "Time, " <<
                     "LB, " << "UB, " << "Iteration "<< endl;
        }  

        Outputs << Instance.substr(0,Instance.size()-4)<<", " << ALPHA<< ", " << NSWCost << ", " << Time << ", "
                << NSWModel.glb<< ", "<< NSWModel.gub<<", "<< NSWModel.iter << ", ";

        Outputs << endl;
        Outputs.close();
    }
    
    
    
    cout << "It's done!"<< endl;
    
    
    //Kill the trash
    for(map<int,std::vector<Route*>>::iterator it = Rv.begin();it!=Rv.end(); ++it){
        for(int i = 0; i<it->second.size(); i++){
            it->second.at(i)->~Route();
        }
    }
    Rv.clear();
    
    delete [] Cij;
    delete [] CPij;
    delete [] Pj;
    delete [] Qj;
    delete [] Qv;
    delete [] Sv;
    delete [] Ev;
    
    return 0;
    
}


bool ReadData()
{
    Inputs.open("Input.txt");
    if (Inputs.fail()){
        std::cout << "The input file 'Input.txt' is not found!!" << std::endl;
        return false;
    }
    
    Inputs >> nD ;
    Inputs >> nV ;
    
    nN = nD + 2*nV + 1;
    
    
    Sv = new int[nV];
    Ev = new int[nV];
    
    Pj = new double[nD] ;
    Qj = new double[nD];
    Qv = new double[nV] ;
    LIMv = new double[nV];
    
    // Prize for each demand point
    for (int i= 0;i<nD; i++){
        Inputs >> Pj[i];
    }
    // Volume for each demand point
    for (int i= 0;i<nD; i++){
        Inputs >> Qj[i];
    }
    // Capacity for each vehicle
    for (int i= 0;i<nV; i++){
        Inputs >> Qv[i];
    }
    
    // Start point for each vehicle
    for (int i= 0;i<nV; i++){
        Inputs >> Sv[i];
    }
    
    // End point for each vehicle
    for (int i= 0;i<nV; i++){
        Inputs >> Ev[i];
    }
    
    // Coordinates
    double Xs[nN];
    double Ys[nN];
    for(int i=0; i<nN; i++)
    {
        Inputs >> Xs[i];
        Inputs >> Ys[i];
    }

    
    Cij = new double[nN*nN];
    CPij = new double[nN*nN];
    
    
    
    for(int i =0; i<nN; i++){
        for(int j=0; j<nN; j++)
        {
            CPij[i*nN+j] = sqrt(pow(Xs[i]-Xs[j],2) +  pow(Ys[i]-Ys[j],2));
        }
    }
    
    for(int i = 0; i<nN; i++)
    {
        for (int j=0; j<nN; j++)
        {
            Cij[i * nN+j] = (1-BETA)*CPij[i*nN+j];
        }
    }
//    for (int i= 0;i<nV; i++){
//        Cij[Sv[i]*nN+Ev[i]] = 0;
//    }
    
    // Time limit for each vehicle
    for (int i= 0;i<nV; i++){
        LIMv[i] = GAMMA*CPij[Sv[i]*nN+Ev[i]];//******************************8
    }
    
    // Prize for each demand point
    for (int i= 0;i<nD; i++){
        Pj[i] = W_LOAD*Pj[i] + W_DIS*CPij[0*nN+i+1];
    }
    
    
    
    ISTAR = new int[nV*nD];
    for(int v=0; v<nV; v++){
        for(int i=1; i<=nD; i++){
            double MinCost = INF;
            for(int j = 1; j<=nD; j++){
                if(i==j){continue;}
                if(MinCost > Pj[j-1] - (CPij[i*nN+j] - Cij[i*nN+j]) - (CPij[j*nN+Ev[v]] - Cij[j*nN+Ev[v]]))
                {
                    MinCost = Pj[j-1]-(CPij[i*nN+j] - Cij[i*nN+j]) - (CPij[j*nN+Ev[v]] - Cij[j*nN+Ev[v]]);
                    ISTAR[v*nD+i-1] = j;      
                }
            }
        }
    }
    

    return true;
    
} 

bool Check_Instance()
{
    for (int i = 1;i<=nD;i++){
        for(int j = 0; j<=nD; j++){
            if(Pj[i-1]<CPij[j*nN+i]){
                std::cout << j <<"-->" << i<< " = "<< CPij[j*nN+i] << "| " << Pj[i-1]; 
                return false;
            }
        }
    }
    return true;
}

void test(){
    
    int v = 0;
    vector <int> Sol1 = {12, 0, 8, 4, 7, 1, 17};
    
    
    cout<< " Cost = " << CostCal(Sol1) << "  length = " << CostCal(Sol1)/(1-BETA) << " <= " << LIMv[v] << " : "
                        << "  load = " << 0<< " <= " << Qv[v]<< endl;

}


double CostCal(vector<int> &Seq)
{
    double cost = 0;
    for (int i=0; i<Seq.size()-1; i++)
    {
//        cout << Seq.at(i) << "-->" << Seq.at(i+1) << "= " << Cij[Seq.at(i)*nN+Seq.at(i+1)]<< endl;
        cost += Cij[Seq.at(i)*nN+Seq.at(i+1)];
    }
    cost -= Cij[Seq.at(0)*nN + *(Seq.end()-1)];
    return cost;
}


