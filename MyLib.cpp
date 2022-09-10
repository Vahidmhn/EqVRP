#include "MyLib.h"


std::string Str(int x){
	std::stringstream S;
	S << x;
	return S.str();
}

void Print(double* Vec, int size)
{
    for (int i=0;i<size;i++){
        std::cout << Vec[i] << ", ";
    }
    std::cout << std::endl;
}

void Print(int* Vec, int size)
{
    for (int i=0;i<size;i++){
        std::cout << Vec[i] << ", ";
    }
    std::cout << std::endl;
}

void Print(bool* Vec, int size)
{
    for (int i=0;i<size;i++){
        std::cout << Vec[i] << ", ";
    }
    std::cout << std::endl;
}


void Print(double* Vec, int m,int n){
    for (int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            std::cout << Vec[i*n+j] << ", ";
        }        
        std::cout << std::endl;
    }
}

void Print(int* Vec, int m,int n){
    for (int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            std::cout << Vec[i*n+j] << ", ";
        }        
        std::cout << std::endl;
    }
}

std::vector<int> Array2Vector(int* Array, int size)
{
    std::vector<int> out(Array,Array+size);
    return out;
}

std::vector<double> Array2Vector(double* Array, int size)
{
    std::vector<double> out(Array,Array+size);
    return out;
}

std::vector<bool> Array2Vector(bool* Array, int size)
{
    std::vector<bool> out(Array,Array+size);
    return out;
}

bool IsEqual(std::vector<double>& Vec1, std::vector<double>& Vec2)
{
    if(Vec1.size() != Vec2.size()){return false;}
    for(int i=0; i<Vec1.size(); i++){
        if (std::abs(Vec1.at(i)-Vec2.at(i)) > 0.000001 ){return false;}
    }
    return true;
}

bool IsEqual(std::vector<int>& Vec1, std::vector<int>& Vec2)
{
    if(Vec1.size() != Vec2.size()){return false;}
    for(int i=0; i<Vec1.size(); i++){
        if (Vec1.at(i)!=Vec2.at(i) ){return false;}
    }
    return true;
}

bool IsEqual(std::vector<bool>& Vec1, std::vector<bool>& Vec2)
{
    if(Vec1.size() != Vec2.size()){return false;}
    for(int i=0; i<Vec1.size(); i++){
        if (Vec1.at(i)!=Vec2.at(i) ){return false;}
    }
    return true;
}

bool IsEqual(bool Vec1[],bool* Vec2, int n)
{
    
    for(int i=0; i<n; i++){
        if (Vec1[i]!=Vec2[i] ){return false;}
    }
    return true;
}



void Print(std::vector<int>& Vec)
{
    for (int i=0;i<Vec.size();i++){
        std::cout << Vec.at(i) << ", ";
    }
    std::cout << std::endl;
}

void Print(std::vector<double>& Vec)
{
    for (int i=0;i<Vec.size();i++){
        std::cout << Vec.at(i) << ", ";
    }
    std::cout << std::endl;
}

void Print(std::vector<int>& Vec, std::fstream & file)
{
    for (int i=0;i<Vec.size();i++){
        file << Vec.at(i) << ", ";
    }
    file << std::endl;
}

void Print(std::vector<double>& Vec, std::fstream & file)
{
    for (int i=0;i<Vec.size();i++){
        file << Vec.at(i) << ", ";
    }
    file << std::endl;
}

int Randi(int Begining, int End)
{
//    srand(time(NULL));
    
    return (rand() % (End-Begining+1)) + Begining;
}


double Rand()
{
//    srand (time(NULL));
    
    return (double)rand()/RAND_MAX;
}

double Sum(double* Array, int n){
    double S = 0;
    for(int i=0; i<n; i++)
    {
        S += Array[i];
    }
    return S;
}

std::vector<int> Vector(int From, int To){
    std::vector<int> Vec(To-From+1);
    for(int i =0; i<To-From+1; i++)
    {
        Vec.at(i) = From+i;
    }
    return Vec;
}

bool IsInSet(int element, std::vector<int>& Set)
{
    for(int i=0; i<Set.size(); i++)
    {
        if (Set.at(i) == element){return true;}
    }
    return false;
}

bool IsSubSet(std::vector<int>& Set1, std::vector<int>& Set2){
    bool Key;
    for(int i = 0;i<Set1.size(); i++)
    {
        Key = false;
        for(int j=0; j<Set2.size(); j++)
        {
            if(Set1.at(i) == Set2.at(j))
            {
                Key = true;
                break;
            }
        }
        if (Key == false){return false;}
    }
    return true;
}

std::vector<int> SetMinus(std::vector<int>& Set1, std::vector<int>& Set2)
{
    std::vector<int> Out;
    for(int i=0; i<Set1.size();i++){
        if(!IsInSet(Set1.at(i), Set2) ){
            Out.push_back(Set1.at(i));
        }
    }
    return Out;
}

int FindIn(int item, std::vector<int>& Vec){
    int ind = -1;
    for(int i=0; i<Vec.size(); i++){
        if(Vec.at(i) == item){
            ind = i;
            break;
        }
    }
    return ind;
}

int VectorProd(bool* Array1, bool* Array2, int size)
{
    int out = 0;
    for (int i=0; i<size; i++){
        out += Array1[i]*Array2[i]; 
    }
    return out;
}

int VectorProd(int* Array1, int* Array2, int size)
{
    int out = 0;
    for (int i=0; i<size; i++){
        out += Array1[i]*Array2[i]; 
    }
    return out;
}
double VectorProd(double* Array1, double* Array2, int size)
{
    double out = 0;
    for (int i=0; i<size; i++){
        out += Array1[i]*Array2[i]; 
    }
    return out;
}


std::vector<int> Sort(double* Vec,int N)
{
    std::vector<int> SortedInd = Vector(0,N-1);
    for(int i=1; i<N; i++){
        double tempVal = Vec[i];
        int tempInd = SortedInd[i];
        for (int j=i-1; j>=0; j--){
            if(tempVal > Vec[j]){
                Vec[j+1] = Vec[j];
                SortedInd[j+1] = SortedInd[j];
                if(j==0){
                    Vec[0] = tempVal;
                    SortedInd[0] = tempInd;
                }
            }else{
                Vec[j+1] = tempVal;
                SortedInd[j+1] = tempInd;
                break;
            }
        }
    }
    return SortedInd;
}
std::vector<int> Sort(int* Vec, int N)
{
    std::vector<int> SortedInd = Vector(0,N-1);
    for(int i=1; i<=N; i++){
        int tempVal = Vec[i];
        int tempInd = SortedInd[i];
        for (int j=i-1; j==0; j--){
            if(tempVal > Vec[j]){
                Vec[j+1] = Vec[j];
                SortedInd[j+1] = SortedInd[j];
                if(j==0){
                    Vec[0] = tempVal;
                    SortedInd[0] = tempInd;
                }
            }else{
                Vec[j+1] = tempVal;
                SortedInd[j+1] = tempInd;
                break;
            }
        }
    }
    return SortedInd;
}


bool IsFileEmpty(std::fstream& file)
{
    return file.tellp() == 0;
}


double get_wall_time()
{
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
