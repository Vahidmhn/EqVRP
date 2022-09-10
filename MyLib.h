#ifndef MYLIB_H
#define MYLIB_H

#include <vector>
#include <string>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>
#include <fstream>
#include "CommonLibs.h"
#include <vector>

std::string Str(int x);

void Print(double* Vec, int size);
void Print(int* Vec, int size);
void Print(bool* Vec, int size);
void Print(double* Vec, int m,int n);
void Print(int* Vec, int m,int n);

void Print(std::vector<int>& Vec);
void Print(std::vector<double>& Vec);
void Print(std::vector<int>& Vec, std::fstream & file);
void Print(std::vector<double>& Vec, std::fstream & file);

bool IsEqual(std::vector<double>& Vec1, std::vector<double>& Vec2);
bool IsEqual(std::vector<int>& Vec1, std::vector<int>& Vec2);
bool IsEqual(std::vector<bool>& Vec1, std::vector<bool>& Vec2);
bool IsEqual(bool Vec1[], bool* Vec2,int n);

std::vector<int> Array2Vector(int* Array, int size);
std::vector<double> Array2Vector(double* Array, int size);
std::vector<bool> Array2Vector(bool* Array, int size);

std::vector<int> Vector(int From, int To);

bool IsInSet(int element, std::vector<int>& Set);
bool IsSubSet(std::vector<int>& Set1, std::vector<int>& Set2);
std::vector<int> SetMinus(std::vector<int>& Set1, std::vector<int>& Set2);
int FindIn(int item, std::vector<int>& Vec);

int VectorProd(bool* Array1, bool* Array2, int size);
int VectorProd(int* Array1, int* Array2, int size);
double VectorProd(double* Array1, double* Array2, int size);

double Sum(double* Array, int n);

int Randi(int Begining, int End);
double Rand();

std::vector<int> Sort(double* Vec, int N);
std::vector<int> Sort(int* Vec, int N);


bool IsFileEmpty(std::fstream& file);

double get_wall_time();

#endif /* MYLIIB_H */

