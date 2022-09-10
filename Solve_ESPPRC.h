/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Solve_ESPPRC.h
 * Author: myubuntu
 *
 * Created on November 16, 2020, 3:53 PM
 */

#ifndef SOLVE_ESPPRC_H
#define SOLVE_ESPPRC_H

#include"CommonLibs.h"
#include "Label.h"
#include "Label_Dict.h"
#include "Route.h"

#include "InitColumns.h"

std::vector<Route*> Solve_ESPPRC(int Veh,double Veh_Dual, double* Demand_Dual,
        std::vector<int>& Tabu , long nRoute, int Min_Or_Max);

std::vector<Route*> Solve_ESPPRC(int Veh,double* Demand_Dual, double Veh_Dual,
        double CompanyCost_Dual, double Hypotenus_Dual_Coeff, 
         std::vector<int>& Tabu , long nRoute, int Min_Or_Max);


#endif /* SOLVE_ESPPRC_H */

