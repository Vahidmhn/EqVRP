#include "Route.h"

Route::Route()
{
    Cost = 0;
    P = 0;
}

Route::~Route()
{
    Seq.clear();
}