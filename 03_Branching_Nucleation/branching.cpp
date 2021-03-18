#include <iostream>
#include <fstream>
#include <random>
#include <cstdlib>
#include <math.h>

#include "branching.h"

using namespace std;

branching :: branching(){}

void branching :: branch_deactivate()
{
    branch_state = 0;
    
    // use value 2 as unbound but inactive when using deactivating nucleators
    // branch_state = 2;
    
    fil_m = 0;
    fil_d = 0;
    
    eps_m1 = 0.0; eps_m2 = 0.0;
    eps_d1 = 0.0; eps_d2 = 0.0;
    
}

void branching :: branch_walk(double fil_length_m,double fil_length_d, double dt, int fil_m_state, int fil_d_state, double v_g, double v_s)
{
    if(fil_m_state == 0)
    {
        eps_m1 = eps_m1 - 0.5*v_g*dt;
        eps_m2 = eps_m2 - 0.5*v_g*dt;
    }
 
    if(fil_m_state == 1 || fil_m_state == 2)
    {
        eps_m1 = eps_m1 + 0.5*v_s*dt;
        eps_m2 = eps_m2 + 0.5*v_s*dt;
    }
 
    if(fil_d_state == 0)
    {
        eps_d1 = eps_d1 - 0.5*v_g*dt;
        eps_d2 = eps_d2 - 0.5*v_g*dt;
    }

    if(fil_d_state == 1 || fil_d_state == 2)
    {
        eps_d1 = eps_d1 + 0.5*v_s*dt;
        eps_d2 = eps_d2 + 0.5*v_s*dt;
    }
    
}

