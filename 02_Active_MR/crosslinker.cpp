#include <iostream>
#include <fstream>
#include <random>
#include <cstdlib>
#include <math.h>

#include "utilities.h"
#include "crosslinker.h"
#include "filament.h"

using namespace std;

crosslinker :: crosslinker(){}

void crosslinker :: cl_init_start(double delta_in, double ksp_in)
{
    cl_delta = delta_in;
    k_spr    = ksp_in;
}

void crosslinker :: cl_init(unsigned int n_dim, unsigned int n_fils, double delta_in, double ksp_in)
{
    epsilon_i = 0.0;
    epsilon_j = 0.0;
    
    cl_delta = delta_in;
    k_spr = ksp_in;
    
    for(unsigned int i = 0; i < n_fils - 1; i++)
    {
        cl_neigh.push_back(0);
        cl_kon_i.push_back(0.0);
    }
}

void crosslinker :: neigh_append(unsigned int fil_id)
{
    cl_neigh.push_back(fil_id);
    cl_kon_i.push_back(0.0);
}

void crosslinker :: nlist_clean()
{
    
    n_neigh = 0;
    
    cl_neigh.resize(0);
    cl_kon_i.resize(0);
}

void crosslinker :: cl_walk_poly_1_type_k(double fil_length_i, double motor_vel, double dt, int fil_state, double v_g, double v_s)
{

    if(fil_state == 0)
    {
        epsilon_i = epsilon_i + (motor_vel - 0.5*v_g)*dt;
    }
    else if(fil_state == 1 || fil_state == 2)
    {
        epsilon_i = epsilon_i + (motor_vel + 0.5*v_s)*dt;
    }
    else if(fil_state == 3)
    {
        epsilon_i = epsilon_i + motor_vel*dt;
    }
    
    if(epsilon_i >= fil_length_i/2.0){state = 0; nlist_clean(); end_state = 0;}

}

void crosslinker :: cl_walk_poly_2_type_k(double fil_length_i,double fil_length_j, double motor_vel, double dt, int fil_i_state, int fil_j_state, double v_g, double v_s)
{

    // Force-velocity relationship. If f_para_i > F_stall: stall
    if(abs(f_para_i) <= F_stall)
    {
        if(fil_i_state == 0)
        {
            epsilon_i = epsilon_i - 0.5*v_g*dt + motor_vel*dt*(1.0-abs(f_para_i/F_stall));
        }
        else if(fil_i_state == 1 || fil_i_state == 2)
        {
            epsilon_i = epsilon_i + (0.5*v_s)*dt + (motor_vel)*dt*(1.0-abs(f_para_i/F_stall));
        }
        else if(fil_i_state == 3)
        {
            epsilon_i = epsilon_i + (motor_vel)*dt*(1.0-abs(f_para_i/F_stall));
        }
    }
    else if(abs(f_para_i) > F_stall)
    {
        if(fil_i_state == 0)
        {
            epsilon_i = epsilon_i - 0.5*v_g*dt;
        }
        else if(fil_i_state == 1 || fil_i_state == 2)
        {
            epsilon_i = epsilon_i + (0.5*v_s)*dt;
        }
    }
    
    if(epsilon_i >= fil_length_i/2.0){epsilon_i = fil_length_i/2.0; end_state_i = 1;}
    
    // Force-velocity relationship. If f_para_j > F_stall: stall
    if(abs(f_para_j) <= F_stall)
    {
        if(fil_j_state == 0)
        {
            epsilon_j = epsilon_j - (0.5*v_g)*dt + (motor_vel)*dt*(1.0-abs(f_para_j/F_stall));
        }
        if(fil_j_state == 1 || fil_j_state == 2)
        {
            epsilon_j = epsilon_j + (0.5*v_s)*dt + (motor_vel)*dt*(1.0-abs(f_para_j/F_stall));
        }
        if(fil_j_state == 3)
        {
            epsilon_j = epsilon_j + (motor_vel)*dt*(1.0-abs(f_para_j/F_stall));
        }
    }
    else if(abs(f_para_j) > F_stall)
    {
        if(fil_j_state == 0)
        {
            epsilon_j = epsilon_j - (0.5*v_g)*dt;
        }
        if(fil_j_state == 1 || fil_j_state == 2)
        {
            epsilon_j = epsilon_j + (0.5*v_s)*dt;
        }
    }
    
    if(epsilon_j >= fil_length_j/2.0){epsilon_j = fil_length_j/2.0; end_state_j = 1;}
}


