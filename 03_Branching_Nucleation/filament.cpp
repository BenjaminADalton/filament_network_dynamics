#include <iostream>
#include <fstream>
#include <random>
#include <cstdlib>
#include <math.h>

#include "utilities.h"
#include "filament.h"

using namespace std;

filament :: filament(){}

void filament :: filament_init(unsigned int n_dim)
{    
    for (unsigned int n = 0; n < n_dim; n++)
    {
        r_com.push_back(0.0);
        u_vec.push_back(0.0);
    }
    
    for (unsigned int n = 0; n < n_dim; n++){r_eqnm.push_back(0.0);}
    for (unsigned int n = 0; n < 3; n++){u_eqnm.push_back(0.0);}
    
    vector<double> buffer_u;
    
    if(n_dim == 2)
    {
        buffer_u.push_back(0.0);
        buffer_u.push_back(0.0);
        
        u_perp.push_back(buffer_u);
    }
    if(n_dim == 3)
    {

        buffer_u.push_back(0.0);
        buffer_u.push_back(0.0);
        buffer_u.push_back(0.0);
        
        u_perp.push_back(buffer_u);
        u_perp.push_back(buffer_u);
    }

    pbc_index[0] = 0.0;
    pbc_index[1] = 0.0;
    pbc_index[2] = 0.0;
    
    fil_id = 0;
    
}

void filament :: filament_drag(double d,double kappa)
{
    
    log_ld = log(fil_length/d);
    
    para_poly = log_ld;
    perp_poly = log_ld;
    rot_poly  = log_ld;

    gamma_para = 2.0*kappa*fil_length/para_poly;
    gamma_perp = 4.0*kappa*fil_length/perp_poly;
    gamma_rot  = (kappa*fil_length*fil_length*fil_length)/(3.0*rot_poly);
    
}

void filament :: filament_dlength(const double dt, const double v_g, const double v_s, int n_dims, double ly)
{
    if(growth_phase == 0)
    {
        double dlength = v_g*dt;
        
        fil_length = fil_length + dlength;
        
        for (int j=0; j < n_dims; j++)
        {
            r_com[j] = r_com[j] + 0.5*dlength*u_vec[j];
        }
	// if(fil_length >= 25.0e-6){growth_phase = 1;}
    }
    
    if(growth_phase == 1 || growth_phase == 2)
    {
        double dlength = -v_s*dt;
        
        fil_length = fil_length + dlength;
        
        for (int j=0; j < n_dims; j++)
        {
            r_com[j] = r_com[j] + 0.5*dlength*u_vec[j];
        }
        
        if(fil_length < fil_min){growth_phase = 2;}
    }
}







