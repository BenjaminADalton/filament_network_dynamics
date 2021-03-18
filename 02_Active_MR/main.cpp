// Date: 14-03-2021
// Author: Benjamin Dalton
// Current Contact: dalton@zedat.fu-berlin.de
//
// Simulate active microrheolgy of filament networks with an embedded probe
// filament and oscillating pinning field.
//
// For description, see the ReadMe.pdf file contained in this directory
// To compile and execute this simulation, run:
//      make clean
//      make
//      ./execute_dynamics
// Since there is no independent input file, any changes made to this main.cpp
// file must be re-compiled.

#include <iostream>
#include <fstream>
#include <chrono>
#include <random>
#include <cstdlib>
#include <ctime>
#include <stdio.h>
#include <math.h>
#include <string>
#include <ctime>
#include <typeinfo>

#include "filament.h"
#include "crosslinker.h"
#include "utilities.h"
#include "binding.h"
#include "pair_struct.h"
#include "pbc.h"
#include "wall.h"
#include "wall_calc.h"

using namespace std;

#define N_DIMS 3            // Number of spatial dimensions (2 0r 3)
#define BC_TYPE 3           // periodic boundar conditions (0,1,2,3 or 4 = see list below)
#define CL 1                // cross-linkers/motors (0 = off, 1 = on, 2 = just forces, skip binding/unbinding stochastics)
#define PIN_FIELD 1         // include fixed position filament structure (0 = off, 1 = on)
#define OSC_SWITCH 1        // turn on/off oscillating crystal

int main()
{
    // length of simulation
    int T_STEP = 100000000;
    
    // real time of each iteration
    const double dt = 2.5e-6;
    
    // number of filaments and crosslinkers
    int N_FILS    = 1;    // initial number of MTs - must be 1 for AMR simulations
    int N_CLS     = 1200; // total number of available CLs
    
    int N_PIN           = 16;      // total number of pins
    int z_stack         = 4;       // must be integer divisable into N_PIN
    double y_cryst      = 0.0e-6;  // place the c.o.m of pinning lattice in y
    double cryst_length = 5.0e-6;  // length of pins
    
    // active micro_rheology oscillating pinning field
    int osc_start     = 5000000;
    double osc_period = 120.0; //  n cycles-per-hour: 2pi/1hr = 2pi/3600s = 0.0017453s^-1
    double osc_amp    = 1.0e-6;
    
    // active micro_rheology force sensor
    double trap_eqm_y = 20.0e-6;
    double k_trap     = 4.0e-5;
    
    // for the filament inital state setup
    double init_cutoff = 2.0;
    double init_length = 2.0e-6;
    double nuc_fil_min = 1.0e-6;

    // box dimensions
    double L_x = 1.0e-6;
    double L_y = 30.0e-6;
    double L_z = 1.0e-6;
    
    // Nucleation and catastrophe rates
    double k_n = 6.0;     // bulk nucleation rate
    double k_c = 0.038;   // catastrophe rate
    
    // velocity for the motor protein (0.0e-6 for passive)
    double motor_vel1 = 0.00e-6; // kinesin - must be either positive or zero
    
    // MT (de)polymerization rates
    const double v_g = 0.25e-6; // 15um/min = 0.25 um/s
    const double v_s = 0.58e-6; // 35um/min = 0.58 um/s
    
    // Kinesin-like or passive
    double k_0_3D    = 250; // (m^-1*s^-1*motor^-1)
    double k_off     = 0.033;
    
    // Spring constants (3.0e-4N/m = 0.3pN/nm)
    double k_spr   = 3.0e-4;

    // Stall force for kinesin and dynein
    double F_stall   = 5.0e-12;
    
    double xz_nuc_frac = 1.1;
    double y_nuc_frac  = 1;

    // (0) uy = 1, (1) random, (2) uy = -1, (3) 3D random,  (4) XWall random : note Cryst uy=1
    int rand_or_switch = 2;
    double rand_orient = 0.0;
    
    // depoly all filaments that turn in channel
    int uvec_stab = 1;
        
    // print frequency and start time for ID trajectories
    const unsigned int t_write = 50000;
    const unsigned int t_print = 500000;
    
    // print frequency and start time for ID trajectories
    int t_write_2   = 10000;
    int start_write = 10000000;
    
    // neighbour list update cycle lengths
    const unsigned int n_list_update       = 750;
    const unsigned int n_list_update_outer = 7500;
    const unsigned int wall_update         = 5000;
    
    // cutoff distance for nestled inner and outer neighbour lists
    const double outer_nlist_cut = 5.0e-07; // 30*sig_min = 4.2e-07
    const double inner_nlist_cut = 1.5e-07; // 15*sig_min = 2.1e-07
    
    // sampling frequencies for processes slower than thermal motion
    const unsigned int cl_neigh_update      = 5000;
    const unsigned int sample_inc_kn        = 10000;
    const unsigned int sample_inc_kcat      = 10000;
    const unsigned int sample_inc_on1       = 10000;
    const unsigned int sample_inc_off1      = 10000;
    const unsigned int sample_inc_pair_on   = 10000;
    const unsigned int sample_inc_pair_off  = 10000;
    
    // check filaments for the depolymerization state
    const unsigned int t_depoly = 5000;
    
    // parameters for fluid interactions - diffusion and drag
    const double k_B     = 1.381e-23;       // Boltzmanns constant
    const double T       = 300;             // temperature
    const double epsilon = k_B*T;           // thermal energy
    const double kTt_2   = 2.0*k_B*T*dt;
    const double pi      = 3.14159265359;
    const double eta     = 0.1;             // solvent viscosity
    const double kappa   = pi*eta;
    const double lambda  = k_B*T/(2.0*pi*eta);

    // filament dimensions and WCA potential parameters
    const double d          = 5.0e-8; // Diameter, not sigma!!!!!
    double sig_shift        = 2.5e-8;
    double beta             = 0.9;
    double sig_shift_wall   = 1.5e-8;
    double beta_wall        = 0.9;
    double sig_shift_min    = sig_shift*pow(2.0,1.0/6.0);
    double sig_min_wall     = sig_shift_wall*pow(2.0,1.0/6.0);
    double beta_sig         = beta*sig_shift;
    double one_on_beta      = sig_shift/beta_sig;
    double beta_sig_wall    = beta_wall*sig_shift_wall;
    double one_on_beta_wall = sig_shift_wall/(beta_sig_wall);

    const double wall_cut  = 15.0; // in units of sig_min_wall
    
    // cross-link/motor parameters
    double eta_sites   = 1.0*1.25e8;    // binding site concentration
    double sig_cl      = 0.3e-6;        // cl search cutoff
    double k_spr_kbT   = k_spr/(k_B*T);
    double cl_delta    = 1.38e-9;
    
    // use for spring force calculation
    double d_ij [3] = {0.0};
    double d_0_vec [3] = {0.0};
    double d_0         = 8.0e-8;
    double d_ij_mag    = 0.0;
    
    // various parameters used and reused throughout
    int n, nt_cl1_1 = 0, nt_cl1_2 = 0, count=0;
    unsigned int i,j,k,l,m;
    double F_ij, F_i_wall, sig_wij, sig_wij_3, sig_wij_6, sig_wij_12;
    double u_mod, umod, u_iw_dot, u_jw_dot, u_id_dot, u_jd_dot;
    double sqrt2 = sqrt(2.0), sqrt2pi = sqrt(2.0*pi), rand1 = 0.0;
    double k_on_dim = 0.0, f_pi = 0.0, f_pj = 0.0, osc_inc = 0.0;
    
    int loop_int = 0, nuc_block = 0, nuc_index = 0;
    double length_acc = 0.0;
    double upper_length = 200.0e-6;
        
    // declare the filament and cross-linker objects
    vector<filament> mt(0);
    vector<crosslinker> cl(N_CLS);
    
    // basis vectors used in the stochastic dynamics
    double *e1_2d = new double[2]; e1_2d[0] = 1; e1_2d[1] = 0;
    double *e1_3d = new double[3]; e1_3d[0] = 1; e1_3d[1] = 0; e1_3d[2] = 0;
    double *e2_3d = new double[3]; e2_3d[0] = 0; e2_3d[1] = 1; e2_3d[2] = 0;
    
    // used for wall boundary conditions (don't need otherwise, but here they are anyway)
    vector<double> uwall1(3); uwall1[0] = 1.0; uwall1[1] = 0.0; uwall1[2] = 0.0;
    vector<double> uwall2(3); uwall2[0] = -1.0; uwall2[1] = 0.0; uwall2[2] = 0.0;
    vector<double> uwall3(3); uwall3[0] = 0.0; uwall3[1] = 0.0; uwall3[2] = 1.0;
    vector<double> uwall4(3); uwall4[0] = 0.0; uwall4[1] = 0.0; uwall4[2] = -1.0;

    // Use for the wall B.Cs 3D
    vector<double> ri_s3d(3);
    vector<double> uwall_3d(3);
    double L_wall_3D = 20.0*L_y;
    
    double tot_fil_length = 0.0;
    double av_tot_fil_length = 0.0;
    
    double gamma_array_3d[3][3];
    double cl_check [3] = {0.0};
    
    // if(BC_TYPE == 1){pbc pbc1; pbc1.pbc_constructor(N_DIMS);}
    pbc pbc;
    pbc.pbc_constructor(N_DIMS, L_x, L_y, L_z);
    
    // declare some objects
    wall wall;
    wall_calc wall_calc;
    wall_calc.wall_calc_init(N_DIMS,epsilon,sig_shift_wall,beta_sig_wall,one_on_beta_wall);

    // initiate binding class for motor binding kinetics
    binding binding;
    binding.binding_init(N_DIMS,eta_sites);

    // constants used for the LJ-WCA Fij calculation
    double f_beta  = -24.0*(epsilon/sig_shift)*(2.0*pow(one_on_beta,13.0) - pow(one_on_beta,7.0));
    double df_beta = 24.0*(epsilon/sig_shift)*(26.0*pow(sig_shift,13.0)/pow(beta_sig,14.0) - 7.0*pow(sig_shift,7.0)/pow(beta_sig,8.0));
        
    int cl_nst1 = 0, cl_nst2 = 0;
    
    // lists used for book-keeping the CL binding kinetics
    vector<int> cl_st1_list; for(i=0; i < N_CLS; i++){cl_st1_list.push_back(0);}
    vector<int> cl_st2_list; for(i=0; i < N_CLS; i++){cl_st2_list.push_back(0);}
    
    int filament_id_count = 0;
    
    // pass motor parameters to CL objects
    for(i=0;i < N_CLS; i++)
    {
        cl[i].type      = 1;
        cl[i].k_spr     = k_spr;
        cl[i].k_spr_kbT = k_spr_kbT;
        cl[i].cl_delta  = cl_delta;
        cl[i].k_on      = k_0_3D;
        cl[i].k_off     = k_off;
        cl[i].F_stall   = F_stall;
    }
    
    // construct some vectors used throughout
    vector<double> box_vec,ri_image,rj_image,rcl_image,ri_pair_image;
    
    for(unsigned int n = 0; n < N_DIMS; n++)
    {
        ri_image.push_back(0.0);
        rj_image.push_back(0.0);
        rcl_image.push_back(0.0);
    }
    
    box_vec.push_back(L_x); box_vec.push_back(L_y); box_vec.push_back(L_z);

    // Set up vectors for Gram Schmidt
    vector<double> u_perp_3D_1(3,0.0);
    vector<double> u_perp_3D_2(3,0.0);
    
    vector<double> network_fils;
    vector<double> network_fils_filtered;
    
    // Set up clock and random number generators
    time_t timev; time(&timev);
    
    // Set up clock properties:
    clock_t tic_full = clock();
    time_t now       = time(0);
    tm*ltm           = localtime(&now);
    int srand_int    = 195;
    // int srand_int = 1 + (ltm->tm_hour)*(ltm->tm_min)*(ltm->tm_sec);
    char* time_date = ctime(&now);
    
    mt19937 generator (srand_int);
    srand             (srand_int);
    
    uniform_real_distribution<double> dis(0.0, 1.0);
    
    // Random number generator seed
    cout << endl << "srand_int (random number seed): " << srand_int << endl;
    
    // ************************* Set up pairs and create the inital state configuration *************************************************

    pair_struct pair;
    pair.pair_init(N_DIMS);
    pair.outer_nlist_init(N_DIMS);
    pair.inner_nlist_init(N_DIMS);
    pair.nlist_init_state(N_DIMS);

    // build the pinning field structure
    if(PIN_FIELD == 1)
    {
        for(k=0; k < z_stack; k++)
        {
            for(i=0; i < N_PIN/z_stack; i++)
            {
                filament mt_new;
                mt_new.filament_init(N_DIMS);
                mt_new.fil_length = cryst_length;
                
                mt_new.r_com[0] = (L_x/(2.0*double(N_PIN/z_stack)) + double(z_stack*i)*L_x/(double(N_PIN)));
                mt_new.r_com[1] = y_cryst;
                mt_new.r_com[2] = (k+1)*L_z/double(z_stack + 1);

                mt_new.u_vec[0] = 0.0;
                mt_new.u_vec[1] = 1.0;
                mt_new.u_vec[2] = 0.0;
                
                mt_new.growth_phase = 3; // no polymerization dynamics
                mt.push_back(mt_new);
            }
        }
    }
    
    
    // use fro loop control
    int N_accum = N_PIN; N_FILS = N_FILS + N_PIN;
    
    loop_int = 1;
    
    // initialize non-pin filaments, check for intercepts and reject
    for(i = N_PIN; i < N_FILS; i++)
    {
        // create the AMR probe filament
        filament mt_new;
        mt_new.filament_init(N_DIMS);
        mt_new.fil_length = init_length;
        
        mt_new.r_com[0] = L_x/2.0 + L_x/8.0;
        mt_new.r_com[1] = trap_eqm_y;
        mt_new.r_com[2] = L_z/2.0;
        mt_new.u_vec[0] = 0.0;
        mt_new.u_vec[1] = -1.0;
        mt_new.u_vec[2] = 0.0;
        
        mt_new.growth_phase = 3;
        
        filament_id_count = filament_id_count + 1; mt_new.fil_id = filament_id_count;
        mt.push_back(mt_new); N_accum = N_accum + 1;
    }
    
    // ******************************************************************************************************************
    
    // prepare the drag coefficients for each filament (update every time when polymerizing)
    for (j=0; j < N_FILS; j++){mt[j].filament_drag(d,kappa);}

    // open text files to be used for data writing
    ofstream trajectory;
    trajectory.open ("trajectory.txt");
    ofstream trajectory_trap;
    trajectory_trap.open ("trajectory_trap.txt");
    ofstream cl_trajectory;
    cl_trajectory.open ("trajectory_cl.txt");
    ofstream id_trajectory;
    id_trajectory.open ("trajectory_id.txt");
    ofstream id_orientation;
    id_orientation.open ("orientation_id.txt");
    ofstream cl_nt;
    cl_nt.open ("trajectory_cl_nt.txt");
    ofstream cl_id;
    cl_id.open ("cl_pairs.txt");
    
    for(j=0; j < N_FILS; j++){tot_fil_length = tot_fil_length + mt[j].fil_length;}
    
    std::cout << "Start: number of CLs:   " << N_CLS << "    " << cl_nst1 << "    " << cl_nst2 << endl;
    
    cout << endl << "Start Main Loop:" << endl;
    
    for(i=0;i<T_STEP;i++)
    {
        //****************************** Mass turn-over ************************************************
        
        // bulk nucleation of new filament
        if(i % sample_inc_kn == 0)
        {
            if(dis(generator) <= sample_inc_kn*k_n*dt)
            {
                filament mt_new;
                mt_new.filament_init(N_DIMS);
                mt_new.fil_length = nuc_fil_min;
                
                loop_int = 0; nuc_block = 0; int nuc_index = 0;
                
                while(nuc_block == 0)
                {
                    nuc_index = nuc_index + 1; nuc_block = 1; u_mod = 0.0;
                    
                    for (unsigned int n = 0; n < N_DIMS; n++)
                    {
                        mt_new.r_com[n] = (box_vec[n]/xz_nuc_frac)*dis(generator);
                    }

                    // shift by 2.5e-6 to fill the pinning region
                    mt_new.r_com[1] = (box_vec[1] + 2.5e-6)*dis(generator) - 2.5e-6;
                    mt_new.u_vec[0] = 0; mt_new.u_vec[2] = 0;
                    
                    // Choosing the orientations as parallel or antiparallel
                    if(rand_or_switch == 0){rand_orient = 0.5;}
                    else if(rand_or_switch == 1){rand_orient = dis(generator) - 0.5;}
                    else if(rand_or_switch == 2){rand_orient = - 0.5;}
                    
                    if(rand_orient <= 0.0){mt_new.u_vec[1] = -1;}
                    else if(rand_orient > 0.0){mt_new.u_vec[1] = 1;}

                    for (unsigned int n = 0; n < N_DIMS; n++){u_mod = u_mod + mt_new.u_vec[n]*mt_new.u_vec[n];}
                    u_mod = sqrt(u_mod);
                    for (unsigned int n = 0; n < N_DIMS; n++){mt_new.u_vec[n] = mt_new.u_vec[n]/u_mod;}

                    // reject nucleation when a new filament intercepts a pre-existing filament
                    for(unsigned int k = 0; k < N_FILS; k++)
                    {
                        
                        for(n=0; n<N_DIMS; n++)
                        {
                            ri_image[n] = mt[k].r_com[n];
                            rj_image[n] = mt_new.r_com[n];
                        }
                        
                        ri_image = pair.image_calc_BC_pair(BC_TYPE, N_DIMS, ri_image, rj_image, box_vec);
                        
                        pair.least_sep_calc_single(ri_image, rj_image, mt[k].u_vec,mt_new.u_vec, mt[k].fil_length,mt_new.fil_length,N_DIMS);
                        
                        if(pair.w_ij_mag_single <= 1.2*sig_shift_min)
                        {
                            nuc_block = 0;
                            break;
                        }
                    }
                    
                    if(nuc_index == 1000){cout << "Nuc LooP 1000!!!!!" << endl;}
                    if(nuc_index > 20000){abort();}
                }
                
                // Add new filament to the mt vector
                filament_id_count = filament_id_count + 1; mt_new.fil_id = filament_id_count; mt_new.fil_time_nuc = i;
                
                mt.push_back(mt_new);
                N_FILS = N_FILS + 1;
                
                int new_mt = N_FILS - 1;
                
                // add new pairs to inner and outer lists
                for(unsigned int k = 0; k < N_FILS-1; k++)
                {
                    for(n=0; n<N_DIMS; n++)
                    {
                        ri_image[n] = mt[new_mt].r_com[n];
                        rj_image[n] = mt[k].r_com[n];
                    }
                    
                    ri_image = pair.image_calc_BC_pair(BC_TYPE, N_DIMS, ri_image, rj_image, box_vec);
                    
                    pair.least_sep_calc_single(ri_image, rj_image, mt[new_mt].u_vec ,mt[k].u_vec, mt[new_mt].fil_length,mt[k].fil_length,N_DIMS);
                    
                    if(pair.w_ij_mag_single <= inner_nlist_cut)
                    {
                        pair.inner_nlist_build(k,new_mt);
                    }
                    if(pair.w_ij_mag_single <= outer_nlist_cut)
                    {
                        pair.outer_nlist_build(k,new_mt);
                    }
                }
                
                // Update the wall neighbor list: add new filament if proximate
                if(BC_TYPE == 2 || BC_TYPE == 3) // herepbc_wall
                {
                    ri_s3d[0]  = 0.0;
                    ri_s3d[1]  = mt[new_mt].r_com[1];
                    ri_s3d[2]  = mt[new_mt].r_com[2];
                    
                    // must normalise:
                    uwall_3d[0] = 0.0;
                    uwall_3d[1] = mt[new_mt].u_vec[1];
                    uwall_3d[2] = mt[new_mt].u_vec[2];
                    
                    umod = sqrt(uwall_3d[1]*uwall_3d[1] + uwall_3d[2]*uwall_3d[2]);
                    
                    uwall_3d[1] = uwall_3d[1]/umod;
                    uwall_3d[2] = uwall_3d[2]/umod;
                    
                    wall_calc.wall_calc_class(mt[new_mt].r_com, ri_s3d, mt[new_mt].u_vec, uwall_3d, mt[new_mt].fil_length, L_wall_3D, N_DIMS,new_mt);
                    
                    if(wall_calc.w_ij_mag_out <= wall_cut*sig_min_wall){wall.wall_x1_add(new_mt);}
                    
                    ri_s3d[0]  = L_x;
                    
                    wall_calc.wall_calc_class(mt[new_mt].r_com, ri_s3d, mt[new_mt].u_vec,uwall_3d,mt[new_mt].fil_length,L_wall_3D,N_DIMS,new_mt);
                    
                    if(wall_calc.w_ij_mag_out <= wall_cut*sig_min_wall){wall.wall_x2_add(new_mt);}
                    
                    ri_s3d[0]  = mt[new_mt].r_com[0];
                    ri_s3d[1]  = mt[new_mt].r_com[1];
                    ri_s3d[2]  = 0.0;
                    
                    // must normalise:
                    uwall_3d[0] = mt[new_mt].u_vec[0];
                    uwall_3d[1] = mt[new_mt].u_vec[1];
                    uwall_3d[2] = 0.0;
                    
                    umod        = sqrt(uwall_3d[0]*uwall_3d[0] + uwall_3d[1]*uwall_3d[1]);
                    uwall_3d[0] = uwall_3d[0]/umod;
                    uwall_3d[1] = uwall_3d[1]/umod;
                    
                    wall_calc.wall_calc_class(mt[new_mt].r_com, ri_s3d, mt[new_mt].u_vec,uwall_3d,mt[new_mt].fil_length,L_wall_3D,N_DIMS,new_mt);
                    
                    if(wall_calc.w_ij_mag_out <= wall_cut*sig_min_wall){wall.wall_z1_add(new_mt);}
                    
                    ri_s3d[2]  = L_z;
                    
                    wall_calc.wall_calc_class(mt[new_mt].r_com, ri_s3d, mt[new_mt].u_vec,uwall_3d,mt[new_mt].fil_length,L_wall_3D,N_DIMS,new_mt);
                    
                    if(wall_calc.w_ij_mag_out <= wall_cut*sig_min_wall){wall.wall_z2_add(new_mt);}
                }
                
                // NOTE: have removed the update of CL list to include new MT.
                // can wait until the next cycle through CL section, there it
                // will be added
            }
        }
        
        // catastrophe events
        if(i % sample_inc_kcat == 0)
        {
            for(n=0; n < N_FILS; n++)
            {
                if(mt[n].growth_phase == 0)
                {
                    if(dis(generator) <= sample_inc_kcat*k_c*dt)
                    {
                        mt[n].growth_phase = 1;
                    }
                }
            }
        }
        
        // Remove filaments that have completely depolymerised and update pair structures
        if (i % t_depoly == 0)
        {
            for (n=N_PIN; n < N_FILS; n++)
            {
                if(mt[n].growth_phase == 2)
                {
                    for (j=0; j < cl_nst1; j++)
                    {
                        if(cl[cl_st1_list[j]].fil_i == n)
                        {
                            cl[cl_st1_list[j]].state = 0;
                            cl[cl_st1_list[j]].nlist_clean();
                            cl[cl_st1_list[j]].end_state = 0;
                        }
                    }
                    
                    for (j=0; j < cl_nst2; j++)
                    {
                        if(cl[cl_st2_list[j]].fil_i == n)
                        {
                            cl[cl_st2_list[j]].state = 1; // NOTE: need to recalculate the neighbour list for this state-1 cl
                            cl[cl_st2_list[j]].nlist_clean();
                            cl[cl_st2_list[j]].epsilon_i = cl[cl_st2_list[j]].epsilon_j;
                            cl[cl_st2_list[j]].fil_i = cl[cl_st2_list[j]].fil_j;
                            cl[cl_st2_list[j]].epsilon_j = 0.0;
                            cl[cl_st2_list[j]].fil_j = 0;
                            cl[cl_st2_list[j]].end_state_i = 0;
                            cl[cl_st2_list[j]].end_state_j = 0;

                            if(CL == 1) // Note (999):
                            {
                            
                                cl_check [0] = 0.0; cl_check [1] = 0.0; cl_check [2] = 0.0;
                                
                                for(k=0; k<N_DIMS; k++)
                                {
                                    cl_check[k] = mt[cl[cl_st2_list[j]].fil_i].r_com[k] + cl[cl_st2_list[j]].epsilon_i*mt[cl[cl_st2_list[j]].fil_i].u_vec[k];
                                }
                                
                                // re - create nlist for the cl returned to state 1
                                int cl_inc_off = 0;
                                
                                for (m = 0; m < N_FILS; m++)
                                {
                                    if(m != cl[cl_st2_list[j]].fil_i)
                                    {
                                        
                                        rcl_image = binding.image_calc_BC_cl(BC_TYPE, N_DIMS, mt[m].r_com, mt[cl[cl_st2_list[j]].fil_i].r_com, box_vec);
                                        
                                        if(binding.nlist_cl(rcl_image,cl_check,mt[m].u_vec,mt[m].fil_length,N_DIMS) <= sig_cl)
                                        {
                                            cl[cl_st2_list[j]].neigh_append(m);
                                            cl_inc_off = cl_inc_off + 1;
                                        }
                                    }
                                }
                                
                                cl[cl_st2_list[j]].n_neigh = cl_inc_off;
                            }
                        }
                        if(cl[cl_st2_list[j]].fil_j == n)
                        {
                            cl[cl_st2_list[j]].state = 1; // NOTE: need to recalculate the neighbour list for this state-1 cl
                            cl[cl_st2_list[j]].nlist_clean();
                            cl[cl_st2_list[j]].epsilon_j = 0.0;
                            cl[cl_st2_list[j]].fil_j = 0;
                            cl[cl_st2_list[j]].end_state_j = 0;
                            
                            if(CL == 1)
                            {
                                cl_check [0] = 0.0; cl_check [1] = 0.0; cl_check [2] = 0.0;
                                
                                for(k=0; k<N_DIMS; k++)
                                {
                                    cl_check[k] = mt[cl[cl_st2_list[j]].fil_i].r_com[k] + cl[cl_st2_list[j]].epsilon_i*mt[cl[cl_st2_list[j]].fil_i].u_vec[k];
                                }
                                
                                // re - create nlist for the cl returned to state 1
                                int cl_inc_off = 0;

                                for (m = 0; m < N_FILS; m++)
                                {
                                    if(m != cl[cl_st2_list[j]].fil_i)
                                    {
                                        rcl_image = binding.image_calc_BC_cl(BC_TYPE, N_DIMS, mt[m].r_com, mt[cl[cl_st2_list[j]].fil_i].r_com, box_vec);
                                        
                                        if(binding.nlist_cl(rcl_image,cl_check,mt[m].u_vec,mt[m].fil_length,N_DIMS) <= sig_cl)
                                        {
                                            cl[cl_st2_list[j]].neigh_append(m);
                                            cl_inc_off = cl_inc_off + 1;
                                        }
                                    }
                                }
                                cl[cl_st2_list[j]].n_neigh = cl_inc_off;
                                
                            }
                        }
                    }
                    
                    // If location n is the final mt[] element
                    // remove from last location and clean up pair list
                    if(n == N_FILS - 1)
                    {
                        mt.pop_back();
                        
                        for(j=0; j < pair.n_inner_pair; j++)
                        {
                            if(pair.inner_nlist[j][0] ==  N_FILS-1 || pair.inner_nlist[j][1] ==  N_FILS-1)
                            {
                                pair.inner_nlist.erase(pair.inner_nlist.begin() + j);
                                pair.n_inner_pair = pair.n_inner_pair - 1;
                                j = j-1;
                            }
                        }

                        for(j=0; j < pair.n_outer_pair; j++)
                        {
                            if(pair.outer_nlist[j][0] ==  N_FILS-1 || pair.outer_nlist[j][1] ==  N_FILS-1)
                            {
                                pair.outer_nlist.erase(pair.outer_nlist.begin() + j);
                                pair.n_outer_pair = pair.n_outer_pair - 1;
                                j = j-1;
                            }
                        }

                    }
                    
                    // location n is NOT the final mt[] element
                    // more complicated. must restructure lists
                    else if(n < N_FILS - 1)
                    {

                        for(j=0; j < N_DIMS; j++)
                        {
                            mt[n].r_com[j] = mt[N_FILS-1].r_com[j];
                            mt[n].u_vec[j] = mt[N_FILS-1].u_vec[j];
                        }

                        mt[n].fil_length   = mt[N_FILS-1].fil_length;
                        mt[n].fil_id       = mt[N_FILS-1].fil_id;
                        mt[n].fil_time_nuc = mt[N_FILS-1].fil_time_nuc;
                        mt[n].growth_phase = mt[N_FILS-1].growth_phase;

                        for(j=0; j < cl_nst1; j++)
                        {
                            if(cl[cl_st1_list[j]].fil_i == N_FILS-1){cl[cl_st1_list[j]].fil_i = n;}
                        }

                        for(j=0; j < cl_nst2; j++)
                        {
                            if(cl[cl_st2_list[j]].fil_i == N_FILS-1){cl[cl_st2_list[j]].fil_i = n;}
                            if(cl[cl_st2_list[j]].fil_j == N_FILS-1){cl[cl_st2_list[j]].fil_j = n;}
                        }

                        // first remove all instances of the depolymerized filament id from the id lists
                        // second substitute all location with the (N_FILS-1)th element with the exchanging
                        // ID of the old one
                        for(j=0; j < pair.n_outer_pair; j++)
                        {
                            if(pair.outer_nlist[j][0] ==  n || pair.outer_nlist[j][1] ==  n)
                            {
                                pair.outer_nlist.erase(pair.outer_nlist.begin() + j);
                                pair.n_outer_pair = pair.n_outer_pair - 1;
                                j=j-1;
                            }
                        }

                        for(j=0; j < pair.n_outer_pair; j++)
                        {
                            if(pair.outer_nlist[j][0] ==  N_FILS - 1 ){pair.outer_nlist[j][0] = n;}
                            if(pair.outer_nlist[j][1] ==  N_FILS - 1 ){pair.outer_nlist[j][1] = n;}
                        }
                        
                        // first erase all elements with id from the removed filament (eg - remove all with id 0)
                        for(j=0; j < pair.n_inner_pair; j++)
                        {
                            if(pair.inner_nlist[j][0] ==  n || pair.inner_nlist[j][1] ==  n)
                            {
                                pair.inner_nlist.erase(pair.inner_nlist.begin() + j);
                                pair.n_inner_pair = pair.n_inner_pair - 1;
                                j=j-1;
                            }
                        }

                        for(j=0; j < pair.n_inner_pair; j++)
                        {
                            if(pair.inner_nlist[j][0] ==  N_FILS - 1 ){pair.inner_nlist[j][0] = n;}
                            if(pair.inner_nlist[j][1] ==  N_FILS - 1 ){pair.inner_nlist[j][1] = n;}
                        }
                        
                        mt.pop_back();
                    }

                    if(BC_TYPE == 2 || BC_TYPE == 3)
                    {
                        wall.wall_3d_remove(n,N_FILS-1);
                    }

                    N_FILS = N_FILS - 1;
                    pair.n_pair = (N_FILS*(N_FILS-1))/2;

                    cl_nst1 = binding.cl_state_count(cl, N_CLS, 1, cl_st1_list);
                    cl_nst2 = binding.cl_state_count(cl, N_CLS, 2, cl_st2_list);
                }
            }

            if(uvec_stab == 1)
            {
                for (n=N_PIN; n < N_FILS; n++)
                {
                    if(mt[n].u_vec[1] > 0)
                    {
                        mt[n].growth_phase = 1;
                        cout << i << " - " << n << " - uvec depoly. Length: " << mt[n].fil_length << endl;
                    }
                }
            }
        }
        
        //*************** cross-linker kinetics - binding, unbinding and neighbour lists ***************
        
        if(CL == 1)
        {
            // ********************  create state 1 bingings ********************
            //
            // We assume un-bound cross-links/motors are uniformly distributed in volume. Here we calculate
            // is an event occures where an unbound cross-linker/motor will be taken out of the environment
            // and located in state-1 somewhere on a filament. Location selection is from a uniform distribution
            // along the total length of all filaments. See cross-linker binding stochastics, 2nd paragraph
            
            // create a state-1 binding
            if(i % sample_inc_on1 == 0)
            {
                // Initial filament binding
                tot_fil_length    = 0.0;
                
                for(j=0; j < N_FILS; j++){tot_fil_length = tot_fil_length + mt[j].fil_length;}
                
                // finding freely unbound motor (0-state) to place new 1-state. Create MT neighbour list
                for(j=0; j < N_CLS; j++)
                {
                    if(cl[j].state == 0)
                    {
                        rand1 = dis(generator);
                        k_on_dim = 0.0;
                
                        // rate of single filament binding
                        k_on_dim = cl[j].k_on*tot_fil_length*dt*sample_inc_on1;
                        
                        // binding event will occur:
                        if(rand1 <= k_on_dim)
                        {
                            cl[j].n_neigh   = 0;
                            cl[j].nlist_clean();
                            cl[j].state     = 1;
                            
                            double rand2 = dis(generator);
                            double l_r2 = rand2*tot_fil_length;
                            double length_acc = 0;
                            
                            // calculate total filament location
                            for(n=0; n < N_FILS; n++)
                            {
                                length_acc = length_acc + mt[n].fil_length;
                                
                                if(l_r2 < length_acc)
                                {
                                    cl[j].fil_i     = n;
                                    cl[j].epsilon_i = mt[n].fil_length/2.0 - (length_acc-l_r2);
                                    break;
                                }
                            }
                            
                            cl_check [0] = 0.0; cl_check [1] = 0.0; cl_check [2] = 0.0;
                            
                            for(n=0; n<N_DIMS; n++)
                            {
                                cl_check[n] = mt[cl[j].fil_i].r_com[n] + cl[j].epsilon_i*mt[cl[j].fil_i].u_vec[n];
                            }
                            
                            // once cl is placed in state-1 we need to create a list of all filaments with which it can potential interact
                            // and form a pair state binding
                            int cl_inc = 0;
                            for (m=0; m < N_FILS; m++)
                            {
                                if(m != cl[j].fil_i)
                                {
                                    // rcl_image = binding.image_calc_BC_cl(BC_TYPE, N_DIMS, ri_com, rj_com, box_vec);
                                    rcl_image = binding.image_calc_BC_cl(BC_TYPE, N_DIMS, mt[m].r_com, mt[cl[j].fil_i].r_com, box_vec);
                                    
                                    if(binding.nlist_cl(rcl_image,cl_check,mt[m].u_vec,mt[m].fil_length,N_DIMS) <= sig_cl)
                                    {
                                        cl[j].neigh_append(m);
                                        cl_inc = cl_inc + 1;
                                    }
                                }
                            }
                            cl[j].n_neigh = cl_inc;
                            cl_nst1 = cl_nst1 + 1;
                        }
                    }
                }
            
                // Update the cl1 and cl2 list and counts
                cl_nst1 = binding.cl_state_count(cl, N_CLS, 1, cl_st1_list);
                cl_nst2 = binding.cl_state_count(cl, N_CLS, 2, cl_st2_list);
                
            }
            
            // remove a state-1 binding
            if(i % sample_inc_off1 == 0)
            {
                // ******************** Un-binding from single filament ********************
                // loop though all state-1 filaments and sample whether it unbinds back into free diffusing state
                for(j=0; j < cl_nst1; j++)
                {
                    if(dis(generator) <= cl[cl_st1_list[j]].k_off*dt*sample_inc_off1)
                    {
                        cl[cl_st1_list[j]].state = 0;
                        cl[cl_st1_list[j]].nlist_clean();
                    }
                }
                cl_nst1 = binding.cl_state_count(cl, N_CLS, 1, cl_st1_list);
            }
            
            // periodically update the cl-neighbour list
            if(i % cl_neigh_update == 0)
            {
                for (j=0; j < cl_nst1; j++)
                {
                    int cl_inc = 0;
                    
                    cl[cl_st1_list[j]].nlist_clean();
                    
                    cl_check [0] = 0.0; cl_check [1] = 0.0; cl_check [2] = 0.0;
                    
                    for(n=0; n<N_DIMS; n++)
                    {
                        cl_check[n] = mt[cl[cl_st1_list[j]].fil_i].r_com[n] + cl[cl_st1_list[j]].epsilon_i*mt[cl[cl_st1_list[j]].fil_i].u_vec[n];
                    }
                    
                    for (m = 0; m < N_FILS; m++)
                    {
                        if(m != cl[cl_st1_list[j]].fil_i)
                        {
                            rcl_image = binding.image_calc_BC_cl(BC_TYPE, N_DIMS, mt[m].r_com, mt[cl[cl_st1_list[j]].fil_i].r_com, box_vec);
                            
                            if(binding.nlist_cl(rcl_image,cl_check,mt[m].u_vec,mt[m].fil_length,N_DIMS) <= sig_cl)
                            {
                                cl[cl_st1_list[j]].neigh_append(m);
                                cl_inc = cl_inc + 1;
                            }
                        }
                    }
                    cl[cl_st1_list[j]].n_neigh = cl_inc;
                }
            }
            
            // create a state-2 binding
            if(i % sample_inc_pair_on == 0)
            {
                for(j=0; j < cl_nst1; j++)
                {
                    double k_on = 0.0;
                    
                    // Make loop over all filaments that are close that the CL can extend to reach it
                    // cl[X].n_neigh list holds these ID's for state-1 cross-linkers.
                    for (m=0; m < cl[cl_st1_list[j]].n_neigh; m++)
                    {
                        rcl_image = binding.image_calc_BC_cl(BC_TYPE, N_DIMS, mt[cl[cl_st1_list[j]].cl_neigh[m]].r_com, mt[cl[cl_st1_list[j]].fil_i].r_com, box_vec);
                                                
                        cl[cl_st1_list[j]].cl_kon_i[m] = binding.kon_calc(mt[cl[cl_st1_list[j]].fil_i], mt[cl[cl_st1_list[j]].cl_neigh[m]], rcl_image, cl[cl_st1_list[j]], N_DIMS, d_0);
                    
                        k_on = k_on + cl[cl_st1_list[j]].cl_kon_i[m];
                    }
                    
                    double rand_2state = dis(generator);
                    
                    // Use kon to determine if pair-binding will occur for cross-linker j (in state-1). See Eq 4 and 6 in document.
                    if(rand_2state < cl[cl_st1_list[j]].k_off*dt*eta_sites*sample_inc_pair_on*k_on)
                    {
                        // event will occur. use stored cl_kon_i for that cross-linker to find which filament form it's neighbor set
                        // will be involved in the pair.
                        double k_rand = dis(generator)*k_on;
                        double k_accum = 0.0;
                        
                        for(m=0; m < cl[cl_st1_list[j]].n_neigh; m++)
                        {
                            k_accum = k_accum + cl[cl_st1_list[j]].cl_kon_i[m];
                            
                            // The pair filament has been selected, now find the location on the filament to place second
                            // cross-linker head
                            if(k_rand < k_accum)
                            {
                                cl[cl_st1_list[j]].state = 2;
                                cl[cl_st1_list[j]].fil_j = cl[cl_st1_list[j]].cl_neigh[m];
                                
                                // create cross-linker position vector (need to change later, inefficient double calculations)
                                
                                cl_check [0] = 0.0; cl_check [1] = 0.0; cl_check [2] = 0.0;
                                
                                for(n=0; n<N_DIMS; n++)
                                {
                                    cl_check[n] = mt[cl[cl_st1_list[j]].fil_i].r_com[n] + cl[cl_st1_list[j]].epsilon_i*mt[cl[cl_st1_list[j]].fil_i].u_vec[n];
                                }
                                
                                rcl_image = binding.image_calc_BC_cl(BC_TYPE, N_DIMS, mt[cl[cl_st1_list[j]].fil_j].r_com, mt[cl[cl_st1_list[j]].fil_i].r_com, box_vec);
                                
                                cl[cl_st1_list[j]].epsilon_j = binding.epsilon_j_calc(mt[cl[cl_st1_list[j]].fil_i], mt[cl[cl_st1_list[j]].fil_j], rcl_image, cl[cl_st1_list[j]], N_DIMS, m, dis(generator),d_0);

                                cl[cl_st1_list[j]].nlist_clean();
                                
                                break;
                            }
                        }
                    }
                    
                }
                cl_nst1 = binding.cl_state_count(cl, N_CLS, 1, cl_st1_list);
                cl_nst2 = binding.cl_state_count(cl, N_CLS, 2, cl_st2_list);
            }
            
            // remove a state-2 binding
            if(i % sample_inc_pair_off == 0)
            {
                if(cl_nst2 > 0)
                {
                    for(j=0; j < cl_nst2; j++)
                    {
                        
                        rcl_image = binding.image_calc_BC_cl(BC_TYPE, N_DIMS, mt[cl[cl_st2_list[j]].fil_i].r_com, mt[cl[cl_st2_list[j]].fil_j].r_com, box_vec);
                        
                        for(n=0; n<N_DIMS; n++)
                        {
                            d_ij[n] = (rcl_image[n] + cl[cl_st2_list[j]].epsilon_i*mt[cl[cl_st2_list[j]].fil_i].u_vec[n]) - (mt[cl[cl_st2_list[j]].fil_j].r_com[n] + cl[cl_st2_list[j]].epsilon_j*mt[cl[cl_st2_list[j]].fil_j].u_vec[n]);
                        }
                        
                        d_ij_mag = 0.0;
                        
                        for(n=0; n<N_DIMS; n++)
                        {
                            d_ij_mag = d_ij_mag + d_ij[n]*d_ij[n];
                        }
                        
                        d_ij_mag = sqrt(d_ij_mag);

                        double rand_2_1 = dis(generator);
                        
                        if(rand_2_1 < dt*sample_inc_pair_off*cl[cl_st2_list[j]].k_off*exp(cl[cl_st2_list[j]].k_spr_kbT*cl[cl_st2_list[j]].cl_delta*(d_ij_mag - d_0)))
                        {
                            cl[cl_st2_list[j]].state = 1;
                            cl[cl_st2_list[j]].nlist_clean();
                            
                            double bin_rand = dis(generator) - 0.5;
                            
                            if(bin_rand <= 0.0)
                            {
                                cl[cl_st2_list[j]].epsilon_j = 0.0;
                                cl[cl_st2_list[j]].fil_j = 0;
                                cl[cl_st2_list[j]].state = 1;
                                cl[cl_st2_list[j]].nlist_clean();
                                cl[cl_st2_list[j]].end_state_j = 0;

                            }
                            else if(bin_rand > 0.0)
                            {
                                cl[cl_st2_list[j]].epsilon_i = cl[cl_st2_list[j]].epsilon_j;
                                cl[cl_st2_list[j]].fil_i = cl[cl_st2_list[j]].fil_j;


                                cl[cl_st2_list[j]].epsilon_j = 0.0;
                                cl[cl_st2_list[j]].fil_j = 0;
                                cl[cl_st2_list[j]].state = 1; // NOTE: need to recalculate the neighbour list for this state-1 cl
                                cl[cl_st2_list[j]].nlist_clean();
                                cl[cl_st2_list[j]].end_state_i = 0;
                                cl[cl_st2_list[j]].end_state_j = 0;
                            }
                            
                            cl_check [0] = 0.0; cl_check [1] = 0.0; cl_check [2] = 0.0;
                            
                            for(n=0; n < N_DIMS; n++)
                            {
                                cl_check[n] = mt[cl[cl_st2_list[j]].fil_i].r_com[n] + cl[cl_st2_list[j]].epsilon_i*mt[cl[cl_st2_list[j]].fil_i].u_vec[n];
                            }
                            
                            // re-create nlist for the cl returned to state 1
                            int cl_inc_off = 0;
                            
                            for (m = 0; m < N_FILS; m++)
                            {
                                if(m != cl[cl_st2_list[j]].fil_i)
                                {
                                    rcl_image = binding.image_calc_BC_cl(BC_TYPE, N_DIMS, mt[m].r_com, mt[cl[cl_st2_list[j]].fil_i].r_com, box_vec);
                                    
                                    if(binding.nlist_cl(rcl_image,cl_check,mt[m].u_vec,mt[m].fil_length,N_DIMS) <= sig_cl)
                                    {
                                        cl[cl_st2_list[j]].neigh_append(m);
                                        cl_inc_off = cl_inc_off + 1;
                                    }
                                }
                            }
                            cl[cl_st2_list[j]].n_neigh = cl_inc_off;
                        }
                    }
                    cl_nst1 = binding.cl_state_count(cl, N_CLS, 1, cl_st1_list);
                    cl_nst2 = binding.cl_state_count(cl, N_CLS, 2, cl_st2_list);
                }
            }
        }

        // MT length dynamics
        for (j=N_PIN; j < N_FILS; j++)
        {
            mt[j].filament_dlength(dt,v_g,v_s,N_DIMS,upper_length);
        }
        
        // state-1 motor walking - kinesin type only
        for (j=0; j < cl_nst1; j++)
        {
            for (j=0; j < cl_nst1; j++)
            {
                if(cl[cl_st1_list[j]].type == 1)
                {
                    cl[cl_st1_list[j]].cl_walk_poly_1_type_k(mt[cl[cl_st1_list[j]].fil_i].fil_length, motor_vel1, dt, mt[cl[cl_st1_list[j]].fil_i].growth_phase, v_g, v_s);
                }
            }
        }
        
        // state-2 motor walking - kinesin type only
        for (j=0; j < cl_nst2; j++)
        {
            if(cl[cl_st2_list[j]].type == 1){cl[cl_st2_list[j]].cl_walk_poly_2_type_k(mt[cl[cl_st2_list[j]].fil_i].fil_length,mt[cl[cl_st2_list[j]].fil_j].fil_length, motor_vel1, dt, mt[cl[cl_st2_list[j]].fil_i].growth_phase, mt[cl[cl_st2_list[j]].fil_j].growth_phase, v_g, v_s);}
        }
        
        // unbinding for motors at the end of filaments
        // need to update properties and n-lists etc.
        for (j=0; j < cl_nst2; j++)
        {
            if(cl[cl_st2_list[j]].end_state_i == 1 && cl[cl_st2_list[j]].end_state_j == 1)
            {
                cout << "WARNING: Both CL2 in end state" << endl;
            }
            
            if(cl[cl_st2_list[j]].end_state_i == 1)
            {
                // only detach connection j, and shuffle properties if need
                cl[cl_st2_list[j]].state = 1;
                cl[cl_st2_list[j]].nlist_clean();
                cl[cl_st2_list[j]].epsilon_i = cl[cl_st2_list[j]].epsilon_j;
                cl[cl_st2_list[j]].fil_i = cl[cl_st2_list[j]].fil_j;
                cl[cl_st2_list[j]].epsilon_j = 0.0;
                cl[cl_st2_list[j]].fil_j = 0;
                cl[cl_st2_list[j]].end_state_i = 0;
                cl[cl_st2_list[j]].end_state_j = 0;

                cl_check [0] = 0.0; cl_check [1] = 0.0; cl_check [2] = 0.0;
                
                for(n=0; n<N_DIMS; n++){cl_check[n] = mt[cl[cl_st2_list[j]].fil_i].r_com[n] + cl[cl_st2_list[j]].epsilon_i*mt[cl[cl_st2_list[j]].fil_i].u_vec[n];}
                
                // re-create nlist for the cl returned to state 1
                int cl_inc_off = 0;
                
                for (m = 0; m < N_FILS; m++)
                {
                    if(m != cl[cl_st2_list[j]].fil_i)
                    {
                        rcl_image = binding.image_calc_BC_cl(BC_TYPE, N_DIMS, mt[m].r_com , mt[cl[cl_st2_list[j]].fil_i].r_com, box_vec);
                        
                        if(binding.nlist_cl(rcl_image,cl_check,mt[m].u_vec,mt[m].fil_length,N_DIMS) <= sig_cl)
                        {
                            cl[cl_st2_list[j]].neigh_append(m); cl_inc_off = cl_inc_off + 1;
                        }
                    }
                }
                cl[cl_st2_list[j]].n_neigh = cl_inc_off;
            }
            if( cl[cl_st2_list[j]].end_state_j == 1)
            {
                cl[cl_st2_list[j]].state = 1; // NOTE: need to recalculate the neighbour list for this state-1 cl
                cl[cl_st2_list[j]].nlist_clean();
                cl[cl_st2_list[j]].epsilon_j = 0.0;
                cl[cl_st2_list[j]].fil_j = 0;
                cl[cl_st2_list[j]].end_state_j = 0;

                cl_check [0] = 0.0; cl_check [1] = 0.0; cl_check [2] = 0.0;
                
                for(n=0; n<N_DIMS; n++){cl_check[n] = mt[cl[cl_st2_list[j]].fil_i].r_com[n] + cl[cl_st2_list[j]].epsilon_i*mt[cl[cl_st2_list[j]].fil_i].u_vec[n];}
                
                // re-create nlist for the cl returned to state 1
                int cl_inc_off = 0;
                
                for (m = 0; m < N_FILS; m++)
                {
                    if(m != cl[cl_st2_list[j]].fil_i)
                    {
                        rcl_image = binding.image_calc_BC_cl(BC_TYPE, N_DIMS, mt[m].r_com , mt[cl[cl_st2_list[j]].fil_i].r_com, box_vec);
                        
                        if(binding.nlist_cl(rcl_image,cl_check,mt[m].u_vec,mt[m].fil_length,N_DIMS) <= sig_cl)
                        {
                            cl[cl_st2_list[j]].neigh_append(m);
                            cl_inc_off = cl_inc_off + 1;
                        }
                    }
                }
                cl[cl_st2_list[j]].n_neigh = cl_inc_off;
            }
        }
         
        // perform state count 
        cl_nst1 = binding.cl_state_count(cl, N_CLS, 1, cl_st1_list);
        cl_nst2 = binding.cl_state_count(cl, N_CLS, 2, cl_st2_list);

        // Calculate drag coefficients for changing filament lengths
        for (j=N_PIN; j < N_FILS; j++){mt[j].filament_drag(d,kappa);}
        
        // oscillate pinning field
        if(i >= osc_start)
        {
            for (j=0; j < N_PIN; j++)
            {
                mt[j].r_com[1] = y_cryst + osc_amp*sin(osc_period*0.0017453*dt*osc_inc);
            }
            osc_inc = osc_inc + 1.0;
        }
        
        // include thermal noise contribution (Brownian-force)
        for (j = N_PIN; j < N_FILS; j++)
        {
            gram_schmidt_3d(mt[j].u_vec, mt[j].u_perp, e1_3d, e2_3d);
            stoch_3d_util(mt[j], kTt_2);
            
            u_mod = 0.0; for (unsigned int n = 0; n < N_DIMS; n++){u_mod = u_mod + mt[j].u_vec[n]*mt[j].u_vec[n];}
            u_mod = sqrt(u_mod); for (unsigned int n = 0; n < N_DIMS; n++){mt[j].u_vec[n] = mt[j].u_vec[n]/u_mod;}
        }
        
        // re-initialise the equations of motion
        for (j=N_PIN; j < N_FILS; j++){for(n=0; n<N_DIMS; n++){mt[j].r_eqnm[n] = 0.0;} mt[j].u_eqnm[0] = 0.0; mt[j].u_eqnm[1] = 0.0; mt[j].u_eqnm[2] = 0.0;}
        
        if(BC_TYPE == 0) // No PBC and no Walls - infinite space
        {
            if (i % n_list_update_outer == 0)
            {
                // Initialise the outer neighbour list
                pair.outer_nlist_init(N_DIMS);
                pair.inner_nlist_init(N_DIMS);
                for (unsigned int i = 0; i < N_FILS-1; i++)
                {
                    for(unsigned int j = i+1; j < N_FILS; j++)
                    {
                        pair.least_sep_calc_single(mt[i].r_com, mt[j].r_com, mt[i].u_vec ,mt[j].u_vec, mt[i].fil_length,mt[j].fil_length,N_DIMS);
                        if(pair.w_ij_mag_single <= outer_nlist_cut){pair.outer_nlist_build(i,j);}
                    }
                }
            }
            if (i % n_list_update == 0)
            {
                pair.inner_nlist_init(N_DIMS);
                
                for(unsigned int k = 0; k < pair.n_outer_pair; k++)
                {
                    
                    pair.least_sep_calc_single(mt[pair.outer_nlist[k][0]].r_com, mt[pair.outer_nlist[k][1]].r_com, mt[pair.outer_nlist[k][0]].u_vec ,mt[pair.outer_nlist[k][1]].u_vec, mt[pair.outer_nlist[k][0]].fil_length,mt[pair.outer_nlist[k][1]].fil_length,N_DIMS);
                    if(pair.w_ij_mag_single <= inner_nlist_cut){pair.inner_nlist_build(pair.outer_nlist[k][0],pair.outer_nlist[k][1]);}
                }
                pair.inner_pair_constructor(N_DIMS);
            }
            {
                
                for(unsigned int k = 0; k < pair.n_inner_pair; k++)
                {
                    pair.least_sep_calc(mt[pair.inner_nlist[k][0]].r_com,mt[pair.inner_nlist[k][1]].r_com, mt[pair.inner_nlist[k][0]].u_vec ,mt[pair.inner_nlist[k][1]].u_vec, mt[pair.inner_nlist[k][0]].fil_length,mt[pair.inner_nlist[k][1]].fil_length,N_DIMS,k);
                }
            }
        }
        if(BC_TYPE == 1) // PBC in (x, y, z) for 3D
        {
            if (i % n_list_update_outer == 0)
            {
                // Initialise the outer neighbour list
                pair.outer_nlist_init(N_DIMS);
                pair.inner_nlist_init(N_DIMS);
                for (unsigned int i = 0; i < N_FILS-1; i++)
                {
                    for(unsigned int j = i+1; j < N_FILS; j++)
                    {
                        for(n=0; n<N_DIMS; n++)
                        {
                            ri_image[n] = mt[i].r_com[n];
                            rj_image[n] = mt[j].r_com[n];
                            
                            if(mt[i].r_com[n] - mt[j].r_com[n] > box_vec[n]/2.0){ri_image[n] = mt[i].r_com[n]-box_vec[n];}
                            else if(mt[i].r_com[n] - mt[j].r_com[n] < -box_vec[n]/2.0){ri_image[n] = mt[i].r_com[n] + box_vec[n];}
                        }
                        pair.least_sep_calc_single(ri_image, rj_image, mt[i].u_vec ,mt[j].u_vec, mt[i].fil_length,mt[j].fil_length,N_DIMS);
                        if(pair.w_ij_mag_single <= outer_nlist_cut){pair.outer_nlist_build(i,j);}
                    }
                }
            }
            if (i % n_list_update == 0)
            {
                pair.inner_nlist_init(N_DIMS);
                for(unsigned int k = 0; k < pair.n_outer_pair; k++)
                {
                    for(n=0;n<N_DIMS;n++)
                    {
                        ri_image[n] = mt[pair.outer_nlist[k][0]].r_com[n]; rj_image[n] = mt[pair.outer_nlist[k][1]].r_com[n];
                    }
                    
                    for(n=0;n<N_DIMS;n++)
                    {
                        if(mt[pair.outer_nlist[k][0]].r_com[n] - mt[pair.outer_nlist[k][1]].r_com[n] > box_vec[n]/2.0){ri_image[n] = mt[pair.outer_nlist[k][0]].r_com[n]-box_vec[n];}
                        else if(mt[pair.outer_nlist[k][0]].r_com[n] - mt[pair.outer_nlist[k][1]].r_com[n] < -box_vec[n]/2.0){ri_image[n] = mt[pair.outer_nlist[k][0]].r_com[n] + box_vec[n];}
                    }
                    
                    pair.least_sep_calc_single(ri_image, rj_image, mt[pair.outer_nlist[k][0]].u_vec ,mt[pair.outer_nlist[k][1]].u_vec, mt[pair.outer_nlist[k][0]].fil_length,mt[pair.outer_nlist[k][1]].fil_length,N_DIMS);
                    
                    if(pair.w_ij_mag_single <= inner_nlist_cut){pair.inner_nlist_build(pair.outer_nlist[k][0],pair.outer_nlist[k][1]);}
                }
                
                // pair.inner_pair_constructor(N_DIMS);
                
            }
            {
                for(unsigned int k = 0; k < pair.n_inner_pair; k++)
                {
                    for(unsigned int n = 0; n < N_DIMS; n++)
                    {
                        ri_image[n] = mt[pair.inner_nlist[k][0]].r_com[n]; rj_image[n] = mt[pair.inner_nlist[k][1]].r_com[n];
                        
                        if(mt[pair.inner_nlist[k][0]].r_com[n] - mt[pair.inner_nlist[k][1]].r_com[n] > box_vec[n]/2.0){ri_image[n] = mt[pair.inner_nlist[k][0]].r_com[n] - box_vec[n];}
                        else if(mt[pair.inner_nlist[k][0]].r_com[n] - mt[pair.inner_nlist[k][1]].r_com[n] < - box_vec[n]/2.0){ri_image[n] = mt[pair.inner_nlist[k][0]].r_com[n] + box_vec[n];}
                    }
                    pair.least_sep_calc(ri_image,rj_image, mt[pair.inner_nlist[k][0]].u_vec ,mt[pair.inner_nlist[k][1]].u_vec, mt[pair.inner_nlist[k][0]].fil_length,mt[pair.inner_nlist[k][1]].fil_length,N_DIMS,k);
                }
            }
        }
        if(BC_TYPE == 2) // Walls in (x, z) with PBC in y for 3D
        {
            if (i % n_list_update_outer == 0)
            {
                // Initialise the outer neighbour list
                pair.outer_nlist_init(N_DIMS); pair.inner_nlist_init(N_DIMS);
                for (unsigned int k = 0; k < N_FILS-1; k++)
                {
                    for(unsigned int j = k+1; j < N_FILS; j++)
                    {
                        for(n=0; n<N_DIMS; n++)
                        {
                            ri_image[n] = mt[k].r_com[n]; rj_image[n] = mt[j].r_com[n];
                        }
                        
                        if(mt[k].r_com[1] - mt[j].r_com[1] > box_vec[1]/2.0){ri_image[1] = mt[k].r_com[1] - box_vec[1];}
                        else if(mt[k].r_com[1] - mt[j].r_com[1] < -box_vec[1]/2.0){ri_image[1] = mt[k].r_com[1] + box_vec[1];}
                        
                        pair.least_sep_calc_single(ri_image, rj_image, mt[k].u_vec ,mt[j].u_vec, mt[k].fil_length,mt[j].fil_length,N_DIMS);
                        if(pair.w_ij_mag_single <= outer_nlist_cut){pair.outer_nlist_build(k,j);}
                    }
                }
            }
            if (i % n_list_update == 0)
            {
                pair.inner_nlist_init(N_DIMS);
                
                for(unsigned int k = 0; k < pair.n_outer_pair; k++)
                {
                    for(n=0;n<N_DIMS;n++)
                    {
                        ri_image[n] = mt[pair.outer_nlist[k][0]].r_com[n]; rj_image[n] = mt[pair.outer_nlist[k][1]].r_com[n];
                    }
                    
                    if(mt[pair.outer_nlist[k][0]].r_com[1] - mt[pair.outer_nlist[k][1]].r_com[1] > box_vec[1]/2.0){ri_image[1] = mt[pair.outer_nlist[k][0]].r_com[1]-box_vec[1];}
                    else if(mt[pair.outer_nlist[k][0]].r_com[1] - mt[pair.outer_nlist[k][1]].r_com[1] < -box_vec[1]/2.0){ri_image[1] = mt[pair.outer_nlist[k][0]].r_com[1] + box_vec[1];}
                    
                    pair.least_sep_calc_single(ri_image, rj_image, mt[pair.outer_nlist[k][0]].u_vec ,mt[pair.outer_nlist[k][1]].u_vec, mt[pair.outer_nlist[k][0]].fil_length,mt[pair.outer_nlist[k][1]].fil_length,N_DIMS);
                    
                    if(pair.w_ij_mag_single <= inner_nlist_cut){pair.inner_nlist_build(pair.outer_nlist[k][0],pair.outer_nlist[k][1]);}
                }
            }
            {
                for(unsigned int k = 0; k < pair.n_inner_pair; k++)
                {
                    for(unsigned int n = 0; n < N_DIMS; n++)
                    {
                        ri_image[n] = mt[pair.inner_nlist[k][0]].r_com[n];rj_image[n] = mt[pair.inner_nlist[k][1]].r_com[n];
                    }
                    
                    if(mt[pair.inner_nlist[k][0]].r_com[1] - mt[pair.inner_nlist[k][1]].r_com[1] > box_vec[1]/2.0){ri_image[1] = mt[pair.inner_nlist[k][0]].r_com[1] - box_vec[1];}
                    else if(mt[pair.inner_nlist[k][0]].r_com[1] - mt[pair.inner_nlist[k][1]].r_com[1] < - box_vec[1]/2.0){ri_image[1] = mt[pair.inner_nlist[k][0]].r_com[1] + box_vec[1];}
                    
                    pair.least_sep_calc(ri_image,rj_image, mt[pair.inner_nlist[k][0]].u_vec ,mt[pair.inner_nlist[k][1]].u_vec, mt[pair.inner_nlist[k][0]].fil_length,mt[pair.inner_nlist[k][1]].fil_length,N_DIMS,k);
                }
            }
            
        }
        if(BC_TYPE == 3) // Walls in (x, z) open BC in y for 3D
        {
            if (i % n_list_update_outer == 0)
            {
                // Initialise the outer neighbour list
                pair.outer_nlist_init(N_DIMS);
                
                for (unsigned int k = 0; k < N_FILS-1; k++)
                {
                    for(unsigned int j = k+1; j < N_FILS; j++)
                    {
                        
                        pair.least_sep_calc_single(mt[k].r_com, mt[j].r_com, mt[k].u_vec ,mt[j].u_vec, mt[k].fil_length,mt[j].fil_length,N_DIMS);
                        
                        if(k < N_PIN && j < N_PIN){pair.w_ij_mag_single = 1.0;}
                        if(pair.w_ij_mag_single <= outer_nlist_cut){pair.outer_nlist_build(k,j);}
                    }
                }
            }
            if (i % n_list_update == 0)
            {
                pair.inner_nlist_init(N_DIMS);
                
                for(unsigned int k = 0; k < pair.n_outer_pair; k++)
                {
                    pair.least_sep_calc_single(mt[pair.outer_nlist[k][0]].r_com, mt[pair.outer_nlist[k][1]].r_com, mt[pair.outer_nlist[k][0]].u_vec ,mt[pair.outer_nlist[k][1]].u_vec, mt[pair.outer_nlist[k][0]].fil_length,mt[pair.outer_nlist[k][1]].fil_length,N_DIMS);
                    
                    if(pair.w_ij_mag_single <= inner_nlist_cut){pair.inner_nlist_build(pair.outer_nlist[k][0],pair.outer_nlist[k][1]);}
                }
                pair.inner_pair_constructor(N_DIMS);
            }
            {
                for(unsigned int k = 0; k < pair.n_inner_pair; k++)
                {
                    pair.least_sep_calc(mt[pair.inner_nlist[k][0]].r_com,mt[pair.inner_nlist[k][1]].r_com, mt[pair.inner_nlist[k][0]].u_vec ,mt[pair.inner_nlist[k][1]].u_vec, mt[pair.inner_nlist[k][0]].fil_length,mt[pair.inner_nlist[k][1]].fil_length,N_DIMS,k);
                }
            }
        }
        
        //************************ Start force section of main loop **************************
        
        // Calculate force from "WCA potential" along all lines of least of interaction for the inner pair list
        for(unsigned int k = 0; k < pair.n_inner_pair; k++)
        {
            F_ij = 0.0;
            
            if(pair.w_ij_mag[k] <= sig_shift_min)
            {
                if(pair.w_ij_mag[k] > beta_sig)
                {
                    F_ij = -24.0*(epsilon/sig_shift)*(2.0*pow(sig_shift/pair.w_ij_mag[k],13.0) - pow(sig_shift/pair.w_ij_mag[k],7.0));
                }
                else if(pair.w_ij_mag[k] <= beta_sig)
                {
                    F_ij = f_beta + df_beta*(pair.w_ij_mag[k] - beta_sig);
                }
                if(pair.w_ij_mag[k] < 1.0e-13)
                {
                    std :: cout << "WARNING - Filaments intercepting: " << i << "    " << pair.w_ij_mag[k] << "\t" << pair.inner_nlist[k][0] << "\t" << pair.inner_nlist[k][1] << endl;
                    std :: cout << "Filament lengths:     " << mt[pair.inner_nlist[k][0]].fil_length << "\t" << mt[pair.inner_nlist[k][1]].fil_length << endl;
                    
                    F_ij = 0;
                }
            }
            pair.F_ij_vec[k] = F_ij;
        }

        // calculate steric interaction force for Eqns of Mot. needed for both position and orientation equations (steric-force)
        for(unsigned int k = 0; k < pair.n_inner_pair; k++)
        {
            for(n=0;n<N_DIMS;n++)
            {
                mt[pair.inner_nlist[k][0]].r_eqnm[n] = mt[pair.inner_nlist[k][0]].r_eqnm[n] - pair.F_ij_vec[k]*pair.uw_ij[k][n];
                mt[pair.inner_nlist[k][1]].r_eqnm[n] = mt[pair.inner_nlist[k][1]].r_eqnm[n] + pair.F_ij_vec[k]*pair.uw_ij[k][n];
            }
            
            u_iw_dot = mt[pair.inner_nlist[k][0]].u_vec[0]*pair.uw_ij[k][0] + mt[pair.inner_nlist[k][0]].u_vec[1]*pair.uw_ij[k][1] + mt[pair.inner_nlist[k][0]].u_vec[2]*pair.uw_ij[k][2];
            u_jw_dot = mt[pair.inner_nlist[k][1]].u_vec[0]*pair.uw_ij[k][0] + mt[pair.inner_nlist[k][1]].u_vec[1]*pair.uw_ij[k][1] + mt[pair.inner_nlist[k][1]].u_vec[2]*pair.uw_ij[k][2];
            
            for(n=0;n<N_DIMS;n++)
            {
                mt[pair.inner_nlist[k][0]].u_eqnm[n] = mt[pair.inner_nlist[k][0]].u_eqnm[n] + pair.F_ij_vec[k]*pair.lambda_array[k][0]*(u_iw_dot*mt[pair.inner_nlist[k][0]].u_vec[n] - pair.uw_ij[k][n]);
                mt[pair.inner_nlist[k][1]].u_eqnm[n] = mt[pair.inner_nlist[k][1]].u_eqnm[n] - pair.F_ij_vec[k]*pair.lambda_array[k][1]*(u_jw_dot*mt[pair.inner_nlist[k][1]].u_vec[n] - pair.uw_ij[k][n]);
            }
        }
        
        // Calculate the wall boundary conditions (wall-force)
        // herepbc_wall
        if(BC_TYPE == 2 || BC_TYPE == 3)
        {
            // Create a wall pair list for neighboring filaments (3D)
            if (i % wall_update == 0)
            {
                wall.wall_init();
                
                for (j=N_PIN; j < N_FILS; j++)
                {
                    
                    ri_s3d[0]  = 0.0;
                    ri_s3d[1]  = mt[j].r_com[1];
                    ri_s3d[2]  = mt[j].r_com[2];
                    
                    // must normalise:
                    uwall_3d[0] = 0.0;
                    uwall_3d[1] = mt[j].u_vec[1];
                    uwall_3d[2] = mt[j].u_vec[2];

                    umod = sqrt(uwall_3d[1]*uwall_3d[1] + uwall_3d[2]*uwall_3d[2]);
                    
                    uwall_3d[1] = uwall_3d[1]/umod;
                    uwall_3d[2] = uwall_3d[2]/umod;
                    
                    wall_calc.wall_calc_class(mt[j].r_com, ri_s3d, mt[j].u_vec, uwall_3d, mt[j].fil_length, L_wall_3D, N_DIMS,j);
                    
                    if(wall_calc.w_ij_mag_out <= wall_cut*sig_min_wall || wall_calc.w_ij_mag_out == 0.0 || isnan(wall_calc.w_ij_mag_out) == 1){wall.wall_x1_add(j);}
                    
                    ri_s3d[0]  = L_x;
                    
                    wall_calc.wall_calc_class(mt[j].r_com, ri_s3d, mt[j].u_vec,uwall_3d,mt[j].fil_length,L_wall_3D,N_DIMS,j);
                    
                    if(wall_calc.w_ij_mag_out <= wall_cut*sig_min_wall || wall_calc.w_ij_mag_out == 0.0 || isnan(wall_calc.w_ij_mag_out) == 1){wall.wall_x2_add(j);}
                    
                    ri_s3d[0]  = mt[j].r_com[0];
                    ri_s3d[1]  = mt[j].r_com[1];
                    ri_s3d[2]  = 0.0;
                    
                    // must normalise:
                    uwall_3d[0] = mt[j].u_vec[0];
                    uwall_3d[1] = mt[j].u_vec[1];
                    uwall_3d[2] = 0.0;

                    umod = sqrt(uwall_3d[0]*uwall_3d[0] + uwall_3d[1]*uwall_3d[1]);
                    
                    uwall_3d[0] = uwall_3d[0]/umod;
                    uwall_3d[1] = uwall_3d[1]/umod;
                    
                    wall_calc.wall_calc_class(mt[j].r_com, ri_s3d, mt[j].u_vec,uwall_3d,mt[j].fil_length,L_wall_3D,N_DIMS,j);
                    
                    if(wall_calc.w_ij_mag_out <= wall_cut*sig_min_wall || wall_calc.w_ij_mag_out == 0.0 || isnan(wall_calc.w_ij_mag_out) == 1){wall.wall_z1_add(j);}
                    
                    ri_s3d[2]  = L_z;
                    
                    wall_calc.wall_calc_class(mt[j].r_com, ri_s3d, mt[j].u_vec,uwall_3d,mt[j].fil_length,L_wall_3D,N_DIMS,j);
                    
                    if(wall_calc.w_ij_mag_out <= wall_cut*sig_min_wall || wall_calc.w_ij_mag_out == 0.0 || isnan(wall_calc.w_ij_mag_out) == 1){wall.wall_z2_add(j);}
                }
            }
            
            // Calculate wall force for each listed wall-filament interaction (do for 2 walls in xy plane and 2 in yz plane)

            for (j=0; j < wall.x1_n; j++)
            {
                // wall_here
                
                ri_s3d[0]  = 0.0;
                ri_s3d[1]  = mt[wall.wl_x1[j]].r_com[1];
                ri_s3d[2]  = mt[wall.wl_x1[j]].r_com[2];
                
                // must normalise:
                uwall_3d[0] = 0.0;
                uwall_3d[1] = mt[wall.wl_x1[j]].u_vec[1];
                uwall_3d[2] = mt[wall.wl_x1[j]].u_vec[2];

                umod = sqrt(uwall_3d[1]*uwall_3d[1] + uwall_3d[2]*uwall_3d[2]);
                
                uwall_3d[1] = uwall_3d[1]/umod;
                uwall_3d[2] = uwall_3d[2]/umod;
                
                F_i_wall    = 0.0;
                
                wall_calc.wall_calc_class(mt[wall.wl_x1[j]].r_com, ri_s3d, mt[wall.wl_x1[j]].u_vec, uwall_3d, mt[wall.wl_x1[j]].fil_length, L_wall_3D, N_DIMS,j);
                
                if(wall_calc.w_ij_mag_out <= sig_min_wall)
                {
                    
                    F_i_wall = wall_calc.wall_wca_calc();
                    
                    for(n=0;n<N_DIMS;n++)
                    {
                        mt[wall.wl_x1[j]].r_eqnm[n] = mt[wall.wl_x1[j]].r_eqnm[n] - F_i_wall*uwall1[n];
                    }
                    
                    u_iw_dot = mt[wall.wl_x1[j]].u_vec[0]*uwall1[0] + mt[wall.wl_x1[j]].u_vec[1]*uwall1[1] + mt[wall.wl_x1[j]].u_vec[2]*uwall1[2];
                    
                    for(n=0;n<N_DIMS;n++)
                    {
                        mt[wall.wl_x1[j]].u_eqnm[n] = mt[wall.wl_x1[j]].u_eqnm[n] + F_i_wall*wall_calc.lambda[0]*(u_iw_dot*mt[wall.wl_x1[j]].u_vec[n] - uwall1[n]);
                    }
                }
            }

            for (j=0; j < wall.x2_n; j++)
            {
                
                ri_s3d[0]  = L_x;
                ri_s3d[1]  = mt[wall.wl_x2[j]].r_com[1];
                ri_s3d[2]  = mt[wall.wl_x2[j]].r_com[2];
                
                // must normalise:
                uwall_3d[0] = 0.0;
                uwall_3d[1] = mt[wall.wl_x2[j]].u_vec[1];
                uwall_3d[2] = mt[wall.wl_x2[j]].u_vec[2];
                
                umod = sqrt(umod = uwall_3d[1]*uwall_3d[1] + uwall_3d[2]*uwall_3d[2]);
                
                uwall_3d[1] = uwall_3d[1]/umod;
                uwall_3d[2] = uwall_3d[2]/umod;
                
                F_i_wall    = 0.0;
                
                wall_calc.wall_calc_class(mt[wall.wl_x2[j]].r_com, ri_s3d, mt[wall.wl_x2[j]].u_vec,uwall_3d,mt[wall.wl_x2[j]].fil_length,L_wall_3D,N_DIMS,j);
                
                
                if(wall_calc.w_ij_mag_out <= sig_min_wall)
                {
                    
                    F_i_wall = wall_calc.wall_wca_calc();
                    
                    for(n=0;n<N_DIMS;n++)
                    {
                        mt[wall.wl_x2[j]].r_eqnm[n] = mt[wall.wl_x2[j]].r_eqnm[n] - F_i_wall*uwall2[n];
                    }
                    
                    u_iw_dot = mt[wall.wl_x2[j]].u_vec[0]*uwall2[0] + mt[wall.wl_x2[j]].u_vec[1]*uwall2[1] + mt[wall.wl_x2[j]].u_vec[2]*uwall2[2];
                    
                    for(n=0;n<N_DIMS;n++)
                    {
                        mt[wall.wl_x2[j]].u_eqnm[n] = mt[wall.wl_x2[j]].u_eqnm[n] + F_i_wall*wall_calc.lambda[0]*(u_iw_dot*mt[wall.wl_x2[j]].u_vec[n] - uwall2[n]);
                    }
                }
            }

            for (j=0; j < wall.z1_n; j++)
            {
                ri_s3d[0]  = mt[wall.wl_z1[j]].r_com[0];
                ri_s3d[1]  = mt[wall.wl_z1[j]].r_com[1];
                ri_s3d[2]  = 0.0;
                
                // must normalise:
                uwall_3d[0] = mt[wall.wl_z1[j]].u_vec[0];
                uwall_3d[1] = mt[wall.wl_z1[j]].u_vec[1];
                uwall_3d[2] = 0.0;
                
                umod = sqrt(uwall_3d[0]*uwall_3d[0] + uwall_3d[1]*uwall_3d[1]);
                
                uwall_3d[0] = uwall_3d[0]/umod;
                uwall_3d[1] = uwall_3d[1]/umod;
                
                F_i_wall    = 0.0;
                
                wall_calc.wall_calc_class(mt[wall.wl_z1[j]].r_com, ri_s3d, mt[wall.wl_z1[j]].u_vec,uwall_3d,mt[wall.wl_z1[j]].fil_length,L_wall_3D,N_DIMS,j);
                
                if(wall_calc.w_ij_mag_out <= sig_min_wall)
                {

                    F_i_wall = wall_calc.wall_wca_calc();
                    
                    for(n=0;n<N_DIMS;n++)
                    {
                        mt[wall.wl_z1[j]].r_eqnm[n] = mt[wall.wl_z1[j]].r_eqnm[n] - F_i_wall*uwall3[n];
                    }
                    
                    u_iw_dot = mt[wall.wl_z1[j]].u_vec[0]*uwall3[0] + mt[wall.wl_z1[j]].u_vec[1]*uwall3[1] + mt[wall.wl_z1[j]].u_vec[2]*uwall3[2];
                    
                    for(n=0;n<N_DIMS;n++)
                    {
                        mt[wall.wl_z1[j]].u_eqnm[n] = mt[wall.wl_z1[j]].u_eqnm[n] + F_i_wall*wall_calc.lambda[0]*(u_iw_dot*mt[wall.wl_z1[j]].u_vec[n] - uwall3[n]);
                    }
                }
            }

            for (j=0; j < wall.z2_n; j++)
            {
                
                ri_s3d[0]  = mt[wall.wl_z2[j]].r_com[0];
                ri_s3d[1]  = mt[wall.wl_z2[j]].r_com[1];
                ri_s3d[2]  = L_z;
                
                // must normalise:
                uwall_3d[0] = mt[wall.wl_z2[j]].u_vec[0];
                uwall_3d[1] = mt[wall.wl_z2[j]].u_vec[1];
                uwall_3d[2] = 0.0;
                
                umod = sqrt(uwall_3d[0]*uwall_3d[0] + uwall_3d[1]*uwall_3d[1]);
                
                uwall_3d[0] = uwall_3d[0]/umod;
                uwall_3d[1] = uwall_3d[1]/umod;
                
                F_i_wall    = 0.0;
                
                wall_calc.wall_calc_class(mt[wall.wl_z2[j]].r_com, ri_s3d, mt[wall.wl_z2[j]].u_vec,uwall_3d,mt[wall.wl_z2[j]].fil_length,L_wall_3D,N_DIMS,j);
                
                if(wall_calc.w_ij_mag_out <= sig_min_wall)
                {
                    F_i_wall = wall_calc.wall_wca_calc();
                    
                    for(n=0;n<N_DIMS;n++)
                    {
                        mt[wall.wl_z2[j]].r_eqnm[n] = mt[wall.wl_z2[j]].r_eqnm[n] - F_i_wall*uwall4[n];
                    }
                    
                    u_iw_dot = mt[wall.wl_z2[j]].u_vec[0]*uwall4[0] + mt[wall.wl_z2[j]].u_vec[1]*uwall4[1] + mt[wall.wl_z2[j]].u_vec[2]*uwall4[2];
                    
                    for(n=0;n<N_DIMS;n++)
                    {
                        mt[wall.wl_z2[j]].u_eqnm[n] = mt[wall.wl_z2[j]].u_eqnm[n] + F_i_wall*wall_calc.lambda[0]*(u_iw_dot*mt[wall.wl_z2[j]].u_vec[n] - uwall4[n]);
                    }
                }
            }
        }
        
        // Add forces due to cross-linker/motor interactions to Eqns of Motion (motor-force)
        if(CL == 1 || CL == 2)
        {
            for(unsigned int k = 0; k < cl_nst2; k++)
            {
                for(n=0;n<N_DIMS;n++)
                {
                    ri_image[n] = mt[cl[cl_st2_list[k]].fil_i].r_com[n];
                    rj_image[n] = mt[cl[cl_st2_list[k]].fil_j].r_com[n];
                }
                
                ri_image = pair.image_calc_BC_pair(BC_TYPE, N_DIMS, ri_image, rj_image, box_vec);
                
                d_ij_mag = 0.0;
                
                for(n=0; n<N_DIMS; n++)
                {
                    d_ij[n] = (ri_image[n] + cl[cl_st2_list[k]].epsilon_i*mt[cl[cl_st2_list[k]].fil_i].u_vec[n]) - (rj_image[n] + cl[cl_st2_list[k]].epsilon_j*mt[cl[cl_st2_list[k]].fil_j].u_vec[n]);
                    d_ij_mag = d_ij_mag + d_ij[n]*d_ij[n];
                }
                
                d_ij_mag = sqrt(d_ij_mag);
                
                for(n=0; n<N_DIMS; n++)
                {
                    d_0_vec[n] = (d_0/d_ij_mag)*d_ij[n];
                }
                
                // update the r_vec equations of motion
                for(n=0; n<N_DIMS; n++)
                {
                    mt[cl[cl_st2_list[k]].fil_i].r_eqnm[n] = mt[cl[cl_st2_list[k]].fil_i].r_eqnm[n] - cl[cl_st2_list[k]].k_spr*(d_ij[n]-d_0_vec[n]);
                    mt[cl[cl_st2_list[k]].fil_j].r_eqnm[n] = mt[cl[cl_st2_list[k]].fil_j].r_eqnm[n] + cl[cl_st2_list[k]].k_spr*(d_ij[n]-d_0_vec[n]);
                }
        
                f_pi = 0.0;
                f_pj = 0.0;
                
                if(cl[cl_st2_list[k]].type == 1 && d_ij_mag > d_0) // only calculate under extension
                {
                    
                    double uidijk_dot = mt[cl[cl_st2_list[k]].fil_i].u_vec[0]*d_ij[0] + mt[cl[cl_st2_list[k]].fil_i].u_vec[1]*d_ij[1] + mt[cl[cl_st2_list[k]].fil_i].u_vec[2]*d_ij[2];
                    double ujdjik_dot = mt[cl[cl_st2_list[k]].fil_j].u_vec[0]*(-d_ij[0]) + mt[cl[cl_st2_list[k]].fil_j].u_vec[1]*(-d_ij[1]) + mt[cl[cl_st2_list[k]].fil_j].u_vec[2]*(-d_ij[2]);
                    
                    if(uidijk_dot > 0.0)
                    {
                        for(n=0; n<N_DIMS; n++)
                        {
                            f_pi = f_pi - cl[cl_st2_list[k]].k_spr*(d_ij[n]-d_0_vec[n])*mt[cl[cl_st2_list[k]].fil_i].u_vec[n];
                        }
                    }
                    if(ujdjik_dot > 0.0)
                    {
                        for(n=0; n<N_DIMS; n++)
                        {
                            f_pj = f_pj + cl[cl_st2_list[k]].k_spr*(d_ij[n]-d_0_vec[n])*mt[cl[cl_st2_list[k]].fil_j].u_vec[n];
                        }
                    }
                }
                if(cl[cl_st2_list[k]].type == 2 && d_ij_mag > d_0) // only calculate under extension
                {
                    
                    double uidijd_dot = mt[cl[cl_st2_list[k]].fil_i].u_vec[0]*d_ij[0] + mt[cl[cl_st2_list[k]].fil_i].u_vec[1]*d_ij[1] + mt[cl[cl_st2_list[k]].fil_i].u_vec[2]*d_ij[2];
                    double ujdjid_dot = mt[cl[cl_st2_list[k]].fil_j].u_vec[0]*(-d_ij[0]) + mt[cl[cl_st2_list[k]].fil_j].u_vec[1]*(-d_ij[1]) + mt[cl[cl_st2_list[k]].fil_j].u_vec[2]*(-d_ij[2]);
                    
                    if(uidijd_dot < 0.0)
                    {
                        for(n=0; n<N_DIMS; n++)
                        {
                            f_pi = f_pi - cl[cl_st2_list[k]].k_spr*(d_ij[n]-d_0_vec[n])*mt[cl[cl_st2_list[k]].fil_i].u_vec[n];
                        }
                    }
                    if(ujdjid_dot < 0.0)
                    {
                        for(n=0; n<N_DIMS; n++)
                        {
                            f_pj = f_pj + cl[cl_st2_list[k]].k_spr*(d_ij[n]-d_0_vec[n])*mt[cl[cl_st2_list[k]].fil_j].u_vec[n];
                        }
                    }
                }
  
                u_id_dot = mt[cl[cl_st2_list[k]].fil_i].u_vec[0]*(d_ij[0]-d_0_vec[0]) + mt[cl[cl_st2_list[k]].fil_i].u_vec[1]*(d_ij[1]-d_0_vec[1]) + mt[cl[cl_st2_list[k]].fil_i].u_vec[2]*(d_ij[2]-d_0_vec[2]);
                u_jd_dot = mt[cl[cl_st2_list[k]].fil_j].u_vec[0]*(d_ij[0]-d_0_vec[0]) + mt[cl[cl_st2_list[k]].fil_j].u_vec[1]*(d_ij[1]-d_0_vec[1]) + mt[cl[cl_st2_list[k]].fil_j].u_vec[2]*(d_ij[2]-d_0_vec[2]);
                
                // update the u_vec equations of motion
                for(n=0;n<N_DIMS;n++)
                {
                    mt[cl[cl_st2_list[k]].fil_i].u_eqnm[n] = mt[cl[cl_st2_list[k]].fil_i].u_eqnm[n] + cl[cl_st2_list[k]].k_spr*cl[cl_st2_list[k]].epsilon_i*(u_id_dot*mt[cl[cl_st2_list[k]].fil_i].u_vec[n] - (d_ij[n]-d_0_vec[n]));
                    mt[cl[cl_st2_list[k]].fil_j].u_eqnm[n] = mt[cl[cl_st2_list[k]].fil_j].u_eqnm[n] - cl[cl_st2_list[k]].k_spr*cl[cl_st2_list[k]].epsilon_j*(u_jd_dot*mt[cl[cl_st2_list[k]].fil_j].u_vec[n] - (d_ij[n]-d_0_vec[n]));
                }

                cl[cl_st2_list[k]].f_para_i = f_pi;
                cl[cl_st2_list[k]].f_para_j = f_pj;
            }
        }


        mt[N_PIN].r_eqnm[0] = mt[N_PIN].r_eqnm[0] + k_trap*(L_x/2.0 - mt[N_PIN].r_com[0]);
        mt[N_PIN].r_eqnm[1] = mt[N_PIN].r_eqnm[1] + k_trap*(trap_eqm_y - mt[N_PIN].r_com[1]);
        mt[N_PIN].r_eqnm[2] = mt[N_PIN].r_eqnm[2] + k_trap*(L_z/2.0 - mt[N_PIN].r_com[2]);
        
        
        // numerically update the Eqns of Motion (2D or 3D).
        for (j = N_PIN; j < N_FILS; j++)
        {
            filament_drag_3d_util(mt[j].u_vec,mt[j].r_eqnm,mt[j].u_eqnm,gamma_array_3d,mt[j].gamma_para,mt[j].gamma_perp,mt[j].gamma_rot);
            
            mt[j].r_com[0] = mt[j].r_com[0] + dt*mt[j].r_eqnm[0];
            mt[j].r_com[1] = mt[j].r_com[1] + dt*mt[j].r_eqnm[1];
            mt[j].r_com[2] = mt[j].r_com[2] + dt*mt[j].r_eqnm[2];
            
            mt[j].u_vec[0] = mt[j].u_vec[0] + dt*mt[j].u_eqnm[0];
            mt[j].u_vec[1] = mt[j].u_vec[1] + dt*mt[j].u_eqnm[1];
            mt[j].u_vec[2] = mt[j].u_vec[2] + dt*mt[j].u_eqnm[2];
            
            u_mod = 0.0; for (unsigned int n = 0; n < N_DIMS; n++){u_mod = u_mod + mt[j].u_vec[n]*mt[j].u_vec[n];}
            u_mod = sqrt(u_mod); for (unsigned int n = 0; n < N_DIMS; n++){mt[j].u_vec[n] = mt[j].u_vec[n]/u_mod;}
        }
        
        // correct for PBCs if required (minimal image calculations)
        if(BC_TYPE == 1 )
        {
            for (j=0; j < N_FILS; j++)
            {
                pbc.pbc_calc(mt[j].r_com,N_DIMS,mt[j].pbc_index);
            }
        }
        if(BC_TYPE == 2)
        {
            for (j=0; j < N_FILS; j++)
            {
                pbc.pbc_y_calc(mt[j].r_com,N_DIMS,mt[j].pbc_index);
            }
        }
        
        // ****************** output file writing section **********************
        
        // files to be written from begin of simulation
        if (i % t_write == 0)
        {
            // just print position of pin c.o.m and trap filament
            trajectory_trap << double(i)*dt << "\t";
            for (j = N_PIN - 1; j <= N_PIN; j++)
            {
                trajectory_trap << mt[j].r_com[0] << "\t" << mt[j].r_com[1] << "\t" << mt[j].r_com[2] << "\t";
            }
            trajectory_trap << "\n";
            
            // main trajectory file for filaments
            trajectory << i << "\t" << N_FILS << "\t";
            for (j=0; j < N_FILS; j++)
            {
                trajectory << mt[j].fil_length << "\t"  << mt[j].r_com[0] << "\t" << mt[j].r_com[1] << "\t" << mt[j].r_com[2] << "\t" << mt[j].u_vec[0] << "\t" << mt[j].u_vec[1] << "\t" << mt[j].u_vec[2]<< "\t";
            }
            trajectory << "\n";
 
            cl_trajectory << i << "\t" << cl_nst1 << "\t" << cl_nst2 << "\t";

            for (j=0; j < N_CLS; j++)
            {
                if(cl[j].state == 1)
                {
                    cl_trajectory << cl[j].state << "\t" << cl[j].type << "\t" << mt[cl[j].fil_i].r_com[0] + cl[j].epsilon_i*mt[cl[j].fil_i].u_vec[0] << "\t" << mt[cl[j].fil_i].r_com[1] + cl[j].epsilon_i*mt[cl[j].fil_i].u_vec[1] << "\t" << mt[cl[j].fil_i].r_com[2] + cl[j].epsilon_i*mt[cl[j].fil_i].u_vec[2] << "\t";
                }
                if(cl[j].state == 2)
                {
                    cl_trajectory << cl[j].state << "\t" << cl[j].type << "\t" << mt[cl[j].fil_i].r_com[0] + cl[j].epsilon_i*mt[cl[j].fil_i].u_vec[0] << "\t" << mt[cl[j].fil_i].r_com[1] + cl[j].epsilon_i*mt[cl[j].fil_i].u_vec[1] << "\t" << mt[cl[j].fil_i].r_com[2] + cl[j].epsilon_i*mt[cl[j].fil_i].u_vec[2] << "\t" << mt[cl[j].fil_j].r_com[0] + cl[j].epsilon_j*mt[cl[j].fil_j].u_vec[0] << "\t" << mt[cl[j].fil_j].r_com[1] + cl[j].epsilon_j*mt[cl[j].fil_j].u_vec[1] << "\t" << mt[cl[j].fil_j].r_com[2] + cl[j].epsilon_j*mt[cl[j].fil_j].u_vec[2] << "\t";
                }
            }
            cl_trajectory << endl;
            
        }
        
        // files to be written with a delay (start at start_write)
        if(i >= start_write)
        {
            if (i % t_write_2 == 0)
            {
                // write all filaments with IDs - used for many things
                id_trajectory << double(i)*dt << "\t" << N_FILS - N_PIN << "\t" << filament_id_count << "\t";
                id_orientation << double(i)*dt << "\t" << N_FILS - N_PIN << "\t" << filament_id_count << "\t";
                if(N_DIMS == 3)
                {
                    for (j = N_PIN; j < N_FILS; j++)
                    {
                        id_trajectory << mt[j].fil_id << "\t" << mt[j].r_com[0] - 0.5*mt[j].fil_length*mt[j].u_vec[0] << "\t" << mt[j].r_com[1] - 0.5*mt[j].fil_length*mt[j].u_vec[1] << "\t" << mt[j].r_com[2] - 0.5*mt[j].fil_length*mt[j].u_vec[2] << "\t";
                        id_orientation << mt[j].fil_id << "\t" << mt[j].fil_length << "\t" << mt[j].u_vec[0] << "\t" << mt[j].u_vec[1] << "\t" << mt[j].u_vec[2] << "\t";
                    }
                }
                id_trajectory  << "\n";
                id_orientation << "\n";
                
                // Write out the filament IDs for each pair-bound CL - used for connectivity plots
                cl_id << i << "\t" << N_FILS - N_PIN << "\t";
                for(unsigned int k = 0; k < cl_nst2; k++)
                {
                    if(cl[cl_st2_list[k]].fil_i >= N_PIN && cl[cl_st2_list[k]].fil_j >= N_PIN)
                    {
                        cl_id << cl[cl_st2_list[k]].fil_i << "\t" << cl[cl_st2_list[k]].fil_j << "\t";
                    }
                }
                cl_id << endl;
                
                // write the CL states out for stress calculation rather then for visualisation
                // only uses the 2-state cls
            }
        }
        
        // print to standard output for monitoring the progress of simulations
        if (i % t_print == 0)
        {
            nt_cl1_1 = 0; nt_cl1_2 = 0;
            
            for (j=0; j < N_CLS; j++)
            {
                if(cl[j].state == 2)
                {
                    if(cl[j].type == 1){nt_cl1_1 = nt_cl1_1 + 1;}
                    if(cl[j].type == 2){nt_cl1_2 = nt_cl1_2 + 1;}
                }
            }
            
            clock_t toc_full = clock();
            
            double tic_toc_full = (double)(toc_full-tic_full);
            cout << endl;
            printf(  "\n%s%i%s%i %s%12.4g%s\n", "time step:    " , i,",    N_FILS:    ", N_FILS,",    Runtime: " ,tic_toc_full/CLOCKS_PER_SEC ,"sec");
            std::cout << "cl bound: " << cl_nst1 + cl_nst2 << ",    state 1: " << cl_nst1 << ",    state 2 (type 1): " << nt_cl1_1 << ",    total length: " << tot_fil_length << endl;
            cout << "inner pairs: " << pair.n_inner_pair << " outer pairs: " << pair.n_outer_pair << endl;
            cout << endl;
        }
    }
    // end of main loop
    
    // write configuration and close files
    std::cout << "Number of filaments:   " << N_FILS << endl;
    std::cout << "Number of CLs:   " << N_CLS << "    " << cl_nst1 << "    " << cl_nst2 << endl;
    
    
        // ##############################################################################################################################

        trajectory.close();
        cl_trajectory.close();
        id_trajectory.close();
        id_orientation.close();
        cl_id.close();
        cl_nt.close();
   
    clock_t toc_full = clock();
        
    double tic_toc_full = (double)(toc_full-tic_full);
    cout << "\nTime taken: " << tic_toc_full/CLOCKS_PER_SEC << " (s)" << " , " << ((tic_toc_full/CLOCKS_PER_SEC))/3600.0 << " (hr)"<< endl << endl;

    return 0;
}
