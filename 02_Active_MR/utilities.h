// Utilities library for useful calculations
// 
// Benjamin Dalton 18/11/2015

#ifndef utilities_H
#define utilities_H

#define PI 3.141592653589793
#define PREC2 5e-6

#include "filament.h"

using namespace std;

float myErfInv2(float x);

void gram_schmidt_2d(vector<double>& u_vec, vector<vector<double>>& u_perp_2d, double e1[2]);
void gram_schmidt_3d(vector<double>& u_vec, vector<vector<double>>& u_perp_3d, double e1[3], double e2[3]);

void gram_schmidt_2d_2(vector<double>& u_vec, vector<double>& u_perp_2d, double e1[2]);
void gram_schmidt_3d_2(vector<double>& u_vec, vector<double>& u_perp_3d_1, vector<double>& u_perp_3d_2, double e1[3], double e2[3]);

void filament_drag_2d_util(vector<double>& u_vec, vector<double>& r_eqnm, vector<double>& u_eqnm, double gamma_array[2][2], double gamma_para, double gamma_perp, double gamma_rot);
void filament_drag_3d_util(vector<double>& u_vec, vector<double>& r_eqnm, vector<double>& u_eqnm, double gamma_array[3][3], double gamma_para, double gamma_perp, double gamma_rot);

void stoch_2d_util(filament& mt, const double kTt_2);
void stoch_3d_util(filament& mt, const double kTt_2);

void fil_dl_util(filament& mt, const double dt, const double v_g, const double v_s, int n_dims, double ly);

double uniform_random(int index, int seed);
double uniform_random_Mersenne();

//double normal_random(double mean, double std_dev, mt19937 generator_util_in);
double normal_random(double mean, double std_dev);

double rand_normal(double mean, double stddev);
double rand_normal_Mersenne(double mean, double stddev);

struct wall_struct
{
    vector<double> uw_ij_out;
    vector<double> lambda_out;
    double w_ij_mag_out;
};

struct wall_struct wall_calc_2d(vector<double>& ri_com,vector<double>& rj_com,vector<double>& ui_vec ,vector<double>& uj_vec, double L_i, double L_j,int N_dims, int id);

double speckle_depsilon(double& speck_eps, int& mt_phase, const double v_g, const double v_s, double dt);

#endif