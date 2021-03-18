#ifndef branching_H
#define branching_H

using namespace std;

class branching
{
    
    private:
    public:
    
    int branch_state = 0;
    
    int fil_m = 0; // mother filament
    int fil_d = 0; // daughter filament
    
    double extension_mag = 0.0; // extension away from eqm for detachment
    
    double eps_m1 = 0.0; // mother epsilon 1
    double eps_m2 = 0.0; // mother epsilon 2
    double eps_d1 = 0.0; // daughter epsilon 1
    double eps_d2 = 0.0; // daughter epsilon 2    

    branching();

    void branch_deactivate();
    void branch_walk(double fil_length_i,double fil_length_j, double dt, int fil_i_state, int fil_j_state, double v_g, double v_s);

};

#endif