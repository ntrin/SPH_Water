#pragma once
#include <vector>
#include <cmath>
#include <iostream>

using namespace std;

class SPH_main;

class SPH_particle{
public:
    double x[4], v[4];                //position and velocity
    double a[2];  //acceleration

    double rho[2], P;                    //density and pressure
    double drho;  //density derivative

    bool fluid = true;              //determine whether it is a fluid/boundary

    static SPH_main *main_data;        //link to SPH_main class so that it can be used in calc_index

    int list_num[2];                //index in neighbour finding array

    void calc_index(int step);
};


class SPH_main{
public:
    SPH_main();

    void set_values(double length, double height, double dx_in);
    
    void initialise_grid(void);

    void place_points(double* min, double* max, bool fluid_in);

    void initialise_domain(double length, double height, double water_depth, double src_length, double src_height, int domain_case);

    void add_point(double x, double y, bool fluid_in);
    
    void check_overlap();
    
    double distance(SPH_particle* part, SPH_particle* other_part);

    void allocate_to_grid(int step);            //allocates all the points to the search grid (assumes that index has been appropriately updated)

    void neighbour_iterate(SPH_particle *part, int step);

    void check_boundary(SPH_particle* part);
    
    void time_stepping(double dt, double final_t, double c_cfl);
    

    double smoothing_func(SPH_particle* part);

    double h;                                //smoothing length
    double h_fac;
    double dx;                                //particle initial spacing

    double gamma, rho_0, c_0, mu, g, m;

    double min_x[2], max_x[2];                //dimensions of simulation region

    int max_list[2];

    vector<SPH_particle> particle_list;                        //list of all the particles

    vector<vector<vector<SPH_particle*>>> search_grid;        //Outer 2 are the grid, inner vector is the list of pointers in each cell

    double find_new_dt(double c_cfl);
private:

    double tait_eq(double rho);
};
