#include "../includes/SPH_2D.h"
#include "../includes/file_writer.h"
#include <sstream>

#define WRITE_FILE

/********************************* SPH_particle *****************************************/
SPH_main* SPH_particle::main_data;

void SPH_particle::calc_index(int step) {
	for (int i = 0; i < 2; i++) // finds where the particle is situated on grid
		list_num[i] = int((x[i + step*2] - main_data->min_x[i]) / (2.0 * main_data->h));
}

/********************************** SPH_main ********************************************/

SPH_main::SPH_main() {
	SPH_particle::main_data = this;
}


// initialise the variables
void SPH_main::set_values(double length, double height, double dx_in){
	min_x[0] = 0.0;
	min_x[1] = 0.0;

	max_x[0] = length;
	max_x[1] = height;

	dx = dx_in;

	h_fac = 1.3;
	h = dx * h_fac;

	gamma = 7.0;
	rho_0 = 1000.0;
	c_0 = 20.0;
	mu = 0.001;
	g = 9.81;
	m = dx * dx * rho_0;
}


//initialise the size of the grid
void SPH_main::initialise_grid(void) {
	for (int i = 0; i < 2; i++) {
		min_x[i] -= 2.0 * h;
		max_x[i] += 2.0 * h;            //add buffer for virtual wall particles

		max_list[i] = int((max_x[i] - min_x[i]) / (2.0 * h) + 1.0); // how many cells in
	}                                                             // in work region

	search_grid.resize(max_list[0]);
	for (int i = 0; i < max_list[0]; i++)
		search_grid[i].resize(max_list[1]); // set grid size corresponding to work region
}


//initialise all the particles in the domain
void SPH_main::initialise_domain(double length, double height, double water_depth, double src_length, double src_height, int domain_case)
{
	// add top boundary
	for (double i = dx; i < max_x[0]; i += dx) {
		for (double j = height; j < max_x[1]; j += dx) {
			add_point(i, j, false);
		}
	}

	// add bottom boundary
	for (double i = dx; i < max_x[0] ; i += dx) {
		for (double j = 0; j > min_x[1]; j -= dx)
			add_point(i, j, false);
	}

	// add left boundary
	for (double i = 0; i > min_x[0]; i -= dx) {
		for (double j = dx; j < max_x[1]; j += dx)
			add_point(i, j, false);
	}

	// add right boundary
	for (double i = length; i < max_x[0]; i += dx) {
		for (double j = dx; j < height - dx; j += dx)
			add_point(i, j, false);
	}

	// add bottom corner left boundary
	for (double i = 0; i > min_x[0]; i -= dx) {
		for (double j = 0; j > min_x[1]; j -= dx)
			add_point(i, j, false);
	}

	// CASE 1 - NO SLOPE (SIMPLE CASE)
	if (domain_case == 1) {

		// add fluid
		for (double j = dx; j < water_depth; j += dx) {
			for (double i = dx; i < length - dx; i += dx)
				add_point(i, j, true);
		}

	}
	// CASE 2 - SINGLE SLOPE
	else if (domain_case == 2) {

		// add slope boundary
		double start_xslope = 10.0, stop_yslope = 5.0;
		double imax, j_max;
		double gradient = stop_yslope / (length - start_xslope);
		double c = -gradient * start_xslope;

		for (double i = start_xslope; i < length - dx; i += dx) {
			j_max = gradient * i + c;
			for (double j = dx; j < j_max; j += dx)
				add_point(i, j, false);
		}

		// add fluid
		for (double j = dx; j < water_depth; j += dx) {
			imax = (j - c) / gradient;
			for (double i = dx; i < imax; i += dx)
				add_point(i, j, true);
		}

	}
	// CASE 3 - DOUBLE SLOPE
	else if (domain_case == 3) {

		// add slope boundary
		double grad1 = 0.4;
		double x1 = 10.0, y1 = 0.0;
		double x2 = 15.0, y2 = grad1 * (x2 - x1);
		double x3 = length, y3 = 3.0;
		double grad2 = (y3 - y2) / (x3 - x2);
		double c1 = -grad1 * x1;
		double c2 = (y2 + y3 - grad2 * (x3 + x2)) / 2.0;

		for (double i = x1; i < x2 - dx; i += dx) {
			double j_max = grad1 * i + c1;
			for (double j = dx; j < j_max; j += dx)
				add_point(i, j, false);
		}

		for (double i = x2; i < x3 - dx; i += dx) {
			double j_max = grad2 * i + c2;
			for (double j = dx; j < j_max; j += dx)
				add_point(i, j, false);
		}

		// add fluid
		for (double j = dx; j < water_depth; j += dx) {

			double imax;

			if (j <= y2) imax = (j - c1) / grad1;
			else if (j <= y3) imax = (j - c2) / grad2;
			else imax = length - dx;

			for (double i = dx; i < imax; i += dx)
				add_point(i, j, true);
		}
	}
	else{

		cerr << "Error: Case not implemented!" << endl;
		exit(0);
	}


	// add source
	for (double i = dx; i <= src_length; i += dx) {
		for (double j = water_depth + dx; j < water_depth + src_height; j += dx)
			add_point(i, j, true);
	}

	check_overlap();  // delete overlapping particles (esp important for CASE 3!)

}


//create a particle (fluid/boundary) at the location x, y
void SPH_main::add_point(double x, double y, bool fluid_in) {

	SPH_particle particle;

	// initial position
	particle.x[0] = x;
	particle.x[1] = y;

	particle.calc_index(0);   // initial grid locations

	for (int i = 0; i < 2; i++)
		particle.v[i] = 0;  // initial velocity

	particle.rho[0] = rho_0;   // initial density

	particle.P = tait_eq(rho_0);  // initial pressure

	particle.fluid = fluid_in; // fluid/boundary

	particle_list.push_back(particle);
}


//check if any particle overlaps (called once)
void SPH_main::check_overlap() {

	int cnt_overlap = 0;
	int temp;

	for (int n = 0; n < particle_list.size(); n++) {
		for (int p = 0; p < particle_list.size(); p++) {
			if (n != p)
			{
				double dist;
				if (distance(&particle_list[p], &particle_list[n]) < dx / 2.)
				{
					if (particle_list[p].fluid) {
						particle_list.erase(particle_list.begin() + p); // first one is fluid
					}
					else if (particle_list[n].fluid) {
						particle_list.erase(particle_list.begin() + n); // second one is fluid
					}
					else {
						particle_list.erase(particle_list.begin() + n); // both are boundary - so can just delete any of them
					}
				}

			}

		}
	}

}


//calculate distance between two particles
double SPH_main::distance(SPH_particle* part, SPH_particle* other_part) {

	double dr[2];   //vector from 1st to 2nd particle

	for (int n = 0; n < 2; n++)
		dr[n] = part->x[n] - other_part->x[n];

	return sqrt(dr[0] * dr[0] + dr[1] * dr[1]);;
}


//allocates all the points to the search grid (assumes that index has been appropriately updated)
//needs to be called each time that all the particles have their positions updated
void SPH_main::allocate_to_grid(int step) {

	// clear the grids
	for (int i = 0; i < max_list[0]; i++)
		for (int j = 0; j < max_list[1]; j++)
			search_grid[i][j].clear();

	// update the grid based on current location
	for (unsigned int cnt = 0; cnt < particle_list.size(); cnt++) {
		particle_list[cnt].calc_index(step);
		search_grid[particle_list[cnt].list_num[0]][particle_list[cnt].list_num[1]].push_back(&particle_list[cnt]);
	}
}

double SPH_main::find_new_dt(double c_cfl)
{
	double t_f, t_a;
	double dt_temp = 0.;

	for (int cnt = 0; cnt<particle_list.size(); cnt++)
	{

		t_f = sqrt(h/sqrt(pow(particle_list[cnt].a[0],2)+pow(particle_list[cnt].a[1], 2)));
		t_a =h/(c_0 * sqrt(pow(particle_list[cnt].rho[0]/rho_0, gamma-1)));

		cout << "rho0" << rho_0 << endl;
//		if (cnt == 0)
//		{
//			return t_f * c_cfl; // return each thing instead of returning at the end.
//		}
//
//		if (dt_temp > t_f)
//		{
//			return t_f * c_cfl;
//
//		}
//
//		if (dt_temp > t_a)
//		{
//			return t_f * c_cfl;
//
//		}
	}
	return t_f * c_cfl;

}

//update a and drho for each particle
void SPH_main::neighbour_iterate(SPH_particle* part, int step) {

	// iterates over all particles within 2h of part - can be made more efficient using
	// a stencil and realising that all interactions are symmetric
	SPH_particle* other_part;

	double dist;            //distance between particles
	double dr[2];            //vector from 1st to 2nd particle
	double dv[2];           //relative velocity
	double e[2];            //unit vector from 1st to 2nd particle

	part->a[0] = 0.0;
	part->a[1] = 0.0;
	part->drho = 0.0;

	//    double M_PI = 3.148;

	for (int i = part->list_num[0] - 1; i <= part->list_num[0] + 1; i++)
	{
		if (i >= 0 && i < max_list[0])
		{
			for (int j = part->list_num[1] - 1; j <= part->list_num[1] + 1; j++)
			{
				if (j >= 0 && j < max_list[1])
				{
					for (unsigned int cnt = 0; cnt < search_grid[i][j].size(); cnt++) // Num of particles in grid square (i,j)
					{
						other_part = search_grid[i][j][cnt];
						if (part != other_part) //stops particle interacting with itself
						{
							//Calculates the distance between potential neighbours
							for (int n = 0; n < 2; n++)
							{
								dr[n] = part->x[n + step*2] - other_part->x[n + step*2];
								dv[n] = part->v[n + step*2] - other_part->v[n + step*2];
							}

							dist = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
							e[0] = dr[0] / dist;
							e[1] = dr[1] / dist;

							double q = dist / h;
							double dWdr;

							if (q <= 1.0)  //only particle within 2h
							{
								dWdr = (10. / (7. * M_PI * h * h * h)) * (-3. * q + 9. * q * q / 4.);
								cout << "dwdr " << dWdr << endl;

								for (int k = 0; k < 2; k++)
								{
									part->a[k] += -m * (part->P / pow(part->rho[step], 2) + other_part->P / pow(other_part->rho[step], 2)) * dWdr * e[k]; // first sum
									part->a[k] += mu * m * (1.0 / pow(part->rho[step], 2) + 1.0 / pow(other_part->rho[step], 2)) * dWdr * dv[k] / dist; // second sum                                }
									part->drho += m * dWdr * (dv[0] * e[0] + dv[1] * e[1]);
								}
								if (1.0 < q && q <= 2.0)
								{
									dWdr = (-30. / (28. * M_PI * h * h * h)) * (2. - q) * (2. - q);
									cout << "dwdr " << dWdr << endl;
									for (int k = 0; k < 2; k++)
									{
										part->a[k] += -m * (part->P / pow(part->rho[step], 2) + other_part->P / pow(other_part->rho[step], 2)) * dWdr * e[k]; // first sum
										part->a[k] += mu * m * (1.0 / pow(part->rho[step], 2) + 1.0 / pow(other_part->rho[step], 2)) * dWdr * dv[k] / dist; // second sum                                }
										part->drho += m * dWdr * (dv[0] * e[0] + dv[1] * e[1]);
									}
								}
								cout << "drho " << part->drho << endl;
								cout << "P " << part->P << endl;
							}
						}
					}
				}
			}
			part->a[1] -= g;
			cout << "a0 " << part->a[0] << endl;
			cout << "a1 " << part->a[1] << endl;
		}
	}
}

//implement forward euler method to update values of x, v, rho and P
void SPH_main::time_stepping(double dt, double final_t, double c_cfl) {

	double t_current = 0;
	int count = 0;

	while (t_current < final_t) {

#ifdef WRITE_FILE
		// write to file for every 10 timesteps
		if (count % 10 == 0) {

			stringstream fname;
			fname << "output_" << count << ".vtp";
			write_file(fname.str().c_str(), &particle_list);
		}
#endif

		// do smoothing function for density
		if (count % 10 == 0 && count != 0) {

			vector<double> updated_rho;
			updated_rho.resize(particle_list.size());

			// store the updated density
			for (int n = 0; n < particle_list.size(); n++)
				updated_rho[n] = smoothing_func(&particle_list[n]);

			// update the density
			for (int n = 0; n < particle_list.size(); n++)
				particle_list[n].rho[0] = updated_rho[n];
		}


		for (int step = 0; step < 2; step++){

			// calculate the values of a and drho for all particles
			for (int n = 0; n < particle_list.size(); n++)
				neighbour_iterate(&particle_list[n], step);

			// do the time stepping
			for (int n = 0; n < particle_list.size(); n++)
			{
				if (particle_list[n].fluid) // only update v and x for fluid (not boundary)
				{
					for (int i = 0; i < 2; i++) {
						particle_list[n].x[2 + i - (step*2)] = particle_list[n].x[i] + 0.5*(step + 1.)*dt*particle_list[n].v[i + step*2]; // combine the equations
						//						cout << n << "particle_list" << particle_list[n].x[2 + 1] << endl;
						particle_list[n].v[2 + i - (step*2)] = particle_list[n].v[i] + 0.5*(step + 1.)*dt*particle_list[n].a[i];
					}

					check_boundary(&particle_list[n]); // check if the particles are still within the domain
				}

				//particle_list[n].full_step(dt);
				particle_list[n].rho[1 - step] = particle_list[n].rho[0] + 0.5*(step + 1.)*dt*particle_list[n].drho;

				// update the pressure
				particle_list[n].P = tait_eq(particle_list[n].rho[1 - step]);
				cout << "x " << particle_list[n].x[1];

			}
			allocate_to_grid(step);

		}


		cout << t_current << endl;

		dt = find_new_dt(c_cfl);
		cout<< "new dt is "<<dt<<endl;
		t_current += dt;
		count++;
	}
}


//check if the particle exceeds the boundaries
void SPH_main::check_boundary(SPH_particle* part){

	// check if exceed vertical boundary
	if (part->x[0] < 0.)
	{
		part->x[0] = -part->x[0]; // bounce back
		part->v[0] = 0; part->v[1] = 0; // velocity becomes zero
	}
	else if (part->x[2] < 0.){
		part->x[2] = -part->x[2]; // bounce back
		part->v[2] = 0; part->v[3] = 0; // velocity becomes zero
	}
	else if (part->x[0] > (max_x[0] - 2. * h))
	{
		part->x[0] = 2. * max_x[0] - 4. * h - part->x[0]; // bounce back
		part->v[0] = 0; part->v[1] = 0;  // velocity becomes zero
	}
	else if (part->x[2] > (max_x[0] - 2. * h)){
		part->x[2] = 2. * max_x[2] - 4. * h - part->x[2]; // bounce back
		part->v[2] = 0; part->v[3] = 0;  // velocity becomes zero
	}

	// check if exceed horizontal boundary
	if (part->x[1] < 0.)
	{
		part->x[1] = -part->x[1]; // bounce back
		part->v[0] = 0; part->v[1] = 0; // velocity becomes zero
	}
	else if (part->x[3] < 0.){
		part->x[3] = -part->x[3]; // bounce back
		part->v[2] = 0; part->v[3] = 0; // velocity becomes zero
	}
	else if (part->x[1] > (max_x[1] - 2. * h))
	{
		part->x[1] = 2. * max_x[1] - 4. * h - part->x[1]; // bounce back
		part->v[0] = 0; part->v[1] = 0;  // velocity becomes zero
	}
	else if (part->x[3] > (max_x[3] - 2. * h)){
		part->x[3] = 2. * max_x[3] - 4. * h - part->x[3]; // bounce back
		part->v[2] = 0; part->v[3] = 0;  // velocity becomes zero
	}
}


//smooth out density values for each particle
double SPH_main::smoothing_func(SPH_particle* part) {

	SPH_particle* other_part;
	double dist;    //distance between particles
	double dr[2];   //vector from 1st to 2nd particle
	//    double M_PI = 3.148;
	double W = 0, W_p = 0;

	for (int i = part->list_num[0] - 1; i <= part->list_num[0] + 1; i++)
	{
		if (i >= 0 && i < max_list[0])
		{
			for (int j = part->list_num[1] - 1; j <= part->list_num[1] + 1; j++)
			{
				if (j >= 0 && j < max_list[1])
				{
					for (unsigned int cnt = 0; cnt < search_grid[i][j].size(); cnt++) // Num of particles in grid square (i,j)
					{
						other_part = search_grid[i][j][cnt];

						//Calculates the distance between potential neighbours
						for (int n = 0; n < 2; n++)
							dr[n] = part->x[n] - other_part->x[n];
						dist = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
						double q = dist / h;

						if (q <= 1)  //only particle within 2h
						{
							double temp;
							temp = (10. / (7. * M_PI * h * h)) * (1 - 3. * q * q / 2. + 3. * q * q * q / 4.);
							W += temp;
							W_p += temp / other_part->rho[0];
						}
						else if (1 < q && q <= 2)
						{
							double temp;
							temp = (10. / (28. * M_PI * h * h)) * (2. - q) * (2. - q) * (2. - q);
							W += temp;
							W_p += temp / other_part->rho[0];
						}

					}
				}
			}
		}
	}
	return (W / W_p);
}


// calculate the pressure from density
double SPH_main::tait_eq(double rho)
{
	double B = rho_0 * pow(c_0, 2) / gamma;
	return B * (pow((rho / rho_0), gamma) - 1.0);
}

