#include "../includes/SPH_2D.h"
#include "../includes/file_writer.h"

SPH_main domain;

int main(void)
{
	// Domain configuration
	double length = 20.0, height = 10.0, water_depth = 2.0;
	double src_length = 3.0, src_height = 3.0, dx = 0.2;
	int domain_case = 1;

	domain.set_values(length, height, dx);	// set simulation parameters

	domain.initialise_grid();	// initialise simulation grid

	// places initial points - include the boundary points and the specifics of where the fluid is in the domain
	// slope will add a sloping boundary near the bottom right corner of the domain (can change steepness in the function implementation)
	domain.initialise_domain(length, height, water_depth, src_length, src_height, domain_case);				


	// Time-stepping configuration
	double final_t = 10.0;
	double dt = 0.1 * domain.h / domain.c_0;

	domain.time_stepping(dt, final_t, 0.1);

	return 0;
}
