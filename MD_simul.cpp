//#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>


struct particle {
	double x, y, z, vx, vy, vz, fx, fy, fz;
};

typedef struct particle particle;
typedef struct particle* ptr_particle;

#define  N_PART 250// Number of particles is declared here

particle init(particle p[], double box_length, double delta_t, double temp);
void lattice_position(particle p[], double box_length);
 particle force(particle p[], double* energy);
 void lennard_jones(double xr, particle * p_1, particle * p_2, double x_d, double y_d, double z_d, double* energy);
 void  periodic_bc(double*x_d, double* y_d, double *z_d, double box_length);
 particle integrate(particle p[], double * energy, double delta_t, double *pot_en, double*kin_en);
void read_write(particle p[], FILE*fp, double *energy, double*pot_en, double*kin_en);
double randomize(double*r);
void r_distr(particle p[], FILE*fp);

/*arrays for storing previous accelerations*/
double prev_ax[N_PART];
double prev_ay[N_PART];
double prev_az[N_PART];


int main()
{
	double en, kin_en, pot_en;

	const double delta_t = 0.001; // Setting time-step 

	particle p[N_PART]; // array of particle structs
	int counter=0;


#define TIME_STEPS 100200   // Number of time-steps
	const double box_length = 10.229; // domain
	const double temp = 0.7876; //reduced temperature(94.4 K*kb/epsilon)
	double time = 0;

	FILE *fp;
	
	fp = fopen("/home/joel/Dokument/SF2568/MD_simulation/energy.txt","w+");
       

	init(p, box_length, delta_t, temp); // initializing particle positions and velocities
	force(p, &en); // Force initialization

	 /*MD-loop, integrates the equations of motion with the verlet algorithm for the particles with a specific number of time steps */
	for (int i = 0; i < TIME_STEPS; i++) {
		integrate(p, &en, delta_t, &pot_en, &kin_en);

		
		/*printing energies to file*/
		fprintf(fp, "   ");
		fprintf(fp, "%f", en);
		fprintf(fp, "   ");
		fprintf(fp, "%f", pot_en);
		fprintf(fp, "   ");
		fprintf(fp, "%f", kin_en);
		fprintf(fp, "\n"); 
		time = time + delta_t; 


	}
	

	fclose(fp);

	return 0;
}

particle init(particle p[], double box_length, double delta_t, double temp) {
	double r;

	lattice_position(p, box_length); // Placing particles on cubic lattice
	double sum_vx = 0;
	double sum_vy = 0;
	double sum_vz = 0;


	double kin_en = 0;
	double kin_enx = 0;
	double kin_eny = 0;
	double kin_enz = 0;


	double scale_factorx;
	double scale_factory;
	double scale_factorz;

	/* setting velocities, calculating total energy and momentum for all components*/
	for (int i = 0; i < N_PART; i++) {
		/*Velocities*/
		p[i].vx = randomize(&r) - 0.5;
		p[i].vy = randomize(&r) - 0.5;
		p[i].vz = randomize(&r) - 0.5;

		/*Momentum*/
		sum_vx = sum_vx + p[i].vx;
		sum_vy = sum_vy + p[i].vy;
		sum_vz = sum_vz + p[i].vz;

		/*Kinetic energy*/
		kin_enx = kin_enx + pow(p[i].vx, 2);
		kin_eny = kin_eny + pow(p[i].vy, 2);
		kin_enz = kin_enz + pow(p[i].vz, 2);

	}
	// mean squared velocity
	kin_enx = kin_enx / N_PART;
	kin_eny = kin_eny / N_PART;
	kin_enz = kin_enz / N_PART;



	scale_factorx = sqrt(temp / kin_enx);
	scale_factory = sqrt(temp / kin_eny);
	scale_factorz = sqrt(temp / kin_enz);

	sum_vx = sum_vx / N_PART;
	sum_vy = sum_vy / N_PART;
	sum_vz = sum_vz / N_PART;

	/*Calculating new velocities and positions*/
	for (int j = 0; j < N_PART; j++) {
		p[j].vx = (p[j].vx - sum_vx)*scale_factorx;
		p[j].vy = (p[j].vy - sum_vy)*scale_factory;
		p[j].vz = (p[j].vz - sum_vz)*scale_factorz;

		p[j].x = p[j].x - p[j].vx*delta_t;		// Implement boundary conditions
		p[j].y = p[j].y - p[j].vy*delta_t;
		p[j].z = p[j].z - p[j].vz*delta_t;

		if (p[j].x > box_length) {				// periodic boundary condition
			p[j].x = p[j].x - box_length;
		}

		else if (p[j].x < 0) {					// periodic boundary condition
			p[j].x = p[j].x + box_length;
		}
		if (p[j].y > box_length) {				// periodic boundary condition
			p[j].y = p[j].y - box_length;
		}

		else if (p[j].y < 0) {					// periodic boundary condition
			p[j].y = p[j].y + box_length;
		}
		if (p[j].z > box_length) {				// periodic boundary condition
			p[j].z = p[j].z - box_length;
		}

		else if (p[j].z < 0) {					// periodic boundary condition
			p[j].z = p[j].z + box_length;
		}

	}
	
	sum_vx = 0;
	sum_vy = 0;
	sum_vz = 0;
	kin_en = 0;

	/*Calculating momentum and energy for the system*/
	for (int n = 0; n < N_PART; n++) {
		sum_vx = sum_vx + p[n].vx;
		sum_vy = sum_vy + p[n].vy;
		sum_vz = sum_vz + p[n].vz;

		kin_en = kin_en + pow(p[n].vx, 2) + pow(p[n].vy, 2) + pow(p[n].vz, 2);

	}

	return *p;
}

double randomize(double*r) {
	*r = ((float)rand()) / RAND_MAX;
	return *r;
}


/*Setting initial particle positions */
void lattice_position(particle p[], double box_length) {
	 // distance between particles, same for all dimensions
	const double spacing = box_length / 10;
	double x, y, z;
	x = 0;
	y = 0;
	z = 0;
	int i = 0;


	while (x< box_length && i < N_PART) {
		p[i].x = x;
		/*Setting y for fixed x*/
		while (y< box_length && i<N_PART) {
			p[i].y = y;
			/*Setting z for fixed x and y*/
			while (z< box_length && i<N_PART) {
				p[i].z = z;
				z = z + spacing;
				if ((i + 1)< N_PART) {
					i++;
				}
				p[i].x = x;
				p[i].y = y;
			}
			y = y + spacing;
			z = 0;
			p[i].x = x;
		}
		x = x + spacing;
		y = 0;
		z = 0;
	}

}



 particle force(particle p[], double* energy) {
	*energy = 0;

	double x_d, y_d, z_d, dt; // distance components
	const double box_length = 10.229;
	const double cutoff_dist = 6.25;



	for (int i = 0; i < N_PART; i++) { //setting forces to zero
		p[i].fx = 0.0;
		p[i].fy = 0.0;
		p[i].fz = 0.0;
	}

	/* Calculating distances between particles*/
	for (int j = 0; j <N_PART; j++) {

		for (int i = j + 1; i < N_PART; i++) {

			x_d = p[j].x - p[i].x;
			y_d = p[j].y - p[i].y;
			z_d = p[j].z - p[i].z;
			periodic_bc(&x_d, &y_d, &z_d, box_length);
			dt = x_d*x_d + y_d*y_d + z_d*z_d; // radial distance between particles squared
			if (dt < cutoff_dist) {           // if distance is less than cutoff distance,calculate forces between particles
				lennard_jones(dt, &p[j], &p[i], x_d, y_d, z_d, energy);

			}
		}
	}
	return *p;

}

 void lennard_jones(double r_2, particle * p_1, particle * p_2, double x_d, double y_d, double z_d, double* energy) {

	double ff;
	const double rc_6inv = 1/pow(2.5, 6);
	const double ecut = 4*rc_6inv*(rc_6inv - 1);

	/*Lennard jones potential*/
	double r2_inv = 1 / r_2;
	double r6_inv = r2_inv*r2_inv*r2_inv;
	ff = 48 * r2_inv*r6_inv*(r6_inv - 0.5);

	/*Updating forces*/
	p_1->fx = p_1->fx + ff*x_d;
	p_1->fy = p_1->fy + ff*y_d;
	p_1->fz = p_1->fz + ff*z_d;


	p_2->fx = p_2->fx - ff*x_d;
	p_2->fy = p_2->fy - ff*y_d;
	p_2->fz = p_2->fz - ff*z_d;


	/*updating potential energy*/
	*energy = *energy + 4 * r6_inv*(r6_inv - 1)-ecut;

}


 particle integrate(particle p[], double * energy, double delta_t, double *pot_en, double*kin_en) {
	/*add constants*/
	double k_energy = 0;
	const double box_length = 10.229;
	double delta_tsq = delta_t*delta_t;

	/*Updating positions for all the particles*/
	for (int i = 0; i < N_PART; i++) {
		/*Calculating position for the next time-step*/
		p[i].x = p[i].x + p[i].vx*delta_t + (p[i].fx*delta_tsq)*0.5;

		if (p[i].x > box_length) {			// periodic boundary condition
			p[i].x = p[i].x - box_length;
		}

		else if (p[i].x < 0) {			   // periodic boundary condition
			p[i].x = p[i].x + box_length;
		}
		p[i].y = p[i].y + p[i].vy*delta_t + (p[i].fy*delta_tsq)*0.5;

		if (p[i].y > box_length) {
			p[i].y = p[i].y - box_length;
		}

		else if (p[i].y < 0) {
			p[i].y = p[i].y + box_length;
		}

		p[i].z = p[i].z + p[i].vz*delta_t + (p[i].fz*delta_tsq)*0.5;

		if (p[i].z > box_length) {
			p[i].z = p[i].z - box_length;
		}

		else if (p[i].z < 0) {
			p[i].z = p[i].z + box_length;
		}
		/*Storing previous force*/
		prev_ax[i] = p[i].fx; 
		prev_ay[i] = p[i].fy;
		prev_az[i] = p[i].fz;
	}
	/*Updating forces for the new positions*/
	force(p, energy);

	for (int i = 0; i < N_PART; i++) {
		/*Calculating new velocities*/
		p[i].vx = p[i].vx + delta_t*(prev_ax[i] + p[i].fx)*0.5;
		p[i].vy = p[i].vy + delta_t*(prev_ay[i] + p[i].fy)*0.5;
		p[i].vz = p[i].vz + delta_t*(prev_az[i] + p[i].fz)*0.5;

		/*Calculation of the kinetic energy and momenta*/
		k_energy = k_energy + (p[i].vx*p[i].vx + p[i].vy*p[i].vy + p[i].vz*p[i].vz)*0.5;
	}

	*pot_en = *energy;
	*kin_en = k_energy;
	*energy = *energy + k_energy;
	return *p;
}

 void  periodic_bc(double*x_d, double* y_d, double *z_d, double box_length) {
	double half_box_length = box_length*0.5;
	if (*x_d > half_box_length) {
		*x_d = *x_d - box_length;
	}
	else if (*x_d<-half_box_length) {
		*x_d = *x_d + box_length;
	}
	if (*y_d>half_box_length) {
		*y_d = *y_d - box_length;
	}
	else if (*y_d<-half_box_length) {
		*y_d = *y_d + box_length;
	}
	if (*z_d> half_box_length) {
		*z_d = *z_d - box_length;
	}
	else if (*z_d<-half_box_length) {
		*z_d = *z_d + box_length;
	}

}







