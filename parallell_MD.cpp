//#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <mpi.h>

using namespace std;

struct particle {
  double x, y, z, vx, vy, vz, fx, fy, fz, fxo, fyo, fzo;
};

typedef struct particle particle;
typedef struct particle* ptr_particle;

#define  N_PART 800// Number of particles is declared here

particle init(particle p[], double box [], double delta_t, double temp,int num_part,int proc,int proc_vect []);
void lattice_position(particle p[], double box [],int proc,int N_proc,int num_part,int proc_vect []);
particle force(particle p[], double* energy,int num_part,double box[]);
 void lennard_jones(double xr, particle*  p_1, particle*  p_2, double x_d, double y_d, double z_d, double* energy);
 void  periodic_bc(double*x_d, double* y_d, double *z_d, double box[]);
particle integrate(particle p[], double * energy, double delta_t, double *pot_en, double*kin_en,int num_part,double box[]);
double randomize(double*r);



int main(int argc, char** argv)
{

  /*---- Initialize the MPI environment------*/
  MPI_Init(NULL, NULL);
  // Find out rank, size
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int P;
   MPI_Comm_size(MPI_COMM_WORLD, &P);
   // int N,I; // total number of elements
   /*-------------------------------------------*/
   int proc_vect [3] ; //number of processors per dimension
   proc_vect[0]=2;
    proc_vect[1]=2;
     proc_vect[2]=2;
     int num_part=N_PART/P; // number of particles per processor
   particle p[num_part];
   double box[3];//dimensions of box
   box[0]=10;
   box[1]=10;
   box[2]=10;
   double sub_box[3]; // dimensions per processor
     sub_box[0]=box[0]/proc_vect[0];
     sub_box[1]=box[1]/proc_vect[1];
     sub_box[2]=box[2]/proc_vect[2];
     double en=0;
     double kin_en,pot_en;
     double delta_t = 0.02;
     double temp = 0.7876;
   if (rank==1){ 
     	FILE *fp;
	
	
     fp = fopen("/home/joel/Dokument/SF2568/MD_simulation/dist_1.txt","w+");
lattice_position(p, sub_box,rank,P,num_part, proc_vect);


 
 init( p,sub_box,delta_t,temp,num_part,rank,proc_vect);
 force(p, &en,num_part,sub_box);
  integrate(p, &en, delta_t, &pot_en, &kin_en,num_part,box);
 
  /* for(int j=0;j<num_part;j++){

   cout<<p[j].fx<< " "<<p[j].fy << " "<< p[j].fz<< '\n';
     
   }*/
 
for(int i=0;i<num_part;i++){
		fprintf(fp, "%f",p[i].x );
		fprintf(fp, "   ");
		fprintf(fp, "%f", p[i].y);
		fprintf(fp, "   ");
		fprintf(fp, "%f", p[i].z);
		fprintf(fp, "\n"); 
     }
     fclose(fp);
     }

     if (rank==0){
     	FILE *fp;
	
     fp = fopen("/home/joel/Dokument/SF2568/MD_simulation/dist_0.txt","w+");
lattice_position(p,  sub_box,rank,P,num_part, proc_vect);
init( p,sub_box,delta_t,temp,num_part,rank,proc_vect);
 force(p, &en,num_part,sub_box);
 integrate(p, &en, delta_t, &pot_en, &kin_en,num_part,box);

for(int i=0;i<num_part;i++){
		fprintf(fp, "%f",p[i].x );
		fprintf(fp, "   ");
		fprintf(fp, "%f", p[i].y);
		fprintf(fp, "   ");
		fprintf(fp, "%f", p[i].z);
		fprintf(fp, "\n"); 
     }
     fclose(fp);
     }

      if (rank==2){
     	FILE *fp;
	
     fp = fopen("/home/joel/Dokument/SF2568/MD_simulation/dist_2.txt","w+");
lattice_position(p,  sub_box,rank,P,num_part, proc_vect);
 init( p,sub_box,delta_t,temp,num_part,rank,proc_vect);
 force(p, &en,num_part,sub_box);
 integrate(p, &en, delta_t, &pot_en, &kin_en,num_part,box);
for(int i=0;i<num_part;i++){
		fprintf(fp, "%f",p[i].x );
		fprintf(fp, "   ");
		fprintf(fp, "%f", p[i].y);
		fprintf(fp, "   ");
		fprintf(fp, "%f", p[i].z);
		fprintf(fp, "\n"); 
     }
     fclose(fp);
     }
 if (rank==3){
     	FILE *fp;
	
     fp = fopen("/home/joel/Dokument/SF2568/MD_simulation/dist_3.txt","w+");
lattice_position(p,  sub_box,rank,P,num_part, proc_vect);
 init( p,sub_box,delta_t,temp,num_part,rank,proc_vect);
 force(p, &en,num_part,sub_box);
integrate(p, &en, delta_t, &pot_en, &kin_en,num_part,box);
for(int i=0;i<num_part;i++){
		fprintf(fp, "%f",p[i].x );
		fprintf(fp, "   ");
		fprintf(fp, "%f", p[i].y);
		fprintf(fp, "   ");
		fprintf(fp, "%f", p[i].z);
		fprintf(fp, "\n"); 
     }
     fclose(fp);
     }
  if (rank==4){
     	FILE *fp;
	
     fp = fopen("/home/joel/Dokument/SF2568/MD_simulation/dist_4.txt","w+");
lattice_position(p,  sub_box,rank,P,num_part, proc_vect);
 init( p,sub_box,delta_t,temp,num_part,rank,proc_vect);
force(p, &en,num_part,sub_box);
 integrate(p, &en, delta_t, &pot_en, &kin_en,num_part,box);
for(int i=0;i<num_part;i++){
		fprintf(fp, "%f",p[i].x );
		fprintf(fp, "   ");
		fprintf(fp, "%f", p[i].y);
		fprintf(fp, "   ");
		fprintf(fp, "%f", p[i].z);
		fprintf(fp, "\n"); 
     }
     fclose(fp);
     }
    if (rank==5){
     	FILE *fp;
	
     fp = fopen("/home/joel/Dokument/SF2568/MD_simulation/dist_5.txt","w+");
lattice_position(p,  sub_box,rank,P,num_part, proc_vect);
 init( p,sub_box,delta_t,temp,num_part,rank,proc_vect);
 force(p, &en,num_part,sub_box);
 integrate(p, &en, delta_t, &pot_en, &kin_en,num_part,box);
for(int i=0;i<num_part;i++){
		fprintf(fp, "%f",p[i].x );
		fprintf(fp, "   ");
		fprintf(fp, "%f", p[i].y);
		fprintf(fp, "   ");
		fprintf(fp, "%f", p[i].z);
		fprintf(fp, "\n"); 
     }
     fclose(fp);
     }
  



  
  //MPI_Send(A, size_1, MPI_DOUBLE,p-1, 0, MPI_COMM_WORLD);
  // MPI_Recv(A, size_1, MPI_DOUBLE, p-1, 0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);


  /*
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
   /*for (int i = 0; i < TIME_STEPS; i++) {
		integrate(p, &en, delta_t, &pot_en, &kin_en);

		
		/*printing energies to file*/
   /*		fprintf(fp, "   ");
		fprintf(fp, "%f", en);
		fprintf(fp, "   ");
		fprintf(fp, "%f", pot_en);
		fprintf(fp, "   ");
		fprintf(fp, "%f", kin_en);
		fprintf(fp, "\n"); 
		time = time + delta_t; 
	}
	

	fclose(fp);
   */	
 MPI_Finalize();
	return 0;
}



particle init(particle p[], double box[], double delta_t, double temp,int num_part,int proc,int proc_vect []) {
	double r;
 double box_min_x,box_min_y,box_min_z;
  double box_max_x,box_max_y,box_max_z;
  int p_x,p_y,p_z;


		/*process ID*/
	p_x=(proc/(proc_vect[1]*proc_vect[2]));
	p_y=((proc/proc_vect[2]) % proc_vect[1]);
	p_z=(proc % proc_vect[2]);
	

	if(p_x==proc_vect[0]){
	   box_min_x = box[0]*(p_x-1);
           box_max_x =box[0]*p_x;
	  
	}
	else if(p_x!=proc_vect[0]){
	   box_min_x = box[0]*p_x;
           box_max_x =box[0]*(p_x+1);
	}

	if(p_y==proc_vect[1]){
	   box_min_y = box[1]*(p_y-1);
           box_max_y =box[1]*p_y;
	}
	else if(p_y!=proc_vect[1]){
	   box_min_y = box[1]*p_y;
           box_max_y =box[1]*(p_y+1);
	}

	if(p_z==proc_vect[2]){
	   box_min_z = box[2]*(p_z-1);
           box_max_z =box[2]*p_z;
	}
	else if(p_z!=proc_vect[2]){
	   box_min_z = box[2]*p_z;
           box_max_z =box[2]*(p_z+1);
	}

	//	lattice_position(p, box_length); // Placing particles on cubic lattice
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
	for (int i = 0; i < num_part; i++) {
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
	kin_enx = kin_enx / num_part;
	kin_eny = kin_eny / num_part;
	kin_enz = kin_enz / num_part;

	scale_factorx = sqrt(temp / kin_enx);
	scale_factory = sqrt(temp / kin_eny);
	scale_factorz = sqrt(temp / kin_enz);

	sum_vx = sum_vx / num_part;
	sum_vy = sum_vy / num_part;
	sum_vz = sum_vz / num_part;

	/*Calculating new velocities and positions*/
	for (int j = 0; j < num_part; j++) {
		p[j].vx = (p[j].vx - sum_vx)*scale_factorx;
		p[j].vy = (p[j].vy - sum_vy)*scale_factory;
		p[j].vz = (p[j].vz - sum_vz)*scale_factorz;

		p[j].x = p[j].x - p[j].vx*delta_t;		// Implement boundary conditions
		p[j].y = p[j].y - p[j].vy*delta_t;
		p[j].z = p[j].z - p[j].vz*delta_t;

		if (p[j].x > box_max_x) {				// periodic boundary condition
			p[j].x = p[j].x - (box_max_x-box_min_x);
		}

			else if (p[j].x < box_min_x) {					// periodic boundary condition
			  p[j].x = p[j].x + (box_max_x-box_min_x); 
			}
		if (p[j].y > box_max_y) {				// periodic boundary condition
			p[j].y = p[j].y - (box_max_y-box_min_y);
		}

		else if (p[j].y < box_min_x) {					// periodic boundary condition
			p[j].y = p[j].y + (box_max_y-box_min_y);
		}
		if (p[j].z > box_max_z) {				// periodic boundary condition
			p[j].z = p[j].z - (box_max_z-box_min_z);
		}

		else if (p[j].z < box_min_z) {					// periodic boundary condition
			p[j].z = p[j].z + (box_max_x-box_min_x);
		}
	}
	
	sum_vx = 0;
	sum_vy = 0;
	sum_vz = 0;
	kin_en = 0;

	/*Calculating momentum and energy for the system*/
	for (int n = 0; n < num_part; n++) {
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
void lattice_position(particle p[], double box [],int proc,int N_proc,int num_part,int proc_vect []) {
  double box_min_x,box_min_y,box_min_z;
  double box_max_x,box_max_y,box_max_z;
  int p_x,p_y,p_z;
	double r;
	/*process ID*/
	p_x=(proc/(proc_vect[1]*proc_vect[2]));
	p_y=((proc/proc_vect[2]) % proc_vect[1]);
	p_z=(proc % proc_vect[2]);
	

	if(p_x==proc_vect[0]){
	   box_min_x = box[0]*(p_x-1);
           box_max_x =box[0]*p_x;
	  
	}
	else if(p_x!=proc_vect[0]){
	   box_min_x = box[0]*p_x;
           box_max_x =box[0]*(p_x+1);
	}

	if(p_y==proc_vect[1]){
	   box_min_y = box[1]*(p_y-1);
           box_max_y =box[1]*p_y;
	}
	else if(p_y!=proc_vect[1]){
	   box_min_y = box[1]*p_y;
           box_max_y =box[1]*(p_y+1);
	}

	if(p_z==proc_vect[2]){
	   box_min_z = box[2]*(p_z-1);
           box_max_z =box[2]*p_z;
	}
	else if(p_z!=proc_vect[2]){
	   box_min_z = box[2]*p_z;
           box_max_z =box[2]*(p_z+1);
	}


	/*	if(proc==1){
	  cout<<proc << "\n"; 
	  cout<<box_min_x << " "<< box_max_x << "\n";
	   cout<<box_min_y << " "<< box_max_y << "\n";
	    cout<<box_min_z << " "<< box_max_z << "\n";
	    }*/
	

	/*----placing particles on random positions within domain----*/
	/*	for(int i=0;i<num_part;i++){
	  p[i].x= box_min_x+randomize(&r)*(box_max_x-box_min_x);
	  p[i].y= box_min_y+randomize(&r)*(box_max_y-box_min_y);
	  p[i].z= box_min_z+randomize(&r)*(box_max_z-box_min_z);

	  }*/
	
	double volume=box[0]*box[1]*box[2];
	int part_x,part_y,part_z;
	double x_dist,y_dist,z_dist;
	part_x=(num_part*box[0])/volume;
	part_y=(num_part*box[1])/volume;
	part_z=(num_part*box[2])/volume;
	x_dist=box[0]/part_x;
	y_dist=box[1]/part_y;
	z_dist=box[2]/part_z;

	/*	double x,y,z;
	x=box_min_x;
	y=box_min_y;
	z=box_min_z;*/
	/*	int i,j,k;
	i=0;
	j=0;
	k=0;

        while(x<box_max_x && k<part_x){
	  while(y<box_max_y && j<part_y){
	    while(z<box_max_z && i<part_z){
	      p[i].x=x;
	      p[i].y=y;
	      p[i].z=z;
	      z=z+z_dist;
	      i++;
	    }
	    i=0;
	    j++;
	    y=y+y_dist;
	    z=box_min_z;
	  }
	  j=0;
	  k++;
	  x=x+x_dist;
	  y=box_min_y;
	  }*/

	/*
	for(int i=0;i<part_x;i++){
	  for(int j=0;j<part_y;j++){
	    for(int k=0;k<part_z;k++){
	      p[part_z*part_y*i+part_z*j+k].x=(double)(i*x_dist);
	      p[part_z*part_y*i+part_z*j+k].y=(double)(j*y_dist);
	      p[part_z*part_y*i+part_z*j+k].z=(double)(k*z_dist);
	    }
	  }
	  
	  }*/
	
	
	double spacing = (double)cbrt((double)((box[0]*box[1]*box[2])/num_part));
	double x, y, z;
	x = box_min_x;
	y = box_min_y;
	z = box_min_z;
	int i = 0;


	while (x< box_max_x && i < num_part) {
		p[i].x = x;
		/*Setting y for fixed x*/
		while (y< box_max_y && i<num_part) {
			p[i].y = y;
			/*Setting z for fixed x and y*/
			while (z< box_max_z && i<num_part) {
				p[i].z = z;
				z = z + spacing;
				if ((i + 1)< num_part) {
					i++;
				}
			       	p[i].x = x;
				p[i].y = y;
			}
			y = y + spacing;
			z = box_min_z;
			p[i].x = x;
		}
		x = x + spacing;
		y = box_min_y;
		z = box_min_z;
	}



}



particle force(particle p[], double* energy,int num_part,double box[]) {
	*energy = 0;

	double x_d, y_d, z_d, dt; // distance components
     
	const double cutoff_dist = 6.25;



	for (int i = 0; i < num_part; i++) { //setting forces to zero
		p[i].fx = 0.0;
		p[i].fy = 0.0;
		p[i].fz = 0.0;
	}

	/* Calculating distances between particles*/
	for (int j = 0; j <num_part; j++) {

		for (int i = j + 1; i < num_part; i++) {

			x_d = p[j].x - p[i].x;
			y_d = p[j].y - p[i].y;
			z_d = p[j].z - p[i].z;
			//periodic_bc(&x_d, &y_d, &z_d, box);
			dt = x_d*x_d + y_d*y_d + z_d*z_d; // radial distance between particles squared
			if (dt < cutoff_dist) {           // if distance is less than cutoff distance,calculate forces between particles
				lennard_jones(dt, &p[j], &p[i], x_d, y_d, z_d, energy);

			}
		}
	}
	return *p;

}


 void lennard_jones(double r_2, particle*  p_1, particle*  p_2, double x_d, double y_d, double z_d, double* energy) {

	double ff;
	//	const double rc_6inv = 1/pow(2.5, 6);
	//	const double ecut = 4*rc_6inv*(rc_6inv - 1);

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
	*energy = *energy + 4 * r6_inv*(r6_inv - 1);//-ecut;

}


particle integrate(particle p[], double * energy, double delta_t, double *pot_en, double*kin_en,int num_part,double box[]) {
	/*add constants*/
	double k_energy = 0;
	double delta_tsq = delta_t*delta_t;
	int num_box_x,num_box_y,num_box_z;

	/*Updating positions for all the particles*/
	for (int i = 0; i < num_part; i++) {
		/*Calculating position for the next time-step*/
		p[i].x = p[i].x + p[i].vx*delta_t + (p[i].fx*delta_tsq)*0.5;

		//num_box_x=p[i].x /box[0];
		if (p[i].x > box[0]) {			// periodic boundary condition
			p[i].x = p[i].x - box[0];
		}

		else if (p[i].x < 0) {			   // periodic boundary condition
			p[i].x = p[i].x + box[0];
		}
		p[i].y = p[i].y + p[i].vy*delta_t + (p[i].fy*delta_tsq)*0.5;

		if (p[i].y > box[1]) {
			p[i].y = p[i].y - box[1];
		}

		else if (p[i].y < 0) {
			p[i].y = p[i].y + box[1];
		}

		p[i].z = p[i].z + p[i].vz*delta_t + (p[i].fz*delta_tsq)*0.5;

		if (p[i].z > box[2]) {
			p[i].z = p[i].z - box[2];
		}

		else if (p[i].z < 0) {
			p[i].z = p[i].z + box[2];
		}
		/*Storing previous force*/
		p[i].fxo = p[i].fx; 
		p[i].fyo = p[i].fy;
		p[i].fzo = p[i].fz;
	}
	/*Updating forces for the new positions*/
	force(p, energy,num_part,box);

	for (int i = 0; i < num_part; i++) {
		/*Calculating new velocities*/
		p[i].vx = p[i].vx + delta_t*(p[i].fxo + p[i].fx)*0.5;
		p[i].vy = p[i].vy + delta_t*(p[i].fyo + p[i].fy)*0.5;
		p[i].vz = p[i].vz + delta_t*(p[i].fzo + p[i].fz)*0.5;

		/*Calculation of the kinetic energy and momenta*/
		k_energy = k_energy + (p[i].vx*p[i].vx + p[i].vy*p[i].vy + p[i].vz*p[i].vz)*0.5;
	}

	*pot_en = *energy;
	*kin_en = k_energy;
	*energy = *energy + k_energy;
	return *p;
}

 void  periodic_bc(double*x_d, double* y_d, double *z_d, double box[]) {
   //	double half_box_length = box_length*0.5;
	if (*x_d > box[0]*0.5) {
		*x_d = *x_d - box[0];
	}
	else if (*x_d<-box[0]*0.5) {
		*x_d = *x_d + box[0];
	}
	if (*y_d>box[1]*0.5) {
		*y_d = *y_d - box[1];
	}
	else if (*y_d<-box[1]*0.5) {
		*y_d = *y_d + box[1];
	}
	if (*z_d> box[2]*0.5) {
		*z_d = *z_d - box[2];
	}
	else if (*z_d<-box[2]*0.5) {
		*z_d = *z_d + box[2];
	}

}







