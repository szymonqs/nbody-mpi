#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define WIDTH  1000
#define HEIGHT 1000
#define BODIES 10

typedef struct Particle{
	double x, y; 	// wspolrzedne
	double m;	// masa
	double xv, yv;  // predkosc
	double fx, fy;  // sila
	double ax, ay;  // przyspieszenie
} Particle;

void ComputeForces(Particle *particles, Particle *recvbuf,int rlen,double *max_f){
int i, j;
double xi, yi, mi, rx, ry, mj, r, fx, fy;
double xnew, ynew, rmin;
/* Compute forces (2D only) */
for (i=0; i<BODIES; i++) {
    rmin = 100.0;
    xi   = particles[i].x;
    yi   = particles[i].y;
    fx   = 0.0;
    fy   = 0.0;
    for (j=0; j<rlen; j++) {
	rx = xi - recvbuf[j].x;
	ry = yi - recvbuf[j].y;
	mj = recvbuf[j].m;
	r  = rx * rx + ry * ry;
	/* ignore overlap and same particle */
	if (r == 0.0) continue;
	if (r < rmin) rmin = r;
	/* compute forces */
	r  = r * sqrt(r);
	fx += mj * rx / r;
	fy += mj * ry / r;
	}
    particles[i].fx -= fx;
    particles[i].fy -= fy;
    /* Compute a rough estimate of (1/m)|df / dx| */
    fx = sqrt(fx*fx + fy*fy)/rmin;
    if (fx > *max_f) *max_f = fx;
    }
}
void PrintParticles(Particle *particles, int rank){
	int i;
	for (i=0; i<BODIES; i++) {
	    printf( "x = %f \t y = %f \t proces = %d\n", particles[i].x, particles[i].y, rank);
	}
}

int main(int argc, char* argv[]) {
    int size, rank, i;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    srand(time(0));	
    Particle particles[BODIES];

	// Pobiera dane
    	if(rank==0){
	   	 for (i = 0; i < BODIES; i++){
			particles[i].x =  random() % WIDTH;
			particles[i].y =  random() % HEIGHT;
			particles[i].m =  random() % 900;
			particles[i].xv = 0;
			particles[i].yv = 0;
			particles[i].fx = 0;
			particles[i].fy = 0;
			particles[i].ax = 0;
			particles[i].ay = 0;
	    	}
   	}
	// Rozsyla dane poczatkowe, czyli polozenie i mase
	for (i = 0; i < BODIES; i++) {
		MPI_Bcast(&particles[i].x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&particles[i].y, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&particles[i].m, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}

PrintParticles(particles, rank);

    
MPI_Finalize();
return 0;
}
