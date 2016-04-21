#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define WIDTH  1000
#define HEIGHT 1000
#define BODIES 10

typedef struct Particle{
	float x, y; 	// wspolrzedne
	float m;	// masa
	float xv, yv;  // predkosc
	float fx, fy;  // sila
	float ax, ay;  // przyspieszenie
} Particle;

void ComputeForces(Particle *particles, Particle *recvbuf,int rlen,float *max_f){
int i, j;
float xi, yi, mi, rx, ry, mj, r, fx, fy;
float xnew, ynew, rmin;

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

void WriteParticles(Particle *particles){
	int i=0;
	FILE *fp;
	if ((fp=fopen("out.txt", "w"))==NULL) {
     		printf ("Nie mogę otworzyć pliku\n");
     		exit(1);
     	}
	while(i != BODIES){
   		fprintf (fp, "%f\t%f\t%f\n", particles[i].x, particles[i].y, particles[i].m);
		i++;
	}
   	fclose (fp);
}


void PrintParticles(Particle *particles, int rank){
	int i;
	for (i=0; i<BODIES; i++) {
	    printf( "%f \t %f \t %f\n", particles[i].x, particles[i].y, particles[i].m);
	}
}

int main(int argc, char* argv[]) {
    int size, rank, i=0;
    float x, y, m;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    srand(time(0));	
    Particle particles[BODIES];

	// Pobiera dane
    	if(rank==0){
		i=0;
		FILE *fp;
		if ((fp=fopen("data.txt", "r"))==NULL) {
     			printf ("Nie mogę otworzyć pliku\n");
     		exit(1);
    		}
		while(fscanf(fp, "%f %f %f", &x, &y, &m) != EOF){
			particles[i].x =x;
			particles[i].y =y;
			particles[i].m =m;
			i++;
		}
		fclose(fp);
   	}
	// Rozsyla dane poczatkowe, czyli polozenie i mase
	for (i = 0; i < BODIES; i++) {
		MPI_Bcast(&particles[i].x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&particles[i].y, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&particles[i].m, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}

WriteParticles(particles);

  
MPI_Finalize();
return 0;
}
