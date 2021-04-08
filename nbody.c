/* Code adapted from Danail Obreschow's Rcpp nbody code (Chapter 15 astrostats) */
// includes options for first order integrator and leapfrog
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ics_sun_earth.h"

void compute_accelerations(particle_t*p, double pos[], double acc[],
                           double mass, double rsmoothsqr) {
      if ((mass>0)||(p->mass>0)) {
        /* compute the acc of particle j on particle i */
        // sun is pointer values, earth is variables
        printf("Sun = %e %e %e \n", p->pos[0], p->pos[1], p->pos[2]);
        printf("Earth = %e %e %e\n", pos[0], pos[1], pos[2]);

        double dx = p->pos[0] - pos[0]; // faster than using a vector dx()
        double dy = p->pos[1] - pos[1];
        double dz = p->pos[2] - pos[2];

        double dsqr = fmax(rsmoothsqr,dx*dx+dy*dy+dz*dz);
        double q = G/pow(dsqr,1.5);

        /* exploit force symmetry */
        acc[0] = (p->mass)*q*dx; // faster than using vector operation on a(_,i)
        acc[1] = (p->mass)*q*dy;
        acc[2] = (p->mass)*q*dz;

        p->acc[0] = mass*q*dx;
        p->acc[1] = mass*q*dy;
        p->acc[2] = mass*q*dz;

        /* TODO:implement adaptive time step calc */
        // double dminsqr = 1e99;
        // dminsqr = fmin(dminsqr,dsqr);
        //
        // /* determine adaptive time step (for this pair of particles) */
        // double zmin = 1e99;
        //
        // double zi = dminsqr/(acc[0]*acc[0]+acc[1]*acc[1]+acc[2]*acc[2]);
        // double zj = dminsqr/(p->acc[0]*p->acc[0]+p->acc[1]*p->acc[1]+p->acc[2]*p->acc[2]);
        // double z = fmin(zi,zj);
        //
        // if (z<zmin) zmin = z;
        // double dtvar = pow(zmin,0.25);

      }
    }


void iteration(double dt, double rsmoothsqr){
    /* below uses leapfrog integration method */
    for(int i=0; i<n; i++) {
    particles[i].vel[0] += particles[i].acc[0]*dt/2.0;
    particles[i].vel[1] += particles[i].acc[1]*dt/2.0;
    particles[i].vel[2] += particles[i].acc[2]*dt/2.0;

    particles[i].pos[0] += particles[i].vel[0]*dt;
    particles[i].pos[1] += particles[i].vel[1]*dt;
    particles[i].pos[2] += particles[i].vel[2]*dt;
    }

    /* compute the acceleration of particle j on particle i */
    for(int i=0; i<(n-1); i++) {
      for(int j=i+1; j<n; j++) {
        particle_t*p = &particles[j];
        compute_accelerations(&particles[i], p->pos, p->acc, p->mass, rsmoothsqr);
      }
    }

    /* update velocities */
    for(int i=0; i<n; i++) {
    particles[i].vel[0] += particles[i].acc[0]*dt/2.0;
    particles[i].vel[1] += particles[i].acc[1]*dt/2.0;
    particles[i].vel[2] += particles[i].acc[2]*dt/2.0;
    }
  }

/* compute the new position/velocity under euler method */
void move_particle(particle_t*p, double dt) {

    p->pos[0] += (p->vel[0])*dt ;
    p->pos[1] += (p->vel[1])*dt ;
    p->pos[2] += (p->vel[2])*dt ;

    p->vel[0] += (p->acc[0])*dt;
    p->vel[1] += (p->acc[1])*dt;
    p->vel[2] += (p->acc[2])*dt;
  }

void iteration_euler(double dt, double rsmoothsqr){
    /* compute the acceleration of particle j on particle i */
    for(int i=0; i<(n-1); i++) {
      for(int j=i+1; j<n; j++) {

        particle_t*p = &particles[j];
        compute_accelerations(&particles[i], p->pos, p->acc, p->mass, rsmoothsqr);
                }
            }
    for(int i=0; i<n; i++) {
        move_particle(&particles[i], dt);
     }
  }

int main(int argc, char const *argv[]) {
  /* output to txt file*/
  FILE* fp;
  fp = fopen("output/particles.txt","w+");
  fprintf(fp,"particle time x y z \n");

  /* initialise sim parameters */
  int n_iterations=0;
  //double *dtvar = malloc(n*sizeof(double));
  double t_max, dt, dt_out, rsmoothsqr;
  t_max = pow(4*M_PI*M_PI, 0.5);
  dt_out = 1*24*3600;
  rsmoothsqr = pow(0.01,2); //

   /* leapfrog method requires acceleration be evaluated first */
    for(int i=0; i<(n-1); i++) {
      for(int j=i+1; j<n; j++) {
        particle_t*p = &particles[j];
        /* compute the acceleration of particle j on particle i */
        compute_accelerations(&particles[i], p->pos, p->acc, p->mass, rsmoothsqr);
      }
    }

  dt = 0.01;
  double t=0.0;
  // run sim
  while (t < t_max) {
    n_iterations += 1;
    printf("Iteration: %i\n", n_iterations);
    iteration(dt, rsmoothsqr);

    for(int i=0; i<n; i++) {
          fprintf(fp,"%i %e %e %e %e \n", i, t,
          particles[i].pos[0], particles[i].pos[1], particles[i].pos[2]);
     }

    t += dt;
  }
  fclose(fp);

  return 0;
}
