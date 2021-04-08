#define ICS_H
#define n 2 /* number of particles to simulate */

struct particle;

/* define particle structure */
typedef struct particle{
  double pos[3];		/* position of the particle */
  double vel[3];		/* velocity of the particle */
  double acc[3];	  /* acceleration this particle experiences */
  double mass;			/* mass of the particle */
}particle_t;

/* define constants */
const double G = 1.0;
const double msun = 1.0, mearth = 3.003314e-06; // mass Sun, mass Earth [Msol]

/* array of particle structure, one element per particle */
/* velocities and position defined in natural units, fine if time takes units of
T = [4*pi^2*AU^3/M]*/
particle_t particles[n] = {{{0, 0, 0},
                            {0, 0, 0},
                            {0, 0, 0}, msun},

                           {{1, 0, 0},
                            {0, 1, 0},
                            {0, 0, 0}, mearth}};
