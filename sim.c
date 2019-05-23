#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <omp.h>

#define PRINT_XY_DEGREES 0
#define PARALLEL_RUN 1


#define BUS_RATIO 10
double LANEWIDTH;
#define VD_VARIANCE 2.0
#define XBUFFER 2.0

#define SAMPLE_UNIFORM(min, max) ((double)min + ((double)random()/RAND_MAX)*(max - min))
#define MAX(a, b) (((a) > (b))?(a):(b))
#define MIN(a, b) (((a) <= (b))?(a):(b))

typedef enum {FCA_0, FCA_1, FCA_2, FCA_3} fca_method_t;

typedef	struct {
	struct {
		double position;
		double distance;
		double velocity;
		long double bound;
		int leader;
		double degree;
	}x;
	struct {
		double up, degree_up;
		double dn, degree_dn;
		int up_j, dn_j;
		long double bound_up, bound_dn;
	}y;
}walls_t;

typedef struct {
	double alpha;
	double frame_timegap;
	int ftsx_disabled;
	int fca_nbors, fca_nbors_nudge;
	int fca_sat;
	double duet, duet_vd, duet_vd_slow, duet_dy, duet_v0;

	int *reg;
	int fca_max, four_lane_init;
	int barrier;
	double single_lane_space, single_lane_front_vd;
	int single_lane;
	int dynamic_y_walls, dynamic_x_walls;
	double vd_throttling_horizon;
	int warmup, vd_throttling;
	double influence_radius_meters;
	double fwd_force_max_x, fwd_force_max_y;
	double bwd_force_max_x, bwd_force_max_y;

	fca_method_t fca_method;
	double time_gap_x, time_gap_y;
	int zero_controls, zero_initial_speed, lanedrop, virtual_lanes;
	int n, K, uxmax_hard, uxmin_hard, uymax_hard;
	double T, vd_meters_per_sec_lo, vd_meters_per_sec_hi;
	double roadlen_meters, roadwid_meters;
	double ftsx_zeta, ftsx_hi, ftsx_lo;
	double ftsy_zeta, ftsy_hi;
	int *truck;
	double **x, **y, **vx, **vy, **ux, **uy, *w, *l, *vd;
	double **wall_y_up, **wall_y_dn, **wall_y_up_j, **wall_y_dn_j;
	walls_t**walls;
	double **fx, **fy;
	int **crash;
	double *flow;
	double *vdx_effective, *vdy_effective;
	double coeff_fcax, coeff_fcay, fcax_zeta, fcay_zeta;
	double coeff_fmdl;
	double safety_level_x, safety_level_y;

	double ***degree_x;
	double ***degree_y;
	int crashes, crashes_on_leaders;
}sim_t;


#include "sim_forces.c"

static void
determine_controls(sim_t *sim, int i, int t) {
	if (sim->zero_controls || (sim->lanedrop && i == 0)) {
		sim->ux[i][t] = sim->uy[i][t] = 0;
		return;
	}
	sim->ux[i][t] = MIN(sim->fx[i][t], sim->uxmax_hard);
	sim->ux[i][t] = MAX(sim->fx[i][t], sim->uxmin_hard);
	sim->uy[i][t] = (sim->fy[i][t] >= 0)? MIN(sim->fy[i][t], sim->uymax_hard): MAX(sim->fy[i][t], -sim->uymax_hard);
}

static void
sim_configure(sim_t *sim) {
	char buf[1024];
	int input_int;
	double input_dbl;
	while (fgets(buf, sizeof(buf), stdin)) {
		if (sscanf(buf, "alpha:%lf", &input_dbl) == 1)
			sim->alpha = input_dbl;
		if (sscanf(buf, "frame_timegap:%lf", &input_dbl) == 1)
			sim->frame_timegap = input_dbl;
		/* duet configuration */
		if (sscanf(buf, "duet:%lf", &input_dbl) == 1)
			sim->duet = input_dbl;
		if (sscanf(buf, "duet_vd:%lf", &input_dbl) == 1)
			sim->duet_vd = input_dbl;
		if (sscanf(buf, "duet_v0:%lf", &input_dbl) == 1)
			sim->duet_v0 = input_dbl;
		if (sscanf(buf, "duet_vd_slow:%lf", &input_dbl) == 1)
			sim->duet_vd_slow = input_dbl;
		if (sscanf(buf, "duet_dy:%lf", &input_dbl) == 1)
			sim->duet_dy = input_dbl;

		if (sscanf(buf, "four_lane_init:%d", &input_int) == 1)
			sim->four_lane_init = input_int;
		if (sscanf(buf, "ftsx_disabled:%d", &input_int) == 1)
			sim->ftsx_disabled = input_int;
		if (sscanf(buf, "fca_nbors_nudge:%d", &input_int) == 1)
			sim->fca_nbors_nudge = input_int;
		if (sscanf(buf, "fca_nbors:%d", &input_int) == 1)
			sim->fca_nbors = input_int;
		if (sscanf(buf, "fca_sat:%d", &input_int) == 1)
			sim->fca_sat = input_int;
		if (sscanf(buf, "fca_max:%d", &input_int) == 1)
			sim->fca_max = input_int;
		if (sscanf(buf, "barrier:%d", &input_int) == 1)
			sim->barrier = input_int;
		if (sscanf(buf, "single_lane:%d", &input_int) == 1)
			sim->single_lane = input_int;
		if (sscanf(buf, "single_lane_space:%lf", &input_dbl) == 1)
			sim->single_lane_space = input_dbl;
		if (sscanf(buf, "single_lane_front_vd:%lf", &input_dbl) == 1)
			sim->single_lane_front_vd = input_dbl;
		if (sscanf(buf, "dynamic_x_walls:%d", &input_int) == 1)
			sim->dynamic_x_walls = input_int;
			
		if (sscanf(buf, "dynamic_y_walls:%d", &input_int) == 1)
			sim->dynamic_y_walls = input_int;

		if (sscanf(buf, "fca_method:%d", &input_int) == 1)
			sim->fca_method = input_int;
		if (sscanf(buf, "lanedrop:%d", &input_int) == 1)
			sim->lanedrop = input_int;
		if (sscanf(buf, "zero_initial_speed:%d", &input_int) == 1)
			sim->zero_initial_speed = input_int;
		if (sscanf(buf, "zero_controls:%d", &input_int) == 1)
			sim->zero_controls = input_int;
		if (sscanf(buf, "virtual_lanes:%d", &input_int) == 1)
			sim->virtual_lanes = input_int;

		if (sscanf(buf, "n:%d", &input_int) == 1)
			sim->n = input_int;
		if (sscanf(buf, "K:%d", &input_int) == 1)
			sim->K = input_int;
		if (sscanf(buf, "T:%lf", &input_dbl) == 1)
			sim->T = input_dbl;

		if (sscanf(buf, "influence_radius_meters:%lf", &input_dbl) == 1)
			sim->influence_radius_meters= input_dbl;

		if (sscanf(buf, "time_gap_x:%lf", &input_dbl) == 1)
			sim->time_gap_x= input_dbl;
		if (sscanf(buf, "time_gap_y:%lf", &input_dbl) == 1)
			sim->time_gap_y= input_dbl;

		if (sscanf(buf, "coeff_fcax:%lf", &input_dbl) == 1)
			sim->coeff_fcax = input_dbl;
		if (sscanf(buf, "coeff_fcay:%lf", &input_dbl) == 1)
			sim->coeff_fcay = input_dbl;

		if (sscanf(buf, "fcax_zeta:%lf", &input_dbl) == 1)
			sim->fcax_zeta = input_dbl;
		if (sscanf(buf, "fcay_zeta:%lf", &input_dbl) == 1)
			sim->fcay_zeta = input_dbl;

		if (sscanf(buf, "safety_level_x:%lf", &input_dbl) == 1)
			sim->safety_level_x = input_dbl;

		if (sscanf(buf, "safety_level_y:%lf", &input_dbl) == 1)
			sim->safety_level_y = input_dbl;

		if (sscanf(buf, "fwd_force_max_x:%lf", &input_dbl) == 1)
			sim->fwd_force_max_x = input_dbl;
		if (sscanf(buf, "fwd_force_max_y:%lf", &input_dbl) == 1)
			sim->fwd_force_max_y = input_dbl;

		if (sscanf(buf, "bwd_force_max_x:%lf", &input_dbl) == 1)
			sim->bwd_force_max_x = input_dbl;
		if (sscanf(buf, "bwd_force_max_y:%lf", &input_dbl) == 1)
			sim->bwd_force_max_y = input_dbl;

		if (sscanf(buf, "coeff_fmdl:%lf", &input_dbl) == 1)
			sim->coeff_fmdl = input_dbl;

		if (sscanf(buf, "vd_meters_per_sec_lo:%lf", &input_dbl) == 1)
			sim->vd_meters_per_sec_lo = input_dbl;
		if (sscanf(buf, "vd_meters_per_sec_hi:%lf", &input_dbl) == 1)
			sim->vd_meters_per_sec_hi = input_dbl;

		if (sscanf(buf, "uxmax_hard:%lf", &input_dbl) == 1)
			sim->uxmax_hard = input_dbl;
		if (sscanf(buf, "uxmin_hard:%lf", &input_dbl) == 1)
			sim->uxmin_hard = input_dbl;
		if (sscanf(buf, "uymax_hard:%lf", &input_dbl) == 1)
			sim->uymax_hard = input_dbl;

		if (sscanf(buf, "roadlen_meters:%lf", &input_dbl) == 1)
			sim->roadlen_meters = input_dbl;
		if (sscanf(buf, "roadwid_meters:%lf", &input_dbl) == 1)
			sim->roadwid_meters = input_dbl;

		if (sscanf(buf, "ftsx_zeta:%lf", &input_dbl) == 1)
			sim->ftsx_zeta = input_dbl;
		if (sscanf(buf, "ftsx_hi:%lf", &input_dbl) == 1)
			sim->ftsx_hi = input_dbl;
		if (sscanf(buf, "ftsx_lo:%lf", &input_dbl) == 1)
			sim->ftsx_lo = input_dbl;

		if (sscanf(buf, "ftsy_zeta:%lf", &input_dbl) == 1)
			sim->ftsy_zeta = input_dbl;
		if (sscanf(buf, "ftsy_hi:%lf", &input_dbl) == 1)
			sim->ftsy_hi = input_dbl;
	}

	if (!sim->alpha || !sim->n || !sim->K  || !sim->T
			||(sim->vd_throttling && !sim->vd_throttling_horizon)
			||!sim->coeff_fcax ||!sim->coeff_fcay ||!sim->fcax_zeta ||!sim->fcay_zeta ||!sim->safety_level_x || !sim->safety_level_y
			||!sim->vd_meters_per_sec_lo ||!sim->vd_meters_per_sec_hi
			||!sim->roadlen_meters ||!sim->roadwid_meters 
			||!sim->uxmax_hard ||!sim->uxmin_hard ||!sim->uymax_hard
			||!sim->ftsx_zeta || !sim->ftsx_hi || !sim->ftsx_lo 
			||!sim->ftsy_zeta || !sim->ftsy_hi) {
		fprintf(stderr, "incomplete configuration\n");
		exit(1);
	}
	

	LANEWIDTH = 3.4;
	if (sim->four_lane_init)
		LANEWIDTH = 2.6;

	if (sim->duet)
		sim->n = 2;

	sim->fca_nbors = MIN(sim->fca_nbors, sim->n);
}

static void sim_run(sim_t *sim, int K);

static void vehicle_lw(sim_t *sim, int i, int isTruck){

	struct {
		double l, aspect_ratio;
	} types [] = {
		{3.20, 2.0},
		{3.90, 2.3},
		{4.25, 2.4},
		{4.55, 2.5},
		{4.60, 2.6},
		{5.15, 2.8},
	};
	int draw = (int)SAMPLE_UNIFORM(0, 5);

	sim->l[i] = types[draw].l;
	sim->w[i] = sim->l[i] / types[draw].aspect_ratio;

	if (isTruck == 1){
		sim->l[i] = SAMPLE_UNIFORM(12, 17);
		sim->w[i] = 2.5;
	}
}

/* initalize dimensions */
static void
sim_initialize_lw(sim_t *sim)
{
	int i;

	for (i=0; i < sim->n; i++) {
		if (sim->lanedrop && i == 0) {
			sim->l[i] = 0.25 * sim->roadlen_meters;
			sim->w[i] = 3.0;
		} else if (sim->barrier && i == sim->n-1) {
			sim->l[i] = 2.5;
			sim->w[i] = sim->roadwid_meters;
		} else {
			if (SAMPLE_UNIFORM(0, 100) < BUS_RATIO)
				sim->truck[i] = 1;
			
			vehicle_lw(sim, i, sim->truck[i]);
		}
	}
}

/* initialize positions */
static void
sim_initialize_xy(sim_t *sim)
{
	const double lanewidth = LANEWIDTH;
	//const int numlanes = MIN(3, round(sim->roadwid_meters / lanewidth));
	const int numlanes = round(sim->roadwid_meters / lanewidth);

	int i, l, j;
	int *lastcar = calloc(sizeof(int), numlanes);
	int *numcars = calloc(sizeof(int), numlanes);

	for (l=0; l < numlanes; l++)
		lastcar[l] = -1;
 	
	/* well-structured placement */
	for (i=0, l=0; i < sim->n; i++, l++) {
		if (l == numlanes)
			l = 0;
		j = lastcar[l];
		sim->x[i][0] = XBUFFER + 0.5*sim->l[i] + ((j == -1)?0:sim->x[j][0]+0.5*sim->l[j]); 
		sim->y[i][0] = l*lanewidth + 0.5 + 0.5*sim->w[i];
 		lastcar[l] = i;
		numcars[l]++;

		for (j=0; j < i; j++) {
			if (i == j)
				continue;
			if (fabs(circular_dx(sim, sim->x[i][0]-sim->x[j][0])) < 0.5*(sim->l[i]+sim->l[j]) 
				&& fabs(sim->y[i][0] - sim->y[j][0]) < 0.5*(sim->w[i]+sim->w[j])) {
				fprintf(stderr, "Overlap while placing %d\n", i);
				exit(1);
			}
		}
	}

	/* determine random distribution of the slack on each lane */
	double **slack = calloc(sizeof(double*), numlanes);
	int v;
	for (l=0; l < numlanes; l++) {
		if ((j = lastcar[l]) == -1)
			continue;
		slack[l] = calloc(sizeof(double), numcars[l]);
		double total = 0;
		for (v=0; v < numcars[l]; v++)
			total += (slack[l][v] = (double)rand()/RAND_MAX);

		for (v=0; v < numcars[l]; v++)
			slack[l][v] /= total;

		total = sim->roadlen_meters - (sim->x[lastcar[l]][0]+0.5*sim->l[lastcar[l]]);
		double s = 0;
		for (v=0; v < numcars[l]; v++) {
			s += slack[l][v];
			slack[l][v] = s * total;
		}
	}

	/* distribute slack over cars on each lane */
	double **slackptr = calloc(sizeof(double*), numlanes);
	for (l=0; l < numlanes; l++)
			slackptr[l] = &slack[l][0];

	for (j=0; j < sim->n; j++) {
		l = j % numlanes;
		sim->x[j][0] += *(slackptr[l]++);
		sim->y[j][0] += 0.5*(double)rand()/RAND_MAX;

	}

	/* convert trucks that are positioned in the 3/4 of the road to cars */
	for (i = 0; i < sim->n; i++){
		if (sim->truck[i] == 1 && sim->y[i][0] < (1.0/3.0)*sim->roadwid_meters){
			vehicle_lw(sim, i, 0);
		}
	}
	

	for (l=0; l < numlanes; l++)
		free(slack[l]);

	free(slackptr);
	free(slack);
	free(lastcar);
	free(numcars);
}

static void
sim_initialize_xy_duet(sim_t *sim)
{
	sim->x[0][0] = 0.5*sim->l[0];
	sim->x[1][0] = sim->l[0] + 0.5*sim->l[1] + sim->duet;

	sim->y[0][0] = sim->roadwid_meters/2 + sim->duet_dy;
	sim->y[1][0] = sim->roadwid_meters/2;
}

/* initalize desired speeds */
static void
sim_initialize_vd(sim_t *sim)
{
	int i;

	for (i=0; i < sim->n; i++) {
		if (sim->lanedrop && i == 0) {
			sim->vx[i][0] = 0;
			sim->vd[i] = 0;
		} else {

			double lo = sim->vd_meters_per_sec_lo, hi = sim->vd_meters_per_sec_hi;
			sim->vd[i] = lo+ (1. - (sim->y[i][0]/sim->roadwid_meters))*(hi-lo);
			
			if (sim->zero_initial_speed)
				sim->vx[i][0] = 0;
			else
				sim->vx[i][0] = 0.5*(sim->vd_meters_per_sec_lo + sim->vd_meters_per_sec_hi);
		}
	}
}

static void
sim_initialize_vd_duet(sim_t *sim)
{
	sim->vx[0][0] = sim->duet_v0;	
	sim->vd[0] = sim->duet_vd;
	sim->vd[1] = sim->vx[1][0] = sim->duet_vd_slow;
}

static void
sim_initialize_xy_single(sim_t *sim)
{
	int i;
	double dx = 0;
	for (i=0; i < sim->n; i++) {
		sim->x[i][0] = 0.5*sim->l[i] + dx;
		dx = sim->x[i][0] + 0.5*sim->l[i] + sim->single_lane_space;

		if (dx > sim->roadlen_meters) {
			fprintf(stderr, "single-lane mode: unable to fit vehicle nr. %d\n", i);
			exit(1);
		}

		sim->y[i][0] = 0.5*sim->roadwid_meters;
	}
}

static void
sim_initialize_vd_single(sim_t *sim)
{
	sim_initialize_vd(sim);
	sim->vdx_effective[sim->n-1] = sim->vd[sim->n-1] = sim->single_lane_front_vd;
}

static void
sim_initialize(sim_t *sim) {
	int i;
	assert((sim->x = calloc(sizeof(double*), sim->n)));
	assert((sim->y = calloc(sizeof(double*), sim->n)));

	assert((sim->vx = calloc(sizeof(double*), sim->n)));
	assert((sim->vy = calloc(sizeof(double*), sim->n)));

	assert((sim->ux = calloc(sizeof(double*), sim->n)));
	assert((sim->uy = calloc(sizeof(double*), sim->n)));

	assert((sim->fx = calloc(sizeof(double*), sim->n)));
	assert((sim->fy = calloc(sizeof(double*), sim->n)));


	assert((sim->walls = calloc(sizeof(walls_t*), sim->n)));


	assert((sim->wall_y_up = calloc(sizeof(double*), sim->n)));
	assert((sim->wall_y_dn = calloc(sizeof(double*), sim->n)));
	assert((sim->wall_y_up_j = calloc(sizeof(double*), sim->n)));
	assert((sim->wall_y_dn_j = calloc(sizeof(double*), sim->n)));

	assert((sim->crash = calloc(sizeof(int*), sim->n)));
	assert((sim->flow = calloc(sizeof(double), sim->K)));
	assert((sim->reg = calloc(sizeof(int), sim->K)));

	assert((sim->l = calloc(sizeof(double), sim->n)));
	assert((sim->w = calloc(sizeof(double), sim->n)));
	assert((sim->vd = calloc(sizeof(double), sim->n)));
	assert((sim->truck = calloc(sizeof(int), sim->n)));
	assert((sim->vdx_effective = calloc(sizeof(double), sim->n)));
	assert((sim->vdy_effective = calloc(sizeof(double), sim->n)));

	assert((sim->degree_x = calloc(sizeof(double**), sim->n)));
	assert((sim->degree_y = calloc(sizeof(double**), sim->n)));

	for (i=0; i < sim->n; i++) {
		assert((sim->x[i] = calloc(sizeof(double), sim->K)));
		assert((sim->y[i] = calloc(sizeof(double), sim->K)));

		assert((sim->vx[i] = calloc(sizeof(double), sim->K)));
		assert((sim->vy[i] = calloc(sizeof(double), sim->K)));

		assert((sim->ux[i] = calloc(sizeof(double), sim->K)));
		assert((sim->uy[i] = calloc(sizeof(double), sim->K)));

		assert((sim->fx[i] = calloc(sizeof(double), sim->K)));
		assert((sim->fy[i] = calloc(sizeof(double), sim->K)));

		assert((sim->walls[i] = calloc(sizeof(walls_t), sim->K)));

		assert((sim->wall_y_up[i] = calloc(sizeof(double), sim->K)));
		assert((sim->wall_y_dn[i] = calloc(sizeof(double), sim->K)));
		assert((sim->wall_y_up_j[i] = calloc(sizeof(double), sim->K)));
		assert((sim->wall_y_dn_j[i] = calloc(sizeof(double), sim->K)));

		assert((sim->crash[i] = calloc(sizeof(int), sim->K)));

#if(PRINT_XY_DEGREES)
		assert((sim->degree_x[i] = calloc(sizeof(double*), sim->K)));
		assert((sim->degree_y[i] = calloc(sizeof(double*), sim->K)));
		int t;
		for (t=0; t < sim->K; t++) {
			assert((sim->degree_x[i][t] = calloc(sizeof(double), sim->n)));
			assert((sim->degree_y[i][t] = calloc(sizeof(double), sim->n)));
		}
#endif

	}

	sim_initialize_lw(sim);
	if (sim->duet > 0) {
		sim_initialize_xy_duet(sim);
		sim_initialize_vd_duet(sim);
	} else if (sim->single_lane) {
		sim_initialize_xy_single(sim);
		sim_initialize_vd_single(sim);
	} else {
		sim_initialize_xy(sim);
		sim_initialize_vd(sim);
	}

#if 0

	// crash resolution period
	do {
		int j, crash_found, crash_resolution_round = 0;
		double ratio_backup = sim->push_repel_ratio;

		sim->warmup = 1;
		sim->push_repel_ratio = 1.1;

 		do {
			crash_found = 0;
			for (i=0; i < sim->n; i++) {
				for (j=0; j < i; j++) {
					if (crash(sim, i, j, 0)) {
						crash_found = 1;
					}
				}
			}
 			sim_run(sim, 2);

			for (i=0; i < sim->n; i++) {
				sim->x[i][0] = sim->x[i][1];
				sim->y[i][0] = sim->y[i][1];
				sim->vx[i][0] = sim->vy[i][0] = sim->ux[i][0] = sim->uy[i][0] = 0;
			}

			crash_resolution_round++;
			fprintf(stderr, "crash resolution rounds: %d\r", crash_resolution_round);
 		} while(crash_found);
		fprintf(stderr, "\n");

		for (i=0; i < sim->n; i++)
			memset(sim->crash[i], 0, sizeof(int)*sim->K);

		sim->warmup = 0;
		sim->push_repel_ratio = ratio_backup;

	}while(0);
#endif
}

static void
sim_run(sim_t *sim, int K) {
	int t, i;
	double sensor_location = sim->roadlen_meters / 2;

 	sim->flow[0] = 0;
 	for (t=0; t < MIN(K, sim->K)-1; t++) {
 		sim->flow[t+1] = sim->flow[t];

		#if(PARALLEL_RUN)
		#pragma omp parallel for
		#endif
		for (i=0; i < sim->n; i++) 
			determine_forces(sim, i, t);
		
		#if(PARALLEL_RUN)
		#pragma omp parallel for
		#endif
		for (i=0; i < sim->n; i++) {
			regulate_forces(sim, i, t);
			determine_controls(sim, i, t);

			double T = sim->alpha * sim->T;
			sim->vx[i][t+1] = sim->vx[i][t] + sim->ux[i][t]*T; 
			sim->x[i][t+1] = sim->x[i][t] + sim->vx[i][t]*T + sim->ux[i][t]*0.5*pow(T, 2);
			sim->x[i][t+1] = fmod(sim->x[i][t+1], sim->roadlen_meters);

			// update flow
			if (sim->x[i][t] <= sensor_location && sim->x[i][t+1] > sensor_location) {
				#pragma omp critical
				sim->flow[t+1]++;
			}

			if (!sim->single_lane) {
				sim->vy[i][t+1] = sim->vy[i][t] + sim->uy[i][t]*T; 
				sim->y[i][t+1] = sim->y[i][t] + sim->vy[i][t]*T + sim->uy[i][t]*0.5*pow(T, 2);
			} else {
				sim->y[i][t+1] = sim->y[i][t];
				sim->vy[i][t+1] = 0;
			}
		}

		if (!sim->warmup)
			fprintf(stderr,"time-step %5d crashes %5d\r", t, sim->crashes);
 	}

	//for (t=0; t < MIN(K, sim->K)-1; t++) 
	//	sim->flow[t] *= 3600/(t * sim->T);
}

static void
sim_print(sim_t *sim) {
	int i, t, K = MIN(10000, sim->K);
	puts("sim = {");
	printf("n:%d,\n", sim->n);
	printf("K:%d,\n", K);
	printf("T:%.1lf,\n", sim->T);
	printf("roadlen_meters:%.1lf,\n", sim->roadlen_meters);
	printf("roadwid_meters:%.1lf,\n", sim->roadwid_meters);
	printf("lanedrop:%d,\n", sim->lanedrop);

	printf("l: [");
	for (i=0; i < sim->n; i++) 
		printf("%.1lf, ", sim->l[i]);
	puts("],");

	printf("w: [");
	for (i=0; i < sim->n; i++) 
		printf("%.1lf, ", sim->w[i]);
	puts("],");

	printf("vd: [");
	for (i=0; i < sim->n; i++) 
		printf("%.1lf, ", sim->vd[i]);
	puts("],");

	printf("flow: [");
	for (t=0; t < sim->K; t++) 
		printf("%.1lf, ", sim->flow[t]);
	puts("],");

	printf("frame: [");
	for (t=0; t < sim->K; t++) {
		double degree_max = 0, degree_avg = 0;
		int j, count = 0;
		for (j=0; j < sim->n; j++) {
			if (sim->walls[j][t].x.leader == -1) 
				continue;
			if (sim->duet > 0 && j != 0)
				continue;
			degree_max = MAX(degree_max, sim->walls[j][t].x.degree);
			degree_avg += sim->walls[j][t].x.degree;
			count++;
		}
		if (count > 0)
			degree_avg /= count;
		printf("{max: %.5lf, avg: %.3lf}, ", degree_max, degree_avg);
	}
	puts("],");

	printf("reg: [");
	for (t=0; t < sim->K; t++) {
		double reg = 100*((double)sim->reg[t] / (double)sim->n);
		printf("%-3.0f, ", reg);
	}
	puts("],");

	printf("x: [\n");
	for (i=0; i < sim->n; i++) {
		printf("\t[");
		for (t=0; t < K; t++)
			printf("%.1lf,", sim->x[i][t]);
		printf("],\n");
	}
	puts("   ],");

	printf("y: [\n");
	for (i=0; i < sim->n; i++) {
		printf("\t[");
		for (t=0; t < K; t++) {
			double y = sim->y[i][t], hw = sim->w[i]*0.5;
			double err = MAX(hw - y, y - sim->roadwid_meters - hw);
			if (err > 0.1)
				fprintf(stderr, "y[%d][%d] out of bounds: %lf\n", i, t, err);
			printf("%.1lf,", y);
		}
		printf("],\n");
	}
	puts("   ],");

	printf("vx: [\n");
	for (i=0; i < sim->n; i++) {
		printf("\t[");
		for (t=0; t < K; t++)
			printf("%.3lf,", sim->vx[i][t]);
		printf("],\n");
	}
	puts("   ],");

#if(PRINT_XY_DEGREES)
	int j;
	printf("degree: [\n");
	for (i=0; i < sim->n; i++) {
		printf("\t[");
		for (t=0; t < K; t++) {
			printf("[");
			for (j=0; j < sim->n; j++) {
				printf("{x:%.0lf, y:%.0lf},", 100*sim->degree_x[i][t][j], 100*sim->degree_y[i][t][j]);
			}
			printf("],");
		}
		printf("],\n");
	}
	puts("   ],");
#endif

#if 0
	printf("vy: [\n");
	for (i=0; i < sim->n; i++) {
		printf("\t[");
		for (t=0; t < K; t++)
			printf("%.1lf,", sim->vy[i][t]);
		printf("],\n");
	}
	puts("   ],");


	printf("fx: [\n");
	for (i=0; i < sim->n; i++) {
		printf("\t[");
		for (t=0; t < K; t++)
			printf("%.1lf,", sim->fx[i][t]);
		printf("],\n");
	}
	puts("   ],");

	printf("fy: [\n");
	for (i=0; i < sim->n; i++) {
		printf("\t[");
		for (t=0; t < K; t++)
			printf("%.1lf,", sim->fy[i][t]);
		printf("],\n");
	}
	puts("   ],");
#endif

	printf("ux: [\n");
	for (i=0; i < sim->n; i++) {
		printf("\t[");
		for (t=0; t < K; t++)
			printf("%.3lf,", sim->ux[i][t]);
		printf("],\n");
	}
	puts("   ],");

	printf("uy: [\n");
	for (i=0; i < sim->n; i++) {
		printf("\t[");
		for (t=0; t < K; t++)
			printf("%.1lf,", sim->uy[i][t]);
		printf("],\n");
	}
	puts("   ],");

	printf("wall_x: [\n");
	for (i=0; i < sim->n; i++) {
		printf("\t[");
		for (t=0; t < K; t++)
			printf("%.1lf,", sim->walls[i][t].x.position);
		printf("],\n");
	}
	puts("   ],");

	printf("wall_y_up: [\n");
	for (i=0; i < sim->n; i++) {
		printf("\t[");
		for (t=0; t < K; t++)
			printf("%.1lf,", sim->wall_y_up[i][t]);
		printf("],\n");
	}
	puts("   ],");

	printf("wall_y_dn: [\n");
	for (i=0; i < sim->n; i++) {
		printf("\t[");
		for (t=0; t < K; t++)
			printf("%.1lf,", sim->wall_y_dn[i][t]);
		printf("],\n");
	}
	puts("   ],");

	printf("leader: [\n");
	for (i=0; i < sim->n; i++) {
		printf("\t[");
		for (t=0; t < K; t++)
			printf("%3d,", sim->walls[i][t].x.leader);
		printf("],\n");
	}
	puts("   ],");

	printf("crash: [\n");
	for (i=0; i < sim->n; i++) {
		printf("\t[");
		for (t=0; t < K; t++)
			printf("%d,", sim->crash[i][t]);
		printf("],\n");
	}
	puts("   ],");
	
	puts("};");
}

int
main(int argc, char **argv) {
	sim_t sim = {0};

	if (argc > 2 && !strncmp(argv[1], "-seed", 5))
		srandom(atoi(argv[2]));

	sim_configure(&sim);
	sim_initialize(&sim);
	sim_run(&sim, sim.K);
	sim_print(&sim);

	double T = sim.alpha * sim.T;
	int k0 = 0.75 * sim.K;
	double duration = fmax(T, (sim.K-2 - k0) * T);
	double flow = (sim.flow[sim.K-2] - sim.flow[k0]) * (3600.0/duration);

	fprintf(stderr, "n %d roadlen %.0lf coeff %.1lf %.1lf fwd %.1lf %.1lf bwd %.1lf %.1lf safety %.1lf %.1lf fcamax %d flow %.1lf crashes %d crashes-on-leaders %d\n", 
		sim.n, 
		sim.roadlen_meters,
		sim.coeff_fcax,
		sim.coeff_fcay,
		sim.fwd_force_max_x, sim.fwd_force_max_y,
		sim.bwd_force_max_x, sim.bwd_force_max_y,
		sim.safety_level_x, sim.safety_level_y,
		sim.fca_max,
		flow, sim.crashes, sim.crashes_on_leaders);

	return 0;
}
