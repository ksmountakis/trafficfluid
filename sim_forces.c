
#define WALL_DN(point, safety, v) U_lemma3(sim, -(safety)*v, MAX(0, sim->y[i][t]-0.5*sim->w[i] - (point)), -sim->uymax_hard)
#define WALL_UP(point, safety, v) U_lemma3(sim, +(safety)*v, MAX(0, (point) - (sim->y[i][t]+0.5*sim->w[i])), -sim->uymax_hard)

static double
circular_dx(sim_t *sim, double dx) {
	double len = sim->roadlen_meters;
	if (fabs(dx) <= 0.5*len)
		return dx;
	else
		return (dx >= 0)?(dx-len):(dx+len);
}

static double
ftsx(sim_t *sim, int i, int t) {
	//return smooth_max(sim->ftsx_lo, sim->ftsx_hi*(1 - exp(-sim->ftsx_zeta*(sim->vdx_effective[i]-sim->vx[i][t]))), 1.0);
	return erfc(sim->vx[i][t]-sim->vd[i]) - 1;
}

static double
ftsy(sim_t *sim, int i, int t) {
	//return sim->ftsy_hi*(2 / (1 + exp(-sim->ftsy_zeta*(sim->vdy_effective[i] - sim->vy[i][t]))) - 1);
	return erfc(sim->vy[i][t]) - 1;
}

static double
fmdl(sim_t *sim, int i, int t) {
	double vd_range = (sim->vd_meters_per_sec_hi - sim->vd_meters_per_sec_lo);
	double point_in_range = 0.5;
	double spread_factor = 1.0;
	double gravity_point;
	if (sim->virtual_lanes) {
		point_in_range = (sim->vd[i] - sim->vd_meters_per_sec_lo)/vd_range;
	}
	gravity_point = sim->w[i]/2 + (0 + point_in_range)*(sim->roadwid_meters - sim->w[i])*spread_factor;
	return sim->coeff_fmdl*(2 / (1 + exp(-sim->ftsy_zeta*(gravity_point - sim->y[i][t]))) - 1);
}

static double
U_lemma3(sim_t *sim, double v, double d, double ubar) {
	double nom, den; 
	double T = sim->T;
	double ret;

	nom = T*ubar- 2*v + sqrt(pow(T,2)*pow(ubar,2)-(4*T*ubar*v)+8*ubar*(T*v-d));
	den = 2*T;

	if (fpclassify((ret = nom/den)) == FP_NAN) {
		ret = ubar;
	}

	return ret;
}

static double
G(double t, double warmup, double length, double cooldown) {
	return fmax(0, fmin(fmin(1. , t/warmup), 1. - (t-warmup-length)/cooldown));
}

static double
H(double t, double warmup, double length, double cooldown) {
	return G(t+warmup+length, warmup, 2*length, cooldown);
}

static long double
fca_mag(sim_t *sim, double dx, double dy, double vxi, double approaching_vx, double approaching_vy, double len, double wid)
{
	switch(sim->fca_method) 
	{
		case FCA_0:
		{
			double magx = pow(dx/len, 2);
			double magy = pow(dy/wid, 4);
			return (1.0/(1.0 + magx + magy));
		}

		case FCA_1:
		{
			double Lx = vxi*sim->time_gap_x;
			double Dx = 0.5*pow(approaching_vx, 2)/sim->coeff_fcax + len;
			double magx = H(dx, 0.5*Lx, 0.5*Lx + Dx + 2.5, 0.5*Lx);

			double Ly = approaching_vy*sim->time_gap_y + 0.5;
			double Dy = wid;
			double magy = pow(H(dy, Ly, Dy + 0.5, Ly), 1);

			//return (fmin(magx, magy));
			return pow(magx * magy, 1);
		}

		default:
			assert(!"unsupported fca_mag method");
	}
	return 0;
}

static double
fca(sim_t *sim, int i, int j, double dx_i_to_j, double dy_i_to_j, int t, double *fcax, double *fcay) {
	double mag, ang;
	
	double len = 0.5*(sim->l[i]+sim->l[j]);
	double wid = 0.5*(sim->w[i]+sim->w[j]);

	*fcax = *fcay = 0;

	double x, x0, y, y0, approaching_speed_x, approaching_speed_y;
	x0 = sim->x[j][t];	
	y0 = sim->y[j][t];
	x = sim->x[i][t]; 
	y = sim->y[i][t];

	if (dx_i_to_j >= 0) {
		// i is approaching j from behind
		// mag depends on where i is located within j's aura
		approaching_speed_x = MAX(0, sim->vx[i][t] - (1. - sim->safety_level_x)*sim->vx[j][t]);
		approaching_speed_y = fabs(sim->vy[i][t] - sim->vy[j][t]);
	} else {
		// j is approaching i from behind
		// mag depends on where i is located within j's aura
		approaching_speed_x = sim->vx[j][t];
		approaching_speed_y = fabs(sim->vy[j][t] - sim->vy[j][t]);
	}

	mag = fca_mag(sim,
				circular_dx(sim, x-x0), y-y0, 
				sim->vx[i][t],
				approaching_speed_x, approaching_speed_y, len, wid);
	
	ang = atan2(dy_i_to_j, dx_i_to_j);
	*fcax = -cos(ang) * mag;
	*fcay = -sin(ang) * mag;

	return mag;
}

static int
crash(sim_t *sim, int i, int j, double dx, double dy, int t) {
	if (i == j)
		return 0;

	if (fabs(dx) < 0.5*(sim->l[i]+sim->l[j]) && fabs(dy) < 0.5*(sim->w[i]+sim->w[j])) {
		if (!sim->warmup) {
			sim->crashes++;
			if (t > 0 && ((sim->walls[i][t].x.leader == j) || (sim->walls[j][t].x.leader == i)))
				sim->crashes_on_leaders++;
		}	
		return 1;
	}

	return 0;
}

static void
walls_init(sim_t *sim, int i, int t, walls_t *walls) {
	walls->x.position = fmod(sim->x[i][t] + 0.5*sim->roadlen_meters , sim->roadlen_meters);
	walls->x.distance = 0.5*sim->roadlen_meters;
	walls->x.velocity = sim->vd[i];
	walls->x.leader = -1;
	walls->x.bound = (sim->vd[i]-sim->vx[i][t])/sim->T;
	walls->x.degree = 0;

	walls->y.dn = 0;
	walls->y.up = sim->roadwid_meters;

	walls->y.up_j = -1;
	walls->y.dn_j = -1;

	walls->y.bound_dn = WALL_DN(0, sim->safety_level_y, sim->vy[i][t]);
	walls->y.bound_up = WALL_UP(sim->roadwid_meters, sim->safety_level_y, sim->vy[i][t]);

	walls->y.degree_up = 0;
	walls->y.degree_dn = 0;
}

static void
walls_update(sim_t *sim, int i, int j, int t, double fcax, double fcay, walls_t *walls) {
	if (j == i) 
		return;
	
	// clip forces within [-1, 1]
	fcax = MIN(1.0, MAX(-1.0, fcax));
 
	// degree_x in [0, 1]: degree to which j is a front obstacle
  	double degree_x = pow(fabs(MIN(fcax, 0)), 0.25);
 
	// determine front bound
 	double dx = circular_dx(sim, sim->x[j][t]-0.5*sim->l[j] - (sim->x[i][t]+0.5*sim->l[i]));
	double vs = MIN(dx/sim->frame_timegap, sim->vx[j][t] - 0.1);
	double uxbound = ((1.-degree_x)*(sim->vd[i]-sim->vx[i][t]) + degree_x*(vs - sim->vx[i][t]))/(1*sim->T);

	if (degree_x > walls->x.degree) {
		walls->x.bound = uxbound;
		walls->x.leader = j;
		walls->x.distance = dx;
		walls->x.position = sim->x[j][t]-0.5*sim->l[j];
		walls->x.velocity = sim->vx[j][t];
		walls->x.degree = degree_x;
	}

#if (PRINT_XY_DEGREES == 1)
	sim->degree_x[i][t][j] = degree_x;
	sim->degree_y[i][t][j] = fabs(fcay);
#endif
	return;
}

typedef struct {
	int j;
	double mag;
	double fcax, fcay;
} nbor_t;

int cmpnbors(const void *n1ptr, const void *n2ptr) {
	nbor_t *n1 = (nbor_t*)n1ptr;
	nbor_t *n2 = (nbor_t*)n2ptr;

	if (n1->mag == n2->mag)
		return 0;
	if (n1->mag < n2->mag)
		return -1;
	
	return +1;
}

static void
determine_forces(sim_t *sim, int i, int t) {
	int j;
	double fcax = 0, fcay = 0;
	nbor_t nbors[sim->n];
	nbor_t nbors_nudge[sim->n];

	memset(nbors, 0, sizeof(nbors));
	memset(nbors_nudge, 0, sizeof(nbors_nudge));

	/* initalizate desired speed */
	sim->vdx_effective[i] = sim->vd[i];
	sim->vdy_effective[i] = 0;
	/* initialize walls */
	walls_init(sim, i, t, &(sim->walls[i][t]));
	/* initalize forces to zero */
 	sim->fx[i][t] = sim->fy[i][t] = 0;

	/* obstacle-related forces */
	double fcax_sum = 0, fcay_sum =0;
	int nbors_added = 0, nbors_added_nudge = 0;
	
	for (j=0; j < sim->n; j++) {
		double dx_i_to_j = circular_dx(sim, sim->x[j][t]-sim->x[i][t]);
		double dy_i_to_j = sim->y[j][t] - sim->y[i][t];
		double fcamag = 0;

		if (j == i)
			continue;

		if (crash(sim, i, j, dx_i_to_j, dy_i_to_j, t))
			sim->crash[i][t] = 1;
		
		if (fabs(dx_i_to_j > sim->influence_radius_meters))
			continue;
		
		fcamag = fca(sim, i, j, dx_i_to_j, dy_i_to_j, t, &fcax, &fcay);
		
#if(PRINT_XY_DEGREES == 2)
		sim->degree_x[i][t][j] = fcax;
		sim->degree_y[i][t][j] = fcay;
#endif

		if (fcax < 0) {
			nbors[nbors_added].j = j;
			nbors[nbors_added].mag = -fcamag;
			nbors[nbors_added].fcax = fcax;
			nbors[nbors_added].fcay = fcay;
			nbors_added++;
		}
		
		if (fcax > 0) {
			nbors_nudge[nbors_added_nudge].j = j;
			nbors_nudge[nbors_added_nudge].mag = -fcamag;
			nbors_nudge[nbors_added_nudge].fcax = fcax;
			nbors_nudge[nbors_added_nudge].fcay = fcay;
			nbors_added_nudge++;
		}
	}

	// sort neighbors and compute fcax_sum from top sim->nbors
	if (sim->fca_nbors) {
		int x;
		qsort(nbors, nbors_added, sizeof(nbor_t), cmpnbors);
		fcax_sum = fcay_sum = 0;
		for (x=0; x < MIN(sim->fca_nbors, nbors_added); x++) {
			fcax_sum += nbors[x].fcax;
			fcay_sum += nbors[x].fcay;

			walls_update(sim, i, nbors[x].j, t, nbors[x].fcax, nbors[x].fcay, &(sim->walls[i][t]));
		}
	}

	if (sim->fca_nbors_nudge) {
		int x;
		qsort(nbors_nudge, nbors_added_nudge, sizeof(nbor_t), cmpnbors);
 		for (x=0; x < MIN(sim->fca_nbors_nudge, nbors_added_nudge); x++) {
			fcax_sum += nbors_nudge[x].fcax * sim->fwd_force_max_x;
			fcay_sum += nbors_nudge[x].fcay * sim->fwd_force_max_y;
		}
	}
	
	if (sim->fca_sat) {
		fcax_sum = MAX(-1, MIN(1, fcax_sum));
		fcay_sum = MAX(-1, MIN(1, fcay_sum));
	}

	sim->fx[i][t] += sim->coeff_fcax * fcax_sum;
	sim->fy[i][t] += sim->coeff_fcay * fcay_sum;

	/* target-speed related forces */
	if (!sim->ftsx_disabled)
		sim->fx[i][t] += ftsx(sim, i, t);
	sim->fy[i][t] += ftsy(sim, i, t);
	sim->fy[i][t] += fmdl(sim, i, t);

	sim->wall_y_up[i][t] = sim->walls[i][t].y.up;
	sim->wall_y_dn[i][t] = sim->walls[i][t].y.dn;
	sim->wall_y_up_j[i][t] = sim->walls[i][t].y.up_j;
	sim->wall_y_dn_j[i][t] = sim->walls[i][t].y.dn_j;
}

static void
regulate_forces(sim_t *sim, int i, int t) {
	/* y wall */
	if (sim->dynamic_y_walls && i != 0) {
		double uy_max = sim->walls[i][t].y.bound_up;
		double uy_min = sim->walls[i][t].y.bound_dn;
		if (uy_max >= -uy_min) {
			sim->fy[i][t] = MAX(sim->fy[i][t], -uy_min);
			sim->fy[i][t] = MIN(sim->fy[i][t], +uy_max);
		}
	}

	/* keeping vehicles within the road */
	sim->fy[i][t] = MAX(sim->fy[i][t], -WALL_DN(0, 1.05, sim->vy[i][t]));
	sim->fy[i][t] = MIN(sim->fy[i][t], +WALL_UP(sim->roadwid_meters, 1.05, sim->vy[i][t]));

	/* x wall */
	if (sim->dynamic_x_walls && sim->walls[i][t].x.leader != -1) {
		if (sim->walls[i][t].x.bound < sim->fx[i][t]) {
			sim->fx[i][t] = sim->walls[i][t].x.bound;
			sim->reg[t]++;
		}
	}

	/* [umin, umax] ranges */
	sim->fx[i][t] = MIN(sim->fx[i][t], sim->uxmax_hard);
	sim->fx[i][t] = MAX(sim->fx[i][t], sim->uxmin_hard);

	sim->fy[i][t] = MIN(sim->fy[i][t], +sim->uymax_hard);
	sim->fy[i][t] = MAX(sim->fy[i][t], -sim->uymax_hard);

	/* non-negative speed */
	sim->fx[i][t] = MAX(sim->fx[i][t], -sim->vx[i][t]/sim->T);

	/* non-excessive speed */
	sim->fx[i][t] = MIN(sim->fx[i][t], (1.2*sim->vd[i] - sim->vx[i][t])/sim->T);

	/* lat. speed never exceeds 0.5 x lon. speed, i.e.
	 * enforcing: |vy[k+1]| <= 0.5 vx[k+1]
	 */
	sim->fy[i][t] = MIN(sim->fy[i][t], +(0.5*sim->vx[i][t] - sim->vy[i][t])/sim->T);
	sim->fy[i][t] = MAX(sim->fy[i][t], -(0.5*sim->vx[i][t] + sim->vy[i][t])/sim->T);

	/* single-lane exceptions */
	if (sim->single_lane)
		sim->fy[i][t] = 0;
}
