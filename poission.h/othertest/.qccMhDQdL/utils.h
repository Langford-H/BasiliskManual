#ifndef BASILISK_HEADER_3
#define BASILISK_HEADER_3
#line 1 "/home/dahuanghhc/basilisk/src/utils.h"
/**
# Various utility functions

## Default parameters and variables

The default maximum timestep and CFL number. */

double DT = 1e10, CFL = 0.5;

/**
Performance statistics are stored in this structure. */

struct {
  // total number of (leaf) cells advanced for this process
  long nc;
  // total number of (leaf) cells advanced for all processes
  long tnc;
  // real time elapsed since the start
  double t;
  // average computational speed (leaves/sec)
  double speed;
  // global timer
  timer gt;
} perf;

/**
Performance statistics are gathered by this function, which is
typically called by the run() loop. */

void update_perf() {
  perf.nc += grid->n;
  perf.tnc += grid->tn;
  perf.t = timer_elapsed (perf.gt);
  perf.speed = perf.tnc/perf.t;
}

/**
## Timing

These functions can be used for profiling. */

typedef struct {
  double cpu;   // CPU time (sec)
  double real;  // Wall-clock time (sec)
  double speed; // Speed (points.steps/sec)
  double min;   // Minimum MPI time (sec)
  double avg;   // Average MPI time (sec)
  double max;   // Maximum MPI time (sec)
  size_t tnc;   // Number of grid points
  long   mem;   // Maximum resident memory (kB)
} timing;

/**
Given a timer, iteration count *i*, total number of cells *tnc* and
array of MPI timings *mpi* (with a size equal to the number of
processes), this function returns the statistics above. */

timing timer_timing (timer t, int i, size_t tnc, double * mpi)
{
  timing s;
@if _MPI
  s.avg = mpi_time - t.tm;
@endif
  clock_t end = clock();
  s.cpu = ((double) (end - t.c))/CLOCKS_PER_SEC;
  s.real = timer_elapsed (t);
  if (tnc == 0) {
    double n = 0;
    foreach(reduction(+:n)) n++;
    s.tnc = n;
    tnc = n*i;
  }
  else
    s.tnc = tnc;
@if _GNU_SOURCE
  struct rusage usage;
  getrusage (RUSAGE_SELF, &usage);
  s.mem = usage.ru_maxrss;
@else
  s.mem = 0;
@endif
@if _MPI
  if (mpi)
    MPI_Allgather (&s.avg, 1, MPI_DOUBLE, mpi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
  s.max = s.min = s.avg;
  mpi_all_reduce (s.max, MPI_DOUBLE, MPI_MAX);
  mpi_all_reduce (s.min, MPI_DOUBLE, MPI_MIN);
  mpi_all_reduce (s.avg, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.real, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.mem, MPI_LONG, MPI_SUM);
  s.real /= npe();
  s.avg /= npe();
  s.mem /= npe();
@else
  s.min = s.max = s.avg = 0.;
@endif
  s.speed = s.real > 0. ? tnc/s.real : -1.;
  return s;
}

/**
This function writes timing statistics on standard output. */

void timer_print (timer t, int i, size_t tnc)
{
  timing s = timer_timing (t, i, tnc, NULL);
  fprintf (fout,
	   "\n# " GRIDNAME 
	   ", %d steps, %g CPU, %.4g real, %.3g points.step/s, %d var\n",
	   i, s.cpu, s.real, s.speed, (int) (datasize/sizeof(double)));
@if _MPI
  fprintf (fout,
	   "# %d procs, MPI: min %.2g (%.2g%%) "
	   "avg %.2g (%.2g%%) max %.2g (%.2g%%)\n",
	   npe(),
	   s.min, 100.*s.min/s.real,
	   s.avg, 100.*s.avg/s.real,
	   s.max, 100.*s.max/s.real);
@endif
}

/**
## Simple field statistics 

The *normf()* function returns the (volume) average, RMS norm, max
norm and volume for field *f*. */

typedef struct {
  double avg, rms, max, volume;
} norm;

norm normf (scalar f)
{
  double avg = 0., rms = 0., max = 0., volume = 0.;
  foreach(reduction(max:max) reduction(+:avg) 
	  reduction(+:rms) reduction(+:volume)) 
    if (f[] != nodata && dv() > 0.) {
      double v = fabs(f[]);
      if (v > max) max = v;
      volume += dv();
      avg    += dv()*v;
      rms    += dv()*sq(v);
    }
  norm n;
  n.avg = volume ? avg/volume : 0.;
  n.rms = volume ? sqrt(rms/volume) : 0.;
  n.max = max;
  n.volume = volume;
  return n;
}

/**
The *statsf()* function returns the minimum, maximum, volume sum,
standard deviation and volume for field *f*. */

typedef struct {
  double min, max, sum, stddev, volume;
} stats;

stats statsf (scalar f)
{
  double min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
  foreach(reduction(+:sum) reduction(+:sum2) reduction(+:volume)
	  reduction(max:max) reduction(min:min)) 
    if (dv() > 0. && f[] != nodata) {
      volume += dv();
      sum    += dv()*f[];
      sum2   += dv()*sq(f[]);
      if (f[] > max) max = f[];
      if (f[] < min) min = f[];
    }
  stats s;
  s.min = min, s.max = max, s.sum = sum, s.volume = volume;
  if (volume > 0.)
    sum2 -= sum*sum/volume;
  s.stddev = sum2 > 0. ? sqrt(sum2/volume) : 0.;
  return s;
}

/**
## Slope limiters 

Given three values, these [slope
limiters](https://en.wikipedia.org/wiki/Flux_limiter#Limiter_functions)
return the corresponding slope-limited gradient. */

static double generic_limiter (double r, double beta)
{
  double v1 = min (r, beta), v2 = min (beta*r, 1.);
  v1 = max (0., v1);
  return max (v1, v2);
}

double minmod (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 1.)*(s1 - s0);
}

double superbee (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 2.)*(s1 - s0);
}

double sweby (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 1.5)*(s1 - s0);
}

/**
This is the [generalised minmod
limiter](https://en.wikipedia.org/wiki/Flux_limiter#Generalised_minmod_limiter).
The $\theta$ global variable can be used to tune the limiting
($\theta=1$ gives minmod, the most dissipative limiter and $\theta=2$
gives superbee, the least dissipative). */
	  
double theta = 1.3;

double minmod2 (double s0, double s1, double s2)
{
  if (s0 < s1 && s1 < s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 < d1) d1 = d2;
    return min(d1, d3);
  }
  if (s0 > s1 && s1 > s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 > d1) d1 = d2;
    return max(d1, d3);
  }
  return 0.;
}

/**
Given a list of scalar fields *f*, this function fills the gradient
fields *g* with the corresponding gradients. If the *gradient*
attribute of a field is set (typically to one of the limiting
functions above), it is used to compute the gradient, otherwise simple
centered differencing is used. */
	  
void gradients (scalar * f, vector * g)
{
  assert (list_len(f) == vectors_len(g));
  foreach() {
    scalar s; vector v;
    for (s,v in f,g) {
      if (s.gradient)
	foreach_dimension() {
#if EMBED
	  if (!fs.x[] || !fs.x[1])
	    v.x[] = 0.;
	  else
#endif
	    v.x[] = s.gradient (s[-1], s[], s[1])/Delta;
	}
      else // centered
	foreach_dimension() {
#if EMBED
	  if (!fs.x[] || !fs.x[1])
	    v.x[] = 0.;
	  else
#endif
	    v.x[] = (s[1] - s[-1])/(2.*Delta);
	}
    }
  }
}

/**
## Other functions

Given a velocity field $\mathbf{u}$, this function fills a scalar
field $\omega$ with the vorticity field
$$
\omega = \partial_x u_y - \partial_y u_x
$$ 
To properly take the metric into account, $\omega$ is computed as an
average over the (surface) element, which is easily obtained as
the circulation of $\mathbf{u}$ along the boundary of the
element. Using the definitions of the [metric
factors](README#general-orthogonal-coordinates) `fm` and `cm`, this
gives the expression below. */

void vorticity (const vector u, scalar omega)
{
  foreach()
    omega[] = ((fm.x[1] - fm.x[])*u.y[] +
	       fm.x[1]*u.y[1] - fm.x[]*u.y[-1] -
	       (fm.y[0,1] - fm.y[])*u.x[] +
	       fm.y[]*u.x[0,-1] - fm.y[0,1]*u.x[0,1])/(2.*cm[]*Delta + SEPS);
}

/**
Given two scalar fields *s* and *sn* this function returns the maximum
of their absolute difference. */

double change (scalar s, scalar sn)
{
  double max = 0.;
  foreach(reduction(max:max)) {
    if (dv() > 0.) {
      double ds = fabs (s[] - sn[]);
      if (ds > max)
	max = ds;
    }
    sn[] = s[];
  }
  return max;
}

/**
These functions return the scalar/vector fields called *name*, or
-1 if they don't exist. */

scalar lookup_field (const char * name)
{
  if (name)
    for (scalar s in all)
      if (!strcmp (s.name, name))
	return s;
  return (scalar){-1};
}

vector lookup_vector (const char * name)
{
  if (name) {
    char component[strlen(name) + 3];
    strcpy (component, name);
    strcat (component, ".x");
    for (scalar s in all)
      if (!strcmp (s.name, component))
	return s.v;
  }
  return (vector){{-1}};
}

/**
The function below traverses the set of sub-segments intersecting the
mesh and spanning the [A:B] segment. The pair of coordinates defining
the sub-segment contained in each cell are defined by `p[0]` and
`p[1]`. */

#if 1 // fixme: foreach_dimension() does not work anymore within macros
@def foreach_segment(_S,_p) {
  coord t = {(_S)[1].x - (_S)[0].x, (_S)[1].y - (_S)[0].y};
  double norm = sqrt(sq(t.x) + sq(t.y));
  assert (norm > 0.);
  t.x = t.x/norm + 1e-6, t.y = t.y/norm - 1.5e-6;
  double alpha = ((_S)[0].x*((_S)[1].y - (_S)[0].y) -
		  (_S)[0].y*((_S)[1].x - (_S)[0].x))/norm;
  foreach()
    if (fabs(t.y*x - t.x*y - alpha) < 0.708*Delta) {
      coord _o = {x,y}, _p[2];
      int _n = 0;
	if (t.x)
	  for (int _i = -1; _i <= 1 && _n < 2; _i += 2) {
	    _p[_n].x = _o.x + _i*Delta/2.;
	    double a = (_p[_n].x - (_S)[0].x)/t.x;
	    _p[_n].y = (_S)[0].y + a*t.y;
	    if (fabs(_p[_n].y - _o.y) <= Delta/2.) {
	      a = clamp (a, 0., norm);
	      _p[_n].x = (_S)[0].x + a*t.x, _p[_n].y = (_S)[0].y + a*t.y;
	      if (fabs(_p[_n].x - _o.x) <= Delta/2. &&
		  fabs(_p[_n].y - _o.y) <= Delta/2.)
		_n++;
	    }
	  }
#if dimension > 1	
	if (t.y)
	  for (int _i = -1; _i <= 1 && _n < 2; _i += 2) {
	    _p[_n].y = _o.y + _i*Delta/2.;
	    double a = (_p[_n].y - (_S)[0].y)/t.y;
	    _p[_n].x = (_S)[0].x + a*t.x;
	    if (fabs(_p[_n].x - _o.x) <= Delta/2.) {
	      a = clamp (a, 0., norm);
	      _p[_n].y = (_S)[0].y + a*t.y, _p[_n].x = (_S)[0].x + a*t.x;
	      if (fabs(_p[_n].y - _o.y) <= Delta/2. &&
		  fabs(_p[_n].x - _o.x) <= Delta/2.)
		_n++;
	    }
	  }
#endif
      if (_n == 2) {
@
#else
@def foreach_segment(_S,_p) {
  coord t = {(_S)[1].x - (_S)[0].x, (_S)[1].y - (_S)[0].y};
  double norm = sqrt(sq(t.x) + sq(t.y));
  assert (norm > 0.);
  t.x = t.x/norm + 1e-6, t.y = t.y/norm - 1.5e-6;
  double alpha = ((_S)[0].x*((_S)[1].y - (_S)[0].y) -
		  (_S)[0].y*((_S)[1].x - (_S)[0].x))/norm;
  foreach()
    if (fabs(t.y*x - t.x*y - alpha) < 0.708*Delta) {
      coord _o = {x,y}, _p[2];
      int _n = 0;
      foreach_dimension()
	if (t.x)
	  for (int _i = -1; _i <= 1 && _n < 2; _i += 2) {
	    _p[_n].x = _o.x + _i*Delta/2.;
	    double a = (_p[_n].x - (_S)[0].x)/t.x;
	    _p[_n].y = (_S)[0].y + a*t.y;
	    if (fabs(_p[_n].y - _o.y) <= Delta/2.) {
	      a = clamp (a, 0., norm);
	      _p[_n].x = (_S)[0].x + a*t.x, _p[_n].y = (_S)[0].y + a*t.y;
	      if (fabs(_p[_n].x - _o.x) <= Delta/2. &&
		  fabs(_p[_n].y - _o.y) <= Delta/2.)
		_n++;
	    }
	  }
      if (_n == 2) {
@
#endif  
@define end_foreach_segment() } } end_foreach(); }

/**
This function returns a summary of the currently-defined fields. */

void fields_stats()
{
  fprintf (ferr, "# t = %g, fields = {", t);
  for (scalar s in all)
    fprintf (ferr, " %s", s.name);
  fputs (" }\n", ferr);
  fprintf (ferr, "# %12s: %12s %12s %12s %12s\n",
	   "name", "min", "avg", "stddev", "max");
  for (scalar s in all) {
    stats ss = statsf (s);
    fprintf (ferr, "# %12s: %12g %12g %12g %12g\n",
	     s.name, ss.min, ss.sum/ss.volume, ss.stddev, ss.max);
  }
}

#include "output.h"

#ifdef DISPLAY
# include "display.h" // nodep
#endif

#endif
