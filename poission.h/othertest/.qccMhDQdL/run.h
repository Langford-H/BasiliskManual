#ifndef BASILISK_HEADER_2
#define BASILISK_HEADER_2
#line 1 "/home/dahuanghhc/basilisk/src/run.h"
/**
# Generic time loop

The `run()` function below implements a generic time loop which
executes events until termination.

The timestep `dt` can be accessed as a global variable. */

double dt = 1.;

#include "utils.h"

trace
void run (void)
{
  iter = 0, t = 0., dt = 1.;
  init_grid (N);

  perf.nc = perf.tnc = 0;
  perf.gt = timer_start();
  while (events (true)) {

    /**
    We store the total number of cells advanced in time for computing
    speed statistics. */

    update_perf();
    iter = inext, t = tnext;
  }

  /**
  Time/speed statistics are written out on standard output. */

  timer_print (perf.gt, iter, perf.tnc);

  free_grid();
}

/**
By default we display the boundaries of the domain. */

event defaults (i = 0) {
  display ("box();");
}

/**
At the end of the run we need to empty the default display, otherwise
we would append multiple copies when re-running. */

event cleanup (t = end, last) {
  display ("", true);
}

#endif
