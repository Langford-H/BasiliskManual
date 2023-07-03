#ifndef BASILISK_HEADER_1
#define BASILISK_HEADER_1
#line 1 "/home/dahuanghhc/basilisk/src/maxruntime.h"
/**
# Controlling the maximum runtime

On parallel machines, runs are often not allowed to exceed a maximum
duration (typically 24 hours). To avoid having the job terminated
brutally by the queueing system, this module adds the option to exit
gracefully after a given runtime.

The maximum runtime is given as a command-line argument (typically in
the running script given to the queueing system) using the standard
format H:M:S (hours, minutes and seconds). This should match the
wall-clock time requested from the queueing system.

When this time minus 5 minutes (to allow for clean termination) is
exceeded the state of the simulation is dumped in the "restart" file
and the program terminates. */

static double _maxruntime = HUGE;

event runtime (i += 10) {
  mpi_all_reduce (perf.t, MPI_DOUBLE, MPI_MAX);
  if (perf.t >= _maxruntime - 300) { // we allow 5 minutes for termination
    dump (file = "restart"); // so that we can restart
    return 1; // exit
  }
}

void maxruntime (int * argc, char * argv[])
{
  for (int i = 0; i < *argc; i++)
    if (!strcmp (argv[i], "--maxruntime") || !strcmp (argv[i], "-m")) {
      if (i + 1 < *argc) {
	char * s = strtok (argv[i + 1], ":");
	int n = 0;
	_maxruntime = 0;
	do {
	  _maxruntime = 60*_maxruntime + atoi(s);
	  n++;
	} while ((s = strtok (NULL, ":")));
	if (n > 3) {
	  fprintf (ferr, "maxruntime: TIME format must be H:M:S\n");
	  exit (1);
	}
      }
      else {
	fprintf (ferr, "usage: %s TIME\n", argv[i]);
	exit (1);
      }
      *argc -= 2;
      for (int j = i; j < *argc; j++)
	argv[j] = argv[j + 2];
    }
}

#endif
