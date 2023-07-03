#ifndef BASILISK_HEADER_27
#define BASILISK_HEADER_27
#line 1 "/home/dahuanghhc/basilisk/src/grid/fpe.h"
// Catching floating-point-exceptions (and other signals)

@include <signal.h>
@include <unistd.h>

static int gdb()
{
  if (last_point.level >= 0) {
    debug (last_point);
    fputc ('\n', stderr);
    fflush (stderr);
  }
  char command[80];
  sprintf (command, "exec xterm -e 'gdb -p %d' & xterm -e 'gnuplot plot -'",
	   getpid());
  return system (command);
}

static void caught_abort (int sig)
{
  fprintf (stderr, "Caught signal %d (Aborted)\n", sig);
  gdb();
}

static void caught_fpe (int sig)
{
  fprintf (stderr, "Caught signal %d (Floating Point Exception)\n", sig);
  gdb();
  exit (1);
}

static void caught_segfault (int sig)
{
  fprintf (stderr, "Caught signal %d (Segmentation Fault)\n", sig);
  gdb();
  exit (2);
}

void catch_fpe (void)
{
  struct sigaction act;
  act.sa_handler = caught_fpe;
  sigemptyset (&act.sa_mask);
  act.sa_flags = 0;
  last_point.level = -1;
  sigaction (8, &act, NULL); // FPE
  act.sa_handler = caught_segfault;
  sigaction (11, &act, NULL); // Segfault
  act.sa_handler = caught_abort;
  act.sa_flags = SA_RESETHAND;
  sigaction (6, &act, NULL); // Abort
}

#endif
