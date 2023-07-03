@if _XOPEN_SOURCE < 700
  @undef _XOPEN_SOURCE
  @define _XOPEN_SOURCE 700
@endif
@if _GNU_SOURCE
@include <stdint.h>
@include <string.h>
@include <fenv.h>
@endif
#define _CATCH
#define dimension 2
#include "common.h"
#ifndef BASILISK_HEADER_0
#define BASILISK_HEADER_0
#line 1 "foreach.c"
#include "grid/quadtree.h"
#include "run.h"
#include "maxruntime.h"

#define MAXTIME 10

scalar * list = {f, f1};
int main(){
	L0 = 16;
	X0 = Y0 = 0;
  init_grid (2); // Initialize a 2 x 2 grid
  origin(X0,Y0);
  run();
}

event initial(t=0)
{
  refine ((x > 8) && (y > 8) && (level < 2)); // Refine to top right corner
  refine ((x > 8) && (x < 12) && (y > 8) && (y < 12) && (level < 5)); // Refine to top right corner
  //unrefine ((x < 8) && (y < 8) && level >= 1); // Coarsen the bottom left corner
}

event test(t = 0)
{
  int i = 0;
  foreach()
  {
    for(s in list)
    foreach_blockf(s)
    {
      s[] = i;
    }
    i++;
  }
}

event end(t = 10)
{}

#endif
