#include "grid/quadtree.h"
#include "run.h"
#include "maxruntime.h"

#define MAXTIME 10

scalar f[];
scalar f1[];
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
  //refine ((x > 8) && (x < 12) && (y > 8) && (y < 12) && (level < 3)); // Refine to top right corner
  //unrefine ((x < 8) && (y < 8) && level >= 1); // Coarsen the bottom left corner
}

event test(t = 0)
{
  int i = 0;
  foreach()
  {
    for(scalar s in list)
    foreach_blockf(s)
    {
      s[] = i;
    }
    i++;
  }
}

event end(t = 10)
{}
