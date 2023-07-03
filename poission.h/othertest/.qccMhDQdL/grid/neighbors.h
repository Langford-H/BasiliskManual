#ifndef BASILISK_HEADER_29
#define BASILISK_HEADER_29
#line 1 "/home/dahuanghhc/basilisk/src/grid/neighbors.h"
#if dimension == 1
@def foreach_neighbor(_s) {
  int _nn = _s + 0 ? _s + 0 : GHOSTS;
  int _i = point.i;
  for (int _k = - _nn; _k <= _nn; _k++) {
    point.i = _i + _k;
    POINT_VARIABLES;
@
@def end_foreach_neighbor()
  }
  point.i = _i;
}
@
@define foreach_neighbor_break() _k = _nn + 1

#elif dimension == 2
@def foreach_neighbor(_s) {
  int _nn = _s + 0 ? _s + 0 : GHOSTS;
  int _i = point.i, _j = point.j;
  for (int _k = - _nn; _k <= _nn; _k++) {
    point.i = _i + _k;
    for (int _l = - _nn; _l <= _nn; _l++) {
      point.j = _j + _l;
      POINT_VARIABLES;
@
@def end_foreach_neighbor()
    }
  }
  point.i = _i; point.j = _j;
}
@
@define foreach_neighbor_break() _k = _l = _nn + 1

#else // dimension == 3
@def foreach_neighbor(_s) {
  int _nn = _s + 0 ? _s + 0 : GHOSTS;
  int _i = point.i, _j = point.j, _k = point.k;
  for (int _l = - _nn; _l <= _nn; _l++) {
    point.i = _i + _l;
    for (int _m = - _nn; _m <= _nn; _m++) {
      point.j = _j + _m;
      for (int _n = - _nn; _n <= _nn; _n++) {
	point.k = _k + _n;
	POINT_VARIABLES;
@
@def end_foreach_neighbor()
      }
    }
  }
  point.i = _i; point.j = _j; point.k = _k;
}
@
@define foreach_neighbor_break() _l = _m = _n = _nn + 1  
#endif // dimension == 3

#endif
