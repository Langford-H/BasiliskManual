#ifndef BASILISK_HEADER_21
#define BASILISK_HEADER_21
#line 1 "/home/dahuanghhc/basilisk/src/grid/tree.h"
#include "mempool.h"

#define TWO_ONE 1 // enforce 2:1 refinement ratio
#define GHOSTS  2

#include "memindex/range.h"

/* By default only one layer of ghost cells is used on the boundary to
   optimise the cost of boundary conditions. */

#ifndef BGHOSTS
@ define BGHOSTS 1
#endif

#define _I     (point.i - GHOSTS)
#if dimension >= 2
# define _J    (point.j - GHOSTS)
#endif
#if dimension >= 3
# define _K    (point.k - GHOSTS)
#endif
#define _DELTA (1./(1 << point.level))

typedef struct {
  unsigned short flags;
  // number of refined neighbors in a 3^dimension neighborhood
  unsigned short neighbors;
  int pid; // process id
} Cell;

enum {
  active    = 1 << 0,
  leaf      = 1 << 1,
  border    = 1 << 2,
  vertex    = 1 << 3,
  user      = 4,

  face_x = 1 << 0
#if dimension >= 2
  , face_y = 1 << 1
#endif
#if dimension >= 3
  , face_z = 1 << 2
#endif
};

@define is_active(cell)  ((cell).flags & active)
@define is_leaf(cell)    ((cell).flags & leaf)
@define is_coarse()      ((cell).neighbors > 0)
@define is_border(cell)  ((cell).flags & border)
@define is_local(cell)   ((cell).pid == pid())
@define is_vertex(cell)  ((cell).flags & vertex)

// Caches

typedef struct {
  int i;
#if dimension >= 2
  int j;
#endif
#if dimension >= 3
  int k;
#endif  
} IndexLevel;

typedef struct {
  IndexLevel * p;
  int n, nm;
} CacheLevel;

typedef struct {
  int i;
#if dimension >= 2
  int j;
#endif
#if dimension >= 3
  int k;
#endif  
  int level, flags;
} Index;

typedef struct {
  Index * p;
  int n, nm;
} Cache;

// Layer

typedef struct {
  Memindex m; // the structure indexing the data
  Mempool * pool; // the memory pool actually holding the data
  long nc;     // the number of allocated elements
  int len;    // the (1D) size of the array
} Layer;

static size_t _size (size_t depth)
{
  return (1 << depth) + 2*GHOSTS;
}

static size_t poolsize (size_t depth, size_t size)
{
  // the maximum amount of data at a given level
#if dimension == 1
  return _size(depth)*size;  
#elif dimension == 2
  return sq(_size(depth))*size;
#else
  return cube(_size(depth))*size;
#endif
}

static Layer * new_layer (int depth)
{
  Layer * l = qmalloc (1, Layer);
  l->len = _size (depth);
  if (depth == 0)
    l->pool = NULL; // the root layer does not use a pool
  else {
    size_t size = sizeof(Cell) + datasize;
    // the block size is 2^dimension*size because we allocate
    // 2^dimension children at a time
    l->pool = mempool_new (poolsize (depth, size), (1 << dimension)*size);
  }
  l->m = mem_new (l->len);
  l->nc = 0;
  return l;
}

static void destroy_layer (Layer * l)
{
  if (l->pool)
    mempool_destroy (l->pool);
  mem_destroy (l->m, l->len);
  free (l);
}

// Tree

typedef struct {
  Grid g;
  Layer ** L; /* the grids at each level */

  Cache        leaves;   /* leaf indices */
  Cache        faces;    /* face indices */
  Cache        vertices; /* vertex indices */
  Cache        refined;  /* refined cells */
  CacheLevel * active;   /* active cells indices for each level */
  CacheLevel * prolongation; /* halo prolongation indices for each level */
  CacheLevel * boundary;  /* boundary indices for each level */
  /* indices of boundary cells with non-boundary parents */
  CacheLevel * restriction;
  
  bool dirty;       /* whether caches should be updated */
} Tree;

#define tree ((Tree *)grid)

struct _Point {
  /* the current cell index and level */
  int i;
#if dimension >= 2
  int j;
#endif
#if dimension >= 3
  int k;
#endif
  int level;
@ifdef foreach_block
  int l;
  @define _BLOCK_INDEX , point.l
@else
  @define _BLOCK_INDEX
@endif
};
static Point last_point;

#define BSIZE 128

static void cache_level_append (CacheLevel * c, Point p)
{
  if (c->n >= c->nm) {
    c->nm += BSIZE;
    qrealloc (c->p, c->nm, IndexLevel);
  }
  c->p[c->n].i = p.i;
#if dimension >= 2
  c->p[c->n].j = p.j;
#endif
#if dimension >= 3
  c->p[c->n].k = p.k;
#endif
  c->n++;
}

static void cache_level_shrink (CacheLevel * c)
{
  if (c->nm > (c->n/BSIZE + 1)*BSIZE) {
    c->nm = (c->n/BSIZE + 1)*BSIZE;
    assert (c->nm > c->n);
    c->p = (IndexLevel *) realloc (c->p, sizeof (Index)*c->nm);
  }
}

static void cache_append (Cache * c, Point p, unsigned short flags)
{
  if (c->n >= c->nm) {
    c->nm += BSIZE;
    qrealloc (c->p, c->nm, Index);
  }
  c->p[c->n].i = p.i;
#if dimension >= 2
  c->p[c->n].j = p.j;
#endif
#if dimension >= 3
  c->p[c->n].k = p.k;
#endif  
  c->p[c->n].level = p.level;
  c->p[c->n].flags = flags;
  c->n++;
}

void cache_shrink (Cache * c)
{
  cache_level_shrink ((CacheLevel *)c);
}

#undef BSIZE

/* low-level memory management */
#if dimension == 1
@define allocated(k,l,n) (mem_allocated (tree->L[point.level]->m, point.i+k))
@define NEIGHBOR(k,l,n)	(mem_data (tree->L[point.level]->m, point.i+k))
@def PARENT(k,l,n) (mem_data (tree->L[point.level-1]->m,
			      (point.i+GHOSTS)/2+k))
@
@def allocated_child(k,l,n) (level < depth() &&
			     mem_allocated(tree->L[point.level+1]->m,
					   2*point.i-GHOSTS+k))
@
@define CHILD(k,l,n)  (mem_data (tree->L[point.level+1]->m, 2*point.i-GHOSTS+k))
#elif dimension == 2
@def allocated(k,l,n) (mem_allocated (tree->L[point.level]->m,
				      point.i+k, point.j+l))
@
@def NEIGHBOR(k,l,n)	(mem_data (tree->L[point.level]->m,
				   point.i+k, point.j+l))
@
@def PARENT(k,l,n) (mem_data (tree->L[point.level-1]->m,
			      (point.i+GHOSTS)/2+k, (point.j+GHOSTS)/2+l))
@
@def allocated_child(k,l,n)  (level < depth() &&
			      mem_allocated (tree->L[point.level+1]->m,
					     2*point.i-GHOSTS+k,
					     2*point.j-GHOSTS+l))
@
@def CHILD(k,l,n)  (mem_data (tree->L[point.level+1]->m,
			      2*point.i-GHOSTS+k, 2*point.j-GHOSTS+l))
@
#else // dimension == 3
@def allocated(a,l,n) (mem_allocated (tree->L[point.level]->m,
				      point.i+a, point.j+l, point.k+n))
@
@def NEIGHBOR(a,l,n)	(mem_data (tree->L[point.level]->m,
				   point.i+a, point.j+l, point.k+n))
@			
@def PARENT(a,l,n) (mem_data (tree->L[point.level-1]->m,
			      (point.i+GHOSTS)/2+a,
			      (point.j+GHOSTS)/2+l,
			      (point.k+GHOSTS)/2+n))
@
@def allocated_child(a,l,n)  (level < depth() &&
			      mem_allocated (tree->L[point.level+1]->m,
					     2*point.i-GHOSTS+a,
					     2*point.j-GHOSTS+l,
					     2*point.k-GHOSTS+n))
@
@def CHILD(a,l,n)  (mem_data (tree->L[point.level+1]->m,
			      2*point.i-GHOSTS+a,
			      2*point.j-GHOSTS+l,
			      2*point.k-GHOSTS+n))
@
#endif // dimension == 3
@define CELL(m) (*((Cell *)(m)))

/***** Multigrid macros *****/
@define depth()        (grid->depth)
@define aparent(k,l,n) CELL(PARENT(k,l,n))
@define child(k,l,n)   CELL(CHILD(k,l,n))

/***** Tree macros ****/
@define cell		 CELL(NEIGHBOR(0,0,0))
@define neighbor(k,l,n)	 CELL(NEIGHBOR(k,l,n))
@def neighborp(l,m,n) (Point) {
    point.i + l,
#if dimension >= 2
    point.j + m,
#endif
#if dimension >= 3
    point.k + n,
#endif
    point.level
    _BLOCK_INDEX
}
@
			
/***** Data macros *****/
@define data(k,l,n)     ((double *) (NEIGHBOR(k,l,n) + sizeof(Cell)))
@define fine(a,k,p,n)   ((double *) (CHILD(k,p,n) + sizeof(Cell)))[_index(a,n)]
@define coarse(a,k,p,n) ((double *) (PARENT(k,p,n) + sizeof(Cell)))[_index(a,n)]

@def POINT_VARIABLES
  VARIABLES
  int level = point.level; NOT_UNUSED(level);
#if dimension == 1
  struct { int x; } child = { 2*((point.i+GHOSTS)%2)-1 };
#elif dimension == 2
  struct { int x, y; } child = {
    2*((point.i+GHOSTS)%2)-1, 2*((point.j+GHOSTS)%2)-1
  };
#else
  struct { int x, y, z; } child = {
    2*((point.i+GHOSTS)%2)-1, 2*((point.j+GHOSTS)%2)-1, 2*((point.k+GHOSTS)%2)-1
  };
#endif
  NOT_UNUSED(child);
  Point parent = point;	NOT_UNUSED(parent);
  parent.level--;
  parent.i = (point.i + GHOSTS)/2;
#if dimension >= 2
  parent.j = (point.j + GHOSTS)/2;
#endif
#if dimension >= 3
  parent.k = (point.k + GHOSTS)/2;
#endif
#if TRASH
  Cell * cellp = point.level <= depth() && allocated(0,0,0) ?
    (Cell *) NEIGHBOR(0,0,0) : NULL;
  NOT_UNUSED(cellp);
#endif
@

#include "foreach_cell.h"

#if dimension == 1
@def foreach_child() {
  int _i = 2*point.i - GHOSTS;
  point.level++;
  for (int _k = 0; _k < 2; _k++) {
    point.i = _i + _k;
    POINT_VARIABLES;
@
@def end_foreach_child()
  }
  point.i = (_i + GHOSTS)/2;
  point.level--;
}
@
@define foreach_child_break() _k = 2
#elif dimension == 2
@def foreach_child() {
  int _i = 2*point.i - GHOSTS, _j = 2*point.j - GHOSTS;
  point.level++;
  for (int _k = 0; _k < 2; _k++) {
    point.i = _i + _k;
    for (int _l = 0; _l < 2; _l++) {
      point.j = _j + _l;
      POINT_VARIABLES;
@
@def end_foreach_child()
    }
  }
  point.i = (_i + GHOSTS)/2; point.j = (_j + GHOSTS)/2;
  point.level--;
}
@
@define foreach_child_break() _k = _l = 2
#else // dimension == 3
@def foreach_child() {
  int _i = 2*point.i - GHOSTS, _j = 2*point.j - GHOSTS, _k = 2*point.k - GHOSTS;
  point.level++;
  for (int _l = 0; _l < 2; _l++) {
    point.i = _i + _l;
    for (int _m = 0; _m < 2; _m++) {
      point.j = _j + _m;
      for (int _n = 0; _n < 2; _n++) {
	point.k = _k + _n;
	POINT_VARIABLES;
@
@def end_foreach_child()
      }
    }
  }
  point.i = (_i + GHOSTS)/2;point.j = (_j + GHOSTS)/2;point.k = (_k + GHOSTS)/2;
  point.level--;
}
@
@define foreach_child_break() _l = _m = _n = 2
#endif // dimension == 3
  
#define update_cache() { if (tree->dirty) update_cache_f(); }

#define is_refined(cell)      (!is_leaf (cell) && cell.neighbors && cell.pid >= 0)
#define is_prolongation(cell) (!is_leaf(cell) && !cell.neighbors && cell.pid >= 0)
#define is_boundary(cell)     (cell.pid < 0)
  
@def is_refined_check() (is_refined(cell) &&
			 point.i > 0 && point.i < (1 << level) + 2*GHOSTS - 1
#if dimension > 1
			 && point.j > 0 && point.j < (1 << level) + 2*GHOSTS - 1
#endif
#if dimension > 2
			 && point.k > 0 && point.k < (1 << level) + 2*GHOSTS - 1
#endif
			 )
@

@def foreach_cache(_cache) {
  OMP_PARALLEL() {
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
  Point point = {0};
  point.i = GHOSTS;
#if dimension > 1
  point.j = GHOSTS;
#endif
#if dimension > 2
  point.k = GHOSTS;
#endif
  int _k; unsigned short _flags; NOT_UNUSED(_flags);
  OMP(omp for schedule(static))
  for (_k = 0; _k < _cache.n; _k++) {
    point.i = _cache.p[_k].i;
#if dimension >= 2
    point.j = _cache.p[_k].j;
#endif
#if dimension >= 3
    point.k = _cache.p[_k].k;
#endif
    point.level = _cache.p[_k].level;
    _flags = _cache.p[_k].flags;
    POINT_VARIABLES;
@
@define end_foreach_cache() } } }

@def foreach_cache_level(_cache,_l) {
  OMP_PARALLEL() {
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
  Point point = {0};
  point.i = GHOSTS;
#if dimension > 1
  point.j = GHOSTS;
#endif
#if dimension > 2
  point.k = GHOSTS;
#endif
  point.level = _l;
  int _k;
  OMP(omp for schedule(static))
  for (_k = 0; _k < _cache.n; _k++) {
    point.i = _cache.p[_k].i;
#if dimension >= 2
    point.j = _cache.p[_k].j;
#endif
#if dimension >= 3
    point.k = _cache.p[_k].k;
#endif
    POINT_VARIABLES;
@
@define end_foreach_cache_level() } } }

@def foreach_boundary_level(_l) {
  if (_l <= depth()) {
    update_cache();
    CacheLevel _boundary = tree->boundary[_l];
    foreach_cache_level (_boundary,_l)
@
@define end_foreach_boundary_level() end_foreach_cache_level(); }}

#define bid(cell) (- cell.pid - 1)

@def foreach_boundary(_b) {
  for (int _l = depth(); _l >= 0; _l--)
    foreach_boundary_level(_l) {
      if (bid(cell) == _b)
	for (int _d = 0; _d < dimension; _d++) {
	  for (int _i = -1; _i <= 1; _i += 2) {
	    if (_d == 0) ig = _i; else if (_d == 1) jg = _i; else kg = _i;
	    if (allocated(-ig,-jg,-kg) &&
		is_leaf (neighbor(-ig,-jg,-kg)) &&
		!is_boundary(neighbor(-ig,-jg,-kg)) &&
		is_local(neighbor(-ig,-jg,-kg))) {
	      point.i -= ig; x -= ig*Delta/2.; 
#if dimension >= 2
	      point.j -= jg; y -= jg*Delta/2.; 
#endif
#if dimension >= 3
	      point.k -= kg; z -= kg*Delta/2.;
#endif
@
@def end_foreach_boundary()
	      point.i += ig; x += ig*Delta/2.;   
#if dimension >= 2
	      point.j += jg; y += jg*Delta/2.; 
#endif
#if dimension >= 3
	      point.k += kg; z += kg*Delta/2.;
#endif
            }
	  }
	  ig = jg = kg = 0;
	}
    } end_foreach_boundary_level(); }
@
  
@def foreach_halo(_name,_l) {
  if (_l <= depth()) {
    update_cache();
    CacheLevel _cache = tree->_name[_l];
    foreach_cache_level (_cache, _l)
@
@define end_foreach_halo() end_foreach_cache_level(); }}

#include "neighbors.h"

static inline bool has_local_children (Point point)
{
  foreach_child()
    if (is_local(cell))
      return true;
  return false;
}

static inline void cache_append_face (Point point, unsigned short flags)
{
  Tree * q = tree;
  cache_append (&q->faces, point, flags);
#if dimension == 2
  if (!is_vertex(cell)) {
    cache_append (&q->vertices, point, 0);
    cell.flags |= vertex;
  }
  foreach_dimension()
    if ((flags & face_y) && !is_vertex(neighbor(1))) {
      cache_append (&q->vertices, neighborp(1), 0);
      neighbor(1).flags |= vertex;
    }
#elif dimension == 3
  foreach_dimension()
    if (flags & face_x)
      for (int i = 0; i <= 1; i++)
	for (int j = 0; j <= 1; j++)
	  if (!is_vertex(neighbor(0,i,j))) {
	    cache_append (&q->vertices, neighborp(0,i,j), 0);
	    neighbor(0,i,j).flags |= vertex;
	  }
#endif
}

#define FBOUNDARY 1 // fixme: this should work with zero

static void update_cache_f (void)
{
  Tree * q = tree;

  foreach_cache (q->vertices)
    if (level <= depth() && allocated(0))
      cell.flags &= ~vertex;
  
  /* empty caches */
  q->leaves.n = q->faces.n = q->vertices.n = 0;
  for (int l = 0; l <= depth(); l++)
    q->active[l].n = q->prolongation[l].n =
      q->boundary[l].n = q->restriction[l].n = 0;
#if FBOUNDARY
  const unsigned short fboundary = 1 << user;
  foreach_cell() {
#else    
  foreach_cell_all() {
#endif
    if (is_local(cell) && is_active(cell)) {
      // active cells
      //      assert (is_active(cell));
      cache_level_append (&q->active[level], point);
    }
#if !FBOUNDARY
    if (is_boundary(cell)) {
      // boundary conditions
      bool has_neighbors = false;
      foreach_neighbor (BGHOSTS)
	if (allocated(0) && !is_boundary(cell)) {
	  has_neighbors = true; break;
	}
      if (has_neighbors)
	cache_level_append (&q->boundary[level], point);
      // restriction for masked cells
      if (level > 0 && is_local(aparent(0)))
	cache_level_append (&q->restriction[level], point);
    }
#else
    // boundaries
    if (!is_boundary(cell)) {
      // look in a 5x5 neighborhood for boundary cells
      foreach_neighbor (BGHOSTS)
	if (allocated(0) && is_boundary(cell) && !(cell.flags & fboundary)) {
	  cache_level_append (&q->boundary[level], point);
	  cell.flags |= fboundary;
	}
    }
    // restriction for masked cells
    else if (level > 0 && is_local(aparent(0)))
      cache_level_append (&q->restriction[level], point);
#endif
    if (is_leaf (cell)) {
      if (is_local(cell)) {
	cache_append (&q->leaves, point, 0);
	// faces
	unsigned short flags = 0;
	foreach_dimension()
	  if (is_boundary(neighbor(-1)) || is_prolongation(neighbor(-1)) ||
	      is_leaf(neighbor(-1)))
	    flags |= face_x;
	if (flags)
	  cache_append (&q->faces, point, flags);
	foreach_dimension()
	  if (is_boundary(neighbor(1)) || is_prolongation(neighbor(1)) ||
	      (!is_local(neighbor(1)) && is_leaf(neighbor(1))))
	    cache_append (&q->faces, neighborp(1), face_x);
	// vertices
	for (int i = 0; i <= 1; i++)
        #if dimension >= 2
	  for (int j = 0; j <= 1; j++)
        #endif
          #if dimension >= 3
	    for (int k = 0; k <= 1; k++)
	  #endif
	      if (!is_vertex(neighbor(i,j,k))) {
		cache_append (&q->vertices, neighborp(i,j,k), 0);
		neighbor(i,j,k).flags |= vertex;
	      }
	// halo prolongation
        if (cell.neighbors > 0)
	  cache_level_append (&q->prolongation[level], point);
      }
      else if (!is_boundary(cell) || is_local(aparent(0))) { // non-local
	// faces
	unsigned short flags = 0;
	foreach_dimension()
	  if (allocated(-1) &&
	      is_local(neighbor(-1)) && is_prolongation(neighbor(-1)))
	    flags |= face_x;
	if (flags)
	  cache_append_face (point, flags);
	foreach_dimension()
	  if (allocated(1) && is_local(neighbor(1)) &&
	      is_prolongation(neighbor(1)))
	    cache_append_face (neighborp(1), face_x);
      }
#if FBOUNDARY // fixme: this should always be included
      continue; 
#endif
    }
  }

  /* optimize caches */
  cache_shrink (&q->leaves);
  cache_shrink (&q->faces);
  cache_shrink (&q->vertices);
  for (int l = 0; l <= depth(); l++) {
    cache_level_shrink (&q->active[l]);
    cache_level_shrink (&q->prolongation[l]);
    cache_level_shrink (&q->boundary[l]);
    cache_level_shrink (&q->restriction[l]);
}
  
  q->dirty = false;

#if FBOUNDARY
  for (int l = depth(); l >= 0; l--)
    foreach_boundary_level (l)
      cell.flags &= ~fboundary;
#endif
  
  // mesh size
  grid->n = q->leaves.n;
  // for MPI the reduction operation over all processes is done by balance()
@if !_MPI
  grid->tn = grid->n;
  grid->maxdepth = grid->depth;
@endif
}

@define foreach() update_cache(); foreach_cache(tree->leaves)
@define end_foreach()   end_foreach_cache()

@def foreach_face_generic()
  update_cache();
  foreach_cache(tree->faces) @
@define end_foreach_face_generic() end_foreach_cache()

@define is_face_x() { int ig = -1; VARIABLES; if (_flags & face_x) {
@define end_is_face_x() }}
      
#if dimension >= 2
@define is_face_y() { int jg = -1; VARIABLES; if (_flags & face_y) {
@define end_is_face_y() }}
#endif
#if dimension >= 3
@define is_face_z() { int kg = -1; VARIABLES; if (_flags & face_z) {
@define end_is_face_z() }}
#endif
    
@def foreach_vertex()
  update_cache();
  foreach_cache(tree->vertices) {
    x -= Delta/2.;
#if dimension >= 2
    y -= Delta/2.;
#endif
#if dimension >= 3
    z -= Delta/2.;
#endif
@
@define end_foreach_vertex() } end_foreach_cache()

#if dimension == 3
# define foreach_edge()				\
    foreach_vertex()				\
      foreach_dimension()			\
        if (is_vertex(neighbor(1)))
#else // dimension < 3
# define foreach_edge() foreach_face(y,x)
#endif

@def foreach_level(l) {
  if (l <= depth()) {
    update_cache();
    CacheLevel _active = tree->active[l];
    foreach_cache_level (_active,l)
@
@define end_foreach_level() end_foreach_cache_level(); }}

@define foreach_coarse_level(l) foreach_level(l) if (!is_leaf(cell)) {
@define end_foreach_coarse_level() } end_foreach_level()

@def foreach_level_or_leaf(l) {
  for (int _l1 = l; _l1 >= 0; _l1--)
    foreach_level(_l1)
      if (_l1 == l || is_leaf (cell)) {
@
@define end_foreach_level_or_leaf() } end_foreach_level(); }

@if TRASH
@ undef trash
@ define trash(list) reset(list, undefined)
@endif

void reset (void * alist, double val)
{
  scalar * list = (scalar *) alist;
  Tree * q = tree;
  /* low-level memory management */
  for (int l = 0; l <= depth(); l++) {
    Layer * L = q->L[l];
    foreach_mem (L->m, L->len, 1) {
      point.level = l;
      for (scalar s in list) {
	if (!is_constant(s))
	  for (int b = 0; b < s.block; b++)
	    data(0,0,0)[s.i + b] = val;
      }
    }
  }
}

static CacheLevel * cache_level_resize (CacheLevel * name, int a)
{
  for (int i = 0; i <= depth() - a; i++)
    free (name[i].p);
  free (name);
  return qcalloc (depth() + 1, CacheLevel);
}
  
static void update_depth (int inc)
{
  Tree * q = tree;
  grid->depth += inc;
  q->L = &(q->L[-1]);
  qrealloc (q->L, grid->depth + 2, Layer *);
  q->L = &(q->L[1]);
  if (inc > 0)
    q->L[grid->depth] = new_layer (grid->depth);
  q->active = cache_level_resize (q->active, inc);
  q->prolongation = cache_level_resize (q->prolongation, inc);
  q->boundary = cache_level_resize (q->boundary, inc);
  q->restriction = cache_level_resize (q->restriction, inc);
}

#if dimension == 1
typedef void (* PeriodicFunction) (Memindex, int, int, void *);

static void periodic_function (Memindex m, int i, int len, void * b,
			       PeriodicFunction f)
{
  f(m, i, len, b);
  if (Period.x) {
    int nl = len - 2*GHOSTS;
    for (int l = - 1; l <= 1; l += 2)
      for (int n = i + l*nl; n >= 0 && n < len; n += l*nl)
	f(m, n, len, b);
  }
}
			       
static void assign_periodic (Memindex m, int i, int len, void * b)
{
  periodic_function (m, i, len, b, mem_assign);
}

static void free_periodic (Memindex m, int i, int len)
{
  periodic_function (m, i, len, NULL, (PeriodicFunction) mem_free);
}
#elif dimension == 2
typedef void (* PeriodicFunction) (Memindex, int, int, int, void *);

static void periodic_function (Memindex m, int i, int j, int len, void * b,
			       PeriodicFunction f)
{
  f(m, i, j, len, b);
  if (Period.x) {
    int nl = len - 2*GHOSTS;
    for (int l = - 1; l <= 1; l += 2)
      for (int n = i + l*nl; n >= 0 && n < len; n += l*nl)
	f(m, n, j, len, b);
    if (Period.y)
      for (int l = - 1; l <= 1; l += 2)
	for (int n = j + l*nl; n >= 0 && n < len; n += l*nl) {
	  f(m, i, n, len, b);
	  for (int o = - 1; o <= 1; o += 2)
	    for (int p = i + o*nl; p >= 0 && p < len; p += o*nl)
	      f(m, p, n, len, b);
	}
  }
  else if (Period.y) {
    int nl = len - 2*GHOSTS;
    for (int l = - 1; l <= 1; l += 2)
      for (int n = j + l*nl; n >= 0 && n < len; n += l*nl)
	f(m, i, n, len, b);
  }
}
			       
static void assign_periodic (Memindex m, int i, int j, int len, void * b)
{
  periodic_function (m, i, j, len, b, mem_assign);
}

static void free_periodic (Memindex m, int i, int j, int len)
{
  periodic_function (m, i, j, len, NULL, (PeriodicFunction) mem_free);
}
#else // dimension == 3
typedef void (* PeriodicFunction) (Memindex, int, int, int, int, void *);

static void periodic_function (Memindex m, int i, int j, int k, int len,
			       void * b, PeriodicFunction f)
{
  f(m, i, j, k, len, b);
  if (Period.x) {
    int nl = len - 2*GHOSTS;
    for (int l = - 1; l <= 1; l += 2)
      for (int n = i + l*nl; n >= 0 && n < len; n += l*nl)
	f(m, n, j, k, len, b);
    if (Period.y) {
      for (int l = - 1; l <= 1; l += 2)
	for (int n = j + l*nl; n >= 0 && n < len; n += l*nl) {
	  f(m, i, n, k, len, b);
	  for (int o = - 1; o <= 1; o += 2)
	    for (int p = i + o*nl; p >= 0 && p < len; p += o*nl)
	      f(m, p, n, k, len, b);
	}
      if (Period.z)
	for (int l = - 1; l <= 1; l += 2)
	  for (int n = k + l*nl; n >= 0 && n < len; n += l*nl) {
	    f(m, i, j, n, len, b);
	    for (int q = - 1; q <= 1; q += 2)
	      for (int r = j + q*nl; r >= 0 && r < len; r += q*nl)
		f(m, i, r, n, len, b);
	    for (int o = - 1; o <= 1; o += 2)
	      for (int p = i + o*nl; p >= 0 && p < len; p += o*nl) {
		f(m, p, j, n, len, b);
		for (int q = - 1; q <= 1; q += 2)
		  for (int r = j + q*nl; r >= 0 && r < len; r += q*nl)
		    f(m, p, r, n, len, b);
	      }
	  }
    }
    else if (Period.z)
      for (int l = - 1; l <= 1; l += 2)
	for (int n = k + l*nl; n >= 0 && n < len; n += l*nl) {
	  f(m, i, j, n, len, b);
	  for (int o = - 1; o <= 1; o += 2)
	    for (int p = i + o*nl; p >= 0 && p < len; p += o*nl)
	      f(m, p, j, n, len, b);
	}
  }
  else if (Period.y) {
    int nl = len - 2*GHOSTS;
    for (int l = - 1; l <= 1; l += 2)
      for (int n = j + l*nl; n >= 0 && n < len; n += l*nl)
	f(m, i, n, k, len, b);
    if (Period.z)
      for (int l = - 1; l <= 1; l += 2)
	for (int n = k + l*nl; n >= 0 && n < len; n += l*nl) {
	  f(m, i, j, n, len, b);
	  for (int o = - 1; o <= 1; o += 2)
	    for (int p = j + o*nl; p >= 0 && p < len; p += o*nl)
	      f(m, i, p, n, len, b);
	}
  }
  else if (Period.z) {
    int nl = len - 2*GHOSTS;
    for (int l = - 1; l <= 1; l += 2)
      for (int n = k + l*nl; n >= 0 && n < len; n += l*nl)
	f(m, i, j, n, len, b);
  }
}

static void assign_periodic (Memindex m, int i, int j, int k, int len, void * b)
{
  periodic_function (m, i, j, k, len, b, mem_assign);
}

static void free_periodic (Memindex m, int i, int j, int k, int len)
{
  periodic_function (m, i, j, k, len, NULL, (PeriodicFunction) mem_free);
}
#endif // dimension == 3

static void alloc_children (Point point)
{
  if (point.level == grid->depth)
    update_depth (+1);
  else if (allocated_child(0,0,0))
    return;
  
  /* low-level memory management */
  Layer * L = tree->L[point.level + 1];
  L->nc++;
  size_t len = sizeof(Cell) + datasize;
  char * b = (char *) mempool_alloc0 (L->pool);
  int i = 2*point.i - GHOSTS;
  for (int k = 0; k < 2; k++, i++) {
#if dimension == 1
    assign_periodic (L->m, i, L->len, b);
    b += len;
#elif dimension == 2
    int j = 2*point.j - GHOSTS;
    for (int l = 0; l < 2; l++, j++) {
      assign_periodic (L->m, i, j, L->len, b);
      b += len;
    }
#else // dimension == 3
    int j = 2*point.j - GHOSTS;
    for (int l = 0; l < 2; l++, j++) {
      int m = 2*point.k - GHOSTS;
      for (int n = 0; n < 2; n++, m++) {
	assign_periodic (L->m, i, j, m, L->len, b);
	b += len;
      }
    }
#endif
  }
  
  int pid = cell.pid;
  foreach_child() {
    cell.pid = pid;
@if TRASH
    for (scalar s in all)
      s[] = undefined;
@endif    
  }
}

#if dimension == 1 
static void free_children (Point point)
{
  /* low-level memory management */
  Layer * L = tree->L[point.level + 1];
  int i = 2*point.i - GHOSTS;
  assert (mem_data (L->m,i));
  mempool_free (L->pool, mem_data (L->m,i));
  for (int k = 0; k < 2; k++, i++)
    free_periodic (L->m, i, L->len);
  if (--L->nc == 0) {
    destroy_layer (L);
    assert (point.level + 1 == grid->depth);
    update_depth (-1);
  }
}
#elif dimension == 2
static void free_children (Point point)
{
  /* low-level memory management */
  Layer * L = tree->L[point.level + 1];
  int i = 2*point.i - GHOSTS, j = 2*point.j - GHOSTS;
  assert (mem_data (L->m,i,j));
  mempool_free (L->pool, mem_data (L->m,i,j));
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++)
      free_periodic (L->m, i + k, j + l, L->len);
  if (--L->nc == 0) {
    destroy_layer (L);
    assert (point.level + 1 == grid->depth);
    update_depth (-1);
  }
}
#else // dimension == 3
static void free_children (Point point)
{
  /* low-level memory management */
  Layer * L = tree->L[point.level + 1];
  int i = 2*point.i - GHOSTS;
  assert (mem_data (L->m,i,2*point.j - GHOSTS,2*point.k - GHOSTS));
  mempool_free (L->pool, mem_data (L->m,
				   i,2*point.j - GHOSTS,2*point.k - GHOSTS));
  for (int k = 0; k < 2; k++, i++) {
    int j = 2*point.j - GHOSTS;
    for (int l = 0; l < 2; l++, j++) {
      int m = 2*point.k - GHOSTS;
      for (int n = 0; n < 2; n++, m++)
	free_periodic (L->m, i, j, m, L->len);
    }
  }
  if (--L->nc == 0) {
    destroy_layer (L);
    assert (point.level + 1 == grid->depth);
    update_depth (-1);
  }
}
#endif // dimension == 3

void increment_neighbors (Point point)
{
  tree->dirty = true;
  if (cell.neighbors++ == 0)
    alloc_children (point);
  foreach_neighbor (GHOSTS/2)
    if (cell.neighbors++ == 0)
      alloc_children (point);
  cell.neighbors--;
}

void decrement_neighbors (Point point)
{
  tree->dirty = true;
  foreach_neighbor (GHOSTS/2)
    if (allocated(0)) {
      cell.neighbors--;
      if (cell.neighbors == 0)
	free_children (point);
    }
  if (cell.neighbors) {
    int pid = cell.pid;
    foreach_child() {
      cell.flags = 0;
      cell.pid = pid;
    }
  }
}

void realloc_scalar (int size)
{
  /* low-level memory management */
  Tree * q = tree;
  size_t oldlen = sizeof(Cell) + datasize;
  size_t newlen = oldlen + size;
  datasize += size;
  /* the root level is allocated differently */
  Layer * L = q->L[0];
  foreach_mem (L->m, L->len, 1) {
#if dimension == 1
    char * p = (char *) realloc (mem_data (L->m, point.i), newlen*sizeof(char));
    assign_periodic (L->m, point.i, L->len, p);
#elif dimension == 2
    char * p = (char *) realloc (mem_data (L->m, point.i, point.j),
				 newlen*sizeof(char));
    assign_periodic (L->m, point.i, point.j, L->len, p);
#else
    char * p = (char *) realloc (mem_data (L->m, point.i, point.j, point.k),
				 newlen*sizeof(char));
    assign_periodic (L->m, point.i, point.j, point.k, L->len, p);
#endif
  }
  /* all other levels */
  for (int l = 1; l <= depth(); l++) {
    Layer * L = q->L[l];
    Mempool * oldpool = L->pool;
    L->pool = mempool_new (poolsize (l, newlen), (1 << dimension)*newlen);
    foreach_mem (L->m, L->len, 2) {
      char * new = (char *) mempool_alloc (L->pool);
#if dimension == 1
      for (int k = 0; k < 2; k++) {
	memcpy (new, mem_data (L->m, point.i + k), oldlen);
	assign_periodic (L->m, point.i + k, L->len, new);
	new += newlen;
      }
#elif dimension == 2
      for (int k = 0; k < 2; k++)
	for (int o = 0; o < 2; o++) {
	  memcpy (new, mem_data (L->m, point.i + k,point.j + o), oldlen);
	  assign_periodic (L->m, point.i + k, point.j + o, L->len, new);
	  new += newlen;
	}
#else // dimension == 3
      for (int l = 0; l < 2; l++)
	for (int m = 0; m < 2; m++)
	  for (int n = 0; n < 2; n++) {
	    memcpy (new, mem_data (L->m, point.i + l, point.j + m, point.k + n),
		    oldlen);
	    assign_periodic (L->m, point.i + l, point.j + m, point.k + n,
				 L->len, new);
	    new += newlen;
	  }
#endif // dimension == 3
    }
    mempool_destroy (oldpool);
  }
}

/* Boundaries */

@define VN v.x
@define VT v.y
@define VR v.z

#define is_neighbor(...) (allocated(__VA_ARGS__) &&		\
			  !is_boundary(neighbor(__VA_ARGS__)))

@if _MPI // fixme
@ define disable_fpe_for_mpi() disable_fpe (FE_DIVBYZERO|FE_INVALID)
@ define enable_fpe_for_mpi()  enable_fpe (FE_DIVBYZERO|FE_INVALID)
@else
@ define disable_fpe_for_mpi()
@ define enable_fpe_for_mpi()
@endif

static inline void no_restriction (Point point, scalar s);

static bool normal_neighbor (Point point, scalar * scalars, vector * vectors)
{
  for (int k = 1; k <= BGHOSTS; k++)
    foreach_dimension()
      for (int i = -k; i <= k; i += 2*k)
	if (is_neighbor(i)) {
	  Point neighbor = neighborp(i);
	  int id = bid(cell);
	  for (scalar s in scalars)
	    foreach_block()
	      s[] = s.boundary[id](neighbor, point, s, NULL);
	  for (vector v in vectors)
	    foreach_block() {
	      scalar vn = VN;
	      v.x[] = vn.boundary[id](neighbor, point, v.x, NULL);
#if dimension >= 2
	      scalar vt = VT;
	      v.y[] = vt.boundary[id](neighbor, point, v.y, NULL);
#endif
#if dimension >= 3
	      scalar vr = VR;
	      v.z[] = vr.boundary[id](neighbor, point, v.z, NULL);
#endif
	    }
	  return true;
	}
  return false;
}

static bool diagonal_neighbor_2D (Point point,
				  scalar * scalars, vector * vectors)
{
#if dimension >= 2
  for (int k = 1; k <= BGHOSTS; k++)
#if dimension == 3
    foreach_dimension()
#endif
      for (int i = -k; i <= k; i += 2*k)
	for (int j = -k; j <= k; j += 2*k)
	  if (allocated(i,j) && is_neighbor(i,j) &&
	      allocated(i,0) && is_boundary(neighbor(i,0)) &&
	      allocated(0,j) && is_boundary(neighbor(0,j))) {
	    Point n = neighborp(i,j),
	      n1 = neighborp(i,0), n2 = neighborp(0,j);
	    int id1 = bid(neighbor(i,0)), id2 = bid(neighbor(0,j));
	    for (scalar s in scalars)
	      foreach_block()
		s[] = (s.boundary[id1](n,n1,s,NULL) +
		       s.boundary[id2](n,n2,s,NULL) -
		       s[i,j]);
	    for (vector v in vectors)
	      foreach_block() {
		scalar vt = VT, vn = VN;
		v.x[] = (vt.boundary[id1](n,n1,v.x,NULL) +
			 vn.boundary[id2](n,n2,v.x,NULL) -
			 v.x[i,j]);
		v.y[] = (vn.boundary[id1](n,n1,v.y,NULL) +
			 vt.boundary[id2](n,n2,v.y,NULL) -
			 v.y[i,j]);
#if dimension == 3
		scalar vr = VR;
		v.z[] = (vr.boundary[id1](n,n1,v.z,NULL) +
			 vr.boundary[id2](n,n2,v.z,NULL) -
			 v.z[i,j]);
#endif
	      }
	    return true;
	  }
#endif // dimension >= 2
  return false;
}

static bool diagonal_neighbor_3D (Point point,
				  scalar * scalars, vector * vectors)
{
#if dimension == 3
  for (int n = 1; n <= BGHOSTS; n++)
    for (int i = -n; i <= n; i += 2*n)
      for (int j = -n; j <= n; j += 2*n)
	for (int k = -n; k <= n; k += 2*n)
	  if (is_neighbor(i,j,k) &&
	      is_boundary(neighbor(i,j,0)) &&
	      is_boundary(neighbor(i,0,k)) &&
	      is_boundary(neighbor(0,j,k))) {
	    Point
	      n0 = neighborp(i,j,k),
	      n1 = neighborp(i,j,0),
	      n2 = neighborp(i,0,k),
	      n3 = neighborp(0,j,k);
	    int
	      id1 = bid(neighbor(i,j,0)),
	      id2 = bid(neighbor(i,0,k)),
	      id3 = bid(neighbor(0,j,k));
	    for (scalar s in scalars)
	      foreach_block()
		s[] = (s.boundary[id1](n0,n1,s,NULL) +
		       s.boundary[id2](n0,n2,s,NULL) +
		       s.boundary[id3](n0,n3,s,NULL) -
		       2.*s[i,j,k]);
	    for (vector v in vectors)
	      foreach_block() {
		scalar vt = VT, vn = VN, vr = VR;
		v.x[] = (vt.boundary[id1](n0,n1,v.x,NULL) +
			 vt.boundary[id2](n0,n2,v.x,NULL) +
			 vn.boundary[id3](n0,n3,v.x,NULL) -
			 2.*v.x[i,j,k]);
		v.y[] = (vt.boundary[id1](n0,n1,v.y,NULL) +
			 vn.boundary[id2](n0,n2,v.y,NULL) +
			 vt.boundary[id3](n0,n3,v.y,NULL) -
			 2.*v.y[i,j,k]);
		v.z[] = (vn.boundary[id1](n0,n1,v.z,NULL) +
			 vr.boundary[id2](n0,n2,v.z,NULL) +
			 vr.boundary[id3](n0,n3,v.z,NULL) -
			 2.*v.z[i,j,k]);
	      }
	    return true;
	  }
#endif
  return false;
}

#if dimension > 1
foreach_dimension()
static Point tangential_neighbor_x (Point point, bool * zn)
{
  for (int k = 1; k <= BGHOSTS; k++)
    for (int j = -k; j <= k; j += 2*k) {
      if (is_neighbor(0,j) || is_neighbor(-1,j)) {
	*zn = false;
	return neighborp(0,j);
      }
#if dimension == 3
      // fixme: what about diagonals?
      if (is_neighbor(0,0,j) || is_neighbor(-1,0,j)) {
	*zn = true;
	return neighborp(0,0,j);
      }
#endif // dimension == 3
    }
  return (Point){.level = -1};
}
#endif // dimension > 1

static inline bool is_boundary_point (Point point) {
  return is_boundary (cell);
}
 
static void box_boundary_level (const Boundary * b, scalar * list, int l)
{
  disable_fpe_for_mpi();
  scalar * scalars = NULL;
  vector * vectors = NULL, * faces = NULL;
  for (scalar s in list)
    if (!is_constant(s) && s.refine != no_restriction) {
      if (s.v.x.i == s.i) {
	if (s.face)
	  faces = vectors_add (faces, s.v);
	else
	  vectors = vectors_add (vectors, s.v);
      }
      else if (s.v.x.i < 0 && s.boundary[0])
	scalars = list_add (scalars, s);
    }
  
  foreach_boundary_level (l) {
    if (!normal_neighbor (point, scalars, vectors) &&
	!diagonal_neighbor_2D (point, scalars, vectors) &&
	!diagonal_neighbor_3D (point, scalars, vectors)) {
      // no neighbors
      for (scalar s in scalars)
	foreach_block()
	  s[] = undefined;
      for (vector v in vectors)
	foreach_block()
	  foreach_dimension()
	    v.x[] = undefined;
    }
    if (faces) {
      int id = bid(cell);
      foreach_dimension()
	for (int i = -1; i <= 1; i += 2) {
	  // normal neighbor for faces
	  if (is_neighbor(i)) {
	    Point neighbor = neighborp(i);
	    for (vector v in faces) {
	      scalar vn = VN;
	      if (vn.boundary[id])
		foreach_block()
		  v.x[(i + 1)/2] = vn.boundary[id](neighbor, point, v.x, NULL);
	    }
	  }
#if dimension > 1
	  else if (i == -1) {
	    // tangential neighbor
	    bool zn;
	    Point neighbor = tangential_neighbor_x (point, &zn);
	    if (neighbor.level >= 0) {
	      int id = is_boundary_point (neighbor) ?
		bid(neighbor(-1)) : bid(cell);
	      for (vector v in faces) {
#if dimension == 2
		scalar vt = VT;
#else // dimension == 3
		scalar vt = zn ? VT : VR;
#endif
		foreach_block()
		  v.x[] = vt.boundary[id](neighbor, point, v.x, NULL);
	      }
	    }
	    else
	      // no neighbor
	      for (vector v in faces)
		foreach_block()
		  v.x[] = 0.;
	  }
#endif // dimension > 1
	}
    }
  }

  free (scalars);
  free (vectors);
  free (faces);
  enable_fpe_for_mpi();
}

#undef is_neighbor
 
@undef VN
@undef VT
@define VN _attribute[s.i].v.x
@define VT _attribute[s.i].v.y

static double masked_average (Point point, scalar s)
{
  double sum = 0., n = 0.;
  foreach_child()
    if (!is_boundary(cell) && s[] != nodata)
      sum += s[], n++;
  return n ? sum/n : nodata;
}

foreach_dimension()
static double masked_average_x (Point point, scalar s)
{
  double sum = 0., n = 0.;
  foreach_child()
    if (child.x < 0 && (!is_boundary(cell) || !is_boundary(neighbor(1))) &&
	s[1] != nodata)
      sum += s[1], n++;
  return n ? sum/n : nodata;
}

static void masked_boundary_restriction (const Boundary * b,
					 scalar * list, int l)
{
  scalar * scalars = NULL;
  vector * faces = NULL;
  for (scalar s in list)
    if (!is_constant(s) && s.refine != no_restriction) {
      if (s.v.x.i == s.i && s.face)
	faces = vectors_add (faces, s.v);
      else
	scalars = list_add (scalars, s);
    }
  
  foreach_halo (restriction, l) {
    for (scalar s in scalars)
      s[] = masked_average (parent, s);
    for (vector v in faces)
      foreach_dimension() {
	double average = masked_average_x (parent, v.x);
	if (is_boundary(neighbor(-1)))
	  v.x[] = average;
	if (is_boundary(neighbor(1)))
	  v.x[1] = average;
      }
  }

  free (scalars);
  free (faces);  
}

#define mask(func) {					\
  foreach_cell_post(!is_leaf(cell)) {			\
    if (is_leaf(cell)) {				\
      int bid = (func);					\
      if (bid >= 0)					\
	cell.pid = - bid - 1;				\
    }							\
    else { /* not a leaf */				\
      int pid = -1;					\
      foreach_child()					\
	if (cell.pid >= 0 || pid < 0)			\
	  pid = cell.pid;				\
      cell.pid = pid;					\
      if (pid < 0) {					\
	/* fixme: call coarsen_cell()? */		\
	cell.flags |= leaf;				\
	decrement_neighbors (point);			\
      }							\
    }							\
  }							\
  tree->dirty = true;					\
}

static void free_cache (CacheLevel * c)
{
  for (int l = 0; l <= depth(); l++)
    free (c[l].p);
  free (c);
}

void free_grid (void)
{
  if (!grid)
    return;
  free_boundaries();
  Tree * q = tree;
  free (q->leaves.p);
  free (q->faces.p);
  free (q->vertices.p);
  free (q->refined.p);
  /* low-level memory management */
  /* the root level is allocated differently */
  Layer * L = q->L[0];
  foreach_mem (L->m, L->len, 1) {
#if dimension == 1
    free (mem_data (L->m, point.i));
#elif dimension == 2
    free (mem_data (L->m, point.i, point.j));
#else // dimension == 3
    free (mem_data (L->m, point.i, point.j, point.k));
#endif // dimension == 3
  }
  for (int l = 0; l <= depth(); l++)
    destroy_layer (q->L[l]);
  q->L = &(q->L[-1]);
  free (q->L);
  free_cache (q->active);
  free_cache (q->prolongation);
  free_cache (q->boundary);
  free_cache (q->restriction);
  free (q);
  grid = NULL;
}

static void refine_level (int depth);

trace
void init_grid (int n)
{
  // check 64 bits structure alignment
  assert (sizeof(Cell) % 8 == 0);
  
  free_grid();
  int depth = 0;
  while (n > 1) {
    if (n % 2) {
      fprintf (stderr, "tree: N must be a power-of-two\n");
      exit (1);
    }
    n /= 2;
    depth++;
  }
  Tree * q = qcalloc (1, Tree);
  grid = (Grid *) q;
  grid->depth = 0;

  /* low-level memory management */
  q->L = qmalloc (2, Layer *);
  /* make sure we don't try to access level -1 */
  q->L[0] = NULL; q->L = &(q->L[1]);
  /* initialise the root cell */
  Layer * L = new_layer (0);
  q->L[0] = L;
#if dimension == 1
  for (int i = Period.x*GHOSTS; i < L->len - Period.x*GHOSTS; i++)
    assign_periodic (L->m, i, L->len, 
		     (char *) calloc (1, sizeof(Cell) + datasize));
  CELL(mem_data (L->m,GHOSTS)).flags |= leaf;
  if (pid() == 0)
    CELL(mem_data (L->m,GHOSTS)).flags |= active;
  for (int k = -GHOSTS*(1 - Period.x); k <= GHOSTS*(1 - Period.x); k++)
    CELL(mem_data (L->m,GHOSTS+k)).pid = (k < 0 ? -1 - left :
					  k > 0 ? -1 - right :
					  0);
  CELL(mem_data (L->m,GHOSTS)).pid = 0;
#elif dimension == 2
  for (int i = Period.x*GHOSTS; i < L->len - Period.x*GHOSTS; i++)
    for (int j = Period.y*GHOSTS; j < L->len - Period.y*GHOSTS; j++)
      assign_periodic (L->m, i, j, L->len,
		       (char *) calloc (1, sizeof(Cell) + datasize));
  CELL(mem_data (L->m,GHOSTS,GHOSTS)).flags |= leaf;
  if (pid() == 0)
    CELL(mem_data (L->m,GHOSTS,GHOSTS)).flags |= active;
  for (int k = - GHOSTS*(1 - Period.x); k <= GHOSTS*(1 - Period.x); k++)
    for (int l = -GHOSTS*(1 - Period.y); l <= GHOSTS*(1 - Period.y); l++)
      CELL(mem_data (L->m,GHOSTS+k,GHOSTS+l)).pid =
	(k < 0 ? -1 - left :
	 k > 0 ? -1 - right :
	 l > 0 ? -1 - top :
	 l < 0 ? -1 - bottom :
	 0);
  CELL(mem_data (L->m,GHOSTS,GHOSTS)).pid = 0;
#else // dimension == 3
  for (int i = Period.x*GHOSTS; i < L->len - Period.x*GHOSTS; i++)
    for (int j = Period.y*GHOSTS; j < L->len - Period.y*GHOSTS; j++)
      for (int k = Period.z*GHOSTS; k < L->len - Period.z*GHOSTS; k++)
	assign_periodic (L->m, i, j, k, L->len,
			 (char *) calloc (1, sizeof(Cell) + datasize));
  CELL(mem_data (L->m,GHOSTS,GHOSTS,GHOSTS)).flags |= leaf;
  if (pid() == 0)
    CELL(mem_data (L->m,GHOSTS,GHOSTS,GHOSTS)).flags |= active;
  for (int k = - GHOSTS*(1 - Period.x); k <= GHOSTS*(1 - Period.x); k++)
    for (int l = -GHOSTS*(1 - Period.y); l <= GHOSTS*(1 - Period.y); l++)
      for (int n = -GHOSTS*(1 - Period.z); n <= GHOSTS*(1 - Period.z); n++)
	CELL(mem_data (L->m,GHOSTS+k,GHOSTS+l,GHOSTS+n)).pid =
	  (k > 0 ? -1 - right :
	   k < 0 ? -1 - left :
	   l > 0 ? -1 - top :
	   l < 0 ? -1 - bottom :
	   n > 0 ? -1 - front :
	   n < 0 ? -1 - back :
	   0);
  CELL(mem_data (L->m,GHOSTS,GHOSTS,GHOSTS)).pid = 0;
#endif // dimension == 3
  q->active = qcalloc (1, CacheLevel);
  q->prolongation = qcalloc (1, CacheLevel);
  q->boundary = qcalloc (1, CacheLevel);
  q->restriction = qcalloc (1, CacheLevel);
  q->dirty = true;
  N = 1 << depth;
@if _MPI
  void mpi_boundary_new();
  mpi_boundary_new();
@endif
  // boundaries
  Boundary * b = qcalloc (1, Boundary);
  b->level = box_boundary_level;
  b->restriction = masked_boundary_restriction;
  add_boundary (b);
  refine_level (depth);
  reset (all, 0.);
  update_cache();
}

#if dimension == 2
void check_two_one (void)
{
  foreach_leaf()
    if (level > 0)
      for (int k = -1; k <= 1; k++)
	for (int l = -1; l <= 1; l++) {
	  /* fixme: all this mess is just to ignore ghost cells */
	  int i = (point.i + GHOSTS)/2 + k;
	  int j = (point.j + GHOSTS)/2 + l;
	  double x = ((i - GHOSTS + 0.5)*_DELTA*2. - 0.5);
	  double y = ((j - GHOSTS + 0.5)*_DELTA*2. - 0.5);
	  if (x > -0.5 && x < 0.5 && y > -0.5 && y < 0.5 && 
	      !(aparent(k,l).flags & active)) {
	    FILE * fp = fopen("check_two_one_loc", "w");
	    fprintf (fp,
		     "# %d %d\n"
		     "%g %g\n%g %g\n",
		     k, l,
		     ((_I + 0.5)*_DELTA - 0.5),
		     ((_J + 0.5)*_DELTA - 0.5),
		     x, y);
	    fclose (fp);
#if 0
	    fp = fopen("check_two_one", "w");
	    output_cells (fp);
	    fclose (fp);
#endif
	    assert (false);
	  }
	}
}
#endif

struct _locate { double x, y, z; };

Point locate (struct _locate p)
{
  for (int l = depth(); l >= 0; l--) {
    Point point = {0};
    point.level = l;
    int n = 1 << point.level;
    point.i = (p.x - X0)/L0*n + GHOSTS;
#if dimension >= 2
    point.j = (p.y - Y0)/L0*n + GHOSTS;
#endif
#if dimension >= 3
    point.k = (p.z - Z0)/L0*n + GHOSTS;
#endif
    if (point.i >= 0 && point.i < n + 2*GHOSTS
#if dimension >= 2
	&& point.j >= 0 && point.j < n + 2*GHOSTS
#endif
#if dimension >= 3
	&& point.k >= 0 && point.k < n + 2*GHOSTS
#endif
	) {
      if (allocated(0) && is_local(cell) && is_leaf(cell))
	return point;
    }
    else
      break;
  }
  Point point = {0};
  point.level = -1;
  return point;
}

// return true if the tree is "full" i.e. all the leaves are at the
// same level
bool tree_is_full()
{
  update_cache();
  return (grid->tn == 1L << grid->maxdepth*dimension);
}

#include "tree-common.h"

// overload the default periodic() function
void tree_periodic (int dir)
{
  int depth = grid ? depth() : -1;
  if (grid)
    free_grid();
  periodic (dir);
  if (depth >= 0)
    init_grid (1 << depth);
}
#define periodic(dir) tree_periodic(dir)
 
@if _MPI
#include "tree-mpi.h"
#include "balance.h"
@else // !_MPI
void mpi_boundary_refine  (scalar * list){}
void mpi_boundary_coarsen (int a, int b){}
void mpi_boundary_update  (scalar * list) {
  for (scalar s in list)
    s.dirty = true;
  boundary (list);
}
@endif // !_MPI

#endif
