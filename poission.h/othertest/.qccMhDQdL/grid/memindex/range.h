#ifndef BASILISK_HEADER_31
#define BASILISK_HEADER_31
#line 1 "/home/dahuanghhc/basilisk/src/grid/memindex/range.h"
/**
# Multi-dimensional arrays with arbitrary bounds

This file implements low-level memory allocation for multi-dimensional
arrays with arbitrary starting and ending indices. This is used to
implement more memory-efficient tree structures.

The `Memalloc` data structure represents the one-dimensional array
which is defined by a generic pointer `p` and the size of its elements
`size` in bytes. 

The range of indices of the array is defined by the `Memrange` data
structure. */

typedef struct {
  void ** p;
  int size;
} Memalloc;

typedef struct {
  int start, end;
} Memrange;

/**
## Memory allocation for one-dimensional arrays

The `memrange_alloc()` function allocates the memory of the arrays
defined by `mem` and updates the corresponding index range `r` for a
(new) element of index `i`.

The `mem` argument is an array of `Memalloc` structures ending with a
NULL pointer for element `p`. */

void memrange_alloc (Memrange * r, Memalloc * mem, int i)
{
  if (r->start == r->end) {
    r->start = i;
    r->end = i + 1;
    for (Memalloc * m = mem; m->p; m++) {
      *m->p = calloc (1, m->size);
      *m->p = (char *)(*m->p) - i*m->size;
    }
  }
  else if (i >= r->end) {
    for (Memalloc * m = mem; m->p; m++) {
      *m->p = realloc ((char *)(*m->p) + r->start*m->size,
		       m->size*(i + 1 - r->start));
      *m->p = (char *)(*m->p) - r->start*m->size;
      memset ((char *)(*m->p) + r->end*m->size, 0, (i - r->end + 1)*m->size);
    }
    r->end = i + 1;
  }
  else if (i < r->start) {
    for (Memalloc * m = mem; m->p; m++) {
      *m->p = realloc ((char *)(*m->p) + r->start*m->size, m->size*(r->end - i));
      memmove ((char *)(*m->p) + (r->start - i)*m->size, *m->p,
	       m->size*(r->end - r->start));
      memset ((char *)(*m->p), 0, (r->start - i)*m->size);
      *m->p = (char *)(*m->p) - i*m->size;
    }
    r->start = i;
  }
}

/**
## Memory deallocation for one-dimensional arrays

The `memrange_free()` function takes the same arguments as
`memrange_alloc()` but performs the inverse operation. It returns
`true` if the arrays are empty after deallocation (i.e. the range of
indices is empty). */

bool memrange_free (Memrange * r, Memalloc * mem, int i)
{
  if (i == r->start) {
    if (i == r->end - 1) {
      for (Memalloc * m = mem; m->p; m++) {
	free ((char *)(*m->p) + r->start*m->size);
	*m->p = NULL;
      }
      r->start = r->end = 0;
      return true;
    }
    else {
      for (i = i + 1; i < r->end &&
	     !*(void **)((char *)(*mem->p) + i*mem->size); i++);
      for (Memalloc * m = mem; m->p; m++) {
	memmove ((char *)(*m->p) + r->start*m->size,
		 (char *)(*m->p) + i*m->size, m->size*(r->end - i));
	*m->p = realloc ((char *)(*m->p) + r->start*m->size,
			 m->size*(r->end - i));
	*m->p = (char *)(*m->p) - i*m->size;
      }
      r->start = i;
    }
  }
  else if (i == r->end - 1) {
    for (i = i - 1; i >= r->start &&
	   !*(void **)((char *)(*mem->p) + i*mem->size); i--);
    r->end = i + 1;
    for (Memalloc * m = mem; m->p; m++) {
      *m->p = realloc ((char *)(*m->p) + r->start*m->size,
		       m->size*(r->end - r->start));
      *m->p = (char *)(*m->p) - r->start*m->size;
    }
  }
  else {
    assert (i > r->start && i < r->end);
    for (Memalloc * m = mem; m->p; m++)
      memset ((char *)(*m->p) + i*m->size, 0, m->size);
  }
  return false;
}

/**
## Multidimensional arrays

The `Memindex` structure defines multi-dimensional arrays with
variable bounds. The interface is that used by [/src/grid/tree.h](). */

struct _Memindex {
  Memrange r1;
#if dimension >= 2
  Memrange * r2;
# if dimension >= 3
  Memrange ** r3;
# endif
#endif
#if dimension == 1
  char ** b;
#elif dimension == 2  
  char *** b;
#else
  char **** b;
#endif
};

#define Memindex struct _Memindex *

/**
The `mem_data()` macros return the data stored at a specific
(multidimensional) index. It assumes that the index is allocated. This
can be checked with `mem_allocated()`. */

#if dimension == 1
# define mem_allocated(m,i)			                        \
  ((i) >= (m)->r1.start && (i) < (m->r1.end) &&				\
   (m)->b[i])
# define mem_data(m,i) ((m)->b[i])
#elif dimension == 2
# define mem_allocated(m,i,j)						\
  ((i) >= (m)->r1.start && (i) < (m->r1.end) &&				\
   (m)->b[i] &&								\
   (j) >= (m)->r2[i].start && (j) < (m)->r2[i].end &&			\
   (m)->b[i][j])
# define mem_data(m,i,j) ((m)->b[i][j])
#else // dimension == 3
# define mem_allocated(m,i,j,k)						\
  ((i) >= (m)->r1.start && (i) < (m->r1.end) &&				\
   (m)->b[i] &&								\
   (j) >= (m)->r2[i].start && (j) < (m)->r2[i].end &&			\
   (m)->b[i][j] &&							\
   (k) >= (m)->r3[i][j].start && (k) < (m)->r3[i][j].end &&		\
   (m)->b[i][j][k])
# define mem_data(m,i,j,k) ((m)->b[i][j][k])
#endif // dimension == 3

/**
The `mem_new()` function returns a new (empty) `Memindex`. */

Memindex mem_new (int len)
{
  Memindex m = calloc (1, sizeof (struct _Memindex));
  return m;
}

/**
The `mem_destroy()` function frees all the memory allocated by a given
`Memindex`. */

void mem_destroy (Memindex m, int len)
{
#if dimension >= 2
  for (int i = m->r1.start; i < m->r1.end; i++)
    if (m->b[i]) {
#if dimension >= 3
      for (int j = m->r2[i].start; j < m->r2[i].end; j++)
	if (m->b[i][j])
	  free (m->b[i][j] + m->r3[i][j].start);
      free (m->r3[i] + m->r2[i].start);      
#endif // dimension >= 3
      free (m->b[i] + m->r2[i].start);
    }
  if (m->b) {
    free (m->r2 + m->r1.start);
#if dimension >= 3
    free (m->r3 + m->r1.start);
#endif
  }
#endif // dimension >= 2
  if (m->b)
    free (m->b + m->r1.start);
  free (m);
}

/**
The `mem_assign()` function assigns a (pointer) value to a given index. */

#if dimension == 1
void mem_assign (Memindex m, int i, int len, void * b)
{
  Memalloc mem[] = {{(void **)&m->b, sizeof(char *)},
		    {NULL}};
  memrange_alloc (&m->r1, mem, i);
  mem_data(m,i) = b;
}
#elif dimension == 2
void mem_assign (Memindex m, int i, int j, int len, void * b)
{
  Memalloc mem[] = {{(void **)&m->b,  sizeof(char **)},
		    {(void **)&m->r2, sizeof(Memrange)},
		    {NULL}};
  memrange_alloc (&m->r1, mem, i);
  Memalloc mem1[] = {{(void **)&m->b[i], sizeof(char *)},
		     {NULL}};
  memrange_alloc (&m->r2[i], mem1, j);
  mem_data(m,i,j) = b;
}
#else // dimension == 3
void mem_assign (Memindex m, int i, int j, int k, int len, void * b)
{
  Memalloc mem[] = {{(void **)&m->b,  sizeof(char ***)},
		    {(void **)&m->r2, sizeof(Memrange)},
		    {(void **)&m->r3, sizeof(Memrange *)},
		    {NULL}};
  memrange_alloc (&m->r1, mem, i);
  Memalloc mem1[] = {{(void **)&m->b[i], sizeof(char **)},
		     {(void **)&m->r3[i], sizeof(Memrange)},
		     {NULL}};
  memrange_alloc (&m->r2[i], mem1, j);
  Memalloc mem2[] = {{(void **)&m->b[i][j], sizeof(char *)},
		     {NULL}};
  memrange_alloc (&m->r3[i][j], mem2, k);
  mem_data(m,i,j,k) = b;
}
#endif // dimension == 3

/**
The `mem_free()` function frees a given index. */

#if dimension == 1
void mem_free (Memindex m, int i, int len)
{
  Memalloc mem[] = {{(void **)&m->b, sizeof(char *)},
		    {NULL}};
  memrange_free (&m->r1, mem, i);
}
#elif dimension == 2
void mem_free (Memindex m, int i, int j, int len)
{
  Memalloc mem[] = {{(void **)&m->b[i], sizeof(char *)},
		    {NULL}};
  if (memrange_free (&m->r2[i], mem, j)) {
    Memalloc mem[] = {{(void **)&m->b, sizeof(char **)},
		      {(void **)&m->r2, sizeof(Memrange)},
		      {NULL}};
    memrange_free (&m->r1, mem, i);
  }
}
#else // dimension == 3
void mem_free (Memindex m, int i, int j, int k, int len)
{
  Memalloc mem[] = {{(void **)&m->b[i][j], sizeof(char *)},		    
		    {NULL}};
  if (memrange_free (&m->r3[i][j], mem, k)) {
    Memalloc mem[] = {{(void **)&m->b[i], sizeof(char **)},
		      {(void **)&m->r3[i], sizeof(Memrange)},
		      {NULL}};
    if (memrange_free (&m->r2[i], mem, j)) {
      Memalloc mem[] = {{(void **)&m->b, sizeof(char ***)},
			{(void **)&m->r2, sizeof(Memrange)},
			{(void **)&m->r3, sizeof(Memrange *)},
			{NULL}};
      memrange_free (&m->r1, mem, i);      
    }
  }
}
#endif // dimension == 3

/**
The `foreach_mem()` macro traverses every `_i` allocated elements of
array `_m` taking into account a periodicity of `_len` (and ghost
cells). */

#if dimension == 1
@def foreach_mem(_m, _len, _i) {
  Point point = {0};
  for (point.i = max(Period.x*GHOSTS, (_m)->r1.start);
       point.i < min(_len - Period.x*GHOSTS, (_m)->r1.end);
       point.i += _i)
    if ((_m)->b[point.i]) {
@
@define end_foreach_mem() }}
#elif dimension == 2
@def foreach_mem(_m, _len, _i) {
  Point point = {0};
  for (point.i = max(Period.x*GHOSTS, (_m)->r1.start);
       point.i < min(_len - Period.x*GHOSTS, (_m)->r1.end);
       point.i += _i)
    if ((_m)->b[point.i])
      for (point.j = max(Period.y*GHOSTS, (_m)->r2[point.i].start);
	   point.j < min(_len - Period.y*GHOSTS, (_m)->r2[point.i].end);
	   point.j += _i)
	if ((_m)->b[point.i][point.j]) {
@
@define end_foreach_mem() }}
#else // dimension == 3
@def foreach_mem(_m, _len, _i) {
  Point point = {0};
  for (point.i = max(Period.x*GHOSTS, (_m)->r1.start);
       point.i < min(_len - Period.x*GHOSTS, (_m)->r1.end);
       point.i += _i)
    if ((_m)->b[point.i])
      for (point.j = max(Period.y*GHOSTS, (_m)->r2[point.i].start);
	   point.j < min(_len - Period.y*GHOSTS, (_m)->r2[point.i].end);
	   point.j += _i)
	if ((_m)->b[point.i][point.j])
	  for (point.k = max(Period.z*GHOSTS, (_m)->r3[point.i][point.j].start);
	       point.k < min(_len - Period.z*GHOSTS, (_m)->r3[point.i][point.j].end);
	       point.k += _i)
	    if ((_m)->b[point.i][point.j][point.k]) {
@
@define end_foreach_mem() }}
#endif // dimension == 3

#endif
