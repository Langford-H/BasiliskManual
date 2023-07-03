#ifndef BASILISK_HEADER_32
#define BASILISK_HEADER_32
#line 1 "/home/dahuanghhc/basilisk/src/grid/mempool.h"
/* A memory pool implementation for fixed-size blocks 
   It uses Simple Segregated Storage, see e.g:
   http://www.boost.org/doc/libs/1_55_0/libs/pool/doc/html/boost_pool/pool/pooling.html
*/

typedef struct _Pool Pool;

struct _Pool {
  Pool * next; // next pool
};

typedef struct {
  char * first, * lastb; // first and last free blocks
  size_t size;           // block size
  size_t poolsize;       // pool size
  Pool * pool, * last;   // first and last pools
} Mempool;

typedef struct {
  char * next;
} FreeBlock;

Mempool * mempool_new (size_t poolsize, size_t size)
{
  // check for 64 bytes alignment
  assert (poolsize % 8 == 0);
  assert (size >= sizeof(FreeBlock));
  // to get the effective pool size, we cap this amount to 2^20 = 1MB
  // i.e. something comparable to the size of a L2 cache
  poolsize = min(1 << 20, poolsize + sizeof(Pool));
  Mempool * m = qcalloc (1, Mempool);
  m->poolsize = poolsize;
  m->size = size;
  return m;
}

void mempool_destroy (Mempool * m)
{
  Pool * p = m->pool;
  while (p) {
    Pool * next = p->next;
    free (p);
    p = next;
  }
  free (m);
}

void * mempool_alloc (Mempool * m)
{
  if (!m->first) {
    // allocate new pool
    Pool * p = (Pool *) malloc (m->poolsize);
    p->next = NULL;
    if (m->last)
      m->last->next = p;
    else
      m->pool = p;
    m->last = p;
    m->first = m->lastb = ((char *)m->last) + sizeof(Pool);
    FreeBlock * b = (FreeBlock *) m->first;
    b->next = NULL;
  }
  void * ret = m->first;
  FreeBlock * b = (FreeBlock *) ret;
  char * next = b->next;
  if (!next) {
    m->lastb += m->size;
    next = m->lastb;
    if (next + m->size > ((char *) m->last) + m->poolsize)
      next = NULL;
    else {
      FreeBlock * b = (FreeBlock *) next;
      b->next = NULL;
    }
  }
  m->first = next;
@if TRASH
  double * v = (double *) ret;
  for (int i = 0; i < m->size/sizeof(double); i++)
    v[i] = undefined;
@endif
  return ret;
}

void * mempool_alloc0 (Mempool * m)
{
  void * ret = mempool_alloc (m);
  memset (ret, 0, m->size);
  return ret;
}

void mempool_free (Mempool * m, void * p)
{
@if TRASH
  double * v = (double *) p;
  for (int i = 0; i < m->size/sizeof(double); i++)
    v[i] = undefined;
@endif
  FreeBlock * b = (FreeBlock *) p;
  b->next = m->first;
  m->first = (char *) p;
}

#endif
