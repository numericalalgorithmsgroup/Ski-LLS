/**
 * @file   Random.cpp
 * @author Xiangrui Meng <mengxr@stanford.edu>
 * @date   Tue Oct  4 23:44:20 2011
 * 
 * @brief  
 * 
 * 
 */ 

#include <cstdlib>
#include <cmath>

#include <pthread.h>

#include "Config.hpp"
#include "Vec.hpp"
#include "LinOp.hpp"
#include "Utils.hpp"
#include "Random.hpp"
#include <time.h>

#define SHR3 (jz=jsr, jsr^=(jsr<<13), jsr^=(jsr>>17), jsr^=(jsr<<5),jz+jsr)
#define UNI (.5 + (signed) SHR3*.2328306e-9)
  
static unsigned int kn[128];
static double       wn[128], fn[128];
  
static const double r = 3.442620;
  
static unsigned int seeds[MAX_NUM_THREADS];

static int table_created = 0;
static int seeds_initialized = 0;

void initialize_seeds()
{
  for( size_t i=0; i<MAX_NUM_THREADS; ++i )
    // seeds[i] = 78436*(i+1) + time(0);
    seeds[i] = 78436*(i+1) + time(0);
  seeds_initialized = 1;
}
  
void create_table()
{
  const double m1 = 2147483648.0;
  double dn = 3.442619855899,tn = dn, vn = 9.91256303526217e-3, q;
  
  /* Set up tables for RNOR */
  q       = vn/exp(-.5*dn*dn);
  kn[0]   = (dn/q)*m1;
  kn[1]   = 0;
  
  wn[0]   = q/m1;
  wn[127] = dn/m1;
  
  fn[0]   = 1.;
  fn[127] = exp(-.5*dn*dn);
  
  for( int i=126; i>=1; i-- )
  {
    dn      = sqrt(-2.*log(vn/dn+exp(-.5*dn*dn)));
    kn[i+1] = (dn/tn)*m1;
    tn      = dn;
    fn[i]   = exp(-.5*dn*dn);
    wn[i]   = dn/m1;
  }
    
  table_created = 1;
}

typedef struct
{
  unsigned int  seed;
  double       *begin;
  double       *end;
} randn_task;

void *randn_thread( void * );

void randn( ptrdiff_t n, double *a ) //probably generating random normal, in parallel
{
  if( !seeds_initialized )
    initialize_seeds();
    
  if( !table_created )
    create_table();

  int n_threads = get_num_cores();
    
  pthread_t  *threads = (pthread_t *)  malloc( n_threads * sizeof(pthread_t) );
  randn_task *tasks   = (randn_task *) malloc( n_threads * sizeof(randn_task) );

  for( int t=0; t<n_threads; ++t )
  {
    tasks[t].seed  = seeds[t];
    tasks[t].begin = a + (long) ceil(1.0*n* t   /n_threads);
    tasks[t].end   = a + (long) ceil(1.0*n*(t+1)/n_threads);
    pthread_create( &threads[t], NULL, &randn_thread, (void *) &tasks[t] );
  }
    
  for( int t=0; t<n_threads; ++t )
  {
    pthread_join( threads[t], NULL );
    seeds[t] = tasks[t].seed;
  }
    
  free(threads);
  free(tasks);
}

void *randn_thread( void *ptr )
{
  randn_task *task = (randn_task *) ptr;

  unsigned int jsr = task->seed;
  unsigned int jz;
  int          hz;
  unsigned int iz;
  double       x, y;

  for( double *cur = task->begin; cur != task->end; cur++ )
  {
    hz = SHR3;
    iz = hz & 127;
    if( fabs(hz) < kn[iz] )
    {
      *cur = (double) (hz*wn[iz]);
    }
    else
    {
      for(;;)
      {
        x=hz*wn[iz];      /* iz==0, handles the base strip */

        if(iz==0)
        {
          do
          {
            x=-log(UNI)*0.2904764; y=-log(UNI);
          }	/* .2904764 is 1/r */
          while (y+y<x*x);
          *cur = (double) ( (hz>0)? r+x : -r-x );
          break;
        }
              
        /* iz>0, handle the wedges of other strips */
        if ( fn[iz]+UNI*(fn[iz-1]-fn[iz]) < exp(-.5*x*x) )
        {
          *cur = (double) x;
          break;
        }
        /* initiate, try to exit for(;;) for loop*/              
        hz = SHR3;
        iz = hz & 127;
        if( fabs(hz) < kn[iz] )
        {
          *cur = (double) (hz*wn[iz]);
          break;
        }
      }
    }
  }

  task->seed = jsr;

  pthread_exit(NULL);
}

void randn( Vec_d v )
{
  if( v.inc() == 1)
    randn( v.n(), v.data() );
  else
  {
    Vec_d tmp = v.copy();
    randn( tmp.n(), tmp.data() );
    copy( tmp, v );
  }
}

void randn( Mat_d A )
{
  if( A.m() == A.ld() )
    randn( A.m()*A.n(), A.data() );
  else
  {
    for( ptrdiff_t j=0; j<A.n(); ++j )
      randn( A.m(), A.data()+j*A.ld() );
  }
}

