
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>

#include "braid.h"
//#include "/mnt/c/Users/mariu/OneDrive/Dokumente/Masterarbeit/Opendihu/opendihu/dependencies/xbraid/src/xbraid-2.3.0/braid/braid.hpp"
// #include "braid_test.h"
#include "PinT_IE_lib.c"
// #include "PinT_IE_lib.cpp"

#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include <petscdraw.h>
#include <petscvec.h>

#include <opendihu.h>

// #include "time_stepping_scheme/implicit_euler.h"

/*--------------------------------------------------------------------------
 * My App and Vector structures
 *--------------------------------------------------------------------------*/
/* can put anything in my app and name it anything as well */
typedef struct _braid_App_struct
{
  MPI_Comm  comm;
  double    tstart;       /* Define the temporal domain */
  double    tstop;
  int       ntime;
  double    xstart;       /* Define the spatial domain */
  double    xstop;
  int       nspace;
  double    matrix[3];    /* the three point spatial discretization stencil */
  double *  g;            /* temporary vector for inversions and mat-vecs */
  double *  sc_info;      /* Runtime information that tracks the space-time grids visited */
  int       print_level;  /* Level of output desired by user (see the -help message below) */

  typedef  TimeSteppingScheme::ImplicitEuler<
    SpatialDiscretization::FiniteElementMethod<
      Mesh::StructuredRegularFixedOfDimension<1>,
      BasisFunction::LagrangeOfOrder<>,
      Quadrature::None,
      Equation::Dynamic::IsotropicDiffusion
    >
  > NestedSolver;

  std::vector<
    std::shared_ptr<
      NestedSolver
    >
  > *implicitEulerSolvers;   //< vector of nested solvers (implicit euler) for solution on different grids

} my_App;

/* Can put anything in my vector and name it anything as well */
typedef struct _braid_Vector_struct
{
   PetscInt     size;
   PetscReal *values;

} my_Vector;

/* create and allocate a vector */
void
create_vector(my_Vector **u,
              int size)
{
   (*u) = (my_Vector *) malloc(sizeof(my_Vector));
   ((*u)->size)   = size;
   ((*u)->values) = (double *) malloc(size*sizeof(double));
}


/*--------------------------------------------------------------------------
 * My integration routines
 *--------------------------------------------------------------------------*/

int my_Step(braid_App        app,
            braid_Vector     ustop,
            braid_Vector     fstop,
            braid_Vector     u,
            braid_StepStatus status)
{
   PetscReal tstart;             /* current time */
   PetscReal tstop;              /* evolve to this time*/
   PetscInt level, i, solver;
   PetscReal deltaX, deltaT;

   // get level and start and stop time from braid status
   braid_StepStatusGetLevel(status, &level);
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);

   // determine, which solver is needed (depends on size)
   solver=log2(u->size - 1);

   // printf("%d\n", level);
   // printf("%d\n", u->size);
   // printf("%d\n", n);

   /* XBraid forcing */

   if(fstop != NULL)
   {
      for(i = 0; i < u->size; i++)
      {
         u->values[i] = u->values[i] + fstop->values[i];
      }
   }

   Vec V; //create vector V with the given values of u
   VecCreateSeqWithArray(PETSC_COMM_SELF, 1, u->size, u->values,&V);
   VecAssemblyBegin(V);
   VecAssemblyEnd(V);

   // Set the initial guess for the solver
   (*app->implicitEulerSolvers)[solver]->data().solution()->valuesGlobal()=V;

   // set time span for the solver, which is calculated by braid status
   (*app->implicitEulerSolvers)[solver]->setTimeSpan(tstart, tstop);

   // run solver
   (*app->implicitEulerSolvers)[solver]->run();

   // functions for debugging:
   // VecView((*app->implicitEulerSolvers)[0]->data().solution()->valuesGlobal(), 	PETSC_VIEWER_STDOUT_SELF);
   // PetscRealView(6, u->values, 0);

   // put the calculated solution into braid vector u
   VecGetArray(&(*((*app->implicitEulerSolvers)[solver]->data().solution()->valuesGlobal())),&(u->values));

   deltaT = tstop - tstart;
   deltaX = (app->xstop - app->xstart) / (ustop->size - 1.0);

   /* Store info on space-time grids visited during the simulation */
   (app->sc_info)[ (2*level) ] = deltaX;
   (app->sc_info)[ (2*level) + 1] = deltaT;

   /* no refinement */
   braid_StepStatusSetRFactor(status, 1);

   return 0;
}


int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
   my_Vector *u;
   int    i;
   int nspace = (app->nspace);
   double deltaX = (app->xstop - app->xstart) / (nspace - 1.0);

   /* Allocate vector */
   create_vector(&u, nspace);

   /* Initialize vector (with correct boundary conditions) */
   if(t == 0.0)
   {
      /* Get the solution at time t=0 */
      get_solution(u->values, u->size, 0.0, app->xstart, deltaX);
   }
   else
   {
      /* Use random values for u(t>0), this measures asymptotic convergence rate */
      for(i=0; i < nspace; i++)
      {
         (u->values)[i] = ((double)braid_Rand())/braid_RAND_MAX;
      }
   }
   // (u->values)[0] =2;
   // (u->values)[1] =2;
   // (u->values)[2] =4;
   // (u->values)[3] =5;
   // (u->values)[4] =2;
   // (u->values)[5] =2;
   *u_ptr = u;

   return 0;
}

int
my_Clone(braid_App     app,
         braid_Vector  u,
         braid_Vector *v_ptr)
{
   my_Vector *v;
   int size = (u->size);
   int i;

   create_vector(&v, size);
   for (i = 0; i < size; i++)
   {
      (v->values)[i] = (u->values)[i];
   }
   *v_ptr = v;

   return 0;
}

int
my_Free(braid_App    app,
        braid_Vector u)
{
   free(u->values);
   free(u);

   return 0;
}

int
my_Sum(braid_App     app,
       double        alpha,
       braid_Vector  x,
       double        beta,
       braid_Vector  y)
{
   int i;
   int size = (y->size);

   for (i = 0; i < size; i++)
   {
      (y->values)[i] = alpha*(x->values)[i] + beta*(y->values)[i];
   }

   return 0;
}


int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{
   int    i;
   int size   = (u->size);
   double dot = 0.0;

   for (i = 0; i < size; i++)
   {
      dot += (u->values)[i]*(u->values)[i];
   }
   *norm_ptr = sqrt(dot);

   return 0;
}

int
my_Access(braid_App          app,
          braid_Vector       u,
          braid_AccessStatus astatus)
{
   int        index, rank, level, done;
   char       filename[255];
   double     t, error;

   braid_AccessStatusGetT(astatus, &t);
   braid_AccessStatusGetTIndex(astatus, &index);
   braid_AccessStatusGetLevel(astatus, &level);
   braid_AccessStatusGetDone(astatus, &done);

   /* Print solution to file if simulation is over */
   if(done)
   {
      MPI_Comm_rank( (app->comm), &rank);
      sprintf(filename, "%s.%07d.%05d", "PinT_diffusion.out", index, rank);
      save_solution(filename, u->values, u->size, app->xstart,
            app->xstop, app->ntime, app->tstart, app->tstop);
   }

   /* IF on the finest level AND print_level is high enough AND at the final time,
    * THEN print out the discretization error */
   if( (level == 0) && ((app->print_level) > 0) && (index == app->ntime) )
   {
      error = compute_error_norm(u->values, app->xstart, app->xstop, u->size, t);
      printf("  Discretization error at final time:  %1.4e\n", error);
      fflush(stdout);
   }

   return 0;
}

int
my_BufSize(braid_App           app,
           int                 *size_ptr,
           braid_BufferStatus  bstatus)
{
   int size = (app->nspace);
   *size_ptr = (size+1)*sizeof(double);
   return 0;
}

int
my_BufPack(braid_App           app,
           braid_Vector        u,
           void               *buffer,
           braid_BufferStatus  bstatus)
{
   double *dbuffer = (double *) buffer;
   int i, size = (u->size);

   dbuffer[0] = size;
   for (i = 0; i < size; i++)
   {
      dbuffer[i+1] = (u->values)[i];
   }

   braid_BufferStatusSetSize( bstatus,  (size+1)*sizeof(double));

   return 0;
}

int
my_BufUnpack(braid_App           app,
             void               *buffer,
             braid_Vector       *u_ptr,
             braid_BufferStatus  bstatus)
{
   my_Vector *u = NULL;
   double    *dbuffer = (double *) buffer;
   int        i, size;

   size = dbuffer[0];
   create_vector(&u, size);

   for (i = 0; i < size; i++)
   {
      (u->values)[i] = dbuffer[i+1];
   }
   *u_ptr = u;

   return 0;
}


/*--------------------------------------------------------------------------
 * my_Residual, my_Coarsen and my_Refine are advanced XBraid options, ignore
 * them until you understand the rest of the driver.
 *--------------------------------------------------------------------------*/

int
my_Residual(braid_App        app,
            braid_Vector     ustop,
            braid_Vector     r,
            braid_StepStatus status)
{
   double tstart;             /* current time */
   double tstop;              /* evolve to this time*/
   int i;
   double x, deltaX, deltaT;

   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   deltaT = tstop - tstart;
   deltaX = (app->xstop - app->xstart) / (ustop->size - 1.0);

   /* Set up matrix stencil for 1D heat equation*/
   compute_stencil(deltaX, deltaT, app->matrix);

   /* Residual r = A*xstop - r - forcing - boundary
    *   note: there are no boundary terms here */
   matvec_tridiag(ustop->values, app->g, ustop->size, app->matrix);
   x = app->xstart;
   for(i = 0; i < r->size; i++)
   {
      r->values[i] = app->g[i] - r->values[i] - deltaT*forcing(tstop, x);
      x = x + deltaX;
   }

   return 0;
}

/* Bilinear Coarsening */
int
my_Coarsen(braid_App              app,
           braid_Vector           fu,
           braid_Vector          *cu_ptr,
           braid_CoarsenRefStatus status)
{

   int csize, level;
   my_Vector *v;

   /* This returns the level for fu */
   braid_CoarsenRefStatusGetLevel(status, &level);

   /* Stop spatial coarsening/interpolation after reaching a grid of 3 points.
    * This is the smallest spatial grid (2 boundary points, one true DOF). */
   if( level < floor(log2(app->nspace)) - 1 )
   {
      csize = ((fu->size) - 1)/2 + 1;
      create_vector(&v, csize);
      coarsen_1D(v->values, fu->values, csize, fu->size);
   }
   else
   {
      /* No coarsening, clone the vector */
      my_Clone(app, fu, &v);
   }

   *cu_ptr = v;

   return 0;
}

/* Bilinear interpolation */
int
my_Interp(braid_App              app,
          braid_Vector           cu,
          braid_Vector          *fu_ptr,
          braid_CoarsenRefStatus status)
{

   int fsize, level;
   my_Vector *v;

   /* This returns the level for fu_ptr */
   braid_CoarsenRefStatusGetLevel(status, &level);

   /* Stop spatial coarsening/interpolation after reaching a grid of 3 points.
    * This is the smallest spatial grid (2 boundary points, one true DOF). */
   if( level < floor(log2(app->nspace)) - 1 )
   {
      fsize = (cu->size - 1)*2 + 1;
      create_vector(&v, fsize);
      interpolate_1D(cu->values, v->values, cu->size, fsize);
   }
   else
   {
      /* No refinement, clone the vector */
      my_Clone(app, cu, &v);
   }

   *fu_ptr = v;

   return 0;
}
