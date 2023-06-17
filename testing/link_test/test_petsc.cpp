
#include <petsc.h>

static char help[] = "Hello world program.\n\n";

int main(int argc, char** argv)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank;

  ierr = PetscInitialize(&argc,&argv,NULL,help);
  CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "Hello world from %d\n", rank);
  CHKERRQ(ierr);
  ierr = PetscFinalize();

  return ierr;
}