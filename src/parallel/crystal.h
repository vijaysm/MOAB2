#ifndef CRYSTAL_H
#define CRYSTAL_H

/* requires <mpi.h>, "types.h", and "errmem.h" */
#if !defined(TYPES_H) || !defined(ERRMEM_H)
#warning "crystal.h" requires "types.h" and "errmem.h"
#endif

#ifdef USE_MPI

typedef struct { uint n; buffer buf; } crystal_buf;

typedef struct crystal_data {
  crystal_buf buffers[3];
  crystal_buf *all, *keep, *send;
  MPI_Comm comm;
  uint num, id;
} crystal_data;

void moab_crystal_init(crystal_data *, MPI_Comm);
void moab_crystal_free(crystal_data *);
void moab_crystal_router(crystal_data *);

#endif

#endif
