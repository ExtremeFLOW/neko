/*
 * p4est_wrap.c
 * Wrapper for p4est and sc libraries.
 *
 *  Created on: Aug 19, 2022
 *  Author: Adam Peplinski
 */

// libsc speciffic definitions
#define SC_ENABLE_MPI
#define SC_ENABLE_MPIIO

// p4est speciffic definitions
// This is dimension rleated; for now 3D only
#define N_DIM 3

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <sc.h>

#if N_DIM == 2
#undef P4_TO_P8
#else
#include <p4est_to_p8est.h>
#endif

#ifndef P4_TO_P8
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_ghost.h>
#include <p4est_nodes.h>
#include <p4est_vtk.h>
#include <p4est_lnodes.h>
#include <p4est_mesh.h>
#include <p4est_iterate.h>
#else
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <p8est_nodes.h>
#include <p8est_vtk.h>
#include <p8est_lnodes.h>
#include <p8est_mesh.h>
#include <p8est_iterate.h>
#endif

#include "./p4est_wrap.h"

/* Global variables
 *  notice I use tree_neko->user_pointer to store ghost quadrants data
 */
static p4est_connectivity_t *connect_neko = NULL; /**< connectivity structure */
static p4est_t *tree_neko = NULL; /**< tree structure */
static p4est_mesh_t *mesh_neko = NULL; /**< mesh structure */
static p4est_ghost_t *ghost_neko = NULL; /**< ghost zone structure */
static p4est_nodes_t *nodes_neko = NULL; /**< vertex numbering structure */
static p4est_lnodes_t *lnodes_neko = NULL; /**< GLL node numbering structure */
static p4est_t *tree_neko_compare = NULL; /**< tree structure to perform comparison */

/* Wrappers */
/* initialise/finalise */
void wp4est_init(MPI_Fint fmpicomm,
		 int catch_signals, int print_backtrace,
		 int log_threshold)
{
  MPI_Comm mpicomm;
  mpicomm = MPI_Comm_f2c(fmpicomm);
  sc_init (mpicomm, catch_signals, print_backtrace,
  	   NULL, log_threshold);
  p4est_init(NULL, log_threshold);
}

void wp4est_finalize(int log_priority)
{
  sc_package_print_summary (log_priority);
  sc_finalize ();
}

/* Connectivity; only required routines */
void wp4est_cnn_del() {
  if (connect_neko) p4est_connectivity_destroy(connect_neko);
  connect_neko =  NULL;
}

void wp4est_cnn_valid(int * is_valid) {
  *is_valid = p4est_connectivity_is_valid(connect_neko);
}

/* tree_ management */
void wp4est_tree_del() {
  if (tree_neko) p4est_destroy(tree_neko);
  tree_neko = NULL;
}

void wp4est_tree_valid(int * is_valid) {
  *is_valid = p4est_is_valid(tree_neko);
}

void wp4est_tree_save(int save_data, char filename[]) {
  p4est_save_ext(filename, tree_neko, save_data,0);
}

void wp4est_tree_load(MPI_Fint fmpicomm, int load_data, char filename[]) {
  MPI_Comm mpicomm;
  mpicomm = MPI_Comm_f2c(fmpicomm);
  tree_neko = p4est_load_ext(filename, mpicomm, sizeof(user_data_t), load_data,1,0,
			     NULL, &connect_neko);
  tree_neko->user_pointer = NULL;
}

/* tree and grid info */
void wp4est_ghost_new() {
  if (ghost_neko) p4est_ghost_destroy(ghost_neko);
  ghost_neko = p4est_ghost_new(tree_neko, P4EST_CONNECT_FULL);
  /* ghost data */
  if (tree_neko->user_pointer) P4EST_FREE(tree_neko->user_pointer);
  tree_neko->user_pointer = (void *) P4EST_ALLOC (user_data_t, ghost_neko->ghosts.elem_count);
  p4est_ghost_exchange_data (tree_neko, ghost_neko, (user_data_t *) tree_neko->user_pointer);
}

void wp4est_ghost_del() {
  if (ghost_neko) p4est_ghost_destroy(ghost_neko);
  ghost_neko = NULL;
  /* ghost data */
  if (tree_neko->user_pointer) P4EST_FREE(tree_neko->user_pointer);
  tree_neko->user_pointer = NULL;
}

void wp4est_mesh_new() {
  if (mesh_neko) p4est_mesh_destroy(mesh_neko);
  mesh_neko = p4est_mesh_new(tree_neko, ghost_neko, P4EST_CONNECT_FULL);
}

void wp4est_mesh_del() {
  if (mesh_neko) p4est_mesh_destroy(mesh_neko);
  mesh_neko = NULL;
}

void wp4est_nodes_new() {
  if (nodes_neko) p4est_nodes_destroy(nodes_neko);
  nodes_neko = p4est_nodes_new(tree_neko, ghost_neko);
}

void wp4est_nodes_del() {
  if (nodes_neko) p4est_nodes_destroy(nodes_neko);
  nodes_neko = NULL;
}

void wp4est_lnodes_new(int degree) {
  if (lnodes_neko) p4est_lnodes_destroy(lnodes_neko);
  lnodes_neko = p4est_lnodes_new(tree_neko, ghost_neko, degree);
}

void wp4est_lnodes_del() {
  if (lnodes_neko) p4est_lnodes_destroy(lnodes_neko);
  lnodes_neko = NULL;
}

/* p4est internal load balance */
void wp4est_part(int partforcoarsen) {
  p4est_partition(tree_neko, partforcoarsen, NULL);
}

/** @brief Iterate over faces to correct Neko boundary conditions
 *
 * @details Required by wp4est_bc_check
 *
 * @param info
 * @param user_data
 */
void iter_bc_chk(p4est_iter_face_info_t * info, void *user_data) {
  user_data_t *data;

  int mpirank = info->p4est->mpirank;
  p4est_gloidx_t gfirst_quad =  info->p4est->global_first_quadrant[mpirank];
  /* ghost data */
  user_data_t *ghost_data = (user_data_t *) info->p4est->user_pointer;

  sc_array_t *sides = &(info->sides);
  p4est_iter_face_side_t *side;
  int nside = (int) sides->elem_count;


  int il, jl, kl; // loop index
  int iwl; // global quad number
  int face; // quad face
  int imsh[2]; // mesh type mark

  if (info->tree_boundary) {
    /* Face is on the outside of the tree; it can be any type of boundary
     * condition. V-type type mesh can have external bc inside the mesh if
     * a neighbor is T-type quad.
     */
    P4EST_ASSERT (nside <= 2);
    if (nside == 1) {
      /* external face; no E (0), P (-1) bc */
      side = p4est_iter_fside_array_index_int (sides, 0);
      face = (int) side->face;
      if (side->is_hanging) {
	/* hanging face */
	for (jl = 0; jl < P4EST_HALF; jl++) {
	  if (!side->is.hanging.is_ghost[jl]) {
	    /* local node */
	    data = (user_data_t *) side->is.hanging.quad[jl]->p.user_data;
	    /* check connection type (should not be E or P) */
	    if (data->bc[face] == 0 || data->bc[face] == -1) {
	      iwl = (int) (gfirst_quad + side->treeid + side->is.hanging.quadid[jl]);
	      printf("Connectivity error: %i %i %i %i\n",
		     tree_neko->mpirank, iwl, face, data->bc[face]);
	      printf("external face marked as internal.\n");
	      SC_ABORT("Aborting: iter_bc_chk\n");
	    }
	  } // is ghost
	} // children loop

      } else {
	if (!side->is.full.is_ghost) {
	  /* local node */
	  data = (user_data_t *) side->is.full.quad->p.user_data;
	  /* check connection type (should not be E or P) */
	  if (data->bc[face] == 0 || data->bc[face] == -1) {
	    iwl = (int) (gfirst_quad + side->treeid + side->is.full.quadid);
	    printf("Connectivity error: %i %i %i %i\n",
		   tree_neko->mpirank, iwl, face, data->bc[face]);
	    printf("external face marked as internal.\n");
	    SC_ABORT("Aborting: iter_bc_chk\n");
	  }
	} // is ghost
      } // hanging

    } else {
      /* internal face; any face type possible */
      /* collect quad type;
       * for different values of imsh V-type elements should point to external bc
       */
      imsh[0] = 0;
      imsh[1] = 0;
      for (il = 0; il < nside; ++il) {
	side = p4est_iter_fside_array_index_int (sides, il);
	if (side->is_hanging) {
	  /* imsh for all children is the same */
	  if (side->is.hanging.is_ghost[0]) {
	    data = (user_data_t *) &ghost_data[side->is.hanging.quadid[0]];
	  } else {
	    data = (user_data_t *) side->is.hanging.quad[0]->p.user_data;
	  }
	} else {
	  if (side->is.full.is_ghost) {
	    data = (user_data_t *) &ghost_data[side->is.full.quadid];
	  } else {
	    data = (user_data_t *) side->is.full.quad->p.user_data;
	  }
	}
	imsh[il] = data->imsh;
      }
      /* test bc */
      for (il = 0; il < nside; ++il) {
	side = p4est_iter_fside_array_index_int (sides, il);
	face = (int) side->face;
	if (side->is_hanging) {
	  /* hanging face */
	  for (jl = 0; jl < P4EST_HALF; jl++) {
	    if (!side->is.hanging.is_ghost[jl]) {
	      /* local node */
	      data = (user_data_t *) side->is.hanging.quad[jl]->p.user_data;
	      /* velocity bc for V-type element neighbor to T-type element
	       * should not be E or P
	       */
	      if (imsh[0] != imsh[1]) {
		if (data->bc[face] == 0 || data->bc[face] == -1) {
		  iwl = (int) (gfirst_quad + side->treeid + side->is.hanging.quadid[jl]);
		  printf("Connectivity error: %i %i %i %i\n",
			 tree_neko->mpirank, iwl, face, data->bc[face]);
		  printf("velocity external face marked as internal.\n");
		  SC_ABORT("Aborting: iter_bc_chk\n");
		}
	      } else {
		/* not velocity or not V-T meshes boundary - all
		 * internal elements
		 */
		if (data->bc[face] != 0 && data->bc[face] != -1) {
		  iwl = (int) (gfirst_quad + side->treeid + side->is.hanging.quadid[jl]);
		  printf("Connectivity error: %i %i %i %i\n",
			 tree_neko->mpirank, iwl, face,data->bc[face]);
		  printf("internal face marked as external.\n");
		  SC_ABORT("Aborting: iter_bc_chk\n");
		}
	      }
	    } // is ghost
	  } // children loop

	} else {
	  if (!side->is.full.is_ghost) {
	    /* local node */
	    data = (user_data_t *) side->is.full.quad->p.user_data;
	    /* velocity bc for V-type element neighbor to T-type element
	     * should not be E or P
	     */
	    if (imsh[0] != imsh[1]) {
	      if (data->bc[face] == 0|| data->bc[face] == -1) {
		iwl = (int) (gfirst_quad + side->treeid + side->is.full.quadid);
		printf("Connectivity error: %i %i %i %i\n",
		       tree_neko->mpirank, iwl, face, data->bc[face]);
		printf("velocity external face marked as internal.\n");
		SC_ABORT("Aborting: iter_bc_chk\n");
	      }
	    } else {
	      /* not velocity or not V-T meshes boundary - all
	       * internal elements
	       */
	      if (data->bc[face] != 0 && data->bc[face] != -1) {
		iwl = (int) (gfirst_quad + side->treeid + side->is.full.quadid);
		printf("Connectivity error: %i %i %i %i\n",
		       tree_neko->mpirank, iwl, face, data->bc[face]);
		printf("internal face marked as external.\n");
		SC_ABORT("Aborting: iter_bc_chk\n");
	      }
	    }
	  } // is ghost
	} // hanging
      } // side loop
    } // forest internal/external face
  } else {
    /* face is on the interior of the tree; all faces are E (0)
     * for V-type mesh external and periodic boundary can be defined on
     * tree faces only
     */
    P4EST_ASSERT (nside == 2);
    for (il = 0; il < nside; ++il) {
      side = p4est_iter_fside_array_index_int (sides, il);
      face = (int) side->face;
      if (side->is_hanging) {
	/* hanging face */
	for (jl = 0; jl < P4EST_HALF; jl++) {
	  if (!side->is.hanging.is_ghost[jl]) {
	    // local node
	    data = (user_data_t *) side->is.hanging.quad[jl]->p.user_data;
	    if (data->bc[face] != 0 ) {
	      iwl = (int) (gfirst_quad + side->treeid + side->is.hanging.quadid[jl]);
	      printf("Connectivity error: %i %i %i %i\n",
		     tree_neko->mpirank, iwl, face, data->bc[face]);
	      printf("tree internal face marked as external.\n");
	      SC_ABORT("Aborting: iter_bc_chk\n");
	    }
	  } // is ghost
	} // children loop

      } else {
	if (!side->is.full.is_ghost) {
	  /* local node */
	  data = (user_data_t *) side->is.full.quad->p.user_data;
	  if (data->bc[face] != 0) {
	    iwl = (int) (gfirst_quad + side->treeid + side->is.full.quadid);
	    printf("Connectivity error: %i %i %i %i\n",
		   tree_neko->mpirank, iwl, face, data->bc[face]);
	    printf("tree internal face marked as external.\n");
	    SC_ABORT("Aborting: iter_bc_chk\n");
	  }
	} // is ghost
      } // hanging
    } // side loop
  } // tree internal/external

}

/* Check boundary conditions */
void wp4est_bc_check() {
#ifdef P4_TO_P8
  p4est_iterate(tree_neko, ghost_neko, NULL, NULL, iter_bc_chk, NULL, NULL);
#else
  p4est_iterate(tree_neko, ghost_neko, NULL, NULL, iter_bc_chk, NULL);
#endif
}

/* routines for data exchange between Neko and p4est */

/** @brief Iterate over elements to count V-mesh elements
 *
 * @details Required by wp4est_msh_get_size
 *
 * @param info
 * @param user_data
 */
void count_mshv(p4est_iter_volume_info_t * info, void *user_data) {
  int loc_level;
  user_data_t *data = (user_data_t *) info->quad->p.user_data;
  int *lmax = (int *) user_data;

  // coult V-type elements
  if (data->imsh == 0) {
    lmax[0] = lmax[0] + 1;
  }
  // find max local level
  loc_level = (int) info->quad->level;
  lmax[1] = (loc_level > lmax[1] ? loc_level : lmax[1]);
}

/* get mesh size */
void wp4est_msh_get_size(int * mdim, int64_t * nelgt, int64_t * nelgto,
			 int32_t * nelt, int * nelv, int * maxl) {
  int lmax[2];
  // mesh dimension
  *mdim = (int) P4EST_DIM;
  // get global number of quadrants
  *nelgt = tree_neko->global_num_quadrants;
  // zero based global position of local quadrants
  *nelgto = tree_neko->global_first_quadrant[tree_neko->mpirank];
  // number of local quadrants
  *nelt = tree_neko->local_num_quadrants;

  // count number of V-mesh elements and find current max level
  lmax[0] = 0;
  lmax[1] = 0;

#ifdef P4_TO_P8
  p4est_iterate(tree_neko, ghost_neko,(void *) &lmax, count_mshv,
		NULL, NULL, NULL);
#else
  p4est_iterate(tree_neko, ghost_neko,(void *) &lmax, count_mshv,
		NULL, NULL);
#endif

  *nelv = lmax[0];
  *maxl = lmax[1];
}

/* get node list size */
void wp4est_nds_get_size(int * nowin, int * nowsh, int * oowin,
			 int * nin, int * nhf, int * nhe) {
  if (nodes_neko) {
    // number of owned independent nodes
    *nowin = (int) nodes_neko->num_owned_indeps;
    // number of owned shared
    *nowsh = (int) nodes_neko->num_owned_shared;
    // position of the first independent owned node
    *oowin = (int) nodes_neko->offset_owned_indeps;
    // numbers of local independent and hanging nodes (face and edge)
    *nin = (int) nodes_neko->indep_nodes.elem_count;
    *nhf = (int) nodes_neko->face_hangings.elem_count;
#ifdef P4_TO_P8
    *nhe = (int) nodes_neko->edge_hangings.elem_count;
#endif
  } else {
    SC_ABORT("nodes_neko not allocated; aborting: wp4est_nds_get_size\n");
  }
}

/* independent node list */
void wp4est_nds_get_ind(int64_t * nglid, int * nown, double * ncoord) {
  int il, jl;
  int64_t id;
  p4est_indep_t *node;
  double vxyz[3];
  if (nodes_neko) {
    sc_array_t  *indep_nodes = &(nodes_neko->indep_nodes);
    // numbers of local independent nodes
    const int ni = (int) nodes_neko->indep_nodes.elem_count;
    // independent node offset
    const int oi = (int) nodes_neko->offset_owned_indeps;
    // loop over nodes owned by other mpi rank
    for(il=0;il<oi;++il){
      // extract node
      node = (p4est_indep_t *) sc_array_index (indep_nodes,il);
      // get global id
      id = (int64_t) node->p.piggy3.local_num + 1; // conversion to fortran numberring
      for(jl=0;jl<nodes_neko->nonlocal_ranks[il];++jl) {
	id += (int64_t) nodes_neko->global_owned_indeps[jl];
      }
      nglid[il] = id;
      // node owner
      nown[il] = (int) nodes_neko->nonlocal_ranks[il];
      // get physical coordinates
      p4est_qcoord_to_vertex (connect_neko, node->p.piggy3.which_tree,
			      node->x, node->y,
#ifdef P4_TO_P8
			      node->z,
#endif
			      vxyz);
      // copy coordinates
      for(jl=0;jl<N_DIM;++jl){
	ncoord[il*N_DIM+jl] = vxyz[jl];
      }
    }
    // loop over nodes owned by this mpi rank
    // local id offset
    id = (int64_t) 1; // conversion to fortran numberring
    for(jl=0;jl<tree_neko->mpirank;++jl) {
      id += (int64_t) nodes_neko->global_owned_indeps[jl];
    }
    for(il=oi;il<ni;++il){
      // extract node
      node = (p4est_indep_t *) sc_array_index (indep_nodes,il);
      // get global id
      nglid[il] = id + (int64_t) node->p.piggy3.local_num;
      // node owner
      nown[il] = (int) tree_neko->mpirank;
      // get physical coordinates
      p4est_qcoord_to_vertex (connect_neko, node->p.piggy3.which_tree,
			      node->x, node->y,
#ifdef P4_TO_P8
			      node->z,
#endif
			      vxyz);
      // copy coordinates
      for(jl=0;jl<N_DIM;++jl){
	ncoord[il*N_DIM+jl] = vxyz[jl];
      }
    }
  } else {
    SC_ABORT("nodes_neko not allocated; aborting: wp4est_nds_get_ind\n");
  }
}

/* face hanging node list */
void wp4est_nds_get_hfc(int * depend, double * ncoord) {
  int il, jl;
#ifdef P4_TO_P8
  const int ndep = 4;
  p8est_hang4_t *node;
#else
  const int ndep = 2;
  p4est_hang2_t *node;
#endif
  double vxyz[3];

  if (nodes_neko) {
    sc_array_t  *nodes = &(nodes_neko->face_hangings);
    // numbers of local face hanging nodes
    const int ni = (int) nodes_neko->face_hangings.elem_count;
    // loop over nodes
    for(il=0;il<ni;++il){
      // extract node
#ifdef P4_TO_P8
      node = (p8est_hang4_t *) sc_array_index (nodes,il);
#else
      node = (p4est_hang2_t *) sc_array_index (nodes,il);
#endif
      // get independent nodes mapping
      for(jl=0;jl<ndep;++jl){
        depend[il*ndep+jl] = (int) node->p.piggy.depends[jl] +1; // conversion to fortran numberring
      }
      // get physical coordinates
      p4est_qcoord_to_vertex (connect_neko, node->p.piggy.which_tree,
                              node->x, node->y,
#ifdef P4_TO_P8
                              node->z,
#endif
                              vxyz);
      // copy coordinates
      for(jl=0;jl<N_DIM;++jl){
        ncoord[il*N_DIM+jl] = vxyz[jl];
      }
    }
  }else {
    SC_ABORT("nodes_neko not allocated; aborting: wp4est_nds_get_hfc\n");
  }
}


/* edge hanging node list */
void wp4est_nds_get_hed(int * depend, double * ncoord) {
  int il, jl;
#ifdef P4_TO_P8
  const int ndep = 2;
  p8est_hang2_t *node;
  double vxyz[3];

  if (nodes_neko) {
    sc_array_t  *nodes = &(nodes_neko->edge_hangings);
    // numbers of local face hanging nodes
    const int ni = (int) nodes_neko->edge_hangings.elem_count;
    // loop over nodes
    for(il=0;il<ni;++il){
      // extract node
      node = (p8est_hang2_t *) sc_array_index (nodes,il);
      // get independent nodes mapping
      for(jl=0;jl<ndep;++jl){
        depend[il*ndep+jl] = (int) node->p.piggy.depends[jl] + 1; // conversion to fortran numberring
      }
      // get physical coordinates
      p4est_qcoord_to_vertex (connect_neko, node->p.piggy.which_tree,
                              node->x, node->y,
                              node->z,
                              vxyz);
      // copy coordinates
      for(jl=0;jl<N_DIM;++jl){
        ncoord[il*N_DIM+jl] = vxyz[jl];
      }
    }
  }else {
    SC_ABORT("nodes_neko not allocated; aborting: wp4est_nds_get_hed\n");
  }
#endif
}

/* get vertex to node mapping */
void wp4est_nds_get_vmap(int * vmap) {
  int il, jl;

  if (nodes_neko) {
    // number of local elements
    const int vi = (int) nodes_neko->num_local_quadrants;
    // quad to vertex local map
    for (il = 0; il < vi; ++il) {
      for (jl = 0; jl < P4EST_CHILDREN; ++jl) {
	vmap[il * P4EST_CHILDREN + jl] = (int) nodes_neko->local_nodes[il * P4EST_CHILDREN + jl] + 1; // conversion to fortran numberring
      }
    }
  }else {
    SC_ABORT("nodes_neko not allocated; aborting: wp4est_nds_get_vnmap\n");
  }
}


/* data type for mesh data transfer between neko and p4est*/
typedef struct transfer_data_s {
  int64_t *gidx; /**< pointer to global element index array */
  int *level; /**< pointer to element level array */
  int *igrp; /**< pointer to element group array */
  int *crv; /**< pointer to face projection array */
  int *bc; /**< pointer to boundary condition array */
  double *coord; /**< pointer to approximate vertex coordinates array */
  int *falg; /**< face alignment */
} transfer_data_t;

/** @brief Iterate over element volumes to transfer approximate coordinates of the element vertices
 *
 * @details Required by wp4est_elm_get_dat
 *
 * @param info
 * @param user_data
 */
void iter_datav(p4est_iter_volume_info_t * info, void *user_data) {
  user_data_t *data = (user_data_t *) info->quad->p.user_data;
  transfer_data_t *trans_data = (transfer_data_t *) user_data;

  // which quad (local and global element number)
  p4est_tree_t *tree;
  p4est_locidx_t iwl;
  p4est_gloidx_t iwlt, iwg;
  int iwli, il, jl;
  p4est_quadrant_t node;
  double vxyz[3];

  // get quad number
  tree = p4est_tree_array_index(info->p4est->trees, info->treeid);
  // local quad number
  iwl = info->quadid + tree->quadrants_offset;
  iwli = (int) iwl;
  // global quad number
  iwg = (p4est_gloidx_t) info->p4est->global_first_quadrant[info->p4est->mpirank] + iwl;

  // global quadrant index
  trans_data->gidx[iwli] = iwg + 1; // conversion to fortran numberring

  // quadrant level
  trans_data->level[iwli] = (int) info->quad->level;

  // element group mark
  trans_data->igrp[iwli] = data->igrp;

  // curvature and boundary data
  for (il = 0; il < P4EST_FACES; il++) {
    trans_data->crv[iwli*P4EST_FACES+il] = data->crv[il];
    trans_data->bc[iwli*P4EST_FACES+il] = data->bc[il];
  }

  // get corner coordinates
  for (il=0;il < P4EST_CHILDREN; ++il){
    p4est_quadrant_corner_node (info->quad, il, &node);
    p4est_qcoord_to_vertex (info->p4est->connectivity, info->treeid,
			    node.x, node.y,
#ifdef P4_TO_P8
			    node.z,
#endif
			    vxyz);
    // copy coordinates
    for(jl=0;jl<N_DIM;++jl){
      trans_data->coord[(iwli*P4EST_CHILDREN+il)*N_DIM+jl] = vxyz[jl];
    }
  }
}

/** @brief Iterate over faces to get alignment
 *
 * @details Required by wp4est_elm_get_dat
 *
 * @param info
 * @param user_data
 */
void iter_algf(p4est_iter_face_info_t * info, void *user_data) {
  // which quad (local element number)
  p4est_tree_t *tree;
  p4est_locidx_t iwl;
  int iwlt;
  // orientation and side info
  int orient = (int) info->orientation;
  sc_array_t *sides = &(info->sides);
  p4est_iter_face_side_t *side;
  int nside = (int) sides->elem_count;
  transfer_data_t *trans_data = (transfer_data_t *) user_data;
  //int *fcs_arr = (int *) user_data;
  int il, jl;
  int iref, pref, pset;
  int8_t face[2], ftmp;

  /* compare orientation of different trees*/
  /* find the reference side; lowest face number */
  P4EST_ASSERT (nside <= 2);
  for (il = 0; il < nside; ++il) {
    side = p4est_iter_fside_array_index_int (sides, il);
    face[il] = side->face;
  }
  if (nside == 1) {
    /* single side no alignment */
    iref = nside +1;
  } else {
    /* 2 sides; find permutation set */
    if (face[0]<face[1]) {
      iref = 0;
#ifdef P4_TO_P8
      pref = p8est_face_permutation_refs[face[0]][face[1]];
      pset = p8est_face_permutation_sets[pref][orient];
#else
      pset = orient;
#endif
    } else {
      iref = 1;
#ifdef P4_TO_P8
      pref = p8est_face_permutation_refs[face[1]][face[0]];
      pset = p8est_face_permutation_sets[pref][orient];
#else
      pset = orient;
#endif
    }
  }

  for (il = 0; il < nside; ++il) {
    side = p4est_iter_fside_array_index_int (sides, il);
    if (il == iref) {
      orient = pset;
    } else {
      orient = 0;
    }
    if (side->is_hanging) {
      /* hanging face */
      for (jl = 0; jl < P4EST_HALF; jl++) {
	if (!side->is.hanging.is_ghost[jl]){
	  // local node
	  tree = p4est_tree_array_index(info->p4est->trees, side->treeid);
	  // local quad number
	  iwl =  side->is.hanging.quadid[jl] + tree->quadrants_offset;
	  iwlt = (int) iwl;
	  iwlt = iwlt*P4EST_FACES + (int) side->face;
	  trans_data->falg[iwlt] = orient;
	}
      }
    } else {
      if (!side->is.full.is_ghost) {
	// local node
	tree = p4est_tree_array_index(info->p4est->trees, side->treeid);
	// local quad number
	iwl =  side->is.full.quadid + tree->quadrants_offset;
	iwlt = (int) iwl;
	iwlt = iwlt*P4EST_FACES + (int) side->face;
	trans_data->falg[iwlt] = orient;
      }
    }
  }
}

/* get element info */
void wp4est_elm_get_dat(int64_t * gidx, int * level, int * igrp, int * crv,
			int * bc,double * coord, int * falg) {
  transfer_data_t transfer_data;

  transfer_data.gidx = gidx;
  transfer_data.igrp = igrp;
  transfer_data.level = level;
  transfer_data.crv = crv;
  transfer_data.bc = bc;
  transfer_data.coord = coord;
  transfer_data.falg = falg;
#ifdef P4_TO_P8
  p4est_iterate(tree_neko, ghost_neko,(void *) &transfer_data, iter_datav,
		iter_algf, NULL, NULL);
#else
  p4est_iterate(tree_neko, ghost_neko,(void *) &transfer_data, iter_datav,
		iter_algf, NULL);
#endif
}

/* degrees of freedom numbering */
void wp4est_elm_get_lnode(int * lnnum, int * lnown, int64_t * lnoff, int * lnodes) {

  int il, jl, kl;
  if (lnodes_neko) {
    const int vnd = (int) lnodes_neko->vnodes;
    const int owned = (int) lnodes_neko->owned_count;
    const int local = (int) lnodes_neko->num_local_nodes;
    const p4est_gloidx_t offset = lnodes_neko->global_offset;
    *lnnum = (int) local;
    *lnown = (int) owned;
    *lnoff = (int64_t) offset;
    for (il = 0; il < lnodes_neko->num_local_elements*vnd; ++il) {
      lnodes[il] = (int) lnodes_neko->element_nodes[il] + 1; // conversion to fortran numberring
    }
  } else {
    SC_ABORT("lnodes_nek not allocated; aborting: wp4est_elem_get_lnode\n");
  }
}

/* get shares list size */
void wp4est_sharers_get_size(int * nrank, int * nshare) {
  int il, jl;
  p4est_lnodes_rank_t *lnode;
  if (lnodes_neko) {
    // number of ranks in a sharers array
    *nrank = (int) lnodes_neko->sharers->elem_count;
    sc_array_t  *sharers = lnodes_neko->sharers;
    jl = 0;
    // loop over mpi rank sharing the node
    for(il=0;il < lnodes_neko->sharers->elem_count;++il){
      lnode = (p4est_lnodes_rank_t *) sc_array_index (sharers,il);
      // number of nodes shared with given rank
      jl += (int) lnode->shared_nodes.elem_count;
      /*      printf("Sharers: %i %i %i %i %i %i %i %i\n",
	     tree_neko->mpirank, il, lnode->rank, (int) lnode->shared_nodes.elem_count,
	     lnode->shared_mine_offset, lnode->shared_mine_count,
	     lnode->owned_offset, lnode->owned_count);*/
    }
    *nshare = jl;
  } else {
    SC_ABORT("lnodes_neko not allocated; aborting: wp4est_sharers_get_size\n");
  }
}

/* independent lnode list */
void wp4est_sharers_get_ind(int64_t * nglid, int * lrank, int * loff, int * lshare) {
  int il, jl;
  p4est_lnodes_rank_t *lnode;
  p4est_locidx_t *snode;
  if (lnodes_neko) {
    const int owned = (int) lnodes_neko->owned_count;
    const int local = (int) lnodes_neko->num_local_nodes;
    const p4est_gloidx_t offset = lnodes_neko->global_offset;
    // get global index list of local nodes
    // local owned
    for (il = 0; il < owned ; ++il) {
      nglid[il] = (p4est_gloidx_t) 1 + offset + il;
    }
    // local not owned
    for (il = owned; il < local ; ++il) {
      nglid[il] = (p4est_gloidx_t) 1 + lnodes_neko->nonlocal_nodes[il-owned]; // conversion to fortran numberring
    }
    // get node sharing info
    sc_array_t  *sharers = lnodes_neko->sharers;
    // loop over mpi rankssharing the node
    loff[0] = 1;
    for(il=0;il < lnodes_neko->sharers->elem_count;++il){
      lnode = (p4est_lnodes_rank_t *) sc_array_index (sharers,il);
      lrank[il] = (int) lnode->rank;
      loff[il+1] = loff[il] + (int) lnode->shared_nodes.elem_count;
      sc_array_t  *node_list = &(lnode->shared_nodes);
      for(jl=0;jl < lnode->shared_nodes.elem_count;++jl){
	snode = (p4est_locidx_t *) sc_array_index(node_list,jl);
	lshare[loff[il]-1+jl] = 1 + *snode; // conversion to fortran numberring
      } 
    }
  } else {
    SC_ABORT("lnodes_neko not allocated; aborting: wp4est_sharers_get_ind\n");
  }
}


/* get hanging element/face/edge information */
void wp4est_hang_get_info(int * hang_elm, int * hang_fsc, int * hang_edg) {

  int il, jl;
  int hanging_face[P4EST_FACES];
#ifdef P4_TO_P8
  int hanging_edge[P8EST_EDGES];
#endif
  if (lnodes_neko) {
    for (il = 0; il < lnodes_neko->num_local_elements; ++il) {
      for (jl = 0; jl < P4EST_FACES; ++jl) {
	hanging_face[jl] = -1;
      }
#ifdef P4_TO_P8
      for (jl = 0; jl < P8EST_EDGES; ++jl) {
	hanging_edge[jl] = -1;
      }
      hang_elm[il] = p4est_lnodes_decode(lnodes_neko->face_code[il],
					 hanging_face, hanging_edge);
      for (jl = 0; jl < P4EST_FACES; ++jl) {
	hang_fsc[il * P4EST_FACES + jl] = hanging_face[jl];
      }
      for (jl = 0; jl < P8EST_EDGES; ++jl) {
	hang_edg[il * P8EST_EDGES + jl] = hanging_edge[jl];
      }
#else
      hang_elm[il] = p4est_lnodes_decode(lnodes_neko->face_code[il],hanging_face);
      for (jl=0;jl<P4EST_FACES;++jl) {
	hang_fsc[il*P4EST_FACES +jl] = hanging_face[jl];
      }
#endif
    }
  } else {
    SC_ABORT("lnodes_neko not allocated; aborting: wp4est_msh_get_tplg\n");
  }
}

/* extract family information */
void wp4est_fml_get_info(int64_t * family, int * nelf) {
  p4est_topidx_t jt;
  p4est_tree_t *tree;
  sc_array_t *tquadrants;
  int isfamily;
  size_t zz, incount, window;
  p4est_quadrant_t *ci[P4EST_CHILDREN];
  int64_t igs, its, iqs;
  int iqf;

  // initialise number of family quads
  iqf = 0;

  // global quadrants shift
  igs = (int64_t) tree_neko->global_first_quadrant[tree_neko->mpirank] + 1; // conversion to fortran numberring

  /* loop over all local trees */
  for (jt = tree_neko->first_local_tree; jt <= tree_neko->last_local_tree; ++jt) {
    tree = p4est_tree_array_index (tree_neko->trees, jt);
    tquadrants = &tree->quadrants;

    window = 0;  // start position of sliding window in array
    incount = tquadrants->elem_count;  // number of quadrants
    its = (int64_t) tree->quadrants_offset; // local quadrant shift
    while (window + P4EST_CHILDREN <= incount) {
      isfamily = 1;
      for (zz = 0; zz < P4EST_CHILDREN; ++zz) {
	ci[zz] = p4est_quadrant_array_index (tquadrants, window + zz);
	if (zz != (size_t) p4est_quadrant_child_id (ci[zz])) {
	  isfamily = 0;
	  break;
	}
      }
      P4EST_ASSERT (!isfamily || p4est_quadrant_is_familypv (ci));
      if (isfamily) {
	// get global element number in the family
	for (zz = 0; zz < P4EST_CHILDREN; ++zz) {
	  iqs = (int64_t) window + zz;
	  family[iqf] = igs + its + iqs;
	  ++iqf;
	}
	window += P4EST_CHILDREN;
      }
      else {
	++window;
      }
    }
  }

  // copy number of family quads
  *nelf = iqf;
}


/* refinement, coarsening, balancing */

/** Mark element for refinement
 *
 * @details Required by p4est_refine_ext
 *
 * @param p4est
 * @param which_tree
 * @param quadrant
 * @return refinement mark
 */
int ref_mark_f (p4est_t * p4est, p4est_topidx_t which_tree,
		p4est_quadrant_t * quadrant)
{
  user_data_t *data = (user_data_t *) quadrant->p.user_data;

  if (data->ref_mark == AMR_RM_H_REF) {
    return 1;
  } else {
    return 0;
  }
}

/** Mark element for coarsening
 *
 *  @details Required by p4est_coarsen_ext
 *
 * @param p4est
 * @param which_tree
 * @param quadrants
 * @return coarsening mark
 */
int crs_mark_f (p4est_t * p4est, p4est_topidx_t which_tree,
		p4est_quadrant_t * quadrants[])
{
  user_data_t *data;
  int iwt, id;

  /* check if all the children are coarsened */
  iwt = 1;
  for(id=0;id<P4EST_CHILDREN;++id){
    /* get refinement mark */
    data = (user_data_t *) quadrants[id]->p.user_data;
    if (data->ref_mark != AMR_RM_H_CRS) {
      iwt = 0;
      break;
    }
  }
  return iwt;
}

/** Refine/coarsen quadrant
 *
 * @details Required by p4est_refine_ext, p4est_coarsen_ext,
 *  p4est_balance_ext. Performs operations on quad related
 *  user data.
 *
 * @param p4est          forest
 * @param which_tree     tree number
 * @param num_outgoing   number of quadrants that are being replaced
 * @param outgoing       replaced quads
 * @param num_incoming   number of quadrants that are being added
 * @param incoming       added quads
 */
void  quad_replace (p4est_t * p4est, p4est_topidx_t which_tree,
		    int num_outgoing, p4est_quadrant_t * outgoing[],
		    int num_incoming, p4est_quadrant_t * incoming[]) {

  int id, il; // loop index
  int ic; // children refinement status
  int id_ch[P4EST_CHILDREN]; // child id
  user_data_t * parent, * child;
  p4est_quadrant_t tmp_q;

  if (num_outgoing == P4EST_CHILDREN && num_incoming == 1) {
    /* this is coarsening; I assume coarsening is not performed at balancing stage
     * and there are no recursive coarsening */
    /* get child id */
    for(id=0;id<P4EST_CHILDREN;id++){
      id_ch[id] = p4est_quadrant_child_id (outgoing[id]);
    }

    /* check consistency of refinement history; no previous actions on children*/
    ic = 0;
    for(id=0;id<P4EST_CHILDREN;id++){
      child = (user_data_t *) outgoing[id]->p.user_data;
      if (child->parent_gln != -1 || child->el_gln == -1) ic = 1;
    }
    if (ic) SC_ABORT("Recursive refinement/coarsening; aborting: quad_replace; coarsen\n");

    /* Merge element data.
     * I'm lazy here and assume refinement just copy curvature and
     * BC data from children at external faces. It means I expect
     * face value for neighbouring elements to be the same, so
     * it is enough to take values of even and odd faces from
     * first and last element respectively.
     * Notice p4est does not give information about ordering of
     * outgoing quads.
     */
    parent = (user_data_t *) incoming[0]->p.user_data;
    /* fill refinement history data */
    /* no parent */
    parent->parent_gln = -1;
    parent->parent_ln = -1;
    parent->parent_nid = -1;
    /* I was coarsened */
    parent->el_gln = -1;
    parent->el_ln = -1;
    parent->el_nid = -1;
    /* reset refine mark */
    parent->ref_mark = AMR_RM_NONE;

    for(id=0;id<P4EST_CHILDREN;id++){
      child = (user_data_t *) outgoing[id]->p.user_data;
      if (id_ch[id] == 0) {
	/* first element; copy element type data */
	parent->imsh = child->imsh;
	parent->igrp = child->igrp;
	/* copy faces 0, 2 and 4 */
	for(il=0;il<P4EST_FACES;il=il+2){
	  /* curvature and boundary flag */
	  parent->crv[il] = child->crv[il];
	  parent->bc[il] = child->bc[il];
	}
      } else if(id_ch[id] == (P4EST_CHILDREN-1)) {
	/* last element; copy faces 1, 3 and 5 */
	for(il=1;il<P4EST_FACES;il=il+2){
	  /* curvature and boundary flag */
	  parent->crv[il] = child->crv[il];
	  parent->bc[il] = child->bc[il];
	}
      }
      /* copy coarsening data */
      parent->children_gln[id_ch[id]] = child->el_gln;
      parent->children_ln[id_ch[id]] = child->el_ln;
      parent->children_nid[id_ch[id]] = child->el_nid;
    }

#ifdef DEBUG
    /* for testing */
    user_data_t * child0, * child1;
    for(id=0;id<P4EST_CHILDREN;++id){
      child = (user_data_t *) outgoing[id]->p.user_data;
      printf("crs chidlren %i %i %i\n",parent->children_gln[id_ch[id]],
	     parent->children_ln[id_ch[id]],parent->children_nid[id_ch[id]]);
      if (id_ch[id] == 0)
	child0 = (user_data_t *) outgoing[id]->p.user_data;
      if (id_ch[id] == (P4EST_CHILDREN-1))
	child1 = (user_data_t *) outgoing[id]->p.user_data;
    }
    printf("crs parent gln %i %i %i %i\n",parent->children_gln[0],parent->children_gln[1],
	   parent->children_gln[2],parent->children_gln[3]);
    printf("crs parent ln %i %i %i %i\n",parent->children_ln[0],parent->children_ln[1],
	   parent->children_ln[2],parent->children_ln[3]);
    printf("crs parent nid %i %i %i %i\n",parent->children_nid[0],parent->children_nid[1],
	   parent->children_nid[2],parent->children_nid[3]);
    for(il=0;il<P4EST_FACES;++il){
      printf("crs CRV BC %i %i %i %i %i %i %i \n",il,
	     parent->crv[il],parent->bc[il],
	     child0->crv[il],child0->bc[il],
	     child1->crv[il],child1->bc[il]);
    }
#endif

  } else if (num_outgoing == 1 && num_incoming == P4EST_CHILDREN) {
    /* this is refinement; this part can be executed at refinement and balancing stage;
     * Even though I do not allow for recursive refinement, there is possibility of refining
     * the previously coarsened element and I have to keep track of refinement history */
    /* get child id */
    for(id=0;id<P4EST_CHILDREN;id++){
      id_ch[id] = p4est_quadrant_child_id (incoming[id]);
    }

    /* check refinement status*/
    parent = (user_data_t *) outgoing[0]->p.user_data;
    if (parent->parent_gln != -1){
      /* refined element; no recursive refinement allowed */
      SC_ABORT("Recursive refinement; aborting: quad_replace; refine\n");
    }
    ic = -1;
    for (il=0; il < P4EST_CHILDREN; il++) {
      if(parent->children_gln[il] != -1) ic = 1;
    }

    for (id=0;id<P4EST_CHILDREN;++id){
      child = (user_data_t *) incoming[id]->p.user_data;
      /* copy data from the parent (filled in refine_fn) */
      child->imsh = parent->imsh;
      child->igrp = parent->igrp;
      //memcpy (child, parent, p4est->data_size);

      /* reset refine mark
       * be careful not to coarsen quads just refined
       */
      child->ref_mark = AMR_RM_NONE;

      /* update refinement history data */
      if (ic == -1 && parent->el_gln != -1){
	/* unchanged element; neither coarsened nor refined;
	 * proceed with refining */

	/* store parent number on neko side */
	child->parent_gln = parent->el_gln;
	child->parent_ln = parent->el_ln;
	child->parent_nid = parent->el_nid;
	/* who am I */
	child->el_gln = id_ch[id];
	child->el_ln = -1;
	child->el_nid = -1;

	/* reset coarsening data */
	for(il=0;il<P4EST_CHILDREN;++il){
	  child->children_gln[il] = -1;
	  child->children_ln[il] = -1;
	  child->children_nid[il] = -1;
	}
      } else if (ic == 1 && parent->el_gln == -1){
	/* coarsened element; put information back so no coarsening/refinement
	 * operation would be performed on neko variables */

	/* store parent number on neko side */
	child->parent_gln = -1;
	child->parent_ln = -1;
	child->parent_nid = -1;
	/* restore old global element number */
	child->el_gln = parent->children_gln[id_ch[id]];
	child->el_ln = parent->children_ln[id_ch[id]];
	child->el_nid = parent->children_nid[id_ch[id]];

	/* reset children data */
	for(il=0;il<P4EST_CHILDREN;++il){
	  child->children_gln[il] = -1;
	  child->children_ln[il] = -1;
	  child->children_nid[il] = -1;
	}
      } else {
	/* wrong option */
	SC_ABORT("Wrong refinement history; aborting: quad_replace; refine\n");
      }

      /* go across all faces */
      for(il=0;il<P4EST_FACES;++il){
	/* test if the face is an external one
	 * (with respect to tree, not the mesh)
	 * find face neighbour
	 */
	p4est_quadrant_face_neighbor (incoming[id], il, &tmp_q);
	/* does neighbour belong to the same tree */
	if(p4est_quadrant_is_inside_root (&tmp_q)){
	  /* internal face*/
	  child->crv[il] = 0;
	  child->bc[il] = 0;
	}
	else{
	  /* external face
	   * copy curvature and bc data */
	  child->crv[il] = parent->crv[il];
	  child->bc[il] = parent->bc[il];
	}

      }
#ifdef DEBUG
      /*for testing */
      printf("ref %i %i\n",child->el_gln,child->parent_gln);
      for(il=0;il<P4EST_FACES;++il){
	printf("BC %i %i %i %i %i \n",il,child->crv[il],child->bc[il],
	       parent->crv[il],parent->bc[il]);
      }
#endif
    }

  } else {
    /* something is wrong */
    SC_ABORT("Wrong in/out quad number; aborting: quad_replace\n");
  }
}

/*tree refinement */
void wp4est_refine(int max_level)
{
  int refine_recursive = 0;
  p4est_refine_ext (tree_neko, refine_recursive, max_level,
		    ref_mark_f, NULL, quad_replace);
}

/* tree coarsening */
void wp4est_coarsen()
{
  int coarsen_recursive = 0;
  int callback_orphans = 0;
  p4est_coarsen_ext (tree_neko, coarsen_recursive, callback_orphans,
		     crs_mark_f, NULL, quad_replace);
}

/* 2:1 tree balancing */
void wp4est_balance()
{
  p4est_balance_ext (tree_neko, P4EST_CONNECT_FULL,
		     NULL, quad_replace);
}

/* Make tree copy for later comparison */
void wp4est_tree_copy(int quad_data) {
  if (tree_neko_compare) p4est_destroy (tree_neko_compare);
  tree_neko_compare = p4est_copy (tree_neko, quad_data);
}

/* Check if three was modified */
void wp4est_tree_check(int * check, int quad_data) {
  if (tree_neko_compare) {
    *check = p4est_is_equal(tree_neko, tree_neko_compare, quad_data);
    p4est_destroy (tree_neko_compare);
    tree_neko_compare = NULL;
  } else {
    SC_ABORT("Tree comparison; aborting: no tree_neko_compare\n");
  }
}

/** @brief Iterate over element volumes to transfer element refinement mark
 *
 * @details Required by wp4est_refm_put
 *
 * @param info
 * @param user_data
 */
void iter_refm(p4est_iter_volume_info_t * info, void *user_data) {
  user_data_t *data = (user_data_t *) info->quad->p.user_data;
  int *ref_mark = (int *) user_data;

  // which quad (local and global element number)
  p4est_tree_t *tree;
  p4est_locidx_t iwl;
  int iwlt;
  int il;// loop index

  // get quad number
  tree = p4est_tree_array_index(info->p4est->trees, info->treeid);
  // local quad number
  iwl = info->quadid + tree->quadrants_offset;
  iwlt = (int) iwl;

  // refinement mark for given quad
  data->ref_mark = ref_mark[iwlt];
}

/* fill ref_mark in p4est block */
void wp4est_refm_put(int * ref_mark) {
#ifdef P4_TO_P8
  p4est_iterate(tree_neko, ghost_neko,(void *) ref_mark, iter_refm,
		NULL, NULL, NULL);
#else
  p4est_iterate(tree_neko, ghost_neko,(void *) ref_mark, iter_refm,
		NULL, NULL);
#endif
}

/* data type for partitioning data transfer between neko and p4est*/
typedef struct transfer_part_s {
  int *gnum; /**< pointer to element level array */
  int *lnum; /**< pointer to element group array */
  int *nid; /**< pointer to mpi rank owning the element */
} transfer_part_t;

/** @brief Iterate over element volumes to transfer element global mapping
 *
 * @details Required by wp4est_refm_put
 *
 * @param info
 * @param user_data
 */
void iter_emap(p4est_iter_volume_info_t * info, void *user_data) {
  user_data_t *data = (user_data_t *) info->quad->p.user_data;
  transfer_part_t *el_map = (transfer_part_t *) user_data;

  // which quad (local and global element number)
  p4est_tree_t *tree;
  p4est_locidx_t iwl;
  int iwlt, iwg;
  int ifc;// loop index

  // get quad number
  tree = p4est_tree_array_index(info->p4est->trees, info->treeid);
  // local quad number
  iwl = info->quadid + tree->quadrants_offset;
  iwlt = (int) iwl;
  // global quad number
  iwg = (int) info->p4est->global_first_quadrant[info->p4est->mpirank] + iwlt;

  // sanity check
  if (el_map->gnum[iwlt] != iwg+1){
    SC_ABORT("Wrong global element number; aborting: iter_emap\n");
  }

  // update element mapping data
  data->el_gln = iwg;
  data->el_ln = el_map->lnum[iwlt];
  data->el_nid = el_map->nid[iwlt];
  data->parent_gln = -1;
  data->parent_ln = -1;
  data->parent_nid = -1;
  for (ifc = 0; ifc < P4EST_CHILDREN; ++ifc) {
    data->children_gln[ifc] = -1;
    data->children_ln[ifc] = -1;
    data->children_nid[ifc] = -1;
  }
}

/* fill element global mapping in p4est block */
void wp4est_egmap_put(int * el_gnum,int * el_lnum,int * el_nid) {
  transfer_part_t el_map;
  el_map.gnum = el_gnum;
  el_map.lnum = el_lnum;
  el_map.nid = el_nid;
#ifdef P4_TO_P8
  p4est_iterate(tree_neko, ghost_neko,(void *) &el_map, iter_emap,
		NULL, NULL, NULL);
#else
  p4est_iterate(tree_neko, ghost_neko,(void *) &el_map, iter_emap,
		NULL, NULL);
#endif
}

/* data type for refinement history data transfer between neko and p4est*/
typedef struct transfer_hst_s {
  int map_nr; /**< local number of unchanged elements */
  int rfn_nr; /**< local number of refined elements */
  int crs_nr; /**< local number of coarsened elements */
  int *elgl_map; /**< element global mapping info for unchanged elements */
  int *elgl_rfn; /**< element global mapping info for refined elements */
  int *elgl_crs; /**< element global mapping info for coarsened elements */
} transfer_hst_t;

#define MIN( a, b ) ( ( a > b) ? b : a )

/** @brief Iterate over element volumes to transfer refinement history data
 *
 * @details Required by wp4est_msh_get_hst
 *
 * @param info
 * @param user_data
 */
void iter_msh_hst(p4est_iter_volume_info_t * info, void *user_data) {
  user_data_t *data = (user_data_t *) info->quad->p.user_data;
  transfer_hst_t *trans_data = (transfer_hst_t *) user_data;

  // which quad (local and global element number)
  p4est_tree_t *tree;
  p4est_locidx_t iwl;
  int iwlt, iwg;
  int ifc, ifl, ic, il;// loop index

  // get quad number
  tree = p4est_tree_array_index(info->p4est->trees, info->treeid);
  // local quad number
  iwl = info->quadid + tree->quadrants_offset;
  iwlt = (int) iwl;
  // global quad number
  iwg = (int) info->p4est->global_first_quadrant[info->p4est->mpirank] + iwlt;

  // check refinement status
  ic = -1;
  for (il=0; il < P4EST_CHILDREN; il++) {
    if(data->children_gln[il] != -1) ic = 1;
  }

  if (data->parent_gln == -1 && ic == -1) {
    // no refinement
    // count elements
    trans_data->map_nr = trans_data->map_nr + 1;
    // set old element position
    trans_data->elgl_map[3*iwlt] = data->el_gln + 1;
    trans_data->elgl_map[3*iwlt+1] = data->el_ln;
    trans_data->elgl_map[3*iwlt+2] = data->el_nid;
  } else if (data->parent_gln != -1) {
    // refinement
    // count elements
    trans_data->rfn_nr = trans_data->rfn_nr + 1;
    // set dummy element map
    trans_data->elgl_map[3*iwlt] = 0;
    trans_data->elgl_map[3*iwlt+1] = 0;
    trans_data->elgl_map[3*iwlt+2] = 0;
    ic = (trans_data->rfn_nr-1)*5;
    // current global element number
    trans_data->elgl_rfn[ic] = iwg + 1;
    // old parent element position
    trans_data->elgl_rfn[ic +1] = data->parent_gln + 1;
    trans_data->elgl_rfn[ic +2] = data->parent_ln;
    trans_data->elgl_rfn[ic +3] = data->parent_nid;
    // child position; numbered 0,..,P4EST_CHILDREN-1
    trans_data->elgl_rfn[ic +4] = data->el_gln;
  } else {
    // coarsening
    // count elements
    trans_data->crs_nr = trans_data->crs_nr + 1;
    // set dummy element map
    trans_data->elgl_map[3*iwlt] = 0;
    trans_data->elgl_map[3*iwlt+1] = 0;
    trans_data->elgl_map[3*iwlt+2] = 0;
    // new global position
    ic =(trans_data->crs_nr - 1)*4*P4EST_CHILDREN;
    trans_data->elgl_crs[ic] = iwg + 1;
    // old global position
    trans_data->elgl_crs[ic+1] = data->children_gln[0] + 1;
    trans_data->elgl_crs[ic+2] = data->children_ln[0];
    trans_data->elgl_crs[ic+3] = data->children_nid[0];
    for (il = 1; il < P4EST_CHILDREN; il++) {
      // new dummy global position
      trans_data->elgl_crs[ic+il*4] = 0;
      // old global position
      trans_data->elgl_crs[ic+il*4+1] = data->children_gln[il] + 1;
      trans_data->elgl_crs[ic+il*4+2] = data->children_ln[il];
      trans_data->elgl_crs[ic+il*4+3] = data->children_nid[il];
    }
  }
}

// get refinement history data to Neko
void wp4est_msh_get_hst(int * map_nr, int * rfn_nr, int * crs_nr, int *elgl_map,
			int * elgl_rfn, int * elgl_crs) {
  transfer_hst_t transfer_data;
  transfer_data.map_nr = 0;
  transfer_data.rfn_nr = 0;
  transfer_data.crs_nr = 0;
  transfer_data.elgl_map = elgl_map;
  transfer_data.elgl_rfn = elgl_rfn;
  transfer_data.elgl_crs = elgl_crs;
#ifdef P4_TO_P8
  p4est_iterate(tree_neko, ghost_neko,(void *) &transfer_data, iter_msh_hst,
		NULL, NULL, NULL);
#else
  p4est_iterate(tree_neko, ghost_neko,(void *) &transfer_data, iter_msh_hst,
		NULL, NULL);
#endif
  *map_nr = transfer_data.map_nr;
  *rfn_nr = transfer_data.rfn_nr;
  *crs_nr = transfer_data.crs_nr;
}

