/** Wrapper interface of p4est and sc libraries
 *
 */
#ifndef NEKO_P4EST_FWRAP_H_
#define NEKO_P4EST_FWRAP_H_

#if !defined(HAVE_P4EST)
#error "p4est_wrapper.h requires p4est to be turned on"
#endif

#if !defined(SC_H)
#error "p4est_wrapper.h requires sc.h"
#endif

#if !defined(P4EST_H) && !defined(P8EST_H)
#error "p4est_wrapper.h requires p4est.h or p8est.h"
#endif

/** Data type for user variables; required by p4est */
typedef struct user_data_s {
  /* Mesh data initialised with p4est */
  int imsh; /**< velocity (0) and temperature (1) mesh indicator */
  int igrp; /**< element group */
  int crv[P4EST_FACES]; /**< curvature data; concerns external faces; requred
			   by face projection */
  int bc[P4EST_FACES]; /**< boundary condition data; 0 internal, -1 periodic,
			  otherwise external */
  /* Data needed for refinement */
  int ref_mark; /**< integer to store refinement mark; definition should be
		   given in .h */
  // to keep track of changes of Neko global element numbering and element
  // distribution element
  int el_gln; /**< element global numbering; neko side */
  int el_ln; /**< element local numbering; neko side */
  int el_nid; /**< mpi rank owning element; neko side */
  // parent
  int parent_gln; /**< parent global numbering; neko side */
  int parent_ln; /**< parent local numbering; neko side */
  int parent_nid; /**< mpi rank owning parent; neko side */
  // children
  int children_gln[P4EST_CHILDREN]; /**< children global numbering; neko side */
  int children_ln[P4EST_CHILDREN]; /**< children local numbering; neko side */
  int children_nid[P4EST_CHILDREN]; /**< mpi rank owning children; neko side */
} user_data_t;

/** Initialise sc and p4est libraries setting log threshold
 *
 * @param[in] fmpicomm           MPI communicator
 * @param[in] catch_signals      If true, signals INT and SEGV are caught.
 * @param[in] print_backtrace    If true, sc_abort prints a backtrace
 * @param[in] log_threshold      Log threshold
 */
void wp4est_init(MPI_Fint fmpicomm, int catch_signals,
		int print_backtrace, int log_threshold);

/** Finalise sc printing package summary
 *
 * @param[in] log_priority       Logging priority
 */
void wp4est_finalize(int log_priority);

/** Destroy mesh connectivity */
void wp4est_cnn_del()
;

/** Check connectivity consistency
 *
 * @param[out] is_valid   non zero for correct connectivity
 */
void wp4est_cnn_valid(int * is_valid)
;

/** Destroy tree */
void wp4est_tree_del()
;

/** Check tree consistency
 *
 * @param[out] is_valid     non zero for correct tree
 */
void wp4est_tree_valid(int * is_valid)
;

/** Save tree to the file
 *
 * @param[in] filename       file name
 */
void wp4est_tree_save(char filename[])
;

/** Load tree from a file
 *
 * @param[in] filename            file name
 */
void wp4est_tree_load(char filename[])
;

/** Build ghost layer */
void wp4est_ghost_new()
;

/** Destroy ghost layer */
void wp4est_ghost_del()
;

/** Generate mesh information */
void wp4est_mesh_new()
;

/** Destroy mesh information */
void wp4est_mesh_del()
;

/** Generate new node information */
void wp4est_nodes_new()
;

/** Destroy node information */
void wp4est_nodes_del()
;

/** Generate new global nodes (GLL points) numbering
 *
 * @param[in] degree   polynomial degree
 */
void wp4est_lnodes_new(int degree)
;

/** Destroy global node numbering */
void wp4est_lnodes_del()
;

/** Forest partitioning for p4est */
void wp4est_part()
;

/** Initialise p4est block with Neko mesh data
 *
 * @param[in] gidx   gobal element index
 * @param[in] imsh   element mesh type
 * @param[in] igrp   element group
 * @param[in] crv    element face curvature flag
 * @param[in] bc     element face boundary condition flag
 */
void wp4est_elm_ini_dat(int64_t * gidx, int * imsh, int * igrp, int * crv,
		      int * bc)
;

/** Check consistency of boundary conditions for V- and T-type mesh
 */
void wp4est_bc_check()
;

/** Get mesh size information to Neko
 *
 * @param[out] mdim    mesh dimension
 * @param[out] nelgt   global element number
 * @param[out] nelgto  element offset (number of elements on lower nid's)
 * @param[out] nelt    number of T-type elements
 * @param[out] nelv    number of V-type elements
 * @param[out] maxl    current max level
 */
void wp4est_msh_get_size(int * mdim, int64_t * nelgt, int64_t * nelgto,
			 int32_t * nelt, int * nelv, int * maxl)
;

/** Get node list size information to Neko
 *
 * @param[out] nowin   number of owned independent nodes
 * @param[out] nowsh   number of owned shared
 * @param[out] oowin   position of the first independent owned node
 * @param[out] nin     numbers of local independent nodes
 * @param[out] nhf     number of local face hanging nodes
 * @param[out] nhe     number of local edge hanging nodes
 */
void wp4est_nds_get_size(int * nowin, int * nowsh, int * oowin,
			 int * nin, int * nhf, int * nhe)
;

/** Get list of independent nodes
 *
 * @details Node coordinates are exact for 0-level mesh or a linear mesh only.
 *
 * @param[out] nglid   node global id
 * @param[out] nown    node owner (mpi rank)
 * @param[out] ncoord  node coordinates
 */
void wp4est_nds_get_ind(int64_t * nglid, int * nown, double * ncoord)
;

/** Get list of face hanging nodes
 *
 * @details Node coordinates are exact for 0-level mesh or a linear mesh only.
 *
 * @param[out] depend  mapping into independent nodes (2 for 2D and 4 for 3D)
 * @param[out] ncoord  node coordinates
 */
void wp4est_nds_get_hfc(int * depend, double * ncoord)
;

/** Get list of edge hanging nodes
 *
 * @details Node coordinates are exact for 0-level mesh or a linear mesh only.
 *
 * @param[out] depend  mapping into independent nodes (2 for 3D)
 * @param[out] ncoord  node coordinates
 */
void wp4est_nds_get_hed(int * depend, double * ncoord)
;

/** Get mapping of element vertices to nodes
 *
 * @param[out] vmap    vertex mapping to nodes (local numberring)
 */
void wp4est_nds_get_vmap(int * vmap)
;

/** Get element info
 *
 * @details It includes approximate (linear elements) physical coordinates of
 *          the p4est block vertices
 *
 * @param[out] gidx  gobal element index
 * @param[out] level element level
 * @param[out] igrp  element group
 * @param[out] crv   face projection flag
 * @param[out] bc    boundary surface flag
 * @param[out] coord physical coordinates
 * @param[out] falg  face alignment
 */
void wp4est_elm_get_dat(int64_t * gidx, int * level, int * igrp, int * crv,
			int * bc, double * coord, int * falg)
;


/** Get global numberring of degrees of freedom in element
 *
 * @param[out] lnnum      number of local nodes
 * @param[out] lnown      number of owned nodes
 * @param[out] lnoff      global node offset
 * @param[out] lnodes     node mapping list for elements
 */
void wp4est_elem_get_lnode(int * lnnum, int * lnown, int64_t * lnoff,
			   int * lnodes);

/** Get sharers list size information to Neko
 *
 * @param[out] nown    number of owned lnodes
 * @param[out] nloc    number of local lnodes
 */
void wp4est_sharers_get_size(int * nrank, int * nshare)
;

/** Get sharers list to Neko
 *
 * @param[out] nglid   node global index
 * @param[out] lrank   MPI rank list
 * @param[out] loff    offset list
 * @param[out] lshare  share node list (local indexes)
 */
void wp4est_sharers_get_ind(int64_t * nglid, int * lrank, int * loff,
			    int * lshare)
;

/** Get hanging objects information
 *
 * @param[out hang_elm   is any of eement's faces/edges hanging
 * @param[out hang_fsc   hanging face list
 * @param[out hang_edg   hanging edge list; 3D mesh only
 */
void wp4est_hang_get_info(int * hang_elm, int * hang_fsc, int * hang_edg)
;

/** Provide information about element families
 *
 * @param[out family array containing global element numbers in the family
 * @param[out] nelf   number of entrances in family array
 */
void wp4est_fml_get_info(int64_t * family, int * nelf)
;

/** Perform tree refinement
 *
 * @param max_level    max refinement level
 */
void wp4est_refine(int max_level)
;

/** Perform tree coarsening
 */
void wp4est_coarsen()
;

/** Perform 2:1 tree balancing
 */
void wp4est_balance()
;

/** Make tree copy for later comparison
 *
 * @param quad_data   do we test quadrant data
 */
void wp4est_tree_copy(int quad_data)
;

/** Check if tree was modified
 *
 * @param check       tree modification marker
 * @param quad_data   do we test quadrant data
 */
void wp4est_tree_check(int * check, int quad_data)
;

/** Fill ref_mark in p4est block
 *
 * @param ref_mark   refinement mark array
 */
void wp4est_refm_put(int * ref_mark)
;

/** Fill element neko element distribution in p4est block
 *
 * @param el_gnum   element global number (neko distribution)
 * @param el_lnum   element local number (neko distribution)
 * @param el_nid   element mpi rank (neko distribution)
 */
void wp4est_egmap_put(int * el_gnum, int * el_lnum, int * el_nid)
;

/** Get refinement/coarsening data to Neko
 *
 * @param map_nr     local number of unchanged elements
 * @param rfn_nr     local number of refined elements
 * @param crs_nr     local number of coarsened elements
 * @param elgl_map   element mapping info for unchanged elements
 * @param elgl_rfn   element mapping info for refined elements
 * @param elgl_crs   element mapping info for coarsened elements
 */
void wp4est_msh_get_hst(int * map_nr, int * rfn_nr, int * crs_nr, int *elgl_map,
			int * elgl_rfn, int * elgl_crs)
;


#endif /* NEKO_P4EST_FWRAP_H_ */
