/*
  This code is based on p4est_lnodes.c from the p4est (2.8.7) source
  code. Its only goal is to provide additional feature not present in
  the original code, which is a separate global numbering of edges.
  To get it I've modified a smallest needed set of routines.
*/

/* p4est speciffic definitions
* This is dimension rleated; for now 3D only */
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
#include <p8est_bits.h>
#include <p8est_communication.h>
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <p8est_lnodes.h>
#include "mesh/mesh_manager/p8est_lnodes_edge.h"

#define P8EST_LN_E_OFFSET 6

static int
p4est_locidx_offset_compare (const void *key, const void *elem)
{
  const p4est_locidx_t *start = (p4est_locidx_t *) elem;
  const p4est_locidx_t *stop = (p4est_locidx_t *) start + 1;
  int                 comp = p4est_locidx_compare (key, start);
  if (comp < 0) {
    return -1;
  }
  comp = p4est_locidx_compare (key, stop);
  if (comp >= 0) {
    return 1;
  }
  return 0;
}

typedef struct p4est_lnodes_dep
{
  p4est_locidx_t      face[P4EST_DIM];
  p4est_locidx_t      edge[P4EST_DIM];
}
p4est_lnodes_dep_t;

typedef struct p4est_lnodes_buf_info
{
  int8_t              type;     /* which nodes it shares */
  int8_t              send_sharers;     /* whether the sharers are included in
					   the message */
  p4est_locidx_t      first_index;      /* inodes array, first node to/from */
  p4est_locidx_t      share_offset;
  int8_t              share_count;
}
p4est_lnodes_buf_info_t;

typedef struct p4est_lnodes_data
{
  p4est_lnodes_dep_t *local_dep;        /* num local quads */
  p4est_lnodes_dep_t *ghost_dep;        /* num ghost quads */
  p4est_locidx_t     *local_elem_nodes; /* num local quads * nodes per q */
  p4est_locidx_t     *poff;     /* mpisize + 1 */
  sc_array_t         *inodes;   /* 2 * p4est_locidx_t */
  sc_array_t         *inode_sharers;    /* int */
  sc_array_t         *send_buf_info;    /* one for each proc: type buf_info_t */
  sc_array_t         *recv_buf_info;    /* one for each proc: type buf_info_t */
  p4est_lnodes_code_t *face_codes;
  int                 nodes_per_elem;
  int                 nodes_per_edge;
  int                *edge_nodes[P8EST_EDGES];
  sc_array_t          send_requests;
  sc_array_t         *send_buf;
  sc_array_t         *touching_procs;
  sc_array_t         *all_procs;
}
p4est_lnodes_data_t;

static inline int
fside_get_fields (p4est_iter_face_side_t * fside, int *is_hanging,
		  p4est_topidx_t * tid, int *f, int8_t ** is_ghost,
		  p4est_locidx_t ** quadid, p4est_quadrant_t *** quad)
{
  int limit;

  *is_hanging = fside->is_hanging;
  *tid = fside->treeid;
  *f = (int) fside->face;
  if (fside->is_hanging) {
    limit = P4EST_HALF;
    *is_ghost = fside->is.hanging.is_ghost;
    *quadid = fside->is.hanging.quadid;
    *quad = fside->is.hanging.quad;
  }
  else {
    limit = 1;
    *is_ghost = &fside->is.full.is_ghost;
    *quadid = &fside->is.full.quadid;
    *quad = &fside->is.full.quad;
  }

  return limit;
}

static void
p4est_lnodes_face_simple_callback (p4est_iter_face_info_t * info, void *Data)
{
  int                 i, f, fdir, limit, cid, xind, *ip;
  sc_array_t         *sides = &(info->sides);
  size_t              zz, count = sides->elem_count;
  p4est_lnodes_data_t *data = (p4est_lnodes_data_t *) Data;
  p4est_iter_face_side_t *fside;
  sc_array_t         *trees = info->p4est->trees;
  p4est_tree_t       *tree;
  p4est_topidx_t      tid;
  sc_array_t         *touching_procs = data->touching_procs;
  sc_array_t          proc_offsets;
  p4est_lnodes_dep_t *local_dep = data->local_dep;
  p4est_lnodes_dep_t *ghost_dep = data->ghost_dep;
  p4est_lnodes_dep_t *dep;
  p4est_locidx_t      quadrants_offset;
  int                 rank = info->p4est->mpirank;
  p4est_lnodes_code_t *face_codes = data->face_codes;
  int8_t             *is_ghost;
  p4est_locidx_t     *quadid;
  p4est_quadrant_t  **quad;
  int                 is_hanging;
  int                 procs[P4EST_HALF];
  p4est_locidx_t      qid[P4EST_HALF];
  int                 j;
  int                 k, c;
  p4est_quadrant_t    tempq;

  P4EST_ASSERT (touching_procs->elem_size == sizeof (int));
  sc_array_truncate (touching_procs);
  /* even though the original is size mpisize+1, proc_offsets uses
   * p4est_locidx_offset_compare, and we don't want to read past the end of the
   * array */
  sc_array_init_data (&proc_offsets, info->ghost_layer->proc_offsets,
		      sizeof (p4est_locidx_t), (size_t) info->p4est->mpisize);

  for (zz = 0; zz < count; zz++) {
    fside = p4est_iter_fside_array_index (sides, zz);
    limit = fside_get_fields (fside, &is_hanging, &tid, &f, &is_ghost,
			      &quadid, &quad);
    fdir = f / 2;
    tree = p4est_tree_array_index (trees, tid);
    quadrants_offset = tree->quadrants_offset;
    k = -1;
    for (i = 0; i < limit; i++) {
      qid[i] = quadid[i];
      if (qid[i] < 0) {
	P4EST_ASSERT (limit == P4EST_HALF);
	if (k < 0) {
	  for (k = 0; k < P4EST_HALF; k++) {
	    if (quadid[k] >= 0) {
	      P4EST_ASSERT (quad[k]);
	      break;
	    }
	  }
	}
	P4EST_ASSERT (k >= 0 && k < P4EST_HALF);
	c = p4est_face_corners[f][i];
	p4est_quadrant_sibling (quad[k], &tempq, c);
	procs[i] = p4est_comm_find_owner (info->p4est, tid,
					  &tempq, info->p4est->mpirank);
	P4EST_ASSERT (procs[i] >= 0 && procs[i] != rank);
	ip = (int *) sc_array_push (touching_procs);
	*ip = procs[i];
      }
      else if (is_ghost[i]) {
	procs[i] = (int) sc_array_bsearch (&proc_offsets, &(qid[i]),
					   p4est_locidx_offset_compare);
	P4EST_ASSERT (procs[i] >= 0 && procs[i] != rank);
	ip = (int *) sc_array_push (touching_procs);
	*ip = procs[i];
      }
      else {
	qid[i] += quadrants_offset;
	procs[i] = rank;
	/* update face code */
	if (is_hanging) {
	  face_codes[qid[i]] |= ((p4est_lnodes_code_t)
				 p4est_face_corners[f][i]);
	  face_codes[qid[i]] |= ((p4est_lnodes_code_t)
				 1 << (P4EST_DIM + f / 2));
	}
      }
    }

    for (i = 0; i < limit; i++) {
      dep = !is_ghost[i] ? &(local_dep[qid[i]]) : &(ghost_dep[qid[i]]);
      if (is_hanging) {
	int                 ndir[2];

	ndir[0] = SC_MIN (((fdir + 1) % 3), ((fdir + 2) % 3));
	ndir[1] = SC_MAX (((fdir + 1) % 3), ((fdir + 2) % 3));

	for (j = 0; j < 2; j++) {
	  xind = i ^ (j + 1);
	  if (is_ghost[xind]) {
	    dep->edge[ndir[j]] = -((p4est_locidx_t) procs[xind] + 3);
	  }
	  else {
	    dep->edge[ndir[j]] = qid[xind];
	  }
	}
	xind = i ^ (P4EST_HALF - 1);
	if (is_ghost[xind]) {
	  dep->face[fdir] = -((p4est_locidx_t) procs[xind] + 3);
	}
	else {
	  dep->face[fdir] = qid[xind];
	}
      }
      else {
	cid = p4est_quadrant_child_id (quad[i]);
	if (p4est_corner_face_corners[cid][f] >= 0) {
	  dep->face[fdir] = -2;
	}
      }
    }
  }
}

static inline int
eside_get_fields (p8est_iter_edge_side_t * eside, int *is_hanging,
		  p4est_topidx_t * tid, int *e, int *o, int8_t ** is_ghost,
		  p4est_locidx_t ** quadid, p4est_quadrant_t *** quad)
{
  int limit;

  *is_hanging = eside->is_hanging;
  *tid = eside->treeid;
  *e = (int) eside->edge;
  *o = (int) eside->orientation;
  if (eside->is_hanging) {
    limit = 2;
    *is_ghost = eside->is.hanging.is_ghost;
    *quadid = eside->is.hanging.quadid;
    *quad = eside->is.hanging.quad;
  }
  else {
    limit = 1;
    *is_ghost = &eside->is.full.is_ghost;
    *quadid = &eside->is.full.quadid;
    *quad = &eside->is.full.quad;
  }

  return limit;
}

static int
p8est_lnodes_edge_simple_callback (p8est_iter_edge_info_t * info, void *Data)
{
  int                 i, limit, e, edir, cid, *ip;
  sc_array_t         *sides = &(info->sides);
  size_t              zz, count = sides->elem_count;
  p4est_lnodes_data_t *data = (p4est_lnodes_data_t *) Data;
  p8est_iter_edge_side_t *eside;
  sc_array_t         *trees = info->p4est->trees;
  p4est_tree_t       *tree;
  p4est_topidx_t      tid;
  sc_array_t          proc_offsets;
  sc_array_t         *touching_procs = data->touching_procs;
  p4est_lnodes_dep_t *local_dep = data->local_dep;
  p4est_lnodes_dep_t *ghost_dep = data->ghost_dep;
  p4est_lnodes_dep_t *dep;
  int8_t             *is_ghost;
  int                 procs[2];
  int                 rank = info->p4est->mpirank;
  p4est_locidx_t     *quadid;
  p4est_quadrant_t  **quad;
  p4est_locidx_t      qid[2];
  p4est_locidx_t      quadrants_offset;
  p4est_lnodes_code_t *face_codes = data->face_codes;
  int                 is_hanging, o, has_local = 0, c;
  p4est_quadrant_t    tempq;

  P4EST_ASSERT (touching_procs->elem_size == sizeof (int));
  sc_array_truncate (touching_procs);
  /* even though the original is size mpisize+1, proc_offsets uses
   * p4est_locidx_offset_compare, and we don't want to read past the end of the
   * array */
  sc_array_init_data (&proc_offsets, info->ghost_layer->proc_offsets,
		      sizeof (p4est_locidx_t), (size_t) info->p4est->mpisize);

  for (zz = 0; zz < count; zz++) {
    eside = p8est_iter_eside_array_index (sides, zz);
    limit = eside_get_fields (eside, &is_hanging, &tid, &e, &o, &is_ghost,
			      &quadid, &quad);
    edir = e / 4;
    tree = p4est_tree_array_index (trees, tid);
    quadrants_offset = tree->quadrants_offset;
    for (i = 0; i < limit; i++) {
      qid[i] = quadid[i];
      if (qid[i] < 0) {
	if (limit == 2 && quadid[i ^ 1] >= 0) {
	  P4EST_ASSERT (quad[i ^ 1]);
	  c = p8est_edge_corners[e][i];
	  p4est_quadrant_sibling (quad[i ^ 1], &tempq, c);
	  procs[i] = p4est_comm_find_owner (info->p4est, tid,
					    &tempq, info->p4est->mpirank);
	  P4EST_ASSERT (procs[i] >= 0 && procs[i] != rank);
	  ip = (int *) sc_array_push (touching_procs);
	  *ip = procs[i];
	}
      }
      else if (is_ghost[i]) {
	procs[i] = (int) sc_array_bsearch (&proc_offsets, &(qid[i]),
					   p4est_locidx_offset_compare);
	P4EST_ASSERT (procs[i] >= 0 && procs[i] != rank);
	ip = (int *) sc_array_push (touching_procs);
	*ip = procs[i];
      }
      else {
	has_local = 1;
	qid[i] += quadrants_offset;
	procs[i] = rank;
	if (is_hanging) {
	  /* update face code */
	  face_codes[qid[i]] |=
	    ((p4est_lnodes_code_t) p8est_edge_corners[e][i]);
	  face_codes[qid[i]] |= ((p4est_lnodes_code_t) 1 << (6 + e / 4));
	}
      }
    }
    for (i = 0; i < limit; i++) {
      if (qid[i] < 0) {
	continue;
      }
      dep = !is_ghost[i] ? &(local_dep[qid[i]]) : &(ghost_dep[qid[i]]);
      if (is_hanging) {
	if (!has_local && qid[i ^ 1] < 0) {
	  dep->edge[edir] = -1;
	}
	else if (is_ghost[i ^ 1]) {
	  dep->edge[edir] = -((p4est_locidx_t) procs[i ^ 1] + 3);
	}
	else {
	  dep->edge[edir] = qid[i ^ 1];
	}
      }
      else {
	cid = p4est_quadrant_child_id (quad[i]);
	if (p8est_edge_corners[e][0] == cid ||
	    p8est_edge_corners[e][1] == cid) {
	  dep->edge[edir] = -2;
	}
      }
    }
  }

  return has_local;
}

static inline void
p4est_lnodes_push_binfo (sc_array_t * touch, sc_array_t * all,
			 sc_array_t * send, sc_array_t * recv,
			 sc_array_t * share, int owner, int rank,
			 int mpisize, int is_remote,
			 int8_t type, p4est_locidx_t nin)
{
  size_t              zz, count = all->elem_count;
  int                *ip, proc;
  p4est_lnodes_buf_info_t *binfo;
  int8_t              scount = -1;
  p4est_locidx_t      offset = (p4est_locidx_t) share->elem_count;

  if (!is_remote) {
    ip = (int *) sc_array_push (share);
    *ip = rank;
    scount = (int8_t) (count + 1);
  }
  for (zz = 0; zz < count; zz++) {
    proc = *((int *) sc_array_index (all, zz));
    if (!is_remote) {
      ip = (int *) sc_array_push (share);
      *ip = proc;
    }
    if (owner == rank) {
      P4EST_ASSERT (proc != rank);
      P4EST_ASSERT (!is_remote);
      P4EST_ASSERT (0 <= proc && proc < mpisize);
      binfo = (p4est_lnodes_buf_info_t *) sc_array_push (&(send[proc]));
      binfo->send_sharers = 1;
      if (touch == NULL ||
	  sc_array_bsearch (touch, &proc, sc_int_compare) >= 0) {
	binfo->send_sharers = 0;
      }
    }
    else if (proc == owner) {
      P4EST_ASSERT (0 <= proc && proc < mpisize);
      binfo = (p4est_lnodes_buf_info_t *) sc_array_push (&(recv[proc]));
      if (!is_remote) {
	binfo->send_sharers = 0;
      }
      else {
	binfo->send_sharers = 1;
      }
    }
    else {
      continue;
    }
    binfo->type = type;
    binfo->first_index = nin;
    if (!is_remote) {
      binfo->share_offset = offset;
      binfo->share_count = scount;
    }
    else {
      binfo->share_offset = -1;
      binfo->share_count = -1;
    }
  }
}

static void
p8est_lnodes_missing_proc_edge (p8est_iter_edge_info_t * info, int side,
				int b, int *mproc)
{
  sc_array_t         *sides = &(info->sides);
  int                 i, nsides = (int) sides->elem_count;
  p8est_iter_edge_side_t *thisside = p8est_iter_eside_array_index_int
    (sides, side);
  p8est_iter_edge_side_t *eside;
  p4est_quadrant_t   *q, tempq, tempr;
  int                 key, test;
  int                 j;
  int                 e = thisside->edge, f;
  int                 edir = e / 4;
  int                 missdir = 3 - edir - b;
  int                 c, c2;

  P4EST_ASSERT (edir != b);
  P4EST_ASSERT (thisside->is_hanging);
  P4EST_ASSERT (p8est_edge_faces[e][b < missdir ? 0 : 1] / 2 == b);
  q = thisside->is.hanging.quad[0];
  if (!q) {
    q = thisside->is.hanging.quad[1];
    P4EST_ASSERT (q);
  }
  key = thisside->faces[b < missdir ? 0 : 1];
  c = p8est_edge_corners[e][0];
  f = p8est_corner_faces[c][b];
  c = p8est_corner_face_corners[c][f];
  c = p8est_face_corners[f][c ^ 3];
  c2 = p8est_edge_corners[e][1];
  c2 = p8est_corner_face_corners[c2][f];
  c2 = p8est_face_corners[f][c2 ^ 3];
  p4est_quadrant_sibling (q, &tempq, c);
  p4est_quadrant_sibling (q, &tempr, c2);
  for (i = 0; i < nsides; i++) {
    if (i == side) {
      continue;
    }
    eside = p8est_iter_eside_array_index_int (sides, i);
    for (j = 0; j < 2; j++) {
      test = eside->faces[j];
      if (test == key) {
	if (!eside->is_hanging && eside->is.full.quad != NULL) {
	  mproc[0] =
	    p4est_comm_find_owner (info->p4est, thisside->treeid, &tempq,
				   info->p4est->mpirank);
	  P4EST_ASSERT (mproc[0] >= 0);
	  mproc[1] =
	    p4est_comm_find_owner (info->p4est, thisside->treeid, &tempr,
				   mproc[0]);
	  P4EST_ASSERT (mproc[1] >= 0);
	  return;
	}
	else {
	  mproc[0] = -1;
	  mproc[1] = -1;
	  return;
	}
      }
    }
  }
  mproc[0] = -1;
  mproc[1] = -1;
}

static void
p8est_lnodes_edge_callback (p8est_iter_edge_info_t * info, void *Data)
{
  int                 i, j, k, xdir[2];
  sc_array_t         *sides = &(info->sides);
  size_t              zz, count = sides->elem_count;
  p4est_lnodes_data_t *data = (p4est_lnodes_data_t *) Data;
  p8est_iter_edge_side_t *eside, *owner_eside;
  sc_array_t         *inodes = data->inodes;
  p4est_locidx_t     *inode;
  sc_array_t         *inode_sharers = data->inode_sharers;
  p4est_lnodes_dep_t *local_dep = data->local_dep;
  p4est_lnodes_dep_t *ghost_dep = data->ghost_dep;
  p4est_lnodes_dep_t *dep;
  p4est_locidx_t     *local_elem_nodes = data->local_elem_nodes;
  sc_array_t         *send_buf_info = data->send_buf_info;
  sc_array_t         *recv_buf_info = data->recv_buf_info;
  sc_array_t         *touching_procs = data->touching_procs;
  sc_array_t         *all_procs = data->all_procs;
  int                *ip;
  p4est_topidx_t      tid, owner_tid;
  sc_array_t         *trees = info->p4est->trees;
  p4est_tree_t       *tree;
  p4est_locidx_t      quadrants_offset;
  p4est_locidx_t     *qids;
  p4est_locidx_t      qid, owner_qid, nqid;
  p4est_locidx_t      num_inodes = (p4est_locidx_t) inodes->elem_count;
  int8_t             *is_ghost, owner_is_ghost;
  int                 e, edir, owner_e, owner_c, o;
  p4est_locidx_t      nid;
  int                 owner_proc, nproc;
  int                 rank = info->p4est->mpirank;
  p4est_quadrant_t   *owner_q = NULL;
  p4est_quadrant_t   *q;
  p4est_quadrant_t  **quad;
  p4est_quadrant_t    tempq, tempr, ownq;
  int                 nodes_per_edge = data->nodes_per_edge;
  int                 nodes_per_elem = data->nodes_per_elem;
  int               **edge_nodes = data->edge_nodes;
  int                 is_hanging;
  int                 limit;
  int                 stride;
  p4est_locidx_t      start_node;
  int8_t              type;
  int                 is_remote, has_local;
  p4est_connectivity_t *conn = info->p4est->connectivity;
  int                 mproc[2][2] = { {-1, -1}, {-1, -1} };

  sc_array_truncate (touching_procs);
  sc_array_truncate (all_procs);
  has_local = p8est_lnodes_edge_simple_callback (info, data);

  owner_eside = p8est_iter_eside_array_index (sides, 0);
  owner_e = owner_eside->edge;
  owner_tid = owner_eside->treeid;
  if (owner_eside->is_hanging) {
    owner_qid = owner_eside->is.hanging.quadid[0];
    owner_is_ghost = owner_eside->is.hanging.is_ghost[0];
    owner_q = owner_eside->is.hanging.quad[0];
  }
  else {
    owner_qid = owner_eside->is.full.quadid;
    owner_is_ghost = owner_eside->is.full.is_ghost;
    owner_q = owner_eside->is.full.quad;
  }
  P4EST_ASSERT (!owner_eside->orientation);
  owner_c = p8est_edge_corners[owner_e][0];
  if (owner_q == NULL) {
    int                 c;
    p4est_qcoord_t      x, y, z, h, l;

    P4EST_ASSERT (count > 1);
    eside = NULL;
    for (zz = 1; zz < count; zz++) {
      eside = p8est_iter_eside_array_index (sides, zz);
      if ((!eside->is_hanging) && eside->is.full.quad) {
	break;
      }
    }
    P4EST_ASSERT (zz < count);
    e = eside->edge;
    tid = eside->treeid;
    o = eside->orientation;
    c = p8est_edge_corners[e][o];
    q = eside->is.full.quad;
    P4EST_ASSERT (q->level < P4EST_QMAXLEVEL);
    p4est_quadrant_corner_descendant (q, &tempr, c, q->level + 1);
    q = &tempr;
    P4EST_ASSERT (c == p4est_quadrant_child_id (q));
    p8est_quadrant_edge_neighbor (q, e, &tempq);
    /* get the coordinates of the edge */
    if (p4est_quadrant_is_inside_root (&tempq)) {
      h = P4EST_QUADRANT_LEN (q->level);
      x = q->x + h * (c & 1);
      y = q->y + h * ((c & 2) >> 1);
      z = q->z + h * ((c & 4) >> 2);
    }
    else if (p8est_quadrant_is_outside_edge (&tempq)) {
      h = P4EST_QUADRANT_LEN (q->level);
      if (e / 4 == 0) {
	l = q->x + h * (c & 1);
      }
      else if (e / 4 == 1) {
	l = q->y + h * ((c & 2) >> 1);
      }
      else {
	l = q->z + h * ((c & 4) >> 2);
      }
      if (o) {
	l = P4EST_ROOT_LEN - l;
      }
      h = P4EST_QUADRANT_LEN (0);
      if (owner_e / 4 == 0) {
	x = l;
	y = h * (owner_e & 1);
	z = h * ((owner_e & 2) >> 1);
      }
      else if (owner_e / 4 == 1) {
	x = h * (owner_e & 1);
	y = l;
	z = h * ((owner_e & 2) >> 1);
      }
      else {
	x = h * (owner_e & 1);
	y = h * ((owner_e & 2) >> 1);
	z = l;
      }
    }
    else {
      /* outside face */
      int                 owner_f, c1, c2, nf;
      p4est_topidx_t      nt;
      /* this uses some knowledge about how iterate orders the sides of a
       * corner that is in the middle of a face */
      P4EST_ASSERT (count == 2 || count == 4);
      eside = p8est_iter_eside_array_index (sides, count / 2);
      P4EST_ASSERT (eside->treeid == owner_tid);
      P4EST_ASSERT (eside->edge != owner_e);

      c1 = p8est_edge_corners[owner_e][0];
      c2 = p8est_edge_corners[eside->edge][1];
      owner_f = p4est_child_corner_faces[c1][c2];
      P4EST_ASSERT (owner_f >= 0);

      nt = conn->tree_to_tree[P4EST_FACES * owner_tid + owner_f];
      nf = conn->tree_to_face[P4EST_FACES * owner_tid + owner_f];

      /* o2 = nf / P4EST_FACES; */
      nf %= P4EST_FACES;

      if ((nt == owner_tid && nf == owner_f) || (zz % 2) == 0) {
	/* q must be on the same side: the corner is in the same coordinates
	 */
	P4EST_ASSERT (!o);
	h = P4EST_QUADRANT_LEN (q->level);
	x = q->x + h * (c & 1);
	y = q->y + h * ((c & 2) >> 1);
	z = q->z + h * ((c & 4) >> 2);
      }
      else {
	P4EST_ASSERT (nt == tid);
	P4EST_ASSERT (count == 4);
	P4EST_EXECUTE_ASSERT_TOPIDX
	  (p4est_quadrant_face_neighbor_extra (q, tid, nf, &tempq, NULL,
					       conn), owner_tid);
	c2 = p4est_quadrant_child_id (&tempq);
	P4EST_ASSERT (p4est_corner_face_corners[c2][owner_f] >= 0);
	h = P4EST_QUADRANT_LEN (tempq.level);
	x = tempq.x + h * (c2 & 1);
	y = tempq.y + h * ((c2 & 2) >> 1);
	z = tempq.z + h * ((c2 & 4) >> 2);
      }
    }
    h = P4EST_QUADRANT_LEN (P4EST_QMAXLEVEL);
    ownq.x = x - h * (owner_c & 1);
    ownq.y = y - h * ((owner_c & 2) >> 1);
    ownq.z = z - h * ((owner_c & 4) >> 2);
    ownq.level = P4EST_QMAXLEVEL;
    owner_proc = p4est_comm_find_owner (info->p4est, owner_tid, &ownq, rank);
  }
  else {
    if (owner_is_ghost) {
      owner_proc = *((int *) sc_array_index (touching_procs, 0));
      P4EST_ASSERT (owner_proc >= 0 && owner_proc != rank);
    }
    else {
      owner_proc = rank;
      tree = p4est_tree_array_index (trees, owner_tid);
      quadrants_offset = tree->quadrants_offset;
      owner_qid += quadrants_offset;
    }
  }
  if (has_local) {
    sc_array_sort (touching_procs, sc_int_compare);
    sc_array_uniq (touching_procs, sc_int_compare);
  }
  /* create nodes */
  for (i = 0; i < nodes_per_edge; i++) {
    inode = (p4est_locidx_t *) sc_array_push (inodes);
    P4EST_ASSERT (inodes->elem_count <= (size_t)
		  (nodes_per_elem * info->p4est->local_num_quadrants));
    inode[0] = owner_proc;
    inode[1] = owner_qid;
  }
  /* point element nodes at created nodes; find all sharing procs */
  is_remote = !has_local;
  if (!is_remote) {
    sc_array_copy (all_procs, touching_procs);
  }
  else {
    ip = (int *) sc_array_push (all_procs);
    *ip = owner_proc;
  }
  for (zz = 0; zz < count; zz++) {
    eside = p8est_iter_eside_array_index (sides, zz);
    limit = eside_get_fields (eside, &is_hanging, &tid, &e, &o, &is_ghost,
			      &qids, &quad);
    tree = p4est_tree_array_index (trees, tid);
    quadrants_offset = tree->quadrants_offset;
    if (!is_hanging && quad[0] == NULL) {
      continue;
    }
    mproc[0][0] = -2;
    mproc[0][1] = -2;
    mproc[1][0] = -2;
    mproc[1][1] = -2;
    for (i = 0; i < limit; i++) {
      qid = qids[i];
      if (qid < 0) {
       continue;
      }
      stride = (o ? -1 : 1);
      if (!is_ghost[i]) {
       qid += quadrants_offset;
       P4EST_ASSERT (qid < info->p4est->local_num_quadrants);
       start_node = num_inodes + (o ? nodes_per_edge - 1 : 0);
       for (k = 0; k < nodes_per_edge; k++, start_node += stride) {
	 nid = qid * nodes_per_elem + edge_nodes[e][k];
	 P4EST_ASSERT (local_elem_nodes[nid] == -1);
	 local_elem_nodes[nid] = start_node;
       }
      }
      if (!is_hanging) {
	continue;
      }
      if (is_remote && qids[i ^ 1] < 0) {
	continue;
      }
      /* get quads that may be dependent because of hanging faces */
      dep = !is_ghost[i] ? &local_dep[qid] : &ghost_dep[qid];
      edir = e / 4;
      for (j = 0; j < 2; j++) {
	xdir[0] = (edir + j + 1) % 3;
	xdir[1] = (edir + 2 - j) % 3;
	if (dep->face[xdir[1]] == -1) {
	  P4EST_ASSERT (is_ghost[i]);
	  if (!is_ghost[i ^ 1]) {
	    P4EST_ASSERT (local_dep[qids[i ^ 1] + quadrants_offset].face
			  [xdir[1]] == -2);
	    dep->face[xdir[1]] = -2;
	  }
	  else if (mproc[j][i] == -2) {
	    p8est_lnodes_missing_proc_edge (info, zz, xdir[1],
					    &(mproc[j][0]));
	    P4EST_ASSERT (mproc[j][0] != -2 && mproc[j][1] != -2);
	    P4EST_ASSERT (mproc[j][0] != rank && mproc[j][1] != rank);
	    for (k = 0; k < 2; k++) {
	      if (mproc[j][k] >= 0) {
		ip = (int *) sc_array_push (all_procs);
		*ip = mproc[j][k];
	      }
	    }
	    dep->face[xdir[1]] = -((p4est_locidx_t) mproc[j][i] + 3);
	  }
	  else {
	    dep->face[xdir[1]] = -((p4est_locidx_t) mproc[j][i] + 3);
	  }
	}
	if (dep->face[xdir[1]] == -2) {
	  continue;
	}
	nqid = dep->edge[xdir[0]];
	if (nqid >= 0) {
	  has_local = 1;
	  start_node = num_inodes + (o ? nodes_per_edge - 1 : 0);
	  for (k = 0; k < nodes_per_edge; k++, start_node += stride) {
	    nid = nqid * nodes_per_elem + edge_nodes[e][k];
	    P4EST_ASSERT (local_elem_nodes[nid] == -1);
	    local_elem_nodes[nid] = start_node;
	  }
	}
	else if (!is_remote) {
	  nproc = nqid;
	  if (nproc == -1) {
	    nproc = mproc[j][i ^ 1];
	    dep->edge[xdir[0]] = -((p4est_locidx_t) nproc + 3);
	  }
	  else {
	    nproc = -(nproc + 3);
	  }
	  P4EST_ASSERT (nproc >= -1);
	  if (nproc >= 0 && nproc != rank) {
	    ip = (int *) sc_array_push (all_procs);
	    *ip = nproc;
	  }
	}
      }
    }
  }
  P4EST_ASSERT (has_local);
  sc_array_sort (all_procs, sc_int_compare);
  sc_array_uniq (all_procs, sc_int_compare);

  count = all_procs->elem_count;
  if (count) {
    type = (int8_t) (P8EST_LN_E_OFFSET + owner_e);
    p4est_lnodes_push_binfo (touching_procs, all_procs, send_buf_info,
			     recv_buf_info, inode_sharers, owner_proc, rank,
			     info->p4est->mpisize, is_remote,
			     type, num_inodes);
  }
  else {
    P4EST_ASSERT (owner_proc == rank);
  }
}

static void
p4est_lnodes_init_data (p4est_lnodes_data_t * data, p4est_t * p4est,
			p4est_ghost_t * ghost_layer, p4est_lnodes_t * lnodes)
{
  int                 i, j, n;

  int                 e;
  int                 npe;
  int                 ecount[P8EST_EDGES];

  p4est_locidx_t      nlq = p4est->local_num_quadrants;
  p4est_locidx_t      ngq = (p4est_locidx_t) ghost_layer->ghosts.elem_count;
  p4est_locidx_t      nldep = nlq;
  p4est_locidx_t      ngdep = ngq;
  int                 mpisize = p4est->mpisize;

  data->nodes_per_elem = P8EST_EDGES;
  npe = data->nodes_per_edge = 1;

  ecount[0] = ecount[1] = ecount[2] = ecount[3] = ecount[4] = ecount[5] = 0;
  ecount[6] = ecount[7] = ecount[8] = ecount[9] = ecount[10] = ecount[11] = 0;

  for (i = 0; i < P8EST_EDGES; i++) {
    data->edge_nodes[i] = P4EST_ALLOC (int, npe);
  }

  int offset = 0;

  for (e = 0; e < P8EST_EDGES; e++) {
    for (i = 0; i < npe; i++) {
      data->edge_nodes[e][ecount[e]++] = offset++;
    }
  }

  data->local_dep = P4EST_ALLOC (p4est_lnodes_dep_t, nldep);
  memset (data->local_dep, -1, nldep * sizeof (p4est_lnodes_dep_t));
  data->ghost_dep = P4EST_ALLOC (p4est_lnodes_dep_t, ngdep);
  memset (data->ghost_dep, -1, ngdep * sizeof (p4est_lnodes_dep_t));

  data->local_elem_nodes = lnodes->element_nodes;

  data->inodes = sc_array_new (2 * sizeof (p4est_locidx_t));
  data->inode_sharers = sc_array_new (sizeof (int));
  data->send_buf_info = P4EST_ALLOC (sc_array_t, mpisize);
  data->recv_buf_info = P4EST_ALLOC (sc_array_t, mpisize);
  for (i = 0; i < mpisize; i++) {
    sc_array_init (&(data->send_buf_info[i]),
		   sizeof (p4est_lnodes_buf_info_t));
    sc_array_init (&(data->recv_buf_info[i]),
		   sizeof (p4est_lnodes_buf_info_t));
  }
  data->face_codes = lnodes->face_code;
  data->poff = P4EST_ALLOC_ZERO (p4est_locidx_t, mpisize + 1);
  data->touching_procs = sc_array_new (sizeof (int));
  data->all_procs = sc_array_new (sizeof (int));
}

static void
p4est_lnodes_reset_data (p4est_lnodes_data_t * data, p4est_t * p4est)
{
  int                 mpisize = p4est->mpisize;
  int                 i;

  sc_array_destroy (data->touching_procs);
  sc_array_destroy (data->all_procs);
  P4EST_FREE (data->poff);

  for (i = 0; i < P8EST_EDGES; i++) {
    P4EST_FREE (data->edge_nodes[i]);
  }

  sc_array_destroy (data->inodes);
  sc_array_destroy (data->inode_sharers);
  for (i = 0; i < mpisize; i++) {
    sc_array_reset (&(data->send_buf_info[i]));
    sc_array_reset (&(data->recv_buf_info[i]));
  }
  P4EST_FREE (data->send_buf_info);
  P4EST_FREE (data->recv_buf_info);
  P4EST_FREE (data->local_dep);
  P4EST_FREE (data->ghost_dep);
  /* do not free face_codes: controlled by lnodes_t */
}

static void
p4est_lnodes_count_send (p4est_lnodes_data_t * data, p4est_t * p4est,
			 p4est_lnodes_t * lnodes)
{
  p4est_locidx_t      nlq = p4est->local_num_quadrants;
  p4est_locidx_t      nlen, nln;
  p4est_locidx_t      li, *lp;
  p4est_locidx_t      inidx;
  sc_array_t         *inodes = data->inodes;
  p4est_locidx_t     *inode;
  p4est_topidx_t     *local_en = data->local_elem_nodes;
  int                 i, j;
  int                 rank = p4est->mpirank;
  int                 mpisize = p4est->mpisize;
  int                 npe = data->nodes_per_elem;
  p4est_locidx_t      count = 0;
  sc_array_t         *send_buf_info = data->send_buf_info;
  sc_array_t         *send_info;
  sc_array_t         *send;
  size_t              zz;
  size_t              zindex;
  p4est_lnodes_buf_info_t *binfo;
  int8_t              type;
  int                 limit;
  int                 nodes_per_edge = data->nodes_per_edge;
  int                 share_count;
  int                 share_proc;
  sc_array_t         *inode_sharers = data->inode_sharers;
  size_t              send_count;
  sc_MPI_Request     *send_request;
  int                 num_send_procs;
  size_t              total_sent;
  int                 mpiret;
  size_t              countz;
  p4est_locidx_t     *poff = data->poff;
  p4est_locidx_t      pcount;

  nlen = ((p4est_locidx_t) npe) * nlq;
  for (li = 0; li < nlen; li++) {
    inidx = local_en[li];
    P4EST_ASSERT (inidx >= 0);
    inode = (p4est_locidx_t *) sc_array_index (inodes, (size_t) inidx);
    /* if this quadrant owns the node */
    if (inode[0] == rank && inode[1] == li / npe) {
      inode[0] = -1;
      inode[1] = count++;
    }
  }
  for (zz = 0; zz < inodes->elem_count; zz++) {
    inode = (p4est_locidx_t *) sc_array_index (inodes, zz);
    if (inode[0] >= 0) {
      P4EST_ASSERT (inode[0] != rank);
      poff[inode[0]]++;
    }
  }

  pcount = 0;
  for (i = 0; i < mpisize; i++) {
    p4est_topidx_t      temp = pcount;

    pcount += poff[i];
    poff[i] = temp;
  }
  poff[mpisize] = pcount;

  lnodes->owned_count = count;
  lnodes->num_local_nodes = nln = (p4est_locidx_t) inodes->elem_count;
  lnodes->nonlocal_nodes = P4EST_ALLOC (p4est_gloidx_t, nln - count);
  memset (lnodes->nonlocal_nodes, -1,
	  (nln - count) * sizeof (p4est_gloidx_t));

  num_send_procs = 0;
  total_sent = 0;
  sc_array_init (&(data->send_requests), sizeof (sc_MPI_Request));
  data->send_buf = P4EST_ALLOC (sc_array_t, mpisize);
  for (i = 0; i < mpisize; i++) {
    sc_array_init (&(data->send_buf[i]), sizeof (p4est_locidx_t));
  }
  for (i = 0; i < mpisize; i++) {
    send_info = &(send_buf_info[i]);
    countz = send_info->elem_count;
    if (countz > 0) {
      P4EST_ASSERT (i != p4est->mpirank);
      send = &(data->send_buf[i]);
      for (zz = 0; zz < countz; zz++) {
	binfo = (p4est_lnodes_buf_info_t *) sc_array_index (send_info, zz);
	zindex = (size_t) binfo->first_index;
	type = binfo->type;
	P4EST_ASSERT (type == P8EST_LN_E_OFFSET);
	limit = nodes_per_edge;

	for (j = 0; j < limit; j++) {
	  lp = (p4est_locidx_t *) sc_array_push (send);
	  inode = (p4est_locidx_t *) sc_array_index (inodes, zindex++);
	  P4EST_ASSERT (inode[0] == -1 && inode[1] >= 0 && inode[1] < count);
	  *lp = inode[1];
	}
	if (binfo->send_sharers) {
	  lp = (p4est_locidx_t *) sc_array_push (send);
	  *lp = (p4est_locidx_t) binfo->share_count;
	  P4EST_ASSERT (binfo->share_count > 0);
	  zindex = (size_t) binfo->share_offset;
	  share_count = (int) binfo->share_count;
	  for (j = 0; j < share_count; j++) {
	    lp = (p4est_locidx_t *) sc_array_push (send);
	    share_proc = *((int *) sc_array_index (inode_sharers, zindex++));
	    *lp = (p4est_locidx_t) share_proc;
	    P4EST_ASSERT (0 <= share_proc && share_proc < mpisize);
	  }
	}
      }
      send_count = send->elem_count;
      send_request = (sc_MPI_Request *) sc_array_push (&data->send_requests);
      mpiret = sc_MPI_Isend (send->array,
			     (int) (send_count * sizeof (p4est_locidx_t)),
			     sc_MPI_BYTE, i, P4EST_COMM_LNODES_PASS,
			     p4est->mpicomm, send_request);
      SC_CHECK_MPI (mpiret);
      num_send_procs++;
      total_sent += (send_count * sizeof (p4est_locidx_t));
    }
  }
  P4EST_VERBOSEF ("Total of %llu bytes sent to %d processes\n",
		  (unsigned long long) total_sent, num_send_procs);
}

static void
p4est_lnodes_recv (p4est_t * p4est, p4est_lnodes_data_t * data,
		   p4est_lnodes_t * lnodes)
{
  int                 mpisize = p4est->mpisize;
  int                 i, j, k;
  int                 limit;
  sc_array_t         *recv, *recv_info;
  sc_array_t         *recv_buf;
  sc_array_t         *recv_buf_info = data->recv_buf_info;
  size_t              count, info_count, zz;
  int                 mpiret;
  sc_MPI_Status       probe_status, recv_status;
  int                 num_recv_procs = 0;
  size_t              total_recv = 0;
  int                *num_recv_expect = P4EST_ALLOC_ZERO (int, mpisize);
  int                 byte_count;
  size_t              elem_count;
  p4est_lnodes_buf_info_t *binfo;
  size_t              zindex;
  int                 nodes_per_edge = data->nodes_per_edge;
  p4est_locidx_t     *lp;
  int                *ip;
  p4est_locidx_t     *inode;
  sc_array_t         *inode_sharers = data->inode_sharers;
  sc_array_t         *inodes = data->inodes;
  int                 share_count;
  sc_array_t         *sorter;
  p4est_gloidx_t     *nonlocal_nodes = lnodes->nonlocal_nodes;
  p4est_locidx_t     *poff = data->poff;

  for (i = 0; i < mpisize; i++) {
    recv_info = &(recv_buf_info[i]);
    count = recv_info->elem_count;
    if (count) {
      P4EST_ASSERT (i != p4est->mpirank);
      P4EST_ASSERT (poff[i + 1] - poff[i] > 0);
      num_recv_procs++;
      num_recv_expect[i]++;
    }
  }

  sorter = sc_array_new (2 * sizeof (p4est_locidx_t));

  recv_buf = P4EST_ALLOC (sc_array_t, mpisize);
  for (i = 0; i < mpisize; i++) {
    sc_array_init (&(recv_buf[i]), sizeof (p4est_locidx_t));
  }
  for (i = 0; i < num_recv_procs; i++) {
    mpiret = sc_MPI_Probe (sc_MPI_ANY_SOURCE, P4EST_COMM_LNODES_PASS,
			   p4est->mpicomm, &probe_status);
    SC_CHECK_MPI (mpiret);
    j = probe_status.MPI_SOURCE;
    P4EST_ASSERT (j != p4est->mpirank && num_recv_expect[j] == 1);
    recv = &(recv_buf[j]);
    recv_info = &(recv_buf_info[j]);
    mpiret = sc_MPI_Get_count (&probe_status, sc_MPI_BYTE, &byte_count);
    SC_CHECK_MPI (mpiret);
    P4EST_ASSERT (byte_count % ((int) sizeof (p4est_locidx_t)) == 0);
    elem_count = ((size_t) byte_count) / sizeof (p4est_locidx_t);
    sc_array_resize (recv, elem_count);
    mpiret = sc_MPI_Recv (recv->array, byte_count, sc_MPI_BYTE, j,
			  P4EST_COMM_LNODES_PASS, p4est->mpicomm,
			  &recv_status);
    SC_CHECK_MPI (mpiret);
    num_recv_expect[j]--;

    info_count = recv_info->elem_count;
    count = 0;
    for (zz = 0; zz < info_count; zz++) {
      binfo = (p4est_lnodes_buf_info_t *) sc_array_index (recv_info, zz);
      P4EST_ASSERT (binfo->type >= P8EST_LN_E_OFFSET);
      limit = nodes_per_edge;
      zindex = (size_t) binfo->first_index;
      for (k = 0; k < limit; k++) {
	inode = (p4est_locidx_t *) sc_array_index (inodes, zindex);
	lp = (p4est_locidx_t *) sc_array_index (recv, count++);
	P4EST_ASSERT (inode[0] == j);
	P4EST_ASSERT (*lp >= 0);
	inode[1] = *lp;
	lp = (p4est_locidx_t *) sc_array_push (sorter);
	lp[0] = (p4est_locidx_t) inode[1];
	lp[1] = (p4est_locidx_t) zindex++;
      }
      if (binfo->send_sharers) {
	lp = (p4est_locidx_t *) sc_array_index (recv, count++);
	share_count = (int) (*lp);
	P4EST_ASSERT (share_count > 0);
	P4EST_ASSERT (binfo->share_count == -1);
	P4EST_ASSERT (binfo->share_offset == -1);
	binfo->share_count = (int8_t) share_count;
	binfo->share_offset = (p4est_locidx_t) inode_sharers->elem_count;
	ip = (int *) sc_array_push_count (inode_sharers, share_count);
	for (k = 0; k < share_count; k++) {
	  lp = (p4est_locidx_t *) sc_array_index (recv, count++);
	  ip[k] = (int) (*lp);
	  P4EST_ASSERT (0 <= *ip && *ip < mpisize);
	}
      }
    }
    P4EST_ASSERT (count == elem_count);
    total_recv += byte_count;
    P4EST_ASSERT ((p4est_locidx_t) sorter->elem_count ==
		  poff[j + 1] - poff[j]);
    sc_array_sort (sorter, p4est_locidx_compare);
    for (zz = 0; zz < sorter->elem_count; zz++) {
      lp = (p4est_locidx_t *) sc_array_index (sorter, zz);
      nonlocal_nodes[poff[j] + zz] = lp[1];
    }
    sc_array_reset (sorter);
  }

  if (data->send_requests.elem_count > 0) {
    mpiret = sc_MPI_Waitall ((int) data->send_requests.elem_count,
			     (sc_MPI_Request *) data->send_requests.array,
			     sc_MPI_STATUSES_IGNORE);
    SC_CHECK_MPI (mpiret);
  }
  sc_array_reset (&data->send_requests);
  for (i = 0; i < mpisize; i++) {
    sc_array_reset (&(data->send_buf[i]));
    sc_array_reset (&(recv_buf[i]));
  }

  P4EST_VERBOSEF ("Total of %llu bytes received from %d processes\n",
		  (unsigned long long) total_recv, num_recv_procs);
  P4EST_FREE (data->send_buf);
  P4EST_FREE (recv_buf);
  P4EST_FREE (num_recv_expect);
  sc_array_destroy (sorter);
}

static              p4est_gloidx_t
p4est_lnodes_global_and_sharers (p4est_lnodes_data_t * data,
				 p4est_lnodes_t * lnodes, p4est_t * p4est)
{
  int                 i, j, k, l;
  int                 mpisize = p4est->mpisize;
  p4est_gloidx_t     *gnodes = lnodes->nonlocal_nodes, gtotal;
  size_t              count, zz;
  p4est_locidx_t     *lp, li, *inode;
  sc_array_t         *inodes = data->inodes;
  p4est_locidx_t     *elnodes = lnodes->element_nodes;
  p4est_locidx_t      nlen = lnodes->num_local_elements * lnodes->vnodes;
#ifdef P4EST_ENABLE_DEBUG
  p4est_locidx_t      num_inodes = (p4est_locidx_t) data->inodes->elem_count;
#endif
  p4est_locidx_t      inidx;
  int                *comm_proc;
  int                 comm_proc_count;
  sc_array_t         *inode_sharers = data->inode_sharers;
  sc_array_t         *sharers;
  p4est_lnodes_rank_t *lrank;
  sc_array_t         *binfo_array;
  p4est_lnodes_buf_info_t *binfo;
  p4est_locidx_t      share_offset, owned_count = lnodes->owned_count;
  int                 share_count;
  int                 limit;
  size_t              zindex;
  int                 nodes_per_edge = data->nodes_per_edge;
  int                 proc;
  int                 shareidx;
  p4est_locidx_t      gidx;
  sc_array_t         *shared_nodes;
  p4est_locidx_t     *global_num_indep;
  p4est_gloidx_t     *global_offsets = P4EST_ALLOC (p4est_gloidx_t,
						    mpisize + 1);
  p4est_locidx_t     *poff = data->poff;

  global_num_indep = lnodes->global_owned_count = P4EST_ALLOC (p4est_locidx_t,
							       mpisize);
  sc_MPI_Allgather (&owned_count, 1, P4EST_MPI_LOCIDX, global_num_indep, 1,
		    P4EST_MPI_LOCIDX, p4est->mpicomm);

  global_offsets[0] = 0;
  for (i = 0; i < mpisize; i++) {
    global_offsets[i + 1] = global_offsets[i] +
      (p4est_gloidx_t) global_num_indep[i];
  }
  lnodes->global_offset = global_offsets[p4est->mpirank];
  gtotal = global_offsets[p4est->mpisize];

  i = p4est->mpirank;
  for (i = 0; i < mpisize; i++) {
    if (i == p4est->mpirank) {
      continue;
    }
    for (j = poff[i]; j < poff[i + 1]; j++) {
      li = gnodes[j];
      inode = (p4est_locidx_t *) sc_array_index (inodes, li);
      P4EST_ASSERT (inode[0] == i);
      gnodes[j] = inode[1] + global_offsets[i];
      inode[1] = j + owned_count;
    }
  }

  for (li = 0; li < nlen; li++) {
    inidx = elnodes[li];
    P4EST_ASSERT (0 <= inidx && inidx < num_inodes);
    inode = (p4est_locidx_t *) sc_array_index (inodes, (size_t) inidx);
    if (inode[0] == -1) {
      P4EST_ASSERT (0 <= inode[1] && inode[1] < lnodes->owned_count);
      elnodes[li] = inode[1];
    }
    else {
      P4EST_ASSERT (inode[0] >= 0 && inode[0] != p4est->mpirank &&
		    inode[0] < mpisize);
      P4EST_ASSERT (inode[1] >= poff[inode[0]] + owned_count &&
		    inode[1] < poff[inode[0] + 1] + owned_count);
      elnodes[li] = inode[1];
    }
  }

  /* figure out all nodes that also share nodes shared by the local process */
  comm_proc = P4EST_ALLOC_ZERO (int, mpisize);
  count = inode_sharers->elem_count;
  for (zz = 0; zz < count; zz++) {
    i = *((int *) sc_array_index (inode_sharers, zz));
    comm_proc[i] = 1;
  }
  /* create an entry in sharers for each such process, providing a map from
   * process id to sharer index */
  comm_proc_count = 0;
  lnodes->sharers = sharers = sc_array_new (sizeof (p4est_lnodes_rank_t));
  for (i = 0; i < mpisize; i++) {
    if (comm_proc[i]) {
      lrank = (p4est_lnodes_rank_t *) sc_array_push (sharers);
      lrank->rank = i;
      sc_array_init (&(lrank->shared_nodes), sizeof (p4est_locidx_t));
      comm_proc[i] = comm_proc_count++;
    }
    else {
      comm_proc[i] = -1;
    }
  }

  /* for every node in a send or receive list, figure out which global node it
   * is, and which processes share it, and add the index in global nodes to that
   * sharer's element_nodes array.
   */
  for (i = 0; i < mpisize; i++) {
    for (j = 0; j < 2; j++) {
      if (j == 0) {
	binfo_array = &(data->send_buf_info[i]);
      }
      else {
	binfo_array = &(data->recv_buf_info[i]);
      }
      count = binfo_array->elem_count;
      for (zz = 0; zz < count; zz++) {
	binfo = (p4est_lnodes_buf_info_t *) sc_array_index (binfo_array, zz);
	P4EST_ASSERT (binfo->type == P8EST_LN_E_OFFSET);
	limit = nodes_per_edge;
	zindex = (size_t) binfo->first_index;
	share_offset = binfo->share_offset;
	share_count = (int) binfo->share_count;
	for (k = 0; k < limit; k++) {
	  inode = (p4est_locidx_t *) sc_array_index (inodes, zindex++);
	  gidx = inode[1];
	  if (j == 0) {
	    P4EST_ASSERT (inode[0] == -1);
	    P4EST_ASSERT (gidx < owned_count);
	    shareidx = comm_proc[i];
	    P4EST_ASSERT (shareidx >= 0);
	    lrank = p4est_lnodes_rank_array_index_int (sharers, shareidx);
	    P4EST_ASSERT (lrank->rank == i);
	    lp = (p4est_locidx_t *) sc_array_push (&(lrank->shared_nodes));
	    *lp = gidx;

	    P4EST_ASSERT (share_count >= 2);
	    proc = *((int *) sc_array_index (inode_sharers,
					     (size_t) share_offset + 1));
	    P4EST_ASSERT (proc != p4est->mpirank);
	    if (proc == i) {
	      shareidx = comm_proc[p4est->mpirank];
	      P4EST_ASSERT (shareidx >= 0);
	      lrank = p4est_lnodes_rank_array_index_int (sharers, shareidx);
	      P4EST_ASSERT (lrank->rank == p4est->mpirank);
	      lp = (p4est_locidx_t *) sc_array_push (&(lrank->shared_nodes));
	      *lp = gidx;
	    }
	  }
	  else {
	    P4EST_ASSERT (inode[0] == i);
	    P4EST_ASSERT (poff[i] <= inode[1] - owned_count &&
			  inode[1] - owned_count < poff[i + 1]);
	    for (l = 0; l < share_count; l++) {
	      proc = *((int *) sc_array_index (inode_sharers,
					       (size_t) share_offset +
					       (size_t) l));
	      shareidx = comm_proc[proc];
	      P4EST_ASSERT (shareidx >= 0);
	      lrank = p4est_lnodes_rank_array_index_int (sharers, shareidx);
	      P4EST_ASSERT (lrank->rank == proc);
	      lp = (p4est_locidx_t *) sc_array_push (&(lrank->shared_nodes));
	      *lp = gidx;
	    }
	  }
	}
      }
    }
  }

  /* for each sharer, figure out which entries in element_nodes are owned by
   * the current process, and which are owned by the sharer's rank */
  for (i = 0; i < comm_proc_count; i++) {
    lrank = p4est_lnodes_rank_array_index_int (sharers, i);
    shared_nodes = &(lrank->shared_nodes);
    count = shared_nodes->elem_count;
    if (count) {
      sc_array_t         *sortshared =
	sc_array_new_size (2 * sizeof (p4est_gloidx_t),
			   count);
      for (zz = 0; zz < count; zz++) {
	p4est_gloidx_t     *gp;

	gidx = *((p4est_locidx_t *) sc_array_index (shared_nodes, zz));
	gp = (p4est_gloidx_t *) sc_array_index (sortshared, zz);
	gp[0] = p4est_lnodes_global_index (lnodes, gidx);
	gp[1] = gidx;
      }
      sc_array_sort (sortshared, p4est_gloidx_compare);
      for (zz = 0; zz < count; zz++) {
	p4est_gloidx_t     *gp;

	gp = (p4est_gloidx_t *) sc_array_index (sortshared, zz);
	*((p4est_locidx_t *) sc_array_index (shared_nodes, zz)) = gp[1];
      }
      sc_array_destroy (sortshared);
    }
    proc = lrank->rank;
    lrank->shared_mine_offset = -1;
    lrank->shared_mine_count = 0;
    for (zz = 0; zz < count; zz++) {
      gidx = *((p4est_locidx_t *) sc_array_index (shared_nodes, zz));
      if (gidx < lnodes->owned_count) {
	if (lrank->shared_mine_count == 0) {
	  lrank->shared_mine_offset = (p4est_locidx_t) zz;
	}
	lrank->shared_mine_count++;
      }
    }
    if (proc == p4est->mpirank) {
      lrank->owned_count = lnodes->owned_count;
      lrank->owned_offset = 0;
    }
    else {
      lrank->owned_offset = poff[proc] + owned_count;
      lrank->owned_count = poff[proc + 1] - poff[proc];
      P4EST_VERBOSEF ("Processor %d shares %llu nodes with processor %d\n",
		      p4est->mpirank, (unsigned long long) count,
		      lrank->rank);
      P4EST_VERBOSEF ("Processor %d owns %d nodes used by processor %d\n",
		      p4est->mpirank, lrank->shared_mine_count, lrank->rank);
      P4EST_VERBOSEF ("Processor %d borrows %d nodes from processor %d\n",
		      p4est->mpirank, lrank->owned_count, lrank->rank);
    }
  }
  P4EST_FREE (comm_proc);
  P4EST_FREE (global_offsets);

  return gtotal;
}

p8est_lnodes_t     *
p8est_lnodes_edge (p8est_t * p8est, p8est_ghost_t * ghost_layer)
{
  p4est_lnodes_data_t data;
  p4est_locidx_t      nel;
  p4est_locidx_t      nlen;
  p4est_lnodes_t     *lnodes = P4EST_ALLOC (p8est_lnodes_t, 1);
  p4est_gloidx_t      gtotal;

  P4EST_GLOBAL_PRODUCTION ("Into " P8EST_STRING "_lnodes_edge\n");
  p4est_log_indent_push ();

  lnodes->mpicomm = p8est->mpicomm;
  lnodes->degree = -4;
  lnodes->num_local_elements = nel = p8est->local_num_quadrants;
  lnodes->vnodes = P8EST_EDGES;

  lnodes->face_code = P4EST_ALLOC_ZERO (p4est_lnodes_code_t, nel);
  nlen = nel * lnodes->vnodes;
  lnodes->element_nodes = P4EST_ALLOC (p4est_locidx_t, nlen);
  memset (lnodes->element_nodes, -1, nlen * sizeof (p4est_locidx_t));

  p4est_lnodes_init_data (&data, p8est, ghost_layer, lnodes);

  p4est_iterate_ext (p8est, ghost_layer, &data, NULL,
		     p4est_lnodes_face_simple_callback,
		     p8est_lnodes_edge_callback, NULL, 1);

  p4est_lnodes_count_send (&data, p8est, lnodes);

  p4est_lnodes_recv (p8est, &data, lnodes);

  gtotal = p4est_lnodes_global_and_sharers (&data, lnodes, p8est);

  p4est_lnodes_reset_data (&data, p8est);

  p4est_log_indent_pop ();
  P4EST_GLOBAL_PRODUCTIONF ("Done " P8EST_STRING "_lnodes_edge with"
			    " %lld global nodes\n",
			    (unsigned long long) gtotal);
  return lnodes;
}
#endif
