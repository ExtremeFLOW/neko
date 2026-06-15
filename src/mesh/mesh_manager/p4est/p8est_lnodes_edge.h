/*
  This code is based on p4est_lnodes.h from the p4est (2.8.7) source
  code. Its only goal is to provide additional feature not present in
  the original code, which is a separate global numbering of edges.
  To get it I've modified a smallest needed set of routines.
*/

#ifndef P4EST_LNODES_EDGE_H
#define P4EST_LNODES_EDGE_H

#include <p8est_ghost.h>
#include <p8est_lnodes.h>

/** Create a global numbering of edges.
 * \param [in] p8est            Valid forest.
 * \param [in] ghost_layer      Valid full ghost layer, i. e. constructed
 *                              by \ref p8est_ghost_new with the same forest
 *                              and argument \ref P8EST_CONNECT_FULL.
 * \param [in] npts             number of points along the edge
 * \return                      Fully initialized nodes structure.
 */
p8est_lnodes_t     *p8est_lnodes_edge (p8est_t * p8est,
				       p8est_ghost_t * ghost_layer, int npts);

#endif
