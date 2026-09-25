/*
 Copyright (c) 2026, The Neko Authors
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.

   * Redistributions in binary form must reproduce the above
     copyright notice, this list of conditions and the following
     disclaimer in the documentation and/or other materials provided
     with the distribution.

   * Neither the name of the authors nor the names of its
     contributors may be used to endorse or promote products derived
     from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 POSSIBILITY OF SUCH DAMAGE.
*/

/**
 * Native Tofu (uTofu) helper routines for the gs_utofu gather-scatter
 * backend (see gs_utofu.F90).
 *
 * The transport is a set of one-sided RDMA puts. Each rank registers its
 * receive buffer once and advertises the (vcq_id, stadd) to its neighbours.
 * A send is then a single utofu_put straight into the neighbour's receive
 * buffer, carrying an `edata` tag in the descriptor; the put's arrival
 * raises a remote MRQ notice carrying that tag, so the *signal is fused
 * into the data transfer* -- no separate atomic/credit message, unlike the
 * coarray backend.
 *
 * Injection is non-blocking and thread-parallel: one injection VCQ is
 * created per OpenMP thread (dealt round-robin over the TNIs), and each
 * thread fires the puts for its share of the peers on its own VCQ -- the
 * uTofu substitute for the per-thread MPI_Isend that Fujitsu MPI's missing
 * MPI_THREAD_MULTIPLE support rules out (cf. the RIKEN-LQCD/jacobi2d
 * reference program). Every VCQ is created with UTOFU_VCQ_FLAG_THREAD_SAFE;
 * without it uTofu takes the unlocked direct-descriptor injection path and
 * any cross-thread touch of the VCQ corrupts its TOQ ("TOQ Direct
 * Descriptor Exception"). With one thread per VCQ the internal lock is
 * uncontended. The number of TNIs to deal the VCQs over is chosen by the
 * caller (NEKO_GS_UTOFU_NTNI, default 1 -- multi-TNI is opt-in until
 * validated on a given machine); NEKO_GS_UTOFU_MASTER_INJECT=1 serialises
 * all injection onto thread 0 as an A/B reference.
 *
 * Build: requires Fujitsu's uTofu headers/runtime (Fugaku). Only compiled
 * when ENABLE_UTOFU. Link with -ltofucom (verify the exact library name in
 * your environment).
 *
 * Items marked "VERIFY" depend on exact utofu.h symbol/field names that
 * vary slightly between Tofu library versions and should be checked when
 * first building on the target system.
 */

#include <stdlib.h>
#include <stdint.h>
#include <utofu.h>

/* Put flags: notice on the local TCQ (so we can reclaim the send buffer)
 * and on the remote MRQ (so the receiver learns of arrival). Cache
 * injection (added at init unless NEKO_GS_UTOFU_CACHE_INJECT=0) drops the
 * received slab straight into the peer's cache -- a win if it is consumed
 * immediately, a loss if it pollutes cache, hence the runtime toggle.
 *
 * UTOFU_ONESIDED_FLAG_STRONG_ORDER (which jacobi2d sets) is deliberately
 * NOT used: it is only required when arrival is detected by watchdog-
 * polling the payload's last word, which needs the data to land in address
 * order. The RMT_PUT MRQ notice used here is raised only after the whole
 * put has been written at the target, so it already carries that
 * guarantee, and strong ordering would just restrict the NIC's packet-
 * level freedom on multi-packet slabs. Revisit if the receive side ever
 * moves to jacobi2d-style payload watchdogs. */
#define GS_UTOFU_PUT_FLAGS_BASE                                          \
  (UTOFU_ONESIDED_FLAG_TCQ_NOTICE | UTOFU_ONESIDED_FLAG_REMOTE_MRQ_NOTICE)

/* --- process-global endpoint state, created once ------------------------ */
static int             gs_utofu_ready = 0;
static size_t          gs_utofu_nvcq  = 0;     /* number of injection VCQs  */
static size_t          gs_utofu_ntni  = 0;     /* TNIs used for injection   */
static utofu_vcq_hdl_t *gs_utofu_vcq    = NULL; /* [nvcq] injection VCQs    */
static utofu_tni_id_t  *gs_utofu_alltni = NULL; /* [nalltni] ALL 1-sided TNIs;
                                                 * first ntni = injection set */
static size_t           gs_utofu_nalltni = 0;
static int              gs_utofu_recv_rr = 0;  /* round-robins recv-VCQ TNIs */
static uint64_t         gs_utofu_edata_max = 255;
static unsigned long    gs_utofu_put_flags = GS_UTOFU_PUT_FLAGS_BASE;
static int              gs_utofu_master_inject = 0;
static size_t           gs_utofu_nrvcq_def = 1;  /* recv VCQs per ctx (req.) */

/* Longs per counter slot: one A64FX L1/L2 cache line (256 bytes), so the
 * per-VCQ completion counters never false-share between threads. */
#define GS_UTOFU_PEND_STRIDE 32

/* --- per-instance context (one per gs_utofu_t) -------------------------- */
typedef struct {
  utofu_stadd_t  *send_stadd;  /* [nvcq] send_buf registered on each inj VCQ */
  /* This instance's OWN receive VCQs, spread over the TNIs. Each receive
   * peer is bound to one of them at init (the receiver advertises that
   * VCQ's id/STADD to the peer), so incoming halo traffic is processed by
   * several TNIs instead of funnelling through one -- the receive-side
   * counterpart of the per-thread injection pool. All are polled by the
   * master; the edata tags are per-instance-unique, so which MRQ a notice
   * appears on carries no protocol meaning. */
  size_t          nrvcq;       /* receive VCQs actually created (>= 1)      */
  utofu_vcq_hdl_t *recv_vcq;   /* [nrvcq]                                   */
  utofu_stadd_t   *recv_stadd; /* [nrvcq] recv_buf registered on each       */
  /* Outstanding (un-reaped) put completions, one padded slot per
   * (injection VCQ, send-buffer half). The put's cbdata points at its
   * slot, so draining a VCQ's TCQ decrements the right counter -- letting
   * nbsend wait only on the half it is about to reuse. The phase barriers
   * in gs_nbsend_utofu guarantee each VCQ (and thus each slot) has a
   * single toucher at any time, so no atomics are needed. */
  long          *send_pending; /* [2 * nvcq * GS_UTOFU_PEND_STRIDE]        */
  /* notices observed early that belong to the *other* parity (next round)  */
  uint64_t      *stash;
  int            stash_n;
  int            stash_cap;
} gs_utofu_ctx_t;

static long *gs_utofu_pend_slot(gs_utofu_ctx_t *ctx, size_t k, int half)
{
  return &ctx->send_pending[(2 * k + (size_t) half) * GS_UTOFU_PEND_STRIDE];
}

/* Drain every completion currently sitting on VCQ k's transmit queue,
 * decrementing the per-(instance, VCQ, half) counter each put carried as
 * cbdata. Must only be called by VCQ k's single toucher of the current
 * phase. */
static void gs_utofu_drain_tcq(size_t k)
{
  void *cbdata;
  while (utofu_poll_tcq(gs_utofu_vcq[k], 0, &cbdata) == UTOFU_SUCCESS) {
    if (cbdata != NULL) (*(long *) cbdata)--;
  }
}

static void gs_utofu_stash_push(gs_utofu_ctx_t *ctx, uint64_t edata)
{
  if (ctx->stash_n == ctx->stash_cap) {
    int cap = ctx->stash_cap ? 2 * ctx->stash_cap : 16;
    ctx->stash = realloc(ctx->stash, cap * sizeof(*ctx->stash));
    ctx->stash_cap = cap;
  }
  ctx->stash[ctx->stash_n++] = edata;
}

/**
 * Create the process-global uTofu endpoint (idempotent).
 *
 * Creates the injection VCQs -- one per OpenMP thread, dealt round-robin
 * over the granted TNIs, every one THREAD_SAFE -- shared by all instances
 * and driven thread-parallel (thread t owns VCQ t). Receive VCQs are
 * created per gs instance in gs_utofu_ctx_create, so each instance owns a
 * private MRQ.
 *
 * @param ntni_req    requested number of TNIs (>=1); clamped to availability
 * @param nthreads    OpenMP thread count = requested injection VCQs (>=1)
 * @param ntni_out    number of TNIs actually in use
 * @param nvcq_out    number of injection VCQs actually created
 * @param nrvcq_out   receive VCQs requested per instance path (the caller
 *                    sizes its id/STADD arrays with this; ctx_create may
 *                    grant fewer)
 * @param edata_max   largest usable edata value (for tag-range checking)
 */
void gs_utofu_init(int ntni_req, int nthreads, int *ntni_out, int *nvcq_out,
                   int *nrvcq_out, uint64_t *edata_max, int *ierr)
{
  *ierr = 0;

  if (!gs_utofu_ready) {
    utofu_tni_id_t *tnis = NULL;
    size_t num_tnis = 0;

    if (utofu_get_onesided_tnis(&tnis, &num_tnis) != UTOFU_SUCCESS
        || num_tnis == 0) {
      *ierr = 1;
      return;
    }

    size_t want_tni = (ntni_req > 0) ? (size_t) ntni_req : 1;
    if (want_tni > num_tnis) want_tni = num_tnis;
    gs_utofu_ntni    = want_tni;
    gs_utofu_nalltni = num_tnis;
    gs_utofu_alltni  = malloc(num_tnis * sizeof(*gs_utofu_alltni));
    for (size_t k = 0; k < num_tnis; k++)
      gs_utofu_alltni[k] = tnis[k];

    /* One injection VCQ per thread, capped by NEKO_GS_UTOFU_NVCQ. The VCQ
     * budget is shared with the per-instance receive VCQs created later
     * (one scalar + one vector per gs instance) -- on Fugaku a TNI grants
     * us ~8 VCQs, so a 12-thread pool can eat a whole TNI; the cap lets a
     * tight budget be rebalanced without rebuilding. THREAD_SAFE is what
     * makes multithreaded operation defined at all: without it uTofu takes
     * the unlocked direct-descriptor injection path and any cross-thread
     * touch corrupts the TOQ (the "TOQ Direct Descriptor Exception" of the
     * first multithreaded attempt). With one thread per VCQ the internal
     * lock is uncontended. */
    size_t want_vcq = (nthreads > 0) ? (size_t) nthreads : 1;
    const char *nv = getenv("NEKO_GS_UTOFU_NVCQ");
    if (nv != NULL) {
      long cap = strtol(nv, NULL, 10);
      if (cap > 0 && (size_t) cap < want_vcq)
        want_vcq = (size_t) cap;
    }

    /* Preferred TNI first, then the rest of the injection set: one
     * exhausted TNI must not end the pool while others still have room.
     * If every TNI in the set is out of VCQs, keep what was granted:
     * threads beyond the count simply do not inject (see
     * gs_utofu_post_puts), and the modulo peer partition shrinks with it. */
    gs_utofu_vcq = malloc(want_vcq * sizeof(*gs_utofu_vcq));
    while (gs_utofu_nvcq < want_vcq) {
      int created = 0;
      for (size_t s = 0; s < want_tni && !created; s++) {
        if (utofu_create_vcq(gs_utofu_alltni[(gs_utofu_nvcq + s) % want_tni],
                             UTOFU_VCQ_FLAG_THREAD_SAFE,
                             &gs_utofu_vcq[gs_utofu_nvcq]) == UTOFU_SUCCESS)
          created = 1;
      }
      if (!created)
        break;
      gs_utofu_nvcq++;
    }
    if (gs_utofu_nvcq == 0) {
      *ierr = 2;
      return;
    }

    /* Start the receive-VCQ round-robin just past the injection TNIs, so
     * receive queues prefer interfaces the pool has not loaded. */
    gs_utofu_recv_rr = (int) (gs_utofu_ntni % gs_utofu_nalltni);

    /* Receive VCQs per instance path: default one per TNI, so an
     * instance's incoming halo traffic is processed by every interface
     * instead of funnelling through one. NEKO_GS_UTOFU_NRVCQ overrides
     * (=1 restores the single-receive-VCQ behaviour). The budget is
     * shared with the injection pool and across instances/paths;
     * ctx_create degrades gracefully to fewer when it runs dry. */
    gs_utofu_nrvcq_def = gs_utofu_nalltni;
    const char *nr = getenv("NEKO_GS_UTOFU_NRVCQ");
    if (nr != NULL) {
      long v = strtol(nr, NULL, 10);
      if (v > 0 && v <= 64) gs_utofu_nrvcq_def = (size_t) v;
    }

    /* Largest edata the hardware will carry. VERIFY: field name and whether
     * the capability is expressed in bytes (assumed here). */
    struct utofu_onesided_caps *caps = NULL;
    if (utofu_query_onesided_caps(tnis[0], &caps) == UTOFU_SUCCESS
        && caps != NULL) {
      size_t eb = caps->max_edata_size;
      /* Capped to a positive signed-64 value: the Fortran side receives this
       * as integer(c_int64_t) and compares against it. */
      if (eb >= 8) gs_utofu_edata_max = (uint64_t) INT64_MAX;
      else if (eb > 0) gs_utofu_edata_max = (1ULL << (8 * eb)) - 1ULL;
    }

    /* Cache injection on by default; NEKO_GS_UTOFU_CACHE_INJECT=0 disables. */
    const char *ci = getenv("NEKO_GS_UTOFU_CACHE_INJECT");
    gs_utofu_put_flags = GS_UTOFU_PUT_FLAGS_BASE;
    if (ci == NULL || ci[0] != '0')
      gs_utofu_put_flags |= UTOFU_ONESIDED_FLAG_CACHE_INJECTION;

    /* A/B switch: serialise all injection onto thread 0 (still spread over
     * every VCQ, which is legal now that they are THREAD_SAFE) to isolate
     * the effect of parallel injection against the same transport. */
    const char *mi = getenv("NEKO_GS_UTOFU_MASTER_INJECT");
    gs_utofu_master_inject = (mi != NULL && mi[0] != '0');

    free(tnis);
    gs_utofu_ready = 1;
  }

  *edata_max   = gs_utofu_edata_max;
  *ntni_out    = (int) gs_utofu_ntni;
  *nvcq_out    = (int) gs_utofu_nvcq;
  *nrvcq_out   = (int) gs_utofu_nrvcq_def;
}

/**
 * Register an instance's send/recv buffers and allocate its context.
 * send_buf is registered on every injection VCQ (the source STADD must
 * belong to the VCQ that issues the put). The instance gets its OWN set of
 * receive VCQs (private MRQs); recv_buf is registered on each, and the
 * caller binds every receive peer to one of them, advertising that VCQ's
 * id/STADD as the peer's put target -- so arrivals land only in this
 * instance's queues, spread over the TNIs.
 *
 * recv_vcq_id / recv_stadd are arrays of (at least) the nrvcq_out value
 * reported by gs_utofu_init; entries [0, *nrvcq_out) are filled. Fewer
 * than requested may be granted when the VCQ budget runs dry (>= 1, or
 * ierr = 1 -- remedy: lower OMP_NUM_THREADS, cap NEKO_GS_UTOFU_NVCQ, or
 * reduce NEKO_GS_UTOFU_NRVCQ).
 */
void gs_utofu_ctx_create(void *send_buf, size_t send_bytes,
                         void *recv_buf, size_t recv_bytes,
                         void **ctx_out, int *nrvcq_out,
                         uint64_t *recv_vcq_id, uint64_t *recv_stadd,
                         int *ierr)
{
  *ierr = 0;
  gs_utofu_ctx_t *ctx = calloc(1, sizeof(*ctx));
  ctx->send_stadd = malloc(gs_utofu_nvcq * sizeof(*ctx->send_stadd));
  ctx->send_pending = calloc(2 * gs_utofu_nvcq * GS_UTOFU_PEND_STRIDE,
                             sizeof(*ctx->send_pending));
  ctx->recv_vcq   = malloc(gs_utofu_nrvcq_def * sizeof(*ctx->recv_vcq));
  ctx->recv_stadd = malloc(gs_utofu_nrvcq_def * sizeof(*ctx->recv_stadd));

  /* Receive VCQs may sit on ANY one-sided TNI: incoming puts are routed
   * by the advertised VCQ id, so spreading the receive processing over all
   * interfaces is free (NEKO_GS_UTOFU_NTNI deliberately confines only the
   * injection side). Round-robin over the full set, sweeping past
   * exhausted TNIs -- the injection pool can eat a whole TNI's VCQ budget
   * (~8 per TNI on Fugaku). Stop early when the budget runs dry; at least
   * one must succeed, since a private MRQ per instance is a correctness
   * requirement. They stay non-THREAD_SAFE (flags 0): they are strictly
   * master-polled, and skipping the internal lock keeps the MRQ sweep --
   * the receive hot path -- cheap. Never touch them off-master. */
  while (ctx->nrvcq < gs_utofu_nrvcq_def) {
    int created = 0;
    for (size_t s = 0; s < gs_utofu_nalltni && !created; s++) {
      utofu_tni_id_t rtni =
        gs_utofu_alltni[((size_t) gs_utofu_recv_rr + s) % gs_utofu_nalltni];
      if (utofu_create_vcq(rtni, 0, &ctx->recv_vcq[ctx->nrvcq])
          == UTOFU_SUCCESS)
        created = 1;
    }
    gs_utofu_recv_rr++;
    if (!created)
      break;
    ctx->nrvcq++;
  }
  if (ctx->nrvcq == 0) {
    *ierr = 1;
    return;
  }

  for (size_t k = 0; k < gs_utofu_nvcq; k++) {
    if (utofu_reg_mem(gs_utofu_vcq[k], send_buf, send_bytes, 0,
                      &ctx->send_stadd[k]) != UTOFU_SUCCESS) {
      *ierr = 2;
      return;
    }
  }
  for (size_t r = 0; r < ctx->nrvcq; r++) {
    if (utofu_reg_mem(ctx->recv_vcq[r], recv_buf, recv_bytes, 0,
                      &ctx->recv_stadd[r]) != UTOFU_SUCCESS) {
      *ierr = 3;
      return;
    }
  }

  for (size_t r = 0; r < ctx->nrvcq; r++) {
    utofu_vcq_id_t vid;
    if (utofu_query_vcq_id(ctx->recv_vcq[r], &vid) != UTOFU_SUCCESS) {
      *ierr = 4;
      return;
    }
    recv_vcq_id[r] = (uint64_t) vid;
    recv_stadd[r]  = (uint64_t) ctx->recv_stadd[r];
  }

  *nrvcq_out = (int) ctx->nrvcq;
  *ctx_out = ctx;
}

void gs_utofu_ctx_free(void *vctx, int *ierr)
{
  *ierr = 0;
  gs_utofu_ctx_t *ctx = vctx;
  if (ctx == NULL) return;

  /* Make sure no put is still reading send_buf before we deregister it.
   * Called from a serial region, so one thread draining every VCQ is safe
   * (the VCQs are THREAD_SAFE, and there is no concurrent toucher). */
  for (int h = 0; h < 2; h++)
    for (size_t k = 0; k < gs_utofu_nvcq; k++)
      while (*gs_utofu_pend_slot(ctx, k, h) > 0)
        gs_utofu_drain_tcq(k);

  for (size_t k = 0; k < gs_utofu_nvcq; k++)
    utofu_dereg_mem(gs_utofu_vcq[k], ctx->send_stadd[k], 0);
  for (size_t r = 0; r < ctx->nrvcq; r++) {
    utofu_dereg_mem(ctx->recv_vcq[r], ctx->recv_stadd[r], 0);
    utofu_free_vcq(ctx->recv_vcq[r]);
  }

  free(ctx->send_stadd);
  free(ctx->recv_vcq);
  free(ctx->recv_stadd);
  free(ctx->send_pending);
  free(ctx->stash);
  free(ctx);
}

/* Post peer i's put on injection VCQ k, charging its completion to the
 * (k, parity) counter slot so the next reuse of this send-buffer half waits
 * on exactly its own puts. The caller must be VCQ k's single toucher of the
 * current phase (the BUSY-retry path drains k's TCQ). */
static int gs_utofu_put_one(gs_utofu_ctx_t *ctx, size_t k, int i,
                            const int64_t *send_off, const int64_t *send_len,
                            const uint64_t *rmt_vcq_id,
                            const uint64_t *rmt_stadd,
                            const int64_t *dst_off0, const int64_t *dst_off1,
                            const int *dst_tag, int parity, int64_t send_base,
                            int elem_size)
{
  const int64_t doff = parity ? dst_off1[i] : dst_off0[i];
  long *pend = gs_utofu_pend_slot(ctx, k, parity);

  utofu_stadd_t lcl = ctx->send_stadd[k]
    + (utofu_stadd_t)((send_base + send_off[i]) * elem_size);
  utofu_stadd_t rmt = (utofu_stadd_t) rmt_stadd[i]
    + (utofu_stadd_t)(doff * elem_size);
  size_t   len   = (size_t)(send_len[i] * elem_size);
  uint64_t edata = (uint64_t) dst_tag[i] * 2u + (uint64_t) parity;

  int rc;
  do {
    rc = utofu_put(gs_utofu_vcq[k], (utofu_vcq_id_t) rmt_vcq_id[i],
                   lcl, rmt, len, edata, gs_utofu_put_flags, pend);
    /* TOQ full (reported as UTOFU_ERR_BUSY): reclaim completions, retry. */
    if (rc == UTOFU_ERR_BUSY)
      gs_utofu_drain_tcq(k);
  } while (rc == UTOFU_ERR_BUSY);

  if (rc == UTOFU_SUCCESS)
    (*pend)++;
  return rc;
}

/**
 * Issue this thread's share of the puts. Called by EVERY thread of the
 * OpenMP team, each passing its own tid: thread tid posts peers
 * i = tid, tid+npost, ... on its own VCQ, where npost = min(nteam, nvcq),
 * so threads beyond the VCQ count post nothing and their peers fall to the
 * others. One thread per VCQ per phase means the completion counters need
 * no atomics; the VCQs are THREAD_SAFE regardless, so a degraded or odd
 * configuration is safe, just slower. Offsets/lengths are in elements;
 * elem_size converts to bytes. The edata tag fuses the slab id (dst_tag)
 * with the parity bit so the receiver can tell this round's notices from a
 * neighbour that has already raced into the next round.
 *
 * With NEKO_GS_UTOFU_MASTER_INJECT=1, thread 0 instead posts every peer,
 * round-robined over all VCQs (the pre-threading behaviour, kept as an A/B
 * reference).
 */
void gs_utofu_post_puts(void *vctx, int tid, int nteam, int nsend,
                        const int64_t *send_off, const int64_t *send_len,
                        const uint64_t *rmt_vcq_id, const uint64_t *rmt_stadd,
                        const int64_t *dst_off0, const int64_t *dst_off1,
                        const int *dst_tag, int parity, int64_t send_base,
                        int elem_size, int *ierr)
{
  *ierr = 0;
  gs_utofu_ctx_t *ctx = vctx;

  if (gs_utofu_master_inject) {
    if (tid != 0) return;
    for (int i = 0; i < nsend; i++) {
      size_t k = (size_t)(i % (int) gs_utofu_nvcq);
      if (gs_utofu_put_one(ctx, k, i, send_off, send_len, rmt_vcq_id,
                           rmt_stadd, dst_off0, dst_off1, dst_tag, parity,
                           send_base, elem_size) != UTOFU_SUCCESS) {
        *ierr = 1;
        return;
      }
    }
    return;
  }

  int npost = (nteam < (int) gs_utofu_nvcq) ? nteam : (int) gs_utofu_nvcq;
  if (tid >= npost) return;

  for (int i = tid; i < nsend; i += npost) {
    if (gs_utofu_put_one(ctx, (size_t) tid, i, send_off, send_len,
                         rmt_vcq_id, rmt_stadd, dst_off0, dst_off1, dst_tag,
                         parity, send_base, elem_size) != UTOFU_SUCCESS) {
      *ierr = 1;
      return;
    }
  }
}

/**
 * Collect the receive slabs that have landed for the current parity.
 *
 * Returns at most `cap` 1-based receive-peer indices in out_idx. A notice
 * whose parity bit does not match the current round is stashed and replayed
 * on the next call (with the flipped parity), keeping the unpack-as-ready
 * model correct in the face of a neighbour running one round ahead. May
 * return 0 (nothing yet) -- the caller polls in a loop.
 */
void gs_utofu_poll_recv(void *vctx, int parity, int cap,
                        int *ncompleted, int *out_idx, int *ierr)
{
  *ierr = 0;
  gs_utofu_ctx_t *ctx = vctx;
  int n = 0;

  /* Replay stashed notices that now match the current parity. */
  int w = 0;
  for (int s = 0; s < ctx->stash_n; s++) {
    uint64_t ed = ctx->stash[s];
    if ((int)(ed & 1u) == parity && n < cap)
      out_idx[n++] = (int)(ed >> 1) + 1;
    else
      ctx->stash[w++] = ed;
  }
  ctx->stash_n = w;

  /* Drain currently-available arrival notices from every one of this
   * instance's MRQs. Each peer's notices always land on its bound VCQ, but
   * the tags are instance-unique, so the queues can be drained in any
   * order into the one result list / stash. */
  struct utofu_mrq_notice notice;
  for (size_t r = 0; r < ctx->nrvcq && n < cap; r++) {
    while (n < cap) {
      int rc = utofu_poll_mrq(ctx->recv_vcq[r], 0, &notice);
      if (rc == UTOFU_ERR_NOT_FOUND) break;
      if (rc != UTOFU_SUCCESS) {
        *ierr = 1;
        return;
      }
      /* VERIFY: RMT_PUT is the receiver-side notice for an incoming put. */
      if (notice.notice_type != UTOFU_MRQ_TYPE_RMT_PUT) continue;

      uint64_t ed = notice.edata;
      if ((int)(ed & 1u) == parity)
        out_idx[n++] = (int)(ed >> 1) + 1;
      else
        gs_utofu_stash_push(ctx, ed);
    }
  }

  *ncompleted = n;
}

/** Block until every put charged to send-buffer half `half` has completed
 *  locally, so that half may be repacked. Called by EVERY thread of the
 *  team: the strided map k = tid, tid+nteam, ... covers each injection VCQ
 *  exactly once, so each VCQ (and its counter slots) has a single toucher
 *  and no atomics are needed. The caller must barrier afterwards -- only
 *  then is the whole half known to be free. With double buffering the half
 *  being reclaimed was posted two rounds ago, so this normally just reaps
 *  already-finished completions. */
void gs_utofu_drain_half(void *vctx, int half, int tid, int nteam, int *ierr)
{
  *ierr = 0;
  gs_utofu_ctx_t *ctx = vctx;
  for (size_t k = (size_t) tid; k < gs_utofu_nvcq; k += (size_t) nteam)
    while (*gs_utofu_pend_slot(ctx, k, half) > 0)
      gs_utofu_drain_tcq(k);
}
