/*
 * External-linkage wrappers for NVTX.
 *
 * NVTX3 (shipped with CUDA 12.9+) is header-only: the entry points are
 * static inline functions and libnvToolsExt is no longer available. These
 * shims provide real symbols so the Fortran bindings work against both
 * NVTX2 (legacy library) and NVTX3 (header-only) installations.
 */

#ifdef HAVE_NVTX

#ifdef HAVE_NVTX3
#include <nvtx3/nvToolsExt.h>
#else
#include <nvToolsExt.h>
#endif

extern "C" {

  void neko_nvtxRangePushA(const char *name)
  {
    nvtxRangePushA(name);
  }

  void neko_nvtxRangePushEx(const nvtxEventAttributes_t *event)
  {
    nvtxRangePushEx(event);
  }

  void neko_nvtxRangePop(void)
  {
    nvtxRangePop();
  }

}

#endif
