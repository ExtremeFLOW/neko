/*
 Copyright (c) 2026, The Neko Authors
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

   * Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
   * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
   * Neither the name of the authors nor the names of its contributors may be
     used to endorse or promote products derived from this software without
     specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES ARE DISCLAIMED.
*/

#include <device/device_config.h>

typedef struct {
  enum { IDP_BUFFER, IDP_INTEGER, IDP_SCALAR } type;
  union {
    void *buffer;
    int integer;
    real scalar;
  };
} idp_arg;

#define RBUF(p) ((idp_arg){.type = IDP_BUFFER, .buffer = (p)})
#define IBUF(p) RBUF(p)
#define RVAL(x) ((idp_arg){.type = IDP_SCALAR, .scalar = (x)})
#define IVAL(x) ((idp_arg){.type = IDP_INTEGER, .integer = (x)})

static void idp_launch(const char *name, int groups, int width, const idp_arg *args,
                       unsigned count);
#define IDP_LAUNCH(kernel, groups, width, ...)                                         \
  do {                                                                                 \
    const idp_arg args[] = {__VA_ARGS__};                                              \
    idp_launch(#kernel, groups, width, args, sizeof(args) / sizeof(args[0]));          \
  } while (0)
#define IDP_ZERO(p, n) IDP_LAUNCH(idp_zero_kernel, idp_blocks(n), 256, RBUF(p), IVAL(n))

enum { validation_summary = 1, observation_summary = 2, limiter_status_summary = 3 };
#include "euler_idp_launch.h"

void gpu_entropy_visc_fraction(void *fraction, void *coeff, void *h, void *wave,
                               real *c, int *n) {
  IDP_LAUNCH(idp_fraction_kernel, idp_blocks(*n), 256, RBUF(fraction), RBUF(coeff),
             RBUF(h), RBUF(wave), RVAL(*c), IVAL(*n));
}
