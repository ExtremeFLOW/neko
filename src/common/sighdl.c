/*
 Copyright (c) 2021, The Neko Authors
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

#include <signal.h>
#include <stdint.h>
#include <time.h>
#include <sys/time.h>

#define SIGHDL_XCPU 0x01
#define SIGHDL_ALRM 0x02

#define SIGHDL_USR1 0x01
#define SIGHDL_USR2 0x02

volatile uint8_t TIMEOUT = 0;
volatile uint8_t USR = 0;

/*!
  Generic signal handler
*/
void sighdl(int sig_code) {
  switch (sig_code) {
  case SIGXCPU:
    TIMEOUT |= SIGHDL_XCPU;
    break;
  case SIGALRM:
    TIMEOUT |= SIGHDL_ALRM;
    break;
  case SIGUSR1:
    USR |= SIGHDL_USR1;
    break;
  case SIGUSR2:
    USR |= SIGHDL_USR2;
    break;
  default:
    break;
  }
}

/*!
  Check if a timeout has occurred
*/
int8_t sighdl_timeout() {
  return TIMEOUT;
}

/*!
  Check if a user signal has been raised
*/
int8_t sighdl_usr() {
  return USR;
}


/*!
  Set a timer
*/
int sighdl_set_timeout(int *sec) {
  struct itimerval itv;
  struct sigaction sig_param;

  TIMEOUT &= ~SIGHDL_ALRM;
  sig_param.sa_handler = sighdl;
  sigemptyset(&sig_param.sa_mask);
  sig_param.sa_flags = SA_RESTART;
  if (sigaction(SIGALRM, &sig_param, 0) < 0)
    return -1;

  itv.it_value.tv_sec = (*sec);
  itv.it_value.tv_usec = 0;
  itv.it_interval.tv_sec = 0;
  itv.it_interval.tv_usec = 0;
  return setitimer(ITIMER_REAL, &itv, 0);
}

/*!
  Trap CPU time limit (SIGXCPU)
*/
int sighdl_trap_cpulimit() {
  struct sigaction sig_param;

  TIMEOUT &= ~SIGHDL_XCPU; 
  sig_param.sa_handler = sighdl;
  sigemptyset(&sig_param.sa_mask);
  sig_param.sa_flags = SA_RESTART;
  return sigaction(SIGXCPU, &sig_param, 0);    
}

/*!
  Trap user signals (SIGUSR1, SIGUSR2)
*/
int sighdl_trap_usr() {
  struct sigaction sig_param;

  USR = 0;  
  sig_param.sa_handler = sighdl;
  sigemptyset(&sig_param.sa_mask);
  sig_param.sa_flags = SA_RESTART;
  if (sigaction(SIGUSR1, &sig_param, 0) < 0)
    return -1;  
  return sigaction(SIGUSR2, &sig_param, 0);    
}
