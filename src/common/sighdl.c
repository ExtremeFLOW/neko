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
