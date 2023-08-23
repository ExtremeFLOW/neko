#if defined(__APPLE__)
#include <sys/ioctl.h>
#include <sys/sysctl.h>
#include <sys/types.h>
#include <unistd.h>
#elif defined(__x86_64__)
#if defined(HAVE_CPUID_H)
#include <cpuid.h>
#endif
#endif

#include <stdio.h>
#include <string.h>
#include <stdint.h>

#define MAXLEN 2048

/**
 * Return the systems' CPU id
 */ 
void system_cpuid(char *name, int len) {  
#if defined(__APPLE__)
  sysctlbyname("machdep.cpu.brand_string", name,(size_t *) &len, NULL, 0);
#elif defined(__x86_64__) && defined(HAVE_CPUID_H)
  uint32_t brand[12];
  __get_cpuid(0x80000002, brand+0x0, brand+0x1, brand+0x2, brand+0x3);
  __get_cpuid(0x80000003, brand+0x4, brand+0x5, brand+0x6, brand+0x7);
  __get_cpuid(0x80000004, brand+0x8, brand+0x9, brand+0xa, brand+0xb);
  strncpy(name, (char *) brand, len);
#elif defined(_ARCH_PPC64) || defined(__aarch64__)
  FILE *fp = fopen("/proc/cpuinfo", "r");
  char buf[MAXLEN];
  const char *delim = ":\n";
#if defined(__aarch64__)
  /* Generic ARM unless we found something known */
  strncpy(name, "ARM", len);
  int cpufj = 0;
#endif
  while (fgets (buf, MAXLEN, fp)) {
#if defined(_ARCH_PPC64)
    if(strstr(buf, "cpu")) {
      char *token = strtok (buf, delim);
      token = strtok (NULL, delim);
      strncpy(name, token+1, len);
      break;
    }
#elif defined(__aarch64__)
    if (strstr(buf, "CPU implementer")) {                   
      char *token = strtok (buf, delim);
      token = strtok (NULL, delim);
      if (strstr(token, "0x46")) {
        cpufj = 1;
        continue;
      }
    }

    if (strstr(buf, "CPU part") && cpufj) {
      char *token = strtok (buf, delim);
      token = strtok (NULL, delim);
      if (strstr(token, "0x001")) {
        strncpy(name, "A64FX", len);
        break;
      }
    }   
#endif
  }
  fclose(fp);
#else
  strncpy(name, "Unknown CPU", len);
#endif
}

