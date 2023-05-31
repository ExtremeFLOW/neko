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

#include <string.h>
#include <stdint.h>

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
#else
  strncpy(name, "Unknown CPU", len);
#endif
}

