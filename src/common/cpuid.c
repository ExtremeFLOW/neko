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

/**
 * Return the systems' CPU id
 */ 
void system_cpuid(char *name, int len) {  
#if defined(__APPLE__)
  sysctlbyname("machdep.cpu.brand_string", name,(size_t *) &len, NULL, 0);
#elif defined(__x86_64__) && defined(HAVE_CPUID_H)
  __get_cpuid(0x80000002, name+0x0, name+0x1, name+0x2, name+0x3);
  __get_cpuid(0x80000003, name+0x4, name+0x5, name+0x6, name+0x7);
  __get_cpuid(0x80000004, name+0x8, name+0x9, name+0xa, name+0xb);
#else
  strncpy(name, "Unknown CPU", len);
#endif
}

