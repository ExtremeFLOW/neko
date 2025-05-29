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
void system_cpuid(char *name, size_t len) {
#if defined(__APPLE__)
  sysctlbyname("machdep.cpu.brand_string", name, &len, NULL, 0);
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
  int cpuarm = 0;
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
      else if(strstr(token, "0x41")) {
	cpuarm = 1;
	continue;
      }
    }

    if (strstr(buf, "CPU part") && cpufj) {
      char *token = strtok (buf, delim);
      token = strtok (NULL, delim);
      if (strstr(token, "0x001")) {
        strncpy(name, "Fujitsu A64FX", len);
        break;
      }
    }

    if (strstr(buf, "CPU part") && cpuarm) {
      char *token = strtok (buf, delim);
      token = strtok (NULL, delim);
      if (strstr(token, "0xd03")) {
        strncpy(name, "ARM Cortex A53", len);
        break;
      }
      else if (strstr(token, "0xd04")) {
        strncpy(name, "ARM Cortex A35", len);
        break;
      }
      else if (strstr(token, "0xd05")) {
        strncpy(name, "ARM Cortex A55", len);
        break;
      }
      else if (strstr(token, "0xd07")) {
        strncpy(name, "ARM Cortex A57", len);
        break;
      }
      else if (strstr(token, "0xd08")) {
        strncpy(name, "ARM Cortex A72", len);
        break;
      }
      else if (strstr(token, "0xd09")) {
        strncpy(name, "ARM Cortex A73", len);
        break;
      }
      else if (strstr(token, "0xd0a")) {
        strncpy(name, "ARM Cortex A75", len);
        break;
      }
      else if (strstr(token, "0xd0b")) {
        strncpy(name, "ARM Cortex A76", len);
        break;
      }
      else if (strstr(token, "0xd0c")) {
        strncpy(name, "ARM Neoverse N1", len);
        break;
      }
      else if (strstr(token, "0xd0d")) {
        strncpy(name, "ARM Cortex A77", len);
        break;
      }
      else if (strstr(token, "0xd0e")) {
        strncpy(name, "ARM Cortex A76AE", len);
        break;
      }
      else if (strstr(token, "0xd40")) {
        strncpy(name, "ARM Neoverse V1", len);
        break;
      }
      else if (strstr(token, "0xd41")) {
        strncpy(name, "ARM Cortex A78", len);
        break;
      }
      else if (strstr(token, "0xd42")) {
        strncpy(name, "ARM Cortex A78AE", len);
        break;
      }
      else if (strstr(token, "0xd44")) {
        strncpy(name, "ARM Cortex X1", len);
        break;
      }
      else if (strstr(token, "0xd46")) {
        strncpy(name, "ARM Cortex A510", len);
        break;
      }
      else if (strstr(token, "0xd4c")) {
        strncpy(name, "ARM Cortex X1C", len);
        break;
      }
      else if (strstr(token, "0xd80")) {
        strncpy(name, "ARM Cortex A520", len);
        break;
      }
      else if (strstr(token, "0xd47")) {
        strncpy(name, "ARM Cortex A520", len);
        break;
      }
      else if (strstr(token, "0xd4d")) {
        strncpy(name, "ARM Cortex A715", len);
        break;
      }
      else if (strstr(token, "0xd48")) {
        strncpy(name, "ARM Cortex X2", len);
        break;
      }
      else if (strstr(token, "0xd49")) {
        strncpy(name, "ARM Neoverse N2", len);
        break;
      }
      else if (strstr(token, "0xd4b")) {
        strncpy(name, "ARM Cortex A78C", len);
        break;
      }
      else if (strstr(token, "0xd4c")) {
        strncpy(name, "ARM Cortex X1C", len);
        break;
      }
      else if (strstr(token, "0xd4e")) {
        strncpy(name, "ARM Cortex X3", len);
        break;
      }
      else if (strstr(token, "0xd4f")) {
        strncpy(name, "ARM Neoverse V2", len);
        break;
      }
      else if (strstr(token, "0xd81")) {
        strncpy(name, "ARM Cortex A720", len);
        break;
      }
      else if (strstr(token, "0xd82")) {
        strncpy(name, "ARM Cortex X4", len);
        break;
      }
      else if (strstr(token, "0xd84")) {
        strncpy(name, "ARM Neoverse V3", len);
        break;
      }
      else if (strstr(token, "0xd85")) {
        strncpy(name, "ARM Cortex X925", len);
        break;
      }
      else if (strstr(token, "0xd87")) {
        strncpy(name, "ARM Cortex A725", len);
        break;
      }
      else if (strstr(token, "0xd8e")) {
        strncpy(name, "ARM Neoverse N3", len);
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
