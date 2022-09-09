from posix import RTLD_LAZY
from ctypes import CDLL, util

#libneko = CDLL(util.find_library("nekointf"), mode=RTLD_LAZY)
libneko = CDLL("libnekointf.so", mode=RTLD_LAZY)
libneko.intf_init.resType = None

def neko_init():
    libneko.intf_init()



