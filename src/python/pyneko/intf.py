from ctypes import CDLL, util, c_char_p, c_int, c_double, POINTER, byref, cast, create_string_buffer
import json

libneko = CDLL(util.find_library("neko"))

libneko.neko_init.resType = None
libneko.neko_finalize.resType = None
libneko.neko_job_info.resType = None

libneko.neko_case_init.argtypes = [POINTER(c_char_p), c_int, POINTER(c_int)]
libneko.neko_case_init.resType = c_int

libneko.neko_case_free.argtypes = [POINTER(c_int)]
libneko.neko_case_free.resType = None

libneko.neko_solve.argtypes = [POINTER(c_int)]
libneko.neko_solve.resType = None

libneko.neko_step.argtypes = [POINTER(c_int), c_double, c_int]

def init():
    libneko.neko_init()

def finalize():
    libneko.neko_finalize()


def job_info():
    libneko.neko_job_info()

def case_init(case_json):
    cp = python_dict_to_fortran(case_json)
    case_descr = c_int()
    libneko.neko_case_init(byref(cp), len(json.dumps(case_json)),
                           byref(case_descr))
    return case_descr

def case_free(case_descr):
    libneko.neko_case_free(byref(case_descr))

def solve(case_descr):
    libneko.neko_solve(byref(case_descr))

def step(case_descr, t, tstep):
    libneko.neko_step(byref(case_descr), t, tstep)


#
# https://degenerateconic.com/fortran-json-python.html
#

def python_str_to_fortran(s):
    return cast(create_string_buffer(s.encode()),c_char_p)

def python_dict_to_fortran(d):
    return python_str_to_fortran(json.dumps(d))
