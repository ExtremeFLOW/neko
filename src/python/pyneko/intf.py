from ctypes import CDLL, util, c_char_p, c_int, POINTER, byref, cast, create_string_buffer
import json

libneko = CDLL(util.find_library("neko"))

libneko.neko_init.resType = None
libneko.neko_job_info.resType = None
#neko.finalize.resType = None

#neko.solve.argtypes = [POINTER(c_char_p), c_int]
#neko.solve.resType = None

def init():
    libneko.neko_init()

def job_info():
    libneko.neko_job_info()

#def finalize():
#    neko.finalize()

#def solve(case_json):
#    cp = python_dict_to_fortran(case_json)
#    neko.solve(byref(cp), len(json.dumps(case_json)))


#
# https://degenerateconic.com/fortran-json-python.html 
#

#def python_str_to_fortran(s):
#    return cast(create_string_buffer(s.encode()),c_char_p)

#def python_dict_to_fortran(d):
#    return python_str_to_fortran(json.dumps(d))
