from ctypes import CDLL, util, c_char_p, c_int, POINTER, byref, cast, create_string_buffer
import json

neko = CDLL(util.find_library("nekointf"))

neko.init.resType = None
neko.finalize.resType = None

neko.solve.argtypes = [POINTER(c_char_p), c_int]
neko.solve.resType = None

def init():
    neko.init()

def finalize():
    neko.finalize()

def solve(case_json):
    cp = python_dict_to_fortran(case_json)
    neko.solve(byref(cp), len(json.dumps(case_json)))


#
# https://degenerateconic.com/fortran-json-python.html 
#

def python_str_to_fortran(s):
    return cast(create_string_buffer(s.encode()),c_char_p)

def python_dict_to_fortran(d):
    return python_str_to_fortran(json.dumps(d))
