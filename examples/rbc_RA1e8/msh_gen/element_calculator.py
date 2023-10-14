import numpy as np


R  = 0.05
r  = 0.7*R
RB = 0.978*R
H = 1


# Case 1
NC = 8
NB = 1
NM = 4

# Case 2
NC = 10
NB = 1
NM = 5

# Case 2
#NC = 13
#NB = 1
#NM = 6

bulk_size = r*2
mid_size = (RB - r)*2


nelv_bulk = (NC-1)**2
nelv_mid = (NC-1)*(NM-1)*4
nelv_wall = (NC-1)*NB*4


print("Elements in the bulk is: " +repr(nelv_bulk))
print("Elements in the mid is: " +repr(nelv_mid))
print("Elements in the wall is: " +repr(nelv_wall))
print("Total: "+repr(nelv_bulk+nelv_mid+nelv_wall))


print("size of the bulk =" +repr(bulk_size))
print("size of the mid =" +repr(mid_size))

elem_size_bulk = bulk_size/(NC-1)
elem_size_mid = (mid_size/2)/((NM-1))

print("elem size in bulk = "+repr(elem_size_bulk))
print("elem size in mid = "+repr(elem_size_mid))

avg = (elem_size_bulk+elem_size_mid)/2

nelv_h = int(H / avg)

print("average size in cross section = "+repr(avg))
print("number of elements required in z = "+repr(nelv_h))

print("Total 3d elements : "+repr((nelv_bulk+nelv_mid+nelv_wall)*nelv_h))
