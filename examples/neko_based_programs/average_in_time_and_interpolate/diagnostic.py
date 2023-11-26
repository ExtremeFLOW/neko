import numpy as np
import math as mt
import csv

file = open("coordinates.csv")
reader = csv.reader(file)
data = list(reader)
xyz_coords = np.array(data,dtype=np.double)
number_of_points = xyz_coords.shape[0]


print(number_of_points)


file = open("interpolated_fields.csv")
reader = csv.reader(file)
data = list(reader)

dat = np.array(data[number_of_points:],dtype=np.double)
t = dat[:,0]

time_samples = int(dat.shape[0]/number_of_points)

print(time_samples)

t = []
field_data = []

for i in range(0,int(time_samples)):
    t.append(dat[i*number_of_points,0])
    field_data.append(dat[i*number_of_points:(i+1)*number_of_points,1:])
        
data=field_data[0]

s=data[:,0]
ss=data[:,1]
s_rms=data[:,2]

print("data from interpolated quantities:")
print(np.min(s)) 
print(np.min(ss)) 
print(np.min(s_rms))

print("data from substracting after interpolating:")
print(np.min(ss-s**2)) 
