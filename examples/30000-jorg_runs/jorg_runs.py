import numpy as np

# Original runs by jorg
Nu = [1e8, 1e9, 1e10,1e11,1e12, 1e13, 1e14, 1e15]
nel = [192000, 192000, 192000, 192000, 537600, 925200, 17145000, 17145000]
lx = [3, 5, 5, 9, 11, 13, 11, 11]


mynel = [500000,500000,500000,500000, 6000000, 6000000, 30000000, 30000000]


for i in range(0,len(Nu)):
    print("Nu = " + repr(Nu[i]))

    number_of_points = nel[i]*(lx[i]+1)**3
    #number_of_points = nel[i]*(lx[i])**3
   

    mylx = (number_of_points/mynel[i])**(1/3)


    print("my nel: "+repr(mynel[i]))
    print("my lx: "+repr(mylx))
    print("my round lx: "+repr(np.ceil(mylx)))
    print("my order: "+repr(np.ceil(mylx)-1))
    print("mynodes: "+repr(np.ceil(mynel[i]/80000)))

    print("--------------------------------------")
