# Import required modules
import numpy as np
from pymech.neksuite import readnek
from pymech.neksuite import writenek
import json

class case_c():
    
    # Contains:
    def __init__(self,
                 params_c,
                 comm):

        self.params = params_c
        self.comm = comm
        self.get_case_info()
 
    def get_case_info(self):
    
        # Associate objects for readability
        params = self.params
        comm   = self.comm

        # Get information from the communicator
        rank     = comm.Get_rank()
        size     = comm.Get_size()
        path     = params.dataPath
        casename = params.casename

        if rank==0:
            buff  = np.zeros((10))
            # Read information about the case (number of elements, etc)
            filename=path+casename+'0.f'+str(0+0).zfill(5)
            data=readnek(filename)
            # simulation parameters
            buff[0]=data.nel #nel
            buff[1]=data.lr1[0] #lx1
            buff[2]=data.lr1[1] #ly1
            buff[3]=data.lr1[2] #lz1
            buff[4]=np.prod(data.lr1) #nxyz
            buff[5]= data.nel*data.lr1[0]*data.lr1[1]*data.lr1[2] #nxyze
            if data.lr1[2]==1:
                buff[6]=2
            else:
                buff[6]=3
        else:
            buff  = np.empty((10))

        comm.Bcast(buff, root=0)

        self.nel=int(buff[0])
        self.lx1=int(buff[1])
        self.ly1=int(buff[2])
        self.lz1=int(buff[3])
        self.nxyz=int(buff[4])
        self.nxyze=int(buff[5])
        self.dim=int(buff[6])
        self.num_fields = int(len(params.fields))
        self.nxyzef = self.nxyze*self.num_fields
        
        # Put the data in json format
        case_data = {
            "Total number of elements": self.nel,
            "lx1": self.lx1,
            "ly1": self.ly1,
            "lz1": self.lz1,
            "nxyz": self.nxyz,
            "nxyze": self.nxyze,
            "number of fields": self.num_fields,
            "nxyzef": self.nxyzef,
            "dimension": self.dim
            }

        case_data_str = json.dumps(case_data, indent=4)
        self.string = case_data_str

        return


    def write_string(self):
        
        print("The case has the following parameters:")
        print(self.string)



