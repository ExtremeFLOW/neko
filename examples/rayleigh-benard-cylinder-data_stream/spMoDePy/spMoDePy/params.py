# Import required modules
import json




class params_c():
    def __init__(self,filename):

        # Open input file
        f = open (filename, "r") 
        params_file = json.loads(f.read())
        params_file_str = json.dumps(params_file, indent=4)
        f.close()

        # Assign variables to object
        if "decomposition" in params_file["case"]:
            self.decomposition = params_file["case"]["decomposition"]
        else:
            print("Specify which decomposition: POD, SPOD or DMD")

        if self.decomposition == "SPOD":
            self.window_size = params_file["case"]["window_size"]
            self.frequency_locations = params_file["case"]["frequency_locations"]
            self.number_of_decompositions = len(self.frequency_locations)

        self.execution_type = params_file["case"]["execution"]["type"]
        
        if self.execution_type == "fromAdios":
            self.adiosPath = params_file["case"]["execution"]["adiosPath"]

        self.number_of_snapshots = params_file["case"]["number_of_snapshots"]
        
        self.batch_size = params_file["case"]["batch_size"]

        self.keep_modes = params_file["case"]["keep_modes"]
        
        self.fields = params_file["case"]["fields"]
        
        self.casename = params_file["case"]["casename"]
        
        self.dataPath = params_file["case"]["dataPath"]
        
        self.update_type = params_file["case"]["update_type"]

        if self.update_type == "global":
            self.ifgbl_update = True
        else:
            self.ifgbl_update = False


        self.autoupdate = params_file["case"]["autoupdate"]["value"]

        if self.autoupdate == "True":
            self.ifautoupdate = True
        else:
            self.ifautoupdate = False
        
        self.restart = params_file["case"]["restart"]["value"]

        if self.restart == "True":
            self.ifrestart = True
            self.restart_file_path = params_file["case"]["restart"]["filePath"]
        else:
            self.ifrestart = False
        
        self.maximun_number_of_modes = params_file["case"]["autoupdate"]["maximun_number_of_modes"]
        
        self.minimun_orthogonality_ratio = params_file["case"]["autoupdate"]["minimun_orthogonality_ratio"]
        
        self.num_fields = len(self.fields)

        self.file = params_file
        self.string = params_file_str


    def write_string(self):
        print("The input file contains:")
        print(self.string)

