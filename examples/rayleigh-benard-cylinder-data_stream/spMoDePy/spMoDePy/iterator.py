import sys
if 'adios2' in sys.modules:
    import adios2

class iterator_c():
    def __init__(self,
                params_c, case_c, logger, comm):

        self.params = params_c
        self.case = case_c
        self.log = logger

        self.iteration = 0
        self.first_iteration = True
        
        self.status = True

        self.maxiterations = params_c.number_of_snapshots

        self.run_update_loop = True
        self.load_buffer1 = True
        self.load_buffer2 = False
        self.update_from_buffer1 = True
        self.update_from_buffer2 = False

    def mark_iteration(self):
        self.iteration += 1
        self.first_iteration = False
        self.run_update_loop = True

    def check_status(self,io_controller):
        if 'adios2' in sys.modules:
            if io_controller.stepStatus == io_controller.okstep:
                self.status = True
            elif io_controller.stepStatus == io_controller.endstream:
                self.status = False

        else:
            if self.iteration == self.maxiterations:
                self.status = False





            


        



