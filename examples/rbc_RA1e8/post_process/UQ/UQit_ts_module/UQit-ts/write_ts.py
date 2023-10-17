#*************************************************************
#     write-related methods for timeAvgUQtool
#*************************************************************
# Saleh Rezaeiravesh, salehr@kth.se
#-------------------------------------------------------------
#

def pw(arg):
    """
    print arg on screen and write the same thing in run.log
    """
    print(arg)
    #fLog=open('tUQ.log','w')
    #fLog.write(arg)
    
def printRepeated(string_to_expand, length):
    """
    Repeat a string 'length' times
    """
    return (string_to_expand * (int(length/len(string_to_expand))+1))[:length]
