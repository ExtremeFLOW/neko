#*********************************************************************
# Algebraic functions to model ACF
#*********************************************************************
#  Saleh Rezaeiravesh, salehr@kth.se
#---------------------------------------------------------------------
import numpy as np
 

class ACF_FUNCS:
    """
    Algebraic functions to model ACF
    Args:
       `k`: numpy array containing integer lags
       `a,b,c`: model parameters
    Returns:
       function value at `k` given parameters (a,b,c)      
    """
    def __init__(self,funType):
        self.funType=funType

        if self.funType=='fun1':
           self.fun=self.fun1
        elif self.funType=='fun2':
           self.fun=self.fun2
        elif self.funType=='fun3':
           self.fun=self.fun3
        elif self.funType=='fun4':
           self.fun=self.fun4

    @classmethod
    def fun1(self,k,a):
        """
        Note parameters in the exp(.) should not become negative => constraint optimization needed
        """
        return(np.exp(-a*k))   

    @classmethod
    def fun2(self,k,a,b,c):
        """
        Note parameters in the exp(.) should not become negative => constraint optimization needed
        """
        ##Multiscale exponential stochastic process
        return(a*np.exp(-b*k)+(1-a)*np.exp(-c*k))   

    @classmethod
    def fun3(self,k,a,b,c):
        """
        Note parameters in the exp(.) should not become negative => constraint optimization needed
        """
        ##Multiscale Gaussian stochastic process
        return(a*np.exp(-(b*k)**2.)+(1-a)*np.exp(-(c*k)**2.))   

    @classmethod
    def fun4(self,k,a,b,c,d,e):
        """
        Note parameters in the exp(.) should not become negative => constraint optimization needed
        """
        ##Multiscale Gaussian stochastic process
        return(a*np.exp(-b*k)+c*np.exp(-d*k)+(1-a-c)*np.exp(-e*k))   
#    
