from __future__ import division

import numpy as np
import numpy.ma as ma
cimport numpy as np
from libc.stdint cimport int32_t
cimport cython
from libc.stdio cimport printf

@cython.embedsignature(True)
@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def get_fit(object theta, object height):
    """
       fits 3 lines to a vertical theta profile
       
       parameters
       ----------

        theta, height: numpy 1d array of floating point numbers
        
       returns:
       --------

       fitvals: numpy 1d array of floating point numbers   
       RSS: numpy 2d array of floating point numbers
       j, k: integers
       
       example
       -------
             
    """   
    theta=np.ascontiguousarray(theta)
    theta=theta.astype(np.float64)
    cdef double* thetaPtr= <double*> np.PyArray_DATA(theta)

    height=np.ascontiguousarray(height)
    height=height.astype(np.float64)
    cdef double* heightPtr= <double*> np.PyArray_DATA(height)
    
    
    cdef np.float64_t[:] fitvals=np.empty([theta.size],dtype=np.float64)
    cdef np.float64_t[:,:] RSS=np.empty([290, 290],dtype=np.float64)
    
    cdef int i, j, k, J, K
    #cdef double num_b_11, num_b_12, num_b_13, dem_b_11, dem_b_12
    #cdef double num_b_21, num_b_22, dem_b_21, dem_b_22, num_a_21, num_a_22
    #cdef double num_b_31, num_b_32, dem_b_31, dem_b_32, num_a_31, num_a_32
    #cdef double b_1, a_1, b_2, a_2, b_3, a_3, num_b, dem_b, num_b2, dem_b2, num_b_3, dem_b_3

#def get_fit(theta, height):
#     """
#        Fitting the local theta profile with three lines
#     
#     """
     
     
     #RSS = np.empty((290, 290))+ np.nan
     #print RSS[0,0]
     for j in range(290):
          if j > 2:
               for k in range(290):
                    if k>j+1 and k<289:
                         b_1 = (np.sum(np.multiply(height[:j], theta[:j])) - 1/j*np.sum(height[:j])*np.sum(theta[:j]))/(np.sum(height[:j]**2) - 1/j*np.sum(height[:j])**2)
                         a_1 = np.sum(np.multiply(height[:j], theta[:j]))/np.sum(height[:j]) - b_1*np.sum(height[:j]**2)/np.sum(height[:j])
                         
                         b_2 = (np.sum(theta[j:k]) - (k-j)*(a_1+b_1*height[j]))/(np.sum(height[j:k]) - (k-j)*height[j])
                         
                         a_2 = np.sum(np.multiply(height[j:k], theta[j:k]))/np.sum(height[j:k]) - b_2*np.sum(height[j:k]**2)/np.sum(height[j:k])

                         b_3 = (np.sum(theta[k:290]) - (290-k)*(a_2+b_2*height[k]))/(np.sum(height[k:290]) - (290-k)*height[k])
                         a_3 = np.sum(np.multiply(height[k:290], theta[k:290]))/np.sum(height[k:290]) - b_3*np.sum(height[k:290]**2)/np.sum(height[k:290])
                         
                         RSS[j, k] = np.sum(np.add(theta[2:j], -(a_1+ b_1*height[2:j]))**2) + np.sum(np.add(theta[j:k], -(a_2+ b_2*height[j:k]))**2) + np.sum(np.add(theta[k:290], -(a_3+ b_3*height[k:290]))**2) 

                         if j==3 and k==5:
                              RSS_min = RSS[j, k]                           
                         
                         if RSS[j, k]<RSS_min: 
                              RSS_min = RSS[j, k]
                              J, K = j, k
     
     return RSS, J, K                                                              
                                                                               


    
    
