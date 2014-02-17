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
    cdef double num_b_11, num_b_12, num_b_13, dem_b_11, dem_b_12
    cdef double num_b_21, num_b_22, dem_b_21, dem_b_22, num_a_21, num_a_22
    cdef double num_b_31, num_b_32, dem_b_31, dem_b_32, num_a_31, num_a_32
    cdef double b_1, a_1, b_2, a_2, b_3, a_3, num_b, dem_b, num_b2, dem_b2, num_b_3, dem_b_3
      
    RSS = np.empty([290, 290]) + np.nan
    for j in range(290):
        if j > 2:            
            for k in range(290):                
                if k>j+1 and k<289:                
                    #b_1, a_1, num_b_11, num_b_12, num_b_13, dem_b_11, dem_b_12 = 0, 0, 0, 0, 0, 0, 0  
                    #for i in range(j):
                    #    num_b_11 = num_b_11 + heightPtr[i]*thetaPtr[i]
                    #    num_b_12 = num_b_12 + heightPtr[i]
                    #    num_b_13 = num_b_13 + thetaPtr[i]
                    #    dem_b_11 = dem_b_11 + heightPtr[i]**2
                    #    dem_b_12 = dem_b_12 + heightPtr[i]

                    #num_b = num_b_11 - 1/j*num_b_12*num_b_13
                    #dem_b = dem_b_11  - 1/j*dem_b_12**2
                    #b_1 = num_b/dem_b
                    #a_1 = num_b_11/num_b_12 - b_1*dem_b_11/num_b_12

                    b_1 = (np.sum(np.multiply(height[:j], theta[:j])) - 1/j*np.sum(height[:j])*np.sum(theta[:j]))/(np.sum(height[:j]**2) - 1/j*np.sum(height[:j])**2)
                    a_1 = np.sum(np.multiply(height[:j], theta[:j]))/np.sum(height[:j]) - b_1*np.sum(height[:j]**2)/np.sum(height[:j])
                         
                    #print 'checking'
                    #num_b_21, num_b_22, dem_b_21, dem_b_22, num_a_21, num_a_22, num_b2, dem_b2, b_2, a_2 = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 
                    #for i in range(k-j):
                    #    num_b_21 = num_b_21 + thetaPtr[j+i]
                    #    num_b_22 = num_b_22 + a_1+b_1*heightPtr[j]
                    #    dem_b_21 = dem_b_21 + heightPtr[j+i]
                    #    dem_b_22 = dem_b_22 + heightPtr[j]
                    #    num_a_21 = num_a_21 + heightPtr[j+i]*thetaPtr[j+i]
                    #    num_a_22 = num_a_22 + heightPtr[j+i]**2

                    #num_b2 = num_b_21 - num_b_22
                    #dem_b2 = dem_b_21 - dem_b_22
                    #b_2 = num_b2/dem_b2
                    #a_2 = num_a_21/dem_b_21 - b_2*num_a_22/dem_b_21
                                   
                    b_2 = (np.sum(theta[j:k]) - (k-j)*(a_1+b_1*height[j]))/(np.sum(height[j:k]) - (k-j)*height[j])                         
                    a_2 = np.sum(np.multiply(height[j:k], theta[j:k]))/np.sum(height[j:k]) - b_2*np.sum(height[j:k]**2)/np.sum(height[j:k])
     
                    #num_b_31, num_b_32, dem_b_31, dem_b_32, num_a_31, num_a_32 = 0, 0, 0, 0, 0, 0
                    #for i in range(298-k):
                    #    num_b_31 = num_b_31 + thetaPtr[k+i]
                    #    num_b_32 = num_b_32 + a_2 + b_2*heightPtr[k]
                    #    dem_b_31 = dem_b_31 + heightPtr[k+i]
                    #    dem_b_32 = dem_b_32 + heightPtr[k]
                    #    num_a_31 = num_a_31 + heightPtr[k+i]*thetaPtr[k+i]
                    #    num_a_32 = num_a_32 + heightPtr[k+i]**2

                    #num_b_3 = num_b_31 - num_b_32
                    #dem_b_3 = dem_b_31 - dem_b_32

                    #b_3 = num_b_3/dem_b_3
                    #a_3 = num_a_31/dem_b_31 - b_3*num_a_32/dem_b_31                                           
                          
                    b_3 = (np.sum(theta[k:290]) - (290-k)*(a_2+b_2*height[k]))/(np.sum(height[k:290]) - (290-k)*height[k])
                    a_3 = np.sum(np.multiply(height[k:290], theta[k:290]))/np.sum(height[k:290]) - b_3*np.sum(height[k:290]**2)/np.sum(height[k:290])
                              
                    #print 'checking'                                                  
                    RSS[j, k] = np.sum(np.add(theta[:j], -(a_1+ b_1*height[:j]))**2) + np.sum(np.add(theta[j:k], -(a_2+ b_2*height[j:k]))**2) + np.sum(np.add(theta[k:290], -(a_3+ b_3*height[k:290]))**2)
                    
                    #RSS_1 = 0
                    #for i in range(j):
                    #    RSS_1 = RSS_1 + (thetaPtr[i] -(a_1 + b_1*heightPtr[i]))**2
                    #print "checking"
                    #print RSS_1, np.sum(np.add(theta[:j], -(a_1+ b_1*height[:j]))**2)
                    #RSS_2 = 0     
                    #for i in range(k-j):
                    #    RSS_2 = RSS_2 + (thetaPtr[j+i] - (a_2 + b_2*heightPtr[j+i]))**2
                    #print RSS_2, np.sum(np.add(theta[j:k], -(a_2+ b_2*height[j:k]))**2)
                    #RSS_3 = 0
                    #for i in range(290-k):
                     #   RSS_3 = RSS_3 + (thetaPtr[k+i] - (a_3 + b_3*heightPtr[k+i]))**2
                    #print RSS_3, np.sum(np.add(theta[k:298], -(a_3+ b_3*height[k:298]))**2)
                    #RSS[j, k] = RSS_1 + RSS_2 + RSS_3

                    #if j==3 and k==5:
                    #    RSS_min = RSS[j, k]                           
                         
                    #if RSS[j, k]<RSS_min: 
                    #    RSS_min = RSS[j, k]
                    #    J, K = j, k
                    #print RSS[j, k], RSS_check

                    #print "checking"

    #RSS = ma.masked_where(np.isnan(RSS), RSS)
    #print RSS.shape
    print np.nanmin(RSS)
    RSSmin = np.nanmin(RSS)
    RSSmin_where = np.where(RSS = RSSmin)
    print RSSmin_where
    J, K = RSSmin_where[0][0], RSSmin_where[0][1]
    print J, K
                     
    
    return RSS, J, K
                              
                                                                                
                                                                              




    
    
