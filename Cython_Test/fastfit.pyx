import numpy as np
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
    cdef double* dataPtr= <double*> np.PyArray_DATA(theta)

    height=np.ascontiguousarray(height)
    height=height.astype(np.float64)
    cdef double* dataPtr1= <double*> np.PyArray_DATA(height)
    
    
    cdef np.float64_t[:] fitvals=np.empty([theta.size],dtype=np.float64)
    cdef np.float64_t[:,:] RSS=np.empty([298, 298],dtype=np.float64)
    
    cdef int i, j, k, J, K
    cdef double num_b_11, num_b_12, num_b_13, dem_b_11, dem_b_12
    cdef double num_b_21, num_b_22, dem_b_21, dem_b_22, num_a_21, num_a_22
    cdef double num_b_31, num_b_32, dem_b_31, dem_b_32, num_a_31, num_a_32
    cdef double b_1, a_1, b_2, a_2, b_3, a_3, num_b, dem_b, num_b2, dem_b2, num_b_3, dem_b_3
    
    for j in range(298):
        if j > 2:
            print j
            for k in range(298):
                
                if k>j+1 and k<279:                
                    b_1, a_1, num_b_11, num_b_12, num_b_13, dem_b_11, dem_b_12 = 0, 0, 0, 0, 0, 0, 0  
                    for i in range(j):
                        num_b_11 = num_b_11 + height[i]*theta[i]
                        num_b_12 = num_b_12 + height[i]
                        num_b_13 = num_b_13 + theta[i]
                        dem_b_11 = dem_b_11 + height[i]**2
                        dem_b_12 = dem_b_12 + height[i]
                        num_b = num_b_11 - 1/j*num_b_12*num_b_13
                        dem_b = dem_b_11  - 1/j*dem_b_12**2
                        b_1 = num_b/dem_b
                        a_1 = num_b_11/num_b_12 - b_1*dem_b_11/num_b_12

                        b_2, a_2, num_b_21, num_b_22, dem_b_21, dem_b_22, num_a_21, num_a_22 = 0, 0, 0, 0, 0, 0, 0, 0
                        for i in range(k-j):
                            num_b_21 = num_b_21 + theta[j+i]
                            num_b_22 = num_b_22 + a_1+b_1*height[j]
                            dem_b_21 = dem_b_21 + height[j+i]
                            dem_b_22 = dem_b_22 + height[j]
                            num_a_21 = num_a_21 + height[j]*theta[j]
                            num_a_22 = num_a_22 + height[j]**2

                        num_b2 = num_b_21 - num_b_22
                        dem_b2 = dem_b_21 - dem_b_22
                        b_2 = num_b2/dem_b2
                        a_2 = num_a_21/dem_b_22 - b_2*num_a_22/dem_b_21
                         
                        b_22 = (np.sum(theta[j:k]) - (k-j)*(a_1+b_1*height[j]))/(np.sum(height[j:k]) - (k-j)*height[j])                         
                        a_22 = np.sum(np.multiply(height[j:k], theta[j:k]))/np.sum(height[j:k]) - b_2*np.sum(height[j:k]**2)/np.sum(height[j:k])
                         
                         #print np.sum(theta[j:k]), num_b_21  NOTE: these are different, np.sum seems to round to 1 decimal place
                         
                        num_b_31, num_b_32, dem_b_31, dem_b_32, num_a_31, num_a_32 = 0, 0, 0, 0, 0, 0
                        for i in range(298-k):
                            num_b_31 = num_b_31 + theta[k+i]
                            num_b_32 = num_b_32 + a_2 + b_2*height[k]
                            dem_b_31 = dem_b_31 + height[k+i]
                            dem_b_32 = dem_b_32 + height[k]
                            num_a_31 = num_a_31 + height[k+i]*theta[k+i]
                            num_a_32 = num_a_32 + height[k+i]**2

                        num_b_3 = num_b_31 - num_b_32
                        dem_b_3 = dem_b_31 - dem_b_32

                        b_3 = num_b_3/dem_b_3
                        a_3 = num_a_31/dem_b_31 - b_2*num_a_32/dem_b_31
                                                  
                         #print np.sum(theta[k:298]), num_b_31 NOTE: sum seems to round up to 1 dec                                   
                         #print np.sum(np.multiply(height[k:298], theta[k:298])), num_a_31 NOTE: again rounding but this time not decimal places!!      
                                                                      
                        RSS[j, k] = np.sum(np.add(theta[2:j], -(a_1+ b_1*height[2:j]))**2) + np.sum(np.add(theta[j:k], -(a_2+ b_2*height[j:k]))**2) + np.sum(np.add(theta[k:298], -(a_2+ b_2*height[k:298]))**2) 
                        RSS_1 = 0
                        for i in range(j):
                            RSS_1 = RSS_1 + (theta[i] -(a_1 + b_1*height[i]))**2
                        RSS_2 = 0     
                        for i in range(k-j):
                            RSS_2 = RSS_2 + (theta[j+i] - (a_2 + b_2*height[j+i]))**2
                            RSS_3 = 0
                        for i in range(298-k):
                              RSS_3 = RSS_3 + (theta[k+i] - (a_3 + b_3*height[k+i]))**2

                        RSS[j, k] = RSS_1 + RSS_2 + RSS_3

                        if j==3 and k==5:
                            RSS_min = RSS[j, k]
                            
                         
                        if RSS[j, k]<RSS_min: 
                            RSS_min = RSS[j, k]
                            J, K = j, k
                        

    return RSS, J, K
                              
                                                                                
                                                                              




    
    
