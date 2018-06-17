def percentile_single_formula(x, p, method=7):
    '''
    Computes the qth percentile of the data. 
    Offers 8 additional methods not available in Numpy.
    
    This method uses a single formula, which is less efficient 
    than the percentile method for continuous methods (4-9).
    Included merely for educational purposes.
    
    Parameters
    ----------
    x : array_like
        Input array or object that can be converted to an array.
    p : float or array_like
        Percentile(s) to compute, must be between 0 and 100 inclusive.
    method : integer
        Optional parameter that specifies one of the nine sampling 
        methods discussed in Hyndman and Fan (1996), must be between
        1 and 9 inclusive. 
    
    Returns:
    ---------------
        Numpy array with the percentile-corresponding values of x.
    
    Methods
    ----------
    Methods 1 to 3 are discontinuous:
    * Method 1: Inverse of empirical distribution function (oldest 
                and most studied method).
    * Method 2: Similar to type 1 but with averaging at discontinuities.
    * Method 3: Nearest even integer in case of a tie.

    Methods 4 to 9 are continuous and equivalent to a linear 
    interpolation between the points (pk,xk) where xk is the kth order 
    statistic. Specific expressions for pk are given below:
    * Method 4: pk=kn. Linear interpolation of the empirical 
                distribution function.
    * Method 5: pk=(k−0.5)/n. Piecewise linear function where the 
                knots are the values midway through the steps of 
                the empirical distribution function 
                (Popular amongst hydrologists, used by Mathematica?).
    * Method 6: pk=k/(n+1), thus pk=E[F(xk)]. Linear interpolation 
                of the expectations for the order statistics for the 
                uniform distribution on [0,1]. 
                (Used by Minitab and SPSS).
    * Method 7: pk=(k−1)/(n−1), thus pk=mode[F(xk)]. Linear 
                interpolation of the modes for the order statistics 
                for the uniform distribution on [0,1]. 
                (Default method of Numpy, R, S, and MS Excell).
    * Method 8: pk=(k−1/3)/(n+1/3), thus pk≈median[F(xk)].Linear 
                interpolation of theapproximate medians for order 
                statistics. The resulting estimates are approximately 
                median-unbiased regardless of the distribution of x 
                (Recommended by Hyndman and Fan (1996)).
    * Method 9: pk=(k−3/8)/(n+1/4), thus pk≈F[E(xk)]if x is normal. 
                The resulting estimates are approximately unbiased for 
                the expected order statistics if x is normally 
                distributed (Used for normal QQ plots).
        
    References
    ----------
    Hyndman, R. J. and Fan, Y. (1996) Sample quantiles in statistical 
        packages, American Statistician 50, 361--365.
    Schoonjans, F., De Bacquer, D., & Schmid, P. (2011). Estimation of 
        population percentiles. Epidemiology (Cambridge, Mass.), 22(5), 
        750.
    '''
    
    # Input preparation
    # -------------------------------------------------------------------------
    if method < 1 or method > 9:
        raise ValueError('Invalid method. Must be an integer beween 1 and 9')
    if np.isscalar(p):
        p = [p]
    
    x = np.asarray(x)
    x.sort()
    p = np.array(p)/100
    n = x.size  
    
    m = (0, 0, -0.5, 0, 0.5, p, 1-p, (p+1)/3, p/4+3/8)[method-1]
    npm = n*p+m
    index = np.floor(npm).astype(np.int)
    index_frac = npm % 1
           
    # Compute indexes and indexes masks / fractionals
    # -------------------------------------------------------------------------
    if method < 4:       
        im0 = (0, 0.5, 0)[method-1] 
        index_mask = np.ones(p.size)
        if method < 3:
            index_mask[(index_frac == 0)] = im0
        else:
            index_mask[(index_frac == 0) & (index%2 == 0)] = im0 
    else:
        index_mask = index_frac

    # Adjust indexes and return output
    # -------------------------------------------------------------------------       
    index_ = index.copy()
    index[index_ <= 0] = 1
    index[index_  > n] = n
    index_[index_ < 0] = 0
    index_[index_ >= n] = n-1
    
    i = x[index - 1]
    j = x[index_]
    
    return (1-index_mask)*i + index_mask*j
