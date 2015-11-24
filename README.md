# Predua source


###### INSTALLATION:

Add the path matlab and matlab/aux to your matlab path. 


###### USAGE:
The file "expect_cond_var_predia.m" is the basic wrapper function that evlaluates the expected conditional variance. 
The data worth of the data is than prior variance - expected conditional variance. 

The file "predia_weight_matrix" is the inner core of the Bayesina likelyhood-based weighting scheme and the cross 
Monte Carlo loop over future observation values.

The file "expect_cond_ent_predia.m" provides a wrapper function for the evaluation expected conditional entropy. 
This is used to evaluate the Mutual Information between observations and prediction quanity. The use fo this function
requires the definition of an probability density estimator.


###### COMMENTS:

Current Matlab version is pretty stable and well tested.

The python version requires more rigorous testing for all use cases
