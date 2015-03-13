from cstruct import *
from predia import *
import numpy as np
import matplotlib.pyplot as plt

prior_data = np.random.random(5000)*20
prior_data.shape = (1,np.size(prior_data))
obs_data = np.array([5,3])
obs_data.shape = (1,2)
obs_err_std = np.ones((1,1)) *2

ctrl = cstruct()
print prior_data
out_ = predia_weight_matrix(ctrl, prior_data,obs_data, obs_err_std)


print prior_data
print out_['weights']

print np.var(prior_data)
print ['post variance: ', weighted_cond_var(ctrl, out_, prior_data)]
print prior_data.shape
print np.shape(out_['weights'])
plt.figure(1)
plt.scatter(prior_data,out_['weights'].T)
plt.grid(True)
plt.show()
