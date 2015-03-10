

import numpy as np
from scipy import special, optimize
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from predia import *

from cstruct import *
    
ctrl = cstruct(); # control structure used in the subfunction of predia            
ctrl.debug = 1
n_mc = 50000

input = np.random.randn(10,n_mc);

for i in xrange(1,9) :
  input[i,:] = input[i-1,:] + np.random.randn(1,n_mc) * 0.5;


input[9,:] = np.random.randn(1,n_mc);

meas_err_std = np.zeros((10,1))
meas_err_std[0:9,0]  = 0.8;
meas_err_std[4:9,0]  = 0.3;
pred_err_std         = 1

output = (input[0,:]+1)**2 - (.5 * input[6,:] + np.random.randn(1,n_mc) * 0.25)**3 + np.random.randn(1,n_mc) * 0.2;

color = 'blue'
scale = 10

plt.figure(1)
plt.hold(True)
n_plot = 50000
#ax1.set_color_cycle(['c', 'm', 'y', 'k'])
for i in xrange(0,9):

    plt.scatter(input[9-i,0:n_plot] , output[0,0:n_plot], c=[(0.1*i),0,1-(0.1*i)], s=scale, label=color, alpha=.3, edgecolors='none')
    print 1
plt.grid(True)
plt.show(2)

#fig = plt.figure(10)
#ax = fig.add_subplot(111, projection='3d')
#n = 10000
#col = output[0:n]
#ax.scatter(input[0,0:n], input[6,0:n], output[0,0:n], c=output[0,0:n], marker='o', alpha=.3)
# ax.scatter(input[0,0:n].T, input[6,0:n].T, output[0,0:n].T, c='red', marker=0, alpha=.3)
#plt.show()

prior_var = np.var(output);
print prior_var
E_cond_var = np.zeros((1,10))
for i in xrange(0,10):
    cond_var = expected_cond_var(ctrl,input[i,:],input[i,0:2000], meas_err_std[i,0],output,pred_err_std)
    print cond_var.shape
    E_cond_var[0,i] = np.mean(cond_var)
    print ['E_var:', E_cond_var]

plt.figure(2)
plt.plot(prior_var- E_cond_var)
plt.grid(True)
plt.show()
