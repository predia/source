import numpy as np
#from scipy import special, optimize
import matplotlib.pyplot as plt
from matplotlib.pyplot import ion
from mpl_toolkits.mplot3d import Axes3D
import predia as pr
import cstruct as cs
from time import time as ti
    
ctrl = cs.cstruct(); # control structure used in the subfunction of predia            
ctrl.debug = 0;
ctrl.generate = 0
ctrl.cal_meth = 2
n_mc   = 20000
n_meas = 2000
ion()
if ctrl.isSetTrue('generate'):
    theta = np.random.randn(10,n_mc);

    for i in xrange(1,9) :
        theta[i,:] = theta[i-1,:] + np.random.randn(1,n_mc) * 0.5;


    theta[9,:] = np.random.randn(1,n_mc);
    output = (theta[0,:]+1)**2 - (.5 * theta[6,:] + np.random.randn(1,n_mc) * 0.25)**3 + np.random.randn(1,n_mc) * 0.2;
else:
    dt     = np.dtype([('time', float)])
    theta  = np.float64(np.fromfile("input.data", dtype=dt))
    output = np.float64(np.fromfile("output.data", dtype=dt))
    theta.shape = (n_mc,10)
    output.shape = (n_mc,1)
    theta  = theta.T
    output = output.T
    print "file loaded"
    
meas_err_std = np.zeros((10,1))
meas_err_std[0:10]  = 0.8;
meas_err_std[4:10]  = 0.3;
pred_err_std         = 1



color = 'blue'
scale = 10

plt.figure(1)
plt.hold(True)
n_plot = 1000
#ax1.set_color_cycle(['c', 'm', 'y', 'k'])
for i in xrange(0,9):
    
    plt.scatter(theta[9-i,0:n_plot] , output[0,0:n_plot], c=[(0.1*i),0,1-(0.1*i)], s=scale, label=color, alpha=.3, edgecolors='none')
    #print 1
    
plt.grid(True)
plt.draw()

#fig = plt.figure(10)
#ax = fig.add_subplot(111, projection='3d')
#n = 10000
#col = output[0:n]
#ax.scatter(theta[0,0:n], theta[6,0:n], output[0,0:n], c=output[0,0:n], marker='o', alpha=.3)
# ax.scatter(theta[0,0:n].T, theta[6,0:n].T, output[0,0:n].T, c='red', marker=0, alpha=.3)
#plt.show()

prior_var = np.var(output);
print prior_var
E_cond_var = np.zeros((1,10))
for i in xrange(0,10):
    tt1 = ti()
    E_cond_var[0,i] = pr.expected_cond_var(ctrl,theta[i,:],theta[i,0:n_meas], meas_err_std[i,0],output,pred_err_std)  
    print 'time:', ti() - tt1, 's'
    print ["E_var:", E_cond_var]

plt.figure(2)
x = np.arange(1,11)
plt.plot(x, prior_var- E_cond_var.T, linestyle='-', color= 'b')
plt.grid(True)
plt.show()

E_cond_var2 = np.zeros((1,10))
for i in xrange(0,10):
    tt1 = ti()
    comb_input = theta[[0, i],:];
    E_cond_var2[0,i] = pr.expected_cond_var(ctrl,comb_input,comb_input[:,0:n_meas], meas_err_std[[0,i],0:1],output,pred_err_std)  
    print 'time:', ti() - tt1, 's'

plt.hold(True)
x = np.arange(1,11)
plt.plot(x, prior_var- E_cond_var2.T, linestyle='-', color= 'r')
plt.grid(True)
plt.show()    


dict_out = pr.predia_weight_matrix(ctrl, comb_input[:,0:10],comb_input[:,0:1], meas_err_std[[0, i],0:1]);
print dict_out['weights']
print pr.weighted_cond_var(ctrl, dict_out, output[0,0:10]);

dict_out = pr.predia_weight_matrix(ctrl, comb_input[:,0:20000],comb_input[:,0:2000], meas_err_std[[0, i],0:1]);
print dict_out['weights']
print np.mean(pr.weighted_cond_var(ctrl, dict_out, output[0,0:20000]))