def expected_cond_var(ctrl,prior_data,obs_data, obs_err_std,pred_data,pred_err_std):

    import numpy as np
    # WRAPPER:
    # Evalulation of the expected conditional variance of the prediction data
    # for prior data given the observation data
    
    if np.size(prior_data.shape) == 1:
        length = np.size(prior_data)
        prior_data.shape = (1,length)
        
    if np.size(obs_data.shape) == 1:
        length = np.size(obs_data)
        obs_data.shape = (1,length)
        
    if np.size(pred_data.shape) == 1:
        length = np.size(pred_data)
        pred_data.shape = (1,length)
    
    ## INIT
    n_data = np.shape(prior_data);
    n_meas = np.shape(obs_data);
    eps    = np.spacing(1)
    if ctrl.isSetTrue('debug'):
        print n_data
        print n_meas
    
    if not(ctrl.isfield('warn_ESS')):
        ctrl.warn_ESS = 50;
    
    if not(n_data[0] == n_meas[0]):
        raise Exception('Dimension of data disagree')
        
    #if not(ctrl.isfield('n_para')
    #    ctrl.n_para = 1; % no paralell computing
    
    split = get_n_splits(ctrl,n_data[1],n_meas[1])

    cond_var = np.zeros(n_meas)
        
    for t in xrange(0,split.n_):
        #print split.part_start[t]
        obs_data_part = obs_data[:,split.part_start[t]:split.part_end[t]];
        
        # calculation of the weight matrix
        dict_out = predia_weight_matrix(ctrl, prior_data,obs_data_part, obs_err_std);
        
        cond_var[split.part_start[t]:split.part_end[t]] = weighted_cond_var(ctrl, dict_out, pred_data);
    
    cond_var = cond_var[~np.isnan(cond_var)];
    E_cond_var = np.mean(cond_var);
    
    if np.min(dict_out['ESS']) < ctrl.warn_ESS:
        n_crit = np.sum(dict_out['ESS'] < ctrl.warn_ESS)
        print n_crit
        print ["shape of n_crit", n_crit]
        if ~ctrl.isSetTrue('no_warning'):
            print "ESS lower than ", np.str(ctrl.warn_ESS), " in ", np.str((np.float(n_crit)/np.float(n_meas[1])*100)), "% of obs. realizations"
        
    return E_cond_var

#########################################################################
# core FUNCITON to calculate the weighting matrxi                       #
#########################################################################
def predia_weight_matrix(ctrl, prior_data,obs_data, obs_err_std):
    
    import numpy as np
    
    if np.rank(prior_data) < 2:
        prior_data = np.array(prior_data)
        prior_data.shape = (1,np.size(prior_data))
        
    if np.rank(obs_data) < 2:
        obs_data = np.array(obs_data)
        obs_data.shape = (1,np.size(obs_data))
    
    if np.rank(obs_err_std) < 2:
        obs_err_std = np.array(obs_err_std)
        obs_err_std.shape = (1,np.size(obs_err_std))
            
    if ctrl.isSetTrue('no_err_marg'):
        marg_factor = 1; # no additional smoothing of the likelood function
    else:
        marg_factor = 2; # additional smoothing of the likelood function by factor of 2
    
    ## INIT
    n_data = np.shape(prior_data);
    n_meas = np.shape(obs_data);
    eps    = np.spacing(1)
    
    if ctrl.isfield('n_mc'):
        n_mc = ctrl.n_mc;
    else:
        n_mc = n_data[1]

    if not(n_data[0] == n_meas[0]):
        raise Exception('Dimension of data disagree')
    else:
        n_dim = n_data[0]

    weights = np.zeros((n_meas[1],n_mc))

########### main calculation ##################
    if ctrl.isSetTrue('debug'):
        print prior_data.shape
        print prior_data.shape
        print [obs_err_std, obs_err_std.shape]
    
    prior_data_norm = np.zeros(np.shape(prior_data))
    obs_data_norm   = np.zeros(np.shape(obs_data))
        
    # normalization with error standart deviation    
    for i in xrange(0,n_dim):
        prior_data_norm[i,:] = prior_data[i,:]  / obs_err_std[i] / np.sqrt(2*marg_factor);
        obs_data_norm [i,:]  = obs_data[i,:]    / obs_err_std[i] / np.sqrt(2*marg_factor);
    
    if ctrl.isSetTrue('no_err_marg'):
        prior_data_norm = prior_data_norm + np.random.randn(np.shapesize(prior_data_norm));
    
    for i_mc in xrange(0,n_mc):
        for i_data in xrange(0,n_dim):
            weights[:,i_mc] = weights[:,i_mc] + ((prior_data_norm[i_data,i_mc]-obs_data_norm[i_data,:])**2);
    
    #print weights
    weights = np.exp(-weights)
    
    #print weights
    ########### post processing  ##################
    
    dict_out = dict()
    
    if ctrl.isfield('bayes_evidence'):
        
        dict_out['bayes_evidence'] = weights
        
    if ctrl.isfield('diag_zero') and ctrl.diag_zero:
    
        weights[[np.arange(0,n_meas[1]), np.arange(0,n_meas[1])]] = 0
            
    # Normalizing of weights    
    sumWeights = np.sum(weights,1)
    sumWeights[sumWeights< 2 * eps] = 1
    sumWeights.shape = (n_meas[1],1)
    weights        = np.tile(1/sumWeights,(1,n_mc))*weights
    
    #avoiding -inf errror when one weight is equal to one
    dict_out['sumSqrWeights'] = np.sum(weights**2,1)
    del_idx = dict_out['sumSqrWeights'] > 0.999;
    n_del = sum(del_idx)
    if n_del > 0:
        dict_out['sumSqrWeights'][del_idx] = np.NaN
        weights[del_idx] = np.NaN
        if ~ctrl.isSetTrue('no_warning'):
            print ["### Warning: ", n_del, " measurement realizations deleted ###"]
    idx_ = np.repeat([True],np.size(del_idx)) 
    idx_[del_idx] = 0
    dict_out['ESS'] = 1 / dict_out['sumSqrWeights']
    #print np.shape(np.tile(1/sumWeights,(1,n_mc)))
    #print np.shape(weights)
    #print ['sum', np.sum(weights,1)]
    #print weights
    dict_out['weights'] = weights
    
    return dict_out

#########################################################################
# FUNCITON to calculate the conditional variance                        #
#########################################################################

def weighted_cond_var(ctrl, dict_weight, targ_data):
    
    import numpy as np
    
    sumSqrWeights = dict_weight['sumSqrWeights']
    sumSqrWeights = sumSqrWeights[~np.isnan(sumSqrWeights)]
    if ctrl.isSetTrue('debug'):
        # targ_data.shape
        print np.shape(dict_weight['weights'])
        print np.shape(dict_weight['sumSqrWeights'])
        print sumSqrWeights.shape
    
    out = np.dot(dict_weight['weights'], targ_data.T**2) - np.dot(dict_weight['weights'], targ_data.T)**2 
    #print['size out', np.shape(out)]
    #print out
    #print sumSqrWeights
    #print ['sum squared weights',(1/(1-sumSqrWeights))]
    out = (1/(1-sumSqrWeights)) * out.T
    print out.shape
    return out
    

#########################################################################
# FUNCITON to calculate how often the predia matrix needs to be devided #
#########################################################################
def get_n_splits(ctrl, n_mc,n_meas):

    class simpl_stuct:
        pass
     
    split = simpl_stuct()
    
    import os, math
    import numpy as np

    mem_bytes = float(os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES'))
    split.n_ = int(math.ceil((8 * n_mc * n_meas *6) / mem_bytes));

    if split.n_ < 1:

        n_tmp = math.floor(n_meas[1]/split.n_)
        split.part_start = np.zeros([split.n_])
        split.part_end   = np.ones([split.n_]) * n_tmp
        split.n_part     = np.ones([split.n_]) * n_tmp
        for t in xrange(1,split.n_-1):
            split.part_start[t] = split.part_end[t-1]+1
            if t == split.n_-1:
                split.part_end[t] = n_meas[1]
            else:
                split.part_end[t] = part_end[t-1] + split

    else:
        split.part_start = np.zeros(1);
        split.part_end   = np.ones(1) * n_meas;
        split.n_part     = np.ones(1) * n_meas;
        
    if ctrl.isSetTrue('debug'):
        print split.part_start
        print split.part_end
        print split.n_part
                
    return split
    

