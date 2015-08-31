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
    
    (n_split,split_start,split_end,s_part) = get_n_splits(ctrl,n_data[1],n_meas[1])
    
    if ctrl.isSetTrue('debug') and n_split > 1:
        print(['splitting in ', n_split, ' pieces'])
        print ['slice start', split_start]
        print ['slice end  ', split_end]
        print ['slice size ', s_part]

    cond_var = np.zeros((1,n_meas[1]))
    
    for t in xrange(0,n_split):
        #print split_start[t]
        obs_data_part = obs_data[:,split_start[t]:split_end[t]];
        # calculation of the weight matrix
        dict_out = predia_weight_matrix(ctrl, prior_data,obs_data_part, obs_err_std);
        cond_var[0,split_start[t]:split_end[t]] = weighted_cond_var(ctrl, dict_out, pred_data);
    
    cond_var = cond_var[~np.isnan(cond_var)];
    E_cond_var = np.mean(cond_var);
    
    if np.min(dict_out['ESS']) < ctrl.warn_ESS:
        n_crit = np.sum(dict_out['ESS'] < ctrl.warn_ESS)
        if ~ctrl.isSetTrue('no_warning'):
            print "ESS lower than ", np.str(ctrl.warn_ESS), " in ", np.str((np.float(n_crit)/np.float(n_meas[1])*100)), "% of obs. realizations"
        
    return E_cond_var

#########################################################################
# core FUNCITON to calculate the weighting matrxi                       #
#########################################################################
def predia_weight_matrix(ctrl, prior_data,obs_data, obs_err_std, aux_data):
    
    import numpy as np
    import numexpr as ne
    
    if np.ndim(prior_data) < 2:
        prior_data = np.array(prior_data)
        prior_data.shape = (1,np.size(prior_data))
        
    if np.ndim(obs_data) < 2:
        obs_data = np.array(obs_data)
        obs_data.shape = (1,np.size(obs_data))
    
    if np.ndim(obs_err_std) < 2:
        obs_err_std = np.array(obs_err_std)
        obs_err_std.shape = (1,np.size(obs_err_std))
            
    if ctrl.isSetTrue('no_err_marg'):
        marg_factor = 1; # no additional smoothing of the likelood function
    else:
        marg_factor = 2; # additional smoothing of the likelood function by factor of 2
    
    if not(ctrl.isfield('cal_meth')):
        ctrl.cal_meth = 1
        print "method set to one"
    
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
    
    # normalization with error standart deviation    
    for i in xrange(0,n_dim):
        prior_data[i,:] /= (obs_err_std[0,i] * np.sqrt(2*marg_factor));
        #print obs_err_std[0,i]
        obs_data [i,:]  /= (obs_err_std[0,i] / np.sqrt(2*marg_factor));
    
    if ctrl.isSetTrue('no_err_marg'):
        prior_data += np.random.randn(np.shapesize(prior_data));
    
    if ctrl.cal_meth == 1:
        for i_mc in xrange(0,n_mc):
            for i_data in xrange(0,n_dim):
                weights[:,i_mc] += ((prior_data[i_data,i_mc]-obs_data[i_data,:])**2)
                
    if ctrl.cal_meth == 2:
        for i_data in xrange(0,n_dim):
            # ALT 1 (fastest)
            aa = np.tile(prior_data[i_data,0:n_mc],(n_meas[1],1))
            bb = np.tile(obs_data[i_data:i_data+1,:].T,(1,n_mc))
            weights += ne.evaluate("(aa - bb)**2")   
            
            # ALT 2
            #weights += (np.tile(prior_data_norm[i_data,0:n_mc],(n_meas[1],1)) - np.tile(obs_data_norm[i_data:i_data+1,:].T,(1,n_mc)))**2
            
    if ctrl.cal_meth == 3:
        for i_mc in xrange(0,n_mc):
            aa = np.subtract(prior_data[:,i_mc],obs_data.T)
            weights[:,i_mc] = ne.evaluate("sum(aa**2,axis=1)")
        
    #print weights
    
    # ALT 1
    #weights = np.exp(-weights)
    
    # ALT 2 (fastest)
    weights = ne.evaluate("exp(-weights)")
    
    # incpoporate prior weights
    if aux_data.has_key('prior_weights'):
        
        # ALT 1
        #weights = weights * np.tile(aux_data['prior_weights'],[n_meas[1],1])
        
        # ALT 2 (fastest)
        for i_meas in xrange(0,n_meas[1]+1):
            weights[i_meas-1:i_meas,:] *= aux_data['prior_weights']
        
        # ALT 3
        #a = np.tile(aux_data['prior_weights'],[n_meas[1],1])
        #weights = ne.evaluate("weights * a")
        
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
    #print ['shape cond_var: ', out.shape]
    return out
    

#########################################################################
# FUNCITON to calculate how often the predia matrix needs to be devided #
#########################################################################
def get_n_splits(ctrl, n_mc,n_meas):
  
    import os, math
    import numpy as np

    mem_bytes = float(os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES'))
    n_split = int(math.ceil((8 * n_mc * n_meas *6) / mem_bytes));

    if n_split > 1:

        n_tmp = math.floor(n_meas/n_split)
        split_start = np.zeros([n_split])
        split_end   = np.ones([n_split]) * n_tmp
        s_part     = np.ones([n_split]) * n_tmp
        for t in xrange(1,n_split):
            split_start[t] = split_end[t-1]
            if t == n_split-1:
                split_end[t] = n_meas
            else:
                split_end[t] = split_end[t-1] + split

    else:
        split_start = np.zeros(1);
        split_end   = np.ones(1) * n_meas;
        s_part     = np.ones(1) * n_meas;
        
    if ctrl.isSetTrue('debug'):
        print split_start
        print split_end
        print s_part
                
    return (n_split,split_start,split_end,s_part)
    

