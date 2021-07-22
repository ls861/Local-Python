""" linmix -- A hierarchical Bayesian approach to linear regression with error in both X and Y.
"""

from __future__ import print_function
from scipy.stats import norm,multivariate_normal
import numpy as np
import sys
import pickle

def task_manager(conn):
    print('task_manager')
    chain = None
    while True:
        message = conn.recv()
        if message['task'] == 'init':
#            print('tm_init')
            chain = Chain(**message['init_args'])
            chain.initial_guess()
        elif message['task'] == 'init_chain':
#            print('tm_init_chain')
            chain.initialize_chain(message['minIter'],message['nBurn'])
        elif message['task'] == 'step':
#            print('tm_step')
            chain.step(message['niter'])
        elif message['task'] == 'extend':
#            print('tm_extend')
            chain.extend(message['niter'])
        elif message['task'] == 'fetch':
#            print('tm_fetch')
            conn.send(chain.__dict__[message['key']])
        elif message['task'] == 'kill':
#            print('tm_kill')
            break
        else:
            print('tm_invalidtask')
            raise ValueError("Invalid task")

class Chain(object):
    def __init__(self, xArr, yArr, xsigArr, ysigArr, xycovArr, delta, nGaussXi, nChains, nGaussBeagle, piBeagle, proposalscale_xi, proposalscale_eta, proposalscale_alpha, proposalscale_beta, proposalscale_sig0, proposalscale_k, rng=None):
        print('__init__ (chain)')
        self.xArr = np.array(xArr, dtype=float)
        self.yArr = np.array(yArr, dtype=float)
        self.xsigArr = np.array(xsigArr, dtype=float)
        self.ysigArr = np.array(ysigArr, dtype=float)
        self.xycovArr = np.array(xycovArr, dtype=float)
        
        # take the first column for all of the above
        self.x = np.array(xArr[:,0], dtype=float)
        self.y = np.array(yArr[:,0], dtype=float)
        self.xsig = np.array(xsigArr[:,0], dtype=float)
        self.ysig = np.array(ysigArr[:,0], dtype=float)
        self.xycov = np.array(xycovArr[:,0], dtype=float)
        
        # lester, these are basically all true anyway..., literally array of TRUE
        self.wxerr = (self.xsig != 0.0)
        self.wyerr = (self.ysig != 0.0)

        self.werrs = werrs = self.wxerr & self.wyerr

#        if len(self.x) == len(self.x[werrs]):
#            print('werrs ALL true')
#        else:
#            print('werrs ERROR')

        self.xvar = self.xsig**2
        self.yvar = self.ysig**2
        
        self.xycorr = np.zeros_like(self.xycov)
        self.xycorr[werrs] = self.xycov[werrs] / (self.xsig[werrs] * self.ysig[werrs])
            
        if delta is None: # for censorship, array of 1s means keep all objects
            self.delta = np.ones((len(self.x)), dtype=bool)
        else:
            self.delta = np.array(delta, dtype=bool)
            
        self.N = len(self.x) # number of objects
        self.nGaussXi = nGaussXi # number of gaussians modelling xi
        self.nChains = nChains # parallel chains
        self.nGaussBeagle = nGaussBeagle
        self.piBeagle = np.array(piBeagle, dtype=float) # weighting of nGaussBeagle
        
        #parameters for the MCMC chain on update_xi 
        self.accept_xi = 0 # these get updated later on
        self.reject_xi = 0 # these get updated later on
        self.proposalscale_xi = proposalscale_xi
        self.proposalscale_eta = proposalscale_eta
        self.proposalscale_alpha = proposalscale_alpha
        self.proposalscale_beta = proposalscale_beta
        self.proposalscale_sig0 = proposalscale_sig0
        self.proposalscale_k = proposalscale_k

        if rng is None:
            self.rng = np.random.RandomState()
        else:
            self.rng = rng

        self.initialized = False # still needs to make initial guesses as below, then set to True

    def initial_guess(self): # Step 1
        print('initial_guess')
        N = self.N
        nGaussXi = self.nGaussXi
        nGaussBeagle = self.nGaussBeagle # 3
        piBeagle = self.piBeagle # the relative amplitudes of the Beagle posterior GMMs are fixed

        if nGaussBeagle == 1: # number of gaussians modelling BEAGLE posterior
            self.GBeagle = np.ones(N, dtype=int)
        else:
            self.GBeagle = np.zeros((N, nGaussBeagle), dtype=int)
            
            # just start by assigning GBeagle to the Gaussian with highest amplitude
            for i in range(N): # for every object
                maxind = np.argmax(piBeagle[i])
                self.GBeagle[i,maxind] = 1
                self.x[i] = self.xArr[i,maxind]
                self.y[i] = self.yArr[i,maxind]
                self.xsig[i] = self.xsigArr[i,maxind]
                self.ysig[i] = self.ysigArr[i,maxind]
                self.xycov[i] = self.xycovArr[i,maxind]
                
                self.xvar[i] = self.xsig[i]**2
                self.yvar[i] = self.ysig[i]**2
                self.xycorr[i] = self.xycov[i] / (self.xsig[i] * self.ysig[i])

        x = self.x
        xsig = self.xsig
        xvar = self.xvar
        
        y = self.y
        ysig = self.ysig

        self.xi_min = 8.5
        self.xi_max = 10.0

        # take a random sample of xs and ys, calculate alpha, beta and sigma
        np.random.seed()
        x_sampled = np.random.normal(x, xsig)
        y_sampled = np.random.normal(y, ysig)
        fit_sample = np.polyfit(x_sampled, y_sampled, 1)
        residual_sample = y_sampled - (fit_sample[0]*x_sampled + fit_sample[1])
        
        self.alpha = fit_sample[1]
        self.beta = fit_sample[0]
        self.sig0 = np.std(residual_sample)
        self.k = max(0.0, np.random.normal(1.0, 0.2))
        print(self.alpha, self.beta, self.sig0, self.k)

        self.mu0 = np.median(x)
        self.wsqr = np.var(x, ddof=1) - np.median(xvar)
        self.wsqr = np.max([self.wsqr, 0.01*np.var(x, ddof=1)])

        # Now get the values for the mixture parameters, first do prior params
        self.mu0min = min(x)
        self.mu0max = max(x)

        mu0g = np.nan
        while not (mu0g > self.mu0min) & (mu0g < self.mu0max):
            mu0g = self.mu0 + (self.rng.normal(scale=np.sqrt(np.var(x, ddof=1) / N)) /
                               np.sqrt(self.nChains/self.rng.chisquare(self.nChains)))
        self.mu0 = mu0g

        # wsqr is the global scale
        self.wsqr *= 0.5 * N / self.rng.chisquare(0.5 * N)
        self.usqrmax = 1.5 * np.var(x, ddof=1)
        self.usqr = 0.5 * np.var(x, ddof=1)
        self.tausqr = 0.5 * self.wsqr * self.nChains / self.rng.chisquare(self.nChains, size=nGaussXi)
        self.mu = self.mu0 + self.rng.normal(scale=np.sqrt(self.wsqr), size=nGaussXi)

        # get initial group proportions and group labels
        pig = np.zeros(self.nGaussXi, dtype=float) # [0. 0. 0.]
        if nGaussXi == 1:
            self.G = np.ones(N, dtype=int)
            self.pi = np.array([1], dtype=float)
        else:
            self.G = np.zeros((N, nGaussXi), dtype=int)
            for i in range(N):
                minind = np.argmin(abs(x[i] - self.mu)) # index of the nearest mean of the xi gaussians
                pig[minind] += 1 # becomes something like [98. 34. 51.], number of items in each column
                self.G[i, minind] = 1
            self.pi = self.rng.dirichlet(pig+1) # eg [0.29385138 0.34768844 0.35846018]
        self.eta = y.copy()
        self.y_ul = y.copy()
        self.xi = x.copy()

        self.cens = np.nonzero(np.logical_not(self.delta))[0]
        self.initialized = True

    # =============================================================================
    # START OF STEP
    # =============================================================================

    def update_cens_y(self):  # Step 2
#        print('update_cens_y')
        todo = self.cens[:]
        while len(todo) > 0:
            self.y[todo] = self.rng.normal(loc=self.eta[todo],
                                           scale=np.sqrt(self.yvar[todo]),
                                           size=len(todo))
            todo = np.nonzero(np.logical_not(self.delta) & (self.y > self.y_ul))[0]

    def update_xi(self):  # Step 3
        print('update_xi')
        wxerr = self.wxerr                                                              # boolean array where xsig != 0

        # P(xi|psi) 
        xi_curr = self.xi                                                               # current xi values
        sigsqr_xi_curr = self.calc_sigsqr(xi_curr)                                      # array of sigsqr values, one per xi (mass dependent)
        xi_prop = self.rng.normal(size=len(xi_curr),scale=self.proposalscale_xi)+xi_curr    # new proposal for each xi
        sigsqr_xi_prop = self.calc_sigsqr(xi_prop)                                      # array of sigsqr values, one per new xi (mass dependent)
        muArr = []                                                                      # mu are means of xi gaussians
        tausqrArr = []                                                                  # tausqr are sigsqr of xi gaussians
        for i in range(len(self.xi)):
            tempIdx = np.where(self.G[i] == 1)[0]                                       # [0], [1] or [2] decides which xi gaussian
            muArr.append(self.mu[tempIdx][0])                                           # chooses mu from [mu1, mu2, mu3]
            tausqrArr.append(self.tausqr[tempIdx][0])                                   # chooses tau from [ta1, ta2, ta3]

        log_p_xi_psi_curr = norm.logpdf(xi_curr,loc=muArr,scale=np.sqrt(tausqrArr))     # gives log( gaussian value at xi_curr ) per object
        log_p_xi_psi_prop = norm.logpdf(xi_prop,loc=muArr,scale=np.sqrt(tausqrArr))     # gives log( gaussian value at xi_prop ) per object
               
        # P(xi|eta,theta) - this is only really true for the proportionality, I'm not calculating the correct normalization
        log_p_xi_eta_theta_curr = norm.logpdf(self.eta, scale=np.sqrt(sigsqr_xi_curr), loc=self.alpha+self.beta*xi_curr)
        log_p_xi_eta_theta_prop = norm.logpdf(self.eta, scale=np.sqrt(sigsqr_xi_prop), loc=self.alpha+self.beta*xi_prop) 

        '''
        # lester test, comparing P(xi|eta,theta) to analytical form
        import matplotlib.pyplot as plt
        log_p_xi_eta_theta_prop_lester = norm.logpdf(self.eta[0],scale=np.sqrt(sigsqr_xi_prop), loc=self.alpha+self.beta*xi_prop)
        plt.scatter(xi_prop, np.e**log_p_xi_eta_theta_prop_lester, marker='x')        
        x_lester = np.linspace(7, 11, 1000)
        sigma_lester = ( self.sig0 * ( ((1.0-self.k)*(x_lester-self.xi_max) / (self.xi_max - self.xi_min)) + 1) )
        fXgY = abs((1/(sigma_lester*((2*np.pi)**0.5))) * np.exp(-0.5 * (( (self.eta[0] - (self.alpha + (self.beta*x_lester))) / sigma_lester ) ** 2) ))
        plt.plot(x_lester, fXgY)
        plt.show()
        '''

        for i in range(len(self.xi)):
            if wxerr[i] == True: # always: boolean array where xsig != 0

#                P(xi|eta,x,y)
                mean, cov = self.get_2d_mean_cov(i)

                log_p_xi_eta_x_y_curr = multivariate_normal.logpdf([xi_curr[i],self.eta[i]],mean,cov)
                log_p_xi_eta_x_y_prop = multivariate_normal.logpdf([xi_prop[i],self.eta[i]],mean,cov)
                
                log_target_curr = log_p_xi_psi_curr[i] + log_p_xi_eta_theta_curr[i] + log_p_xi_eta_x_y_curr
                log_target_prop = log_p_xi_psi_prop[i] + log_p_xi_eta_theta_prop[i] + log_p_xi_eta_x_y_prop

                try:
                  acceptanceProb = min(1,np.exp(log_target_prop-log_target_curr)) # 1 if prop more likely than curr
                except FloatingPointError as exception:
                    print('FLOATING POINT ERROR', exception)
                    if log_target_prop < log_target_curr:
                        acceptanceProb = 0.
                        print('LESTER001')
                    else:
                        acceptanceProb = 1.
                        print('LESTER004')
                
                u = self.rng.uniform() # random between 0 and 1
                if u <= acceptanceProb: # accept proposal
                    self.xi[i] = xi_prop[i]
                    
                    if i == 0:
                        self.accept_xi = self.accept_xi + 1
                else:
                    if i == 0:
                        self.reject_xi = self.reject_xi + 1

                test = self.accept_xi+self.reject_xi
                if self.ichain >= self.nBurn and test > 0 and test % self.testIter == 0:
                    if np.float(self.accept_xi)/np.float(self.accept_xi+self.reject_xi) > 0.5:
                        self.proposalscale_xi = self.proposalscale_xi*1.1
                    elif np.float(self.accept_xi)/np.float(self.accept_xi+self.reject_xi) < 0.4:
                        self.proposalscale_xi = self.proposalscale_xi*0.9

                    self.accept_xi = 0
                    self.reject_xi = 0

            else:
                print('wxerr[i] == FALSE')
    
    def update_eta(self):  # Step 4
        print('update_eta')
        wxerr = self.wxerr
        wyerr = self.wyerr

        etaxyvar = self.yvar * (1.0 - self.xycorr**2)
        etaxy = self.y.copy()
        etaxy[wxerr] += (self.xycov / self.xvar * (self.xi - self.x))[wxerr]

        # Eqn (68)
        sigsqr = self.calc_sigsqr(self.xi)
        sigma_etahat_i_sqr = 1.0/(1.0/etaxyvar + 1.0/sigsqr)
        # Eqn (67)
        etahat_i = (sigma_etahat_i_sqr * (etaxy / etaxyvar
                    + (self.alpha + self.beta * self.xi) / sigsqr))
        # Eqn (66)
        self.eta[wyerr] = self.rng.normal(loc=etahat_i[wyerr],
                                          scale=np.sqrt(sigma_etahat_i_sqr[wyerr]))

    def update_G(self):  # Step 5
#        print('update_G')
        # Eqn (74)
        piNp = self.pi * (1.0/np.sqrt(2.0*np.pi*self.tausqr)
                          * np.exp(-0.5 * (self.xi[:, np.newaxis] - self.mu)**2 / self.tausqr))
        q_ki = piNp / np.sum(piNp, axis=1)[:, np.newaxis]
        # Eqn (73)
        for i in range(self.N):
            self.G[i] = self.rng.multinomial(1, q_ki[i])
     
    def update_alpha_beta_sigma(self):  # Step 6
        print('update_alpha_beta_sigma')
        alpha_curr    = self.alpha
        beta_curr    = self.beta
        sig0_curr      = self.sig0
        k_curr          = self.k

        if self.ichain <= self.nBurn:
            
            alpha_prop, \
            beta_prop, \
            sig0_prop,k_prop = self.rng.normal(loc=[alpha_curr, \
                                                     beta_curr, \
                                                     sig0_curr,k_curr], \
                                                scale = [self.proposalscale_alpha, \
                                                         self.proposalscale_beta, \
                                                         self.proposalscale_sig0,self.proposalscale_k])       
           
        else:
            
            covProposal = np.cov(np.array((self.chain['alpha'][:self.ichain], \
                                           self.chain['beta'][:self.ichain],  \
                                           self.chain['sig0'][:self.ichain],   \
                                           self.chain['k'][:self.ichain])))

            alpha_prop, \
            beta_prop, \
            sig0_prop,k_prop = self.rng.multivariate_normal([alpha_curr, \
                                              beta_curr, \
                                              sig0_curr,k_curr], covProposal)
        
        sigsqr_xi_curr = self.calc_sigsqr(self.xi, sig0=sig0_curr, k=k_curr)
        sigsqr_xi_prop = self.calc_sigsqr(self.xi, sig0=sig0_prop, k=k_prop)     
        
        prior_curr = self.prior_all(alpha_curr, beta_curr, sig0_curr, k_curr)
        prior_prop = self.prior_all(alpha_prop, beta_prop, sig0_prop, k_prop)

        if prior_prop == 0:
            acceptanceProb = 0
            
        else:
            log_target_curr = np.sum(norm.logpdf(self.eta, loc=alpha_curr+beta_curr*self.xi, \
                                                 scale=np.sqrt(sigsqr_xi_curr))) \
                                                 + np.log(prior_curr)
                                                 
            log_target_prop = np.sum(norm.logpdf(self.eta, loc=alpha_prop+beta_prop*self.xi, \
                                                 scale=np.sqrt(sigsqr_xi_prop))) \
                                                 + np.log(prior_prop)

            try:
                acceptanceProb = min(1,np.exp(log_target_prop-log_target_curr))
            except FloatingPointError as exception:
                print('FLOATING POINT ERROR', exception)
                if log_target_prop < log_target_curr:
                    acceptanceProb = 0.
                else:
                    acceptanceProb = 1.

            u = self.rng.uniform()
            if u <= acceptanceProb:
                self.alpha = alpha_prop
                self.beta = beta_prop
                self.sig0 = sig0_prop
                self.k = k_prop
   
    def update_pi(self):  # Step 8
#        print('update_pi')
        # Eqn (82)
        self.nk = np.sum(self.G, axis=0)
#        print(self.nk)
        # Eqn (81)
        self.pi = self.rng.dirichlet(self.nk+1)

    def update_mu(self):  # Step 9
#        print('update_mu')
        Gsum = np.sum(self.G * self.xi[:, np.newaxis], axis=0)
        for k in range(self.nGaussXi):
            if self.nk[k] != 0:
                # Eqn (86)
                Sigma_muhat_k = 1.0/(1.0/self.usqr + self.nk[k]/self.tausqr[k])
                # Eqn (85)
                xibar_k = 1.0/self.nk[k] * Gsum[k]
                # Eqn (84)
                muhat_k = Sigma_muhat_k * (self.mu0/self.usqr + self.nk[k]/self.tausqr[k]*xibar_k)
                # Eqn (83)
                self.mu[k] = self.rng.normal(loc=muhat_k, scale=np.sqrt(Sigma_muhat_k))
            else:
                self.mu[k] = self.rng.normal(loc=self.mu0, scale=np.sqrt(self.usqr))

    def update_tausqr(self):  # Step 10
#        print('update_tausqr')
        # Eqn (88)
        nu_k = self.nk + 1
        # Eqn (89)
        tk_sqr = 1.0/nu_k * (self.wsqr + np.sum(self.G*(self.xi[:, np.newaxis]-self.mu)**2, axis=0))
        # Eqn (87)
        self.tausqr = tk_sqr * nu_k / self.rng.chisquare(nu_k, size=self.nGaussXi)

    def update_mu0(self):  # Step 11
#        print('update_mu0')
        # Eqn (94)
        mubar = np.mean(self.mu)
        # Eqn (93)
        self.mu0 = self.rng.normal(loc=mubar, scale=np.sqrt(self.usqr/self.nGaussXi))

    def update_usqr(self):  # Step 12
#        print('update_usqr')
        # Eqn (96)
        nu_u = self.nGaussXi + 1
        # Eqn (97)
        usqrhat = 1.0/nu_u * (self.wsqr + np.sum((self.mu - self.mu0)**2))
        usqr = np.inf
        while not usqr <= self.usqrmax:
            usqr = usqrhat * nu_u / self.rng.chisquare(nu_u)
        self.usqr = usqr

    def update_wsqr(self):  # Step 13
#        print('update_wsqr')
        # Eqn (102)
        a = 0.5 * (self.nGaussXi + 3)
        # Eqn (103)
        b = 0.5 * (1.0/self.usqr + np.sum(1.0/self.tausqr))
        # Eqn (101)
        self.wsqr = self.rng.gamma(a, 1.0/b)

    def update_GBeagle(self):
#        print('update_GBeagle')
        for i in range(self.N):

            prob_arr = np.zeros(self.nGaussBeagle)
            for j in range(self.nGaussBeagle):
                if self.xsigArr[i][j] == 0 and self.ysigArr[i][j] == 0:
                    prob_arr[j] = 0.
                else:
                    cov = [[self.xsigArr[i][j]**2,self.xycovArr[i][j]],[self.xycovArr[i][j],self.ysigArr[i][j]**2]]
                    prob_arr[j] = np.log(self.piBeagle[i,j])+multivariate_normal.logpdf([self.xi[i],self.eta[i]],[self.xArr[i,j],self.yArr[i,j]],cov)
            if (np.sum(prob_arr) == 0):
                for j in range(self.nGaussBeagle):
                    prob_arr[j] = self.piBeagle[i,j]
                sys.exit()
                #It looks like sometimes with small scatter at the position of some points, the probability is very close to zero
                #In these instances I'm drawing the next gaussian based on input weights... not great, but not sure I can do much else
            else:
                #sometimes the probabilities are really small so we want to be careful - if there is one that is so much smaller we should
                #set the probability to zero
                deltaProb = prob_arr-np.max(prob_arr)
                tempIdx = np.where(deltaProb < -50)[0] #remember this is a difference of 50 in log so it's really a huge difference!
                if len(tempIdx) > 0:
                  prob_arr[tempIdx] = 0.
                  tempIdx = np.where(deltaProb >= -50)[0]
                  prob_arr[tempIdx] = np.exp(prob_arr[tempIdx]-np.min(prob_arr[tempIdx]))
                else:
                  prob_arr = np.exp(prob_arr - np.min(prob_arr))
            prob_arr = prob_arr/np.sum(prob_arr)
            self.GBeagle[i] = self.rng.multinomial(1, prob_arr)
            tempIdx = np.where(self.GBeagle[i] == 1)[0]
            self.x[i] = self.xArr[i,tempIdx]
            self.y[i] = self.yArr[i,tempIdx]
            self.xsig[i] = self.xsigArr[i,tempIdx]
            self.ysig[i] = self.ysigArr[i,tempIdx]
            self.xycov[i] = self.xycovArr[i,tempIdx]
            self.xycorr[i] = self.xycov[i] / (self.xsig[i] * self.ysig[i])
            self.xvar[i] = self.xsig[i]**2
            self.yvar[i] = self.ysig[i]**2        

    def update_chain(self):
#        print('update_chain')
        ximean = np.sum(self.pi * self.mu)
        xisig = np.sqrt(np.sum(self.pi * (self.tausqr + self.mu**2)) - ximean**2)
        
        self.chain['alpha'][self.ichain] = self.alpha
        self.chain['beta'][self.ichain] = self.beta
        self.chain['sig0'][self.ichain] = self.sig0
        self.chain['k'][self.ichain] = self.k
        self.chain['xi_min'][self.ichain] = self.xi_min        
        self.chain['xi_max'][self.ichain] = self.xi_max              
        self.chain['xi'][self.ichain] = self.xi
        self.chain['eta'][self.ichain] = self.eta
        self.chain['pi'][self.ichain] = self.pi
        self.chain['mu'][self.ichain] = self.mu
        self.chain['tausqr'][self.ichain] = self.tausqr
        self.chain['mu0'][self.ichain] = self.mu0
        self.chain['usqr'][self.ichain] = self.usqr
        self.chain['wsqr'][self.ichain] = self.wsqr
        self.chain['ximean'][self.ichain] = ximean
        self.chain['xisig'][self.ichain] = xisig

        self.ichain += 1

    # =============================================================================
    # END OF STEP
    # =============================================================================

    # =============================================================================
    # INITIALISE, STEP and EXTEND
    # =============================================================================
    
    def initialize_chain(self, chain_length, nBurn):
        print('initialize_chain')
        self.chain_dtype = [('alpha', float),
                            ('beta', float),
                            ('sig0', float),
                            ('k', float),
                            ('xi_min', float),
                            ('xi_max', float),                            
                            ('xi', (float, self.N)),
                            ('eta', (float, self.N)),
                            ('pi', (float, self.nGaussXi)),
                            ('mu', (float, self.nGaussXi)),
                            ('tausqr', (float, self.nGaussXi)),
                            ('mu0', float),
                            ('usqr', float),
                            ('wsqr', float),
                            ('ximean', float),
                            ('xisig', float)]
        
        self.chain = np.empty((chain_length,), dtype=self.chain_dtype)
        self.ichain = 0
        self.nBurn = nBurn
        self.testIter = testIter

    def step(self, niter):
        print('step')
        for i in range(niter):
#            print(i)
            old_settings = np.seterr(divide='ignore', invalid='ignore')
            np.seterr(**old_settings)
            
            self.update_cens_y()            
            self.update_xi()
            self.update_eta()
            self.update_G()    
            self.update_alpha_beta_sigma()    
            self.update_pi()
            self.update_mu()
            self.update_tausqr()    
            self.update_mu0()
            self.update_usqr()
            self.update_wsqr()    
            self.update_GBeagle() # THIS IS ADDED FOR GMM    
            self.update_chain()
            
            if (i+1)%(niter/10.0) == 0:
                print(i+1, niter)

    def extend(self, length):
        print('extend')
        extension = np.empty((length), dtype=self.chain_dtype)
        self.chain = np.hstack((self.chain, extension))
        
    # =============================================================================
    # ADDED FOR MH STEP
    # =============================================================================

    def dunif(self, val, min, max):
#        print('dunif')
        if val > max or val < min:
           return 0.
        else:
           return 1./(max-min)

    def prior_alpha(self,alpha):
#        print('prior_alpha')
        min_alpha = -20
        max_alpha = 20
        return self.dunif(alpha, min_alpha, max_alpha)

    def prior_beta(self,beta):
#        print('prior_beta')
        min_beta = -10
        max_beta = 10
        return self.dunif(beta, min_beta, max_beta)

    def prior_sig0(self, sig0):
#        print('prior_sig0')
        min_sig0 = 0.001
        max_sig0 = 5.0
        return self.dunif(sig0, min_sig0, max_sig0)

    def prior_k(self, k):
#        print('prior_k')
        min_k = 0.001
        max_k = 5.0
        return self.dunif(k, min_k, max_k)

    def prior_all(self, alpha, beta, sig0, k):
#        print('prior_all')
        return self.prior_alpha(alpha)*self.prior_beta(beta)*self.prior_sig0(sig0)*self.prior_k(k)
     
    def calc_sigsqr(self, xi, sig0=None, k=None, xi_min=None, xi_max=None):
#        print('calc_sigsqr')
        if sig0 == None:
            sig0 = self.sig0
        if k == None:
            k = self.k
        if xi_min == None:
            xi_min = self.xi_min
        if xi_max == None:
            xi_max = self.xi_max        
            
        sigsqr = ( sig0 * ( ((1.0-k)*(xi-xi_max)/(xi_max-xi_min)) + 1.0 ) ) ** 2.0
        return sigsqr  

    def get_2d_mean_cov(self, i): # 2d means no redshift dependence
#        print('get_2d_mean_cov')
        x = self.x[i]
        y = self.y[i]
        
        xvar = self.xvar[i]
        yvar = self.yvar[i]
        
        xycov = self.xycov[i]

        mean = np.array([x,y])
        cov = np.array([[xvar, xycov],[xycov, yvar]])
        return mean, cov


class LinMix(object):
    """ A class to perform linear regression of `y` on `x` when there are measurement errors in
    both variables.  The regression assumes:

    eta = alpha + beta * xi + epsilon

    x = xi + xerr

    y = eta + yerr

    Here, `alpha` and `beta` are the regression coefficients, `epsilon` is the intrinsic random
    scatter about the regression, `xerr` is the measurement error in `x`, and `yerr` is the
    measurement error in `y`.  `epsilon` is assumed to be normally-distributed with mean zero and
    variance `sigsqr`.  `xerr` and `yerr` are assumed to be normally-distributed with means equal
    to zero, variances `xsig`^2 and `ysig`^2, respectively, and covariance `xycov`. The
    distribution of `xi` is modelled as a mixture of normals, with group proportions `pi`, means
    `mu`, and variances `tausqr`.

    Args:
        x(array_like): The observed independent variable.
        y(array_like): The observed dependent variable.
        xsig(array_like): 1-sigma measurement errors in x.
        ysig(array_like): 1-sigma measurement errors in y.
        xycov(array_like): Covariance between the measurement errors in x and y.
        delta(array_like): Array indicating whether a data point is censored (i.e., not detected),
            or not.  If delta[i] == 1, then the ith source is detected.  If delta[i] == 0, then
            the ith source is not detected and y[i] will be interpreted as an upper limit.  Note
            that if there are censored data points, then the maximum-likelihood estimate
            (alpha, beta, sigsqr) is not valid.  By default, all data points are assumed to be
            detected.
        nGaussXi(int): The number of Gaussians to use in the mixture model for the distribution of xi.
        nChains(int): The number of Monte Carlo Markov Chains to instantiate.
        parallelize(bool): Use a separate thread for each chain.  Only makes sense for nChains > 1.
        seed(int): Random seed.  If `None`, then get seed from np.random.randint().

    Attributes:
        nChains(int): The number of instantiated MCMCs.
        chain(numpy recarray): The concatenated MCMCs themselves.  Actually, only the concatenation
            of the last half of each chain is stored here after convergence is reached.  The
            recarray has the following columns:
                - alpha(float): The regression intercept.
                - beta(float): The regression slope.
                - sigsqr(float): The regression intrinsic scatter.
                - pi(array_like): The mixture model component fractions.
                - mu(array_like): The mixture model component means.
                - tausqr(array_like): The mixture model component variances.
                - mu0(float): The hyperparameter describing the prior variance of the distribution
                    of mixture means.
                - usqr(float): The hyperparameter describing the prior variance of the distribution
                    of mixture variances.
                - wsqr(float): The hyperparameter describing the typical scale for the prior on
                    `usqr` and `tausqr`.
                - ximean(float): The mean of the distribution for the independent latent variable
                    `xi`.
                - xisig(float): The standard deviation of the distribution for the independent
                    latent variable `xi`.
                - corr(float): The linear correlation coefficient between the latent dependent and
                    independent variables `xi` and `eta`.
    """
    def __init__(self, xArr, yArr, xsigArr, ysigArr, xycovArr, proposalscale_xi, proposalscale_eta, proposalscale_alpha, proposalscale_beta, proposalscale_sig0, proposalscale_k, delta=None, nGaussXi=3, nChains=4, parallelize=True, seed=None, nGaussBeagle=1, piBeagle=None):
        print('__init__ linmix')
        self.nChains = nChains
        self.parallelize = parallelize

        if seed is None:
            seed = np.random.randint(2**32-1)

        if piBeagle is None:
            piBeagle = np.ones(len(xArr)) # lester changed x to xArr
        if self.parallelize:
            # Will place 1 chain in 1 thread.
            from multiprocessing import Process, Pipe
            # Create a pipe for each thread.
            self.pipes = []
            slave_pipes = []
            for i in range(self.nChains):
                master_pipe, slave_pipe = Pipe()
                self.pipes.append(master_pipe)
                slave_pipes.append(slave_pipe)

            # Create chain pool.
            self.pool = []
            for sp in slave_pipes:
                self.pool.append(Process(target=task_manager, args=(sp,)))
                self.pool[-1].start()

            init_kwargs0 = {'xArr':xArr,
                            'yArr':yArr,
                            'xsigArr':xsigArr,
                            'ysigArr':ysigArr,
                            'xycovArr':xycovArr,
                            'proposalscale_xi':proposalscale_xi,
                            'proposalscale_eta':proposalscale_eta,
                            'proposalscale_alpha':proposalscale_alpha,
                            'proposalscale_beta':proposalscale_beta,
                            'proposalscale_sig0':proposalscale_sig0,
                            'proposalscale_k':proposalscale_k,
                            'delta':delta,
                            'nGaussXi':nGaussXi,
                            'nChains':self.nChains,
                            'nGaussBeagle':nGaussBeagle,
                            'piBeagle':piBeagle}
            for i, p in enumerate(self.pipes):
                init_kwargs = init_kwargs0.copy()
                init_kwargs['rng'] = np.random.RandomState(seed+i)
                p.send({'task':'init',
                        'init_args':init_kwargs})
        else:
            self._chains = []
            for i in range(self.nChains):
                self._chains.append(Chain(xArr, yArr, xsigArr, ysigArr, xycovArr, delta, nGaussXi, nChains, nGaussBeagle, piBeagle, proposalscale_xi, proposalscale_eta, proposalscale_alpha, proposalscale_beta, proposalscale_sig0, proposalscale_k, rng=None))
                
                self._chains[-1].initial_guess()

    def _get_psi(self):
        print('_get_psi')
        if self.parallelize:
            for p in self.pipes:
                p.send({'task':'fetch',
                        'key':'chain'})
            chains = [p.recv() for p in self.pipes]
            self.pipes[0].send({'task':'fetch',
                                'key':'ichain'})
            ndraw = int(self.pipes[0].recv()/2) # half length of total chain so far
        else:
            chains = [c.chain for c in self._chains]
            ndraw = int(self._chains[0].ichain/2) # half length of total chain so far
        psi={}
        keys=['alpha','beta','sig0','k']
        for key in keys:
            psi[key] = np.vstack([c[key][ndraw:2*ndraw] for c in chains]).T # c['alpha'] is minIter long
        
        # this multiplies [pi1, pi2, pi3] by [mu1, mu2, mu3] per iteration -> [pi1mu1, pi2mu2, pi3mu3]
        # then adds them up, to get pi1mu1 + pi2mu2 + pi3mu3 per iteration
        # then creates [ch1, ch2, ch3, ch4] x number of iterations
        psi['ximean'] = np.vstack([np.sum(c['pi'][ndraw:2*ndraw] * c['mu'][ndraw:2*ndraw], axis=1) for c in chains]).T
        
        # similar to above, but gets [pi1*(tausqr1+mu1), pi2*(tausqr2+mu2), pi3*(tausqr3+mu3)]
        # which becomes [pi1*(tausqr1+mu1) + pi2*(tausqr2+mu2) + pi3*(tausqr3+mu3)] x number of iterations
        psi['xivar'] = np.vstack([np.sum(c['pi'][ndraw:2*ndraw] * (c['tausqr'][ndraw:2*ndraw] + c['mu'][ndraw:2*ndraw]**2), axis=1) for c in chains]).T - psi['ximean']**2

        keys=['alpha','beta','sig0','k','ximean','xivar']
        for key in keys:
            psi[key] = np.hstack((psi[key][:ndraw/2],psi[key][-ndraw/2:]))
        return psi

    def _get_Rhat(self):
        print('_get_Rhat')
        psi = self._get_psi()
        keys=['alpha','beta','sig0','k','ximean','xivar']
        n = len(psi[keys[0]]) # quarter of initial iterations, half discarded, then split in half
        m = len(psi[keys[0]][0]) # 2x number of chains, usually = 8
        psi_bar_dot_j = {}
        psi_bar_dot_dot = {}
        B = {}
        s_j_sqr = {}
        W = {}
        var_plus = {}
        Rhat = {}
        
        for key in keys:
            psi_bar_dot_j[key] = np.empty(m, dtype=float)
            psi_bar_dot_j[key] = np.mean(psi[key], axis=0) # 8x avg values of each chain
            psi_bar_dot_dot[key] = np.mean(psi_bar_dot_j[key]) # 1x average of averages
            B[key] = (m*n/(m-1.0)) * np.mean((psi_bar_dot_j[key]-psi_bar_dot_dot[key])**2) # between sequence   
            s_j_sqr[key] = np.empty(m, dtype=float)
            s_j_sqr[key] = np.sum((psi[key]-psi_bar_dot_j[key])**2, axis=0) / (n-1.0)
            W[key] = np.mean(s_j_sqr[key]) # within sequence
            var_plus[key] = (((n-1.0)/n)*W[key]) + (B[key]/n)
            Rhat[key] = np.sqrt(var_plus[key]/W[key])
        return Rhat

    def _initialize_chains(self, minIter, nBurn):
        print('_initialize_chains')
        if self.parallelize:
            for p in self.pipes:
                p.send({'task':'init_chain',
                        'minIter':minIter,'nBurn':nBurn})
        else:
            for c in self._chains:
                c.initialize_chain(minIter,nBurn)

    def _step(self, niter):
        print('_step')
        if self.parallelize:
            for p in self.pipes:
                p.send({'task':'step',
                        'niter':niter})
        else:
            for c in self._chains:
                c.step(niter)

    def _extend(self, niter):
        print('_extend')
        if self.parallelize:
            for p in self.pipes:
                p.send({'task':'extend',
                        'niter':niter})
        else:
            for c in self._chains:
                c.extend(niter)

    def _build_chain(self, ikeep):
        print('_build_chain')
        if self.parallelize:
            for p in self.pipes:
                p.send({'task':'fetch',
                        'key':'chain'})
            self.chain = np.hstack([p.recv()[ikeep:] for p in self.pipes])
        else:
            self.chain = np.hstack([c.chain[ikeep:] for c in self._chains])

    def run_mcmc(self, minIter=5000, maxIter=100000, silent=False, nBurn=500, checkIter=5000, RhatLimit=5000, testIter=200):
        print('run_mcmc')
        """ Run the Markov Chain Monte Carlo for the LinMix object.

        Bayesian inference is employed, and a Markov chain containing random draws from the
        posterior is developed.  Convergence of the MCMC to the posterior is monitored using the
        potential scale reduction factor (RHAT, Gelman et al. 2004). In general, when RHAT < 1.1
        then approximate convergence is reached.  After convergence is reached, the second halves
        of all chains are concatenated and stored in the `.chain` attribute as a numpy recarray.

        Args:
            minIter(int): The minimum number of iterations to use.
            maxIter(int): The maximum number of iterations to use.
            silent(bool): If true, then suppress updates during sampling.
        """
        
        self._initialize_chains(minIter,nBurn)
        self.nBurn = nBurn # this is just for when MH starts to "learn" from previous
        self.testIter = testIter
        
        Rhat_array = []
        for i in range(0, minIter, checkIter): # 0 to minIter, whilst checking for convergence every "checkIter"
            self._step(checkIter)   
            Rhat = self._get_Rhat()
            Rhat_array.append(Rhat)
            print()
            print('Iteration: ', i+checkIter)
            print('Rhat values:')
            print('alpha:', Rhat['alpha'])
            print('beta:', Rhat['beta'])
            print('sig0:', Rhat['sig0'])
            print('k:', Rhat['k'])
            print('ximean:', Rhat['ximean'])
            print('xivar:', Rhat['xivar'])
            print()
            
        i += checkIter

        while (max(np.array(Rhat.values())[~np.isnan(Rhat.values())]) > RhatLimit) and (i < maxIter):
            self._extend(checkIter)
            self._step(checkIter)
            Rhat = self._get_Rhat()
            Rhat_array.append(Rhat)
            print()
            print('Extended Iteration: ', i+checkIter)
            print('Rhat values:')
            print('alpha:', Rhat['alpha'])
            print('beta:', Rhat['beta'])
            print('sig0:', Rhat['sig0'])
            print('k:', Rhat['k'])
            print('ximean:', Rhat['ximean'])
            print('xivar:', Rhat['xivar'])
            print()
            i += checkIter
        
        self._build_chain(0) # this means build chain from 0th iteration opposed to int(i/2) etc.
       
        # Clean up threads
        if self.parallelize:
            for p in self.pipes:
                p.send({'task':'kill'})

# =============================================================================
# NEW SIMPLIFIED ATTEMPT - LOPPING REDSHIFT
# =============================================================================

#options = ['z1p3', 'z2p0', 'z3p0', 'z4p0', 'z5p0']
#options = ['z1p3']
options = ['z5p0']

for i in range(len(options)):

#    data = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/kelly/data/scenario_5_7_data_{}.p'.format(options[i]),'r'))
    data = pickle.load(open('/Users/lester/Documents/GitHub/Local-Python/Year2/Project1/kelly/data/scenario_6_7_data_{}.p'.format(options[i]),'r'))
    
    piBeagle                = data['amp_GMM']       # 3x probability of each posterior gaussian
    GMMx                    = data['x_GMM']         # 3x posterior means per mass
    GMMy                    = data['y_GMM']         # 3x posterior means per sfr
    GMMxsig                 = data['xsig_GMM']      # 3x posterior sigmas per mass
    GMMysig                 = data['ysig_GMM']      # 3x posterior sigmas per sfr
    GMMxycov                = data['xycov_GMM']     # 3x posterior covar per mass-sfr pair

    nGaussXi                = 3                     # 3 #gaussians modelling xi
    nGaussBeagle            = 3                     # 3 #gaussians modelling BEAGLE posterior

    nChains                 = 1
    minIter                 = 80                   # if i%(niter/10.0) == 0:, TENTH OF THIS NEEDS TO DIVIDE BY 4, see yellow below
    maxIter                 = minIter
    checkIter               = minIter/10           # must be divisible by 4
    nBurn                   = minIter/4           # this is when the MH proposal distributions start adapting (xi, alpha, beta, sig0 and k)
    testIter                = minIter/10                   # xi MH proposal checks every "testIter" iterations and adjusts proposal if needed
    RhatLimit               = 1.01
    
    proposalscale_xi        = 0.25
    proposalscale_eta       = 0.25
    proposalscale_alpha     = 0.5                   # usually 0.005
    proposalscale_beta      = 0.05                  # usually 0.005
    proposalscale_sig0      = 0.05                  # usually 0.005
    proposalscale_k         = 0.05                  # usually 0.01
    
    parallelize             = False
    
    lm = LinMix(GMMx, GMMy, GMMxsig, GMMysig, xycovArr=GMMxycov, nGaussXi=nGaussXi, nGaussBeagle=nGaussBeagle, piBeagle=piBeagle, nChains=nChains, parallelize=parallelize, proposalscale_xi=proposalscale_xi, proposalscale_eta=proposalscale_eta, proposalscale_alpha=proposalscale_alpha, proposalscale_beta=proposalscale_beta, proposalscale_sig0=proposalscale_sig0, proposalscale_k=proposalscale_k)
    
    lm.run_mcmc(minIter=minIter, maxIter=maxIter, nBurn=nBurn, checkIter=checkIter, testIter=testIter, RhatLimit=RhatLimit)
    
    


#print(len(lm.chain['alpha']))
#print(lm.chain['alpha'])
#    pickle.dump(lm.chain, open('/Users/lester/Documents/linmix_files/lm_chain_new_001.p','w'))
#    pickle.dump(lm.chain, open('/Users/lester/Documents/linmix_files/lm_chain_scenario_5_7_data_{}_k0.p'.format(options[i]),'w'))
#    pickle.dump(lm.chain, open('/Users/lester/Documents/linmix_files/lm_chain_scenario_5_7_data_{}_k0_002.p'.format(options[i]),'w'))

#    pickle.dump(lm.chain, open('/Users/lester/Documents/linmix_files/lm_chain_scenario_6_7_data_{}_k0_convergence_test_001.p'.format(options[i]),'w')) # 12000 1000 4


#    pickle.dump(lm.chain, open('/Users/lester/Documents/linmix_files/lm_chain_scenario_6_7_data_{}_k0_convergence_test_003.p'.format(options[i]),'w'))

#    pickle.dump(lm.chain, open('/Users/lester/Documents/linmix_files/lm_chain_scenario_6_7_data_{}_k0_convergence_test_004.p'.format(options[i]),'w'))

#    pickle.dump(lm.chain, open('/Users/lester/Documents/linmix_files/lm_chain_scenario_6_7_data_{}_k_convergence_test_004.p'.format(options[i]),'w'))
    
#    pickle.dump(lm.chain, open('/Users/lester/Documents/linmix_files/lm_chain_scenario_6_7_data_{}_k_convergence_test_005.p'.format(options[i]),'w'))
    
    
    
#    pickle.dump(lm.chain, open('/Users/lester/Documents/linmix_files/lm_chain_scenario_6_7_data_{}_k_convergence_test_006.p'.format(options[i]),'w'))
    
#    pickle.dump(lm.chain, open('/Users/lester/Documents/linmix_files/lm_chain_scenario_6_7_data_{}_k_convergence_test_007.p'.format(options[i]),'w')) 
    
#    pickle.dump(lm.chain, open('/Users/lester/Documents/linmix_files/lm_chain_scenario_6_7_data_{}_k_convergence_test_008.p'.format(options[i]),'w')) 
    
    
    
    
    
    
    
    
