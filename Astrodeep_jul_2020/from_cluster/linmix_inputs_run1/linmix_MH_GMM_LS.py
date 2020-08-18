""" linmix -- A hierarchical Bayesian approach to linear regression with error in both X and Y.
"""

from __future__ import print_function
from scipy.stats import norm,multivariate_normal,chi2,invgamma,truncnorm
import numpy as np
import matplotlib.pyplot as plt
import sys

def task_manager(conn):
    chain = None
    while True:
        message = conn.recv()
        if message['task'] == 'init':
            print('tm_init')
            chain = Chain(**message['init_args'])
            chain.initial_guess()
        elif message['task'] == 'init_chain':
            print('tm_init_chain')
            chain.initialize_chain(message['miniter'],message['n_burn'])
        elif message['task'] == 'step':
            print('tm_step')
            chain.step(message['niter'])
        elif message['task'] == 'extend':
            print('tm_extend')
            chain.extend(message['niter'])
        elif message['task'] == 'fetch':
            print('tm_fetch')
            conn.send(chain.__dict__[message['key']])
        elif message['task'] == 'kill':
            print('tm_kill')
            break
        else:
            print('tm_invalidtask')
            raise ValueError("Invalid task")

class Chain(object):
    def __init__(self, xArr, yArr, xsigArr, ysigArr, xycovArr, delta, K, nchains, nGMM_err, pi_err, rng=None):
        self.xArr = np.array(xArr, dtype=float)
        self.yArr = np.array(yArr, dtype=float)
        self.xsigArr = np.array(xsigArr, dtype=float)
        self.ysigArr = np.array(ysigArr, dtype=float)
        self.xycovArr = np.array(xycovArr, dtype=float)
        self.x = np.array(xArr[:,0], dtype=float)
        self.y = np.array(yArr[:,0], dtype=float)
        self.nGMM_err = nGMM_err
        self.pi_err = np.array(pi_err, dtype=float)
        
        if xsigArr is None:
            print('error, have not coded for this circumstance')
            sys.exit
            self.xsig = np.zeros_like(self.x)
            xycov = np.zeros_like(self.x)
        else:
            self.xsig = np.array(xsigArr[:,0], dtype=float)
        if ysigArr is None:
            print('error, have not coded for this circumstance')
            sys.exit
            self.ysig = np.zeros_like(self.y)
            xycov = np.zeros_like(self.y)
        else:
            self.ysig = np.array(ysigArr[:,0], dtype=float)
        self.wxerr = (self.xsig != 0.0)
        self.wyerr = (self.ysig != 0.0)
        self.werrs = werrs = self.wxerr & self.wyerr
        
        if xycovArr is None:
            print('error, have not coded for this circumstance')
            sys.exit
            self.xycov = np.zeros_like(self.x)
        else:
            self.xycov = np.array(xycovArr[:,0], dtype=float)
        
        self.xycorr = np.zeros_like(self.xycov)
        self.xycorr[werrs] = self.xycov[werrs] / (self.xsig[werrs] * self.ysig[werrs])
        
        self.N = len(self.x)
        self.K = K
        self.nchains = nchains
        
        self.xvar = self.xsig**2
        self.yvar = self.ysig**2
        
        #parameters for the MCMC chain on update_xi
        self.accept_xi = 0
        self.reject_xi = 0
        self.proposal_xi = 0.25
        self.proposal_alpha = 0.001
        self.proposal_beta = 0.01
        #lester
        self.proposal_sig_0 = 0.005
        self.proposal_k = 0.01
        
        proposal_xi = self.proposal_xi
        proposal_alpha = self.proposal_alpha
        proposal_beta = self.proposal_beta
        proposal_sig_0 = self.proposal_sig_0
        proposal_k = self.proposal_k
        print('proposals', proposal_xi, proposal_alpha, proposal_beta, proposal_sig_0, proposal_k)          
        
        if delta is None:
            self.delta = np.ones((self.N), dtype=bool)
        else:
            self.delta = np.array(delta, dtype=bool)

        if rng is None:
            rng = np.random.RandomState()
        self.rng = rng
        self.initialized = False

    def initial_guess(self):  # Step 1
        # For convenience
        #        x = self.x
        #        y = self.y
        #        xycov = self.xycov
        #        xvar = self.xvar
        #        yvar = self.yvar
        N = self.N
        K = self.K
        #The following now get assigned when G_err is assigned
        x = np.zeros(N)
        y = np.zeros(N)
        xycov = np.zeros(N)
        xvar = np.zeros(N)
        yvar = np.zeros(N)
        nGMM_err = self.nGMM_err
        
        pi_err = self.pi_err #the relative amplitudes of the models are fixed
#        print(nGMM_err, 'LESTER')
        if nGMM_err == 1:
            self.G_err = np.ones(N, dtype=int)
        else:
            self.G_err = np.zeros((N, nGMM_err), dtype=int)
#            print(self.G_err)
            #just start by assigning G_err to the Gaussian with highest amplitude
            for i in range(N):
                maxind = np.argmax(pi_err[i])
                self.G_err[i,maxind] = 1
                self.x[i] = self.xArr[i,maxind]
                self.y[i] = self.yArr[i,maxind]
                self.xsig[i] = self.xsigArr[i,maxind]
                self.ysig[i] = self.ysigArr[i,maxind]
                self.xycov[i] = self.xycovArr[i,maxind]
                self.xvar[i] = self.xsig[i]**2
                self.yvar[i] = self.ysig[i]**2
                self.xycorr[i] = self.xycov[i] / (self.xsig[i] * self.ysig[i])
#        print(self.G_err)
        x = self.x
        y = self.y
        xsig = self.xsig
        ysig = self.ysig
        xycov = self.xycov
        xvar = self.xvar
        yvar = self.yvar
        xycorr = self.xycorr


        # Use BCES estimator for initial guess of theta = {alpha, beta, sigsqr}
        self.beta = ((np.cov(x, y, ddof=1)[1, 0] - np.mean(xycov))
                     / (np.var(x, ddof=1) - np.mean(xvar)))
        self.alpha = np.mean(y) - self.beta * np.mean(x)
        #print(self.beta, self.alpha)
        #lester
#        print('h', self.alpha, self.beta)
#        self.alpha = -6.2
#        self.beta = 0.8
        self.sigsqr = np.var(y, ddof=1) - np.mean(yvar) - self.beta * (np.cov(x, y, ddof=1)[1, 0]
                                                                       - np.mean(xycov))
        self.sigsqr = np.max([self.sigsqr,
                              0.05 * np.var(y - self.alpha - self.beta * x, ddof=1)])


        self.mu0 = np.median(x)
        self.wsqr = np.var(x, ddof=1) - np.median(xvar)
        self.wsqr = np.max([self.wsqr, 0.01*np.var(x, ddof=1)])

        # Now get an MCMC value dispersed around above values
        X = np.ones((N, 2), dtype=float)
        X[:, 1] = x
        Sigma = np.linalg.inv(np.dot(X.T, X)) * self.sigsqr
        coef = self.rng.multivariate_normal([0, 0], Sigma)
        chisqr = self.rng.chisquare(self.nchains)
        self.alpha += coef[0] * np.sqrt(1.0/chisqr)
        self.beta += coef[1] * np.sqrt(1.0/chisqr)
        self.sigsqr *= 0.5 * N / self.rng.chisquare(0.5*N)
        
        #The parameters defining the xi-dependence of sigma
        # lester
        self.sig_0 = np.sqrt(self.sigsqr)
        self.k = 1.0                                                            # need to add a random element to this
        self.xi_min = 6.0
        self.xi_max = 12.0
#        print('h', self.sig_0, self.k)

        # Now get the values for the mixture parameters, first do prior params
        self.mu0min = min(x)
        self.mu0max = max(x)

        mu0g = np.nan
        while not (mu0g > self.mu0min) & (mu0g < self.mu0max):
            mu0g = self.mu0 + (self.rng.normal(scale=np.sqrt(np.var(x, ddof=1) / N)) /
                               np.sqrt(self.nchains/self.rng.chisquare(self.nchains)))
        self.mu0 = mu0g

        # wsqr is the global scale
        self.wsqr *= 0.5 * N / self.rng.chisquare(0.5 * N)

        self.usqrmax = 1.5 * np.var(x, ddof=1)
        self.usqr = 0.5 * np.var(x, ddof=1)

        self.tausqr = 0.5 * self.wsqr * self.nchains / self.rng.chisquare(self.nchains, size=K)

        self.mu = self.mu0 + self.rng.normal(scale=np.sqrt(self.wsqr), size=K)

        # get initial group proportions and group labels

        pig = np.zeros(self.K, dtype=float)
        if K == 1:
            self.G = np.ones(N, dtype=int)
            self.pi = np.array([1], dtype=float)
        else:
            self.G = np.zeros((N, K), dtype=int)
            for i in range(N):
                minind = np.argmin(abs(x[i] - self.mu))
                pig[minind] += 1
                self.G[i, minind] = 1
            self.pi = self.rng.dirichlet(pig+1)

        self.eta = y.copy()
        self.y_ul = y.copy()
        self.xi = x.copy()

        self.cens = np.nonzero(np.logical_not(self.delta))[0]

        self.initialized = True

    def update_cens_y(self):  # Step 2
        todo = self.cens[:]
#        print('todo', todo)
        while len(todo) > 0:
            self.y[todo] = self.rng.normal(loc=self.eta[todo],
                                           scale=np.sqrt(self.yvar[todo]),
                                           size=len(todo))
            todo = np.nonzero(np.logical_not(self.delta) & (self.y > self.y_ul))[0]

    def update_xi(self):  # Step 3
#        print('STEP 3')
        wxerr = self.wxerr                                                              # boolean array where xsig != 0

        #update each object in turn
        xi_curr = self.xi                                                               # current xi values
        sigsqr_xi_curr = self.calc_sigsqr(xi_curr)                                      # array of sigsqr values, one per xi
        xi_prop = self.rng.normal(size=len(xi_curr), scale=self.proposal_xi)+xi_curr    # new proposal for each xi
        sigsqr_xi_prop = self.calc_sigsqr(xi_prop)                                      # array of sigsqr values, one per new xi
        muArr = []                                                                      # mu are means of xi gaussians
        tausqrArr = []                                                                  # tausqr are sigsqr of xi gaussians
        for i in range(len(self.xi)):
            tempIdx = np.where(self.G[i] == 1)[0]                                       # [0], [1] or [2] decides which xi gaussian
            muArr.append(self.mu[tempIdx][0])                                           # chooses mu from [mu1, mu2, mu3]
            tausqrArr.append(self.tausqr[tempIdx][0])                                   # chooses tau from [ta1, ta2, ta3]

        
        # P(xi|psi)
        log_p_xi_psi_curr = norm.logpdf(xi_curr,loc=muArr,scale=np.sqrt(tausqrArr))     # gives log( gaussian value at xi_curr ) per object
        log_p_xi_psi_prop = norm.logpdf(xi_prop,loc=muArr,scale=np.sqrt(tausqrArr))     # gives log( gaussian value at xi_prop ) per object
               
        
        # P(xi|eta,theta) - this is only really true for the proportionality, I'm not calculating the correct normalization
        log_p_xi_eta_theta_curr = norm.logpdf(self.eta,scale=np.sqrt(sigsqr_xi_curr),\
                                              loc=self.alpha+self.beta*xi_curr)
        log_p_xi_eta_theta_prop = norm.logpdf(self.eta,scale=np.sqrt(sigsqr_xi_prop),\
                                              loc=self.alpha+self.beta*xi_prop)

# =============================================================================
#         # lester test, comparing P(xi|eta,theta) to analytical form
#         log_p_xi_eta_theta_prop_lester = norm.logpdf(self.eta[0],scale=np.sqrt(sigsqr_xi_prop),\
#                                               loc=self.alpha+self.beta*xi_prop)
#         plt.scatter(xi_prop, np.e**log_p_xi_eta_theta_prop_lester, marker='x')
#         
#         x_lester = np.linspace(7, 11, 1000)
#         sigma_lester = ( self.sig_0 * ( ((1.0-self.k)*(x_lester-self.xi_max) / (self.xi_max - self.xi_min)) + 1) )
#         fXgY = abs((1/(sigma_lester*((2*np.pi)**0.5))) * np.exp(-0.5 * (( (self.eta[0] - (self.alpha + (self.beta*x_lester))) / sigma_lester ) ** 2) ))
#         plt.plot(x_lester, fXgY)
# 
#         plt.show()
# =============================================================================


        for i in range(len(self.xi)):
            if wxerr[i] == True:
                # P(xi|x,y,eta)
                cov = [[self.xsig[i]**2,self.xycov[i]],[self.xycov[i],self.ysig[i]**2]]
                log_p_xi_x_y_eta_curr = multivariate_normal.logpdf([xi_curr[i],self.eta[i]],[self.x[i],self.y[i]],cov)

                # P(xi|x,y,eta)
                log_p_xi_x_y_eta_prop = multivariate_normal.logpdf([xi_prop[i],self.eta[i]],[self.x[i],self.y[i]],cov)

                #print(log_p_xi_psi_curr, log_p_xi_x_y_eta_curr, log_p_xi_eta_theta_curr)
                log_target_current = log_p_xi_psi_curr[i] + log_p_xi_x_y_eta_curr + log_p_xi_eta_theta_curr[i]
                log_target_prop = log_p_xi_psi_prop[i] + log_p_xi_x_y_eta_prop + log_p_xi_eta_theta_prop[i]
                #print("log_target_current: ", log_target_current)
                #print("log_target_prop: ",log_target_prop)

                try:
                  acceptanceProb = min(1,np.exp(log_target_prop-log_target_current))
                except FloatingPointError as exception:
                  #print("log_target_current: ", log_target_current)
                  #print("log_target_prop: ",log_target_prop)
                  #print(log_target_prop-log_target_current)
                  if log_target_prop < log_target_current:
                    acceptanceProb = 0.
                  else:
                    acceptanceProb = 1.
                
                u = self.rng.uniform()
#                if i== 0:
#                    print(xi_curr[i], xi_prop, log_target_current, log_target_prop, np.exp(log_target_prop-log_target_current), u)
                if u <= acceptanceProb:
#                    if i == 0:
#                        print("accept", log_target_current, log_target_prop)
                    self.xi[i] = xi_prop[i]
                    if i == 0:
                        self.accept_xi = self.accept_xi + 1
                else:
#                    if i == 0:
#                        print("reject", log_target_current, log_target_prop)
                    if i == 0:
                        self.reject_xi = self.reject_xi + 1

                test = self.accept_xi+self.reject_xi
                #print (self.accept_xi+self.reject_xi, test % 100)
                if self.ichain >= self.n_burn and test > 0 and test % 200 == 0:
                    if np.float(self.accept_xi)/np.float(self.accept_xi+self.reject_xi) > 0.5:
                        self.proposal_xi = self.proposal_xi*1.1
                    elif np.float(self.accept_xi)/np.float(self.accept_xi+self.reject_xi) < 0.4:
                        self.proposal_xi = self.proposal_xi*0.9
#                    print("self.proposal_xi: ", self.proposal_xi, np.float(self.accept_xi)/np.float(self.accept_xi+self.reject_xi))
                    self.accept_xi = 0
                    self.reject_xi = 0

    #print("self.xi: ", self.xi)


    def update_eta(self):  # Step 4
        
        wxerr = self.wxerr
        wyerr = self.wyerr

        etaxyvar = self.yvar * (1.0 - self.xycorr**2)
        etaxy = self.y.copy()
        #print('etaxy[71] before: ', etaxy[71])
        etaxy[wxerr] += (self.xycov / self.xvar * (self.xi - self.x))[wxerr]
        #print('etaxy[71] after: ', etaxy[71])

        # Eqn (68)
        sigsqr = self.calc_sigsqr(self.xi)
        sigma_etahat_i_sqr = 1.0/(1.0/etaxyvar + 1.0/sigsqr)
        # Eqn (67)
        etahat_i = (sigma_etahat_i_sqr * (etaxy / etaxyvar
                    + (self.alpha + self.beta * self.xi) / sigsqr))
        # Eqn (66)
        self.eta[wyerr] = self.rng.normal(loc=etahat_i[wyerr],
                                          scale=np.sqrt(sigma_etahat_i_sqr[wyerr]))
        #for i in range(len(etahat_i)):
        #    print('update_eta: ', i, etahat_i[i], sigma_etahat_i_sqr[i], self.eta[i], self.y[i], etaxy[i], etaxyvar[i], sigsqr[i], self.alpha+self.beta*self.xi[i], self.alpha, self.beta, self.x[i], self.yvar[i], self.xvar[i])




    def update_G(self):  # Step 5
        # Eqn (74)
        piNp = self.pi * (1.0/np.sqrt(2.0*np.pi*self.tausqr)
                          * np.exp(-0.5 * (self.xi[:, np.newaxis] - self.mu)**2 / self.tausqr))
        q_ki = piNp / np.sum(piNp, axis=1)[:, np.newaxis]
        # Eqn (73)
        for i in range(self.N):
            self.G[i] = self.rng.multinomial(1, q_ki[i])
     
    def update_alpha_beta(self):  # Step 6 
        # - change to MH update
#        X = np.ones((self.N, 2), dtype=float)
#        X[:, 1] = self.xi
#        # Eqn (77)
#        XTXinv = np.linalg.inv(np.dot(X.T, X))
#        Sigma_chat = XTXinv * self.sigsqr
#        # Eqn (76)
#        chat = np.dot(np.dot(XTXinv, X.T), self.eta)
#        # Eqn (75)
#        self.alpha, self.beta = self.rng.multivariate_normal(chat, Sigma_chat)
        alpha_curr = self.alpha
        beta_curr = self.beta
        #update sigsqr values at the same time
        
        # lester
        sig_0_curr = self.sig_0
        k_curr = self.k
        
        # lester
        if self.ichain <= self.n_burn:
#            print(self.ichain, self.n_burn)
            [alpha_prop,beta_prop,sig_0_prop,k_prop] = \
                self.rng.normal(loc=[alpha_curr,beta_curr,sig_0_curr,k_curr], \
                                scale = [self.proposal_alpha,self.proposal_beta,\
                                         self.proposal_sig_0,self.proposal_k])
        else:
            covProposal = np.cov(np.transpose(np.column_stack((self.chain['alpha'][:self.ichain],\
                          self.chain['beta'][:self.ichain],self.chain['sig_0'][:self.ichain],self.chain['k'][:self.ichain]))))
            alpha_prop,beta_prop,sig_0_prop,k_prop = \
                self.rng.multivariate_normal([alpha_curr,beta_curr,sig_0_curr, k_curr], \
                                    covProposal)

#        min_sigma_prop = truncnorm.rvs((0.-min_sigma_curr)/self.proposal_min_sigma,(1E3-min_sigma_curr)/self.proposal_min_sigma,\
#                                       scale=self.proposal_min_sigma,loc=min_sigma_curr)
#
#        sig_alpha_curr = calc_sig_alpha(self, sig_beta=sig_beta_curr)
#        sig_alpha_prop = calc_sig_alpha(self, sig_beta=sig_beta_prop)
        # lester 
        sigsqr_xi_curr = self.calc_sigsqr(self.xi, sig_0=sig_0_curr, k=k_curr)
        sigsqr_xi_prop = self.calc_sigsqr(self.xi, sig_0=sig_0_prop, k=k_prop)

        log_targetProp = np.sum(norm.logpdf(self.eta, loc=alpha_prop+beta_prop*self.xi,scale=np.sqrt(sigsqr_xi_prop)))
        #print('min max sig_sqrt_xi_curr: ', min(sigsqr_xi_curr), max(sigsqr_xi_curr))
        #print('min max sig_sqrt_xi_prop: ', min(sigsqr_xi_prop), max(sigsqr_xi_prop))
        prior_prop = self.prior_all(alpha_prop, beta_prop, sig_0_prop, k_prop)
        #print('********: ', alpha_prop, beta_prop, sig_alpha_prop, sig_beta_prop, min_sigma_prop, prior_prop)
        if prior_prop == 0:
            acceptanceProb = 0
        else:
          log_targetProp = log_targetProp + np.log(prior_prop)
       
          log_targetCurrent = np.sum(norm.logpdf(self.eta, loc=alpha_curr+beta_curr*self.xi,scale=np.sqrt(sigsqr_xi_curr)))
          prior_curr = self.prior_all(alpha_curr, beta_curr, sig_0_curr, k_prop)
          #print('prior_curr: ', alpha_curr, beta_curr, sig_alpha_curr, sig_beta_curr, prior_curr)
          log_targetCurrent = log_targetCurrent + np.log(prior_curr)
        
          #to maintain detailed balance when not using a symmetric proposal disribution on min_sigma
          log_accNum = log_targetProp
          log_accDenom = log_targetCurrent
          #print('log_targetProp, log_targetCurrent: ', log_targetProp, log_targetCurrent)
    
          #print(self.ichain, log_accNum, log_accDenom, log_accNum-log_accDenom)
          try:
            acceptanceProb = min(1,np.exp(log_accNum-log_accDenom))
          except FloatingPointError as exception:
            #print('log_accNum: ', log_accNum)
            #print('log_accDenom: ', log_accDenom)
            #print(log_accNum-log_accDenom)
            if log_accNum < log_accDenom:
              acceptanceProb = 0.
            else:
              acceptanceProb = 1.
        
        #print('acceptance prob: ', acceptanceProb)        
        u = self.rng.uniform()
        if u <= acceptanceProb:
            #print("accept: ", self.ichain, alpha_prop, beta_prop, sig_alpha_prop, sig_beta_prop, acceptanceProb)
            self.alpha = alpha_prop
            self.beta = beta_prop
            # lester
#            self.sig_alpha = sig_alpha_prop
#            self.sig_beta = sig_beta_prop
#            self.min_sigma = min_sigma_prop
            self.sig_0 = sig_0_prop
            self.k = k_prop
            
        #else:
            #print("reject: ", self.ichain, alpha_prop, beta_prop, sig_alpha_prop, sig_beta_prop, acceptanceProb)

    def update_pi(self):  # Step 8
        # Eqn (82)
        self.nk = np.sum(self.G, axis=0)
        # Eqn (81)
        self.pi = self.rng.dirichlet(self.nk+1)

    def update_mu(self):  # Step 9
        Gsum = np.sum(self.G * self.xi[:, np.newaxis], axis=0)
        for k in range(self.K):
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
        # Eqn (88)
        nu_k = self.nk + 1
        # Eqn (89)
        tk_sqr = 1.0/nu_k * (self.wsqr + np.sum(self.G*(self.xi[:, np.newaxis]-self.mu)**2, axis=0))
        # Eqn (87)
        self.tausqr = tk_sqr * nu_k / self.rng.chisquare(nu_k, size=self.K)

    def update_mu0(self):  # Step 11
        # Eqn (94)
        mubar = np.mean(self.mu)
        # Eqn (93)
        self.mu0 = self.rng.normal(loc=mubar, scale=np.sqrt(self.usqr/self.K))

    def update_usqr(self):  # Step 12
        # Eqn (96)
        nu_u = self.K + 1
        # Eqn (97)
        usqrhat = 1.0/nu_u * (self.wsqr + np.sum((self.mu - self.mu0)**2))
        usqr = np.inf
        while not usqr <= self.usqrmax:
            usqr = usqrhat * nu_u / self.rng.chisquare(nu_u)
        self.usqr = usqr

    def update_wsqr(self):  # Step 13
        # Eqn (102)
        a = 0.5 * (self.K + 3)
        # Eqn (103)
        b = 0.5 * (1.0/self.usqr + np.sum(1.0/self.tausqr))
        # Eqn (101)
        self.wsqr = self.rng.gamma(a, 1.0/b)
        
            
    # =============================================================================
    #   
    # =============================================================================

    def initialize_chain(self, chain_length, n_burn):
        self.chain_dtype = [('alpha', float),
                            ('beta', float),
                            ('sig_0', float),
                            ('k', float),
                            ('xi_min', float),
                            ('xi_max', float),                            
                            ('xi', (float, self.N)),
                            ('eta', (float, self.N)),
                            ('pi', (float, self.K)),
                            ('mu', (float, self.K)),
                            ('tausqr', (float, self.K)),
                            ('mu0', float),
                            ('usqr', float),
                            ('wsqr', float),
                            ('ximean', float),
                            ('xisig', float),
                            ('corr', float)]
        
        self.chain = np.empty((chain_length,), dtype=self.chain_dtype)
        self.ichain = 0
        self.n_burn = n_burn

    def extend(self, length):
        extension = np.empty((length), dtype=self.chain_dtype)
        self.chain = np.hstack((self.chain, extension))

    def update_chain(self):
        self.chain['alpha'][self.ichain] = self.alpha
        self.chain['beta'][self.ichain] = self.beta
        
        # lester
        self.chain['sig_0'][self.ichain] = self.sig_0
        self.chain['k'][self.ichain] = self.k
        self.chain['xi_min'][self.ichain] = self.xi_min        
        self.chain['xi_max'][self.ichain] = self.xi_max              
        self.chain['xi'][self.ichain] = self.xi
        self.chain['pi'][self.ichain] = self.pi
        self.chain['mu'][self.ichain] = self.mu
        self.chain['tausqr'][self.ichain] = self.tausqr
        self.chain['mu0'][self.ichain] = self.mu0
        self.chain['usqr'][self.ichain] = self.usqr
        self.chain['wsqr'][self.ichain] = self.wsqr
        ximean = np.sum(self.pi * self.mu)
        self.chain['ximean'][self.ichain] = ximean
        xisig = np.sqrt(np.sum(self.pi * (self.tausqr + self.mu**2)) - ximean**2)
        self.chain['xisig'][self.ichain] = xisig
#        self.chain['corr'][self.ichain] = self.beta * xisig / np.sqrt(self.beta**2 * xisig**2
#                                                                      + self.sig_0) # lester changed sigsqr to sig_0 as sigsqr is only set once in step1 and remains the same
#        if self.ichain % 100 == 0 or self.ichain == 199:
#            print('**', self.ichain, self.chain['alpha'])
        self.ichain += 1

    def step(self, niter):
        for i in range(niter):
            self.update_cens_y()
            old_settings = np.seterr(divide='ignore', invalid='ignore')
            self.update_xi()
            self.update_eta()
            np.seterr(**old_settings)
            self.update_G()
            self.update_alpha_beta()
            self.update_pi()
            self.update_mu()
            self.update_tausqr()
            self.update_mu0()
            self.update_usqr()
            self.update_wsqr()
            self.update_Gerr()
            self.update_chain()
#            print('lester', self.alpha, self.beta, self.sig_0, self.k)

    # =============================================================================
    # ADDED FOR MH STEP
    # =============================================================================

    def update_Gerr(self):
        for i in range(self.N):
            #            #print 'self.rng.multinomial(1, self.pi_err[i]): ', self.rng.multinomial(1, self.pi_err[i])
            #            self.G_err[i] = self.rng.multinomial(1, self.pi_err[i])
            #            tempIdx = np.where(self.G_err[i] == 1)[0]
            #            self.x[i] = self.xArr[i,tempIdx]
            #            self.y[i] = self.yArr[i,tempIdx]
            #            self.xsig[i] = self.xsigArr[i,tempIdx]
            #            self.ysig[i] = self.ysigArr[i,tempIdx]
            #            self.xycov[i] = self.xycovArr[i,tempIdx]
            #            self.xycorr[i] = self.xycov[i] / (self.xsig[i] * self.ysig[i])
            #            self.xvar[i] = self.xsig[i]**2
            #            self.yvar[i] = self.ysig[i]**2
            #print 'self.rng.multinomial(1, self.pi_err[i]): ', self.rng.multinomial(1, self.pi_err[i])
            prob_arr = np.zeros(self.nGMM_err)
            for j in range(self.nGMM_err):
                if self.xsigArr[i][j] == 0 and self.ysigArr[i][j] == 0:
                    prob_arr[j] = 0.
                else:
                    cov = [[self.xsigArr[i][j]**2,self.xycovArr[i][j]],[self.xycovArr[i][j],self.ysigArr[i][j]**2]]
                    prob_arr[j] = np.log(self.pi_err[i,j])+multivariate_normal.logpdf([self.xi[i],self.eta[i]],[self.xArr[i,j],self.yArr[i,j]],cov)
            if (np.sum(prob_arr) == 0):
                for j in range(self.nGMM_err):
#                    print('cov, prob_arr[j]: ', i, cov, prob_arr[j])
#                    print('self.xi[i], self.xArr[i,j], self.eta[i], self.yArr[i,j]: ', self.xi[i], self.xArr[i,j], self.eta[i], self.yArr[i,j])
#                    print('prob_arr: ', self.xsigArr[i], self.ysigArr[i], prob_arr)
                    prob_arr[j] = self.pi_err[i,j]
                sys.exit()
                #It looks like sometimes with small scatter at the position of some points, the probability is very close to zero
                #In these instances I'm drawing the next gaussian based on input weights... not great, but not sure I can do much else
            else:
                #sometimes the probabilities are really small so we want to be careful - if there is one that is so much smaller we should
                #set the probability to zero
                #print('prob_arr', prob_arr)
                deltaProb = prob_arr-np.max(prob_arr)
                #print('deltaProb: ', deltaProb)
                tempIdx = np.where(deltaProb < -50)[0] #remember this is a difference of 50 in log so it's really a huge difference!
                if len(tempIdx) > 0:
                  prob_arr[tempIdx] = 0.
                  tempIdx = np.where(deltaProb >= -50)[0]
                  prob_arr[tempIdx] = np.exp(prob_arr[tempIdx]-np.min(prob_arr[tempIdx]))
                else:
                  prob_arr = np.exp(prob_arr - np.min(prob_arr))
                #print(prob_arr)
            prob_arr = prob_arr/np.sum(prob_arr)
            self.G_err[i] = self.rng.multinomial(1, prob_arr)
            tempIdx = np.where(self.G_err[i] == 1)[0]
            self.x[i] = self.xArr[i,tempIdx]
            self.y[i] = self.yArr[i,tempIdx]
            self.xsig[i] = self.xsigArr[i,tempIdx]
            self.ysig[i] = self.ysigArr[i,tempIdx]
            self.xycov[i] = self.xycovArr[i,tempIdx]
            self.xycorr[i] = self.xycov[i] / (self.xsig[i] * self.ysig[i])
            self.xvar[i] = self.xsig[i]**2
            self.yvar[i] = self.ysig[i]**2

    def dunif(self, val, min, max):
        if val > max or val < min:
            return 0.
        else:
            return 1./(max-min)

    def prior_alpha(self,alpha):
         min_alpha = -20
         max_alpha = 20
         return self.dunif(alpha, min_alpha, max_alpha)

    def prior_beta(self,beta):
         min_beta = -10
         max_beta = 10
         return self.dunif(beta, min_beta, max_beta)

    def prior_sig_0(self, sig_0):
         min_sig_0 = 0.001
         max_sig_0 = 5.0
         return self.dunif(sig_0, min_sig_0, max_sig_0)

    def prior_k(self, k):
         min_k = 0.2
         max_k = 5.0
         return self.dunif(k, min_k, max_k)

    def prior_all(self, alpha, beta, sig_0, k):
         return self.prior_alpha(alpha)*self.prior_beta(beta)*self.prior_sig_0(sig_0)*self.prior_k(k)
     
    def calc_sigsqr(self, xi, sig_0=None, k=None, xi_min=None, xi_max=None):
#        print('REALLY')
        # lester
        if sig_0 == None:
            sig_0 = self.sig_0
        if k == None:
            k = self.k
        if xi_min == None:
            xi_min = self.xi_min
        if xi_max == None:
            xi_max = self.xi_max        
            
        sigsqr = ( sig_0 * ( ((1.0-k)*(xi-xi_max)/(xi_max-xi_min)) + 1.0 ) ) ** 2.0
    #    print('sigsqr values',sig_0, k, xi_max, xi_min)
        return sigsqr  
    
    
    
    
    
    

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
        K(int): The number of Gaussians to use in the mixture model for the distribution of xi.
        nchains(int): The number of Monte Carlo Markov Chains to instantiate.
        parallelize(bool): Use a separate thread for each chain.  Only makes sense for nchains > 1.
        seed(int): Random seed.  If `None`, then get seed from np.random.randint().

    Attributes:
        nchains(int): The number of instantiated MCMCs.
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
    def __init__(self, xArr, yArr, xsigArr=None, ysigArr=None, xycovArr=None, delta=None, K=3,
                 nchains=4, parallelize=True, seed=None, nGMM_err=1, pi_err=None):
        self.nchains = nchains
        # lester
#        parallelize = False
        self.parallelize = parallelize

        if seed is None:
            seed = np.random.randint(2**32-1)

        if pi_err is None:
            pi_err = np.ones(len(xArr)) # lester changed x to xArr
        if self.parallelize:
            # Will place 1 chain in 1 thread.
            from multiprocessing import Process, Pipe
            # Create a pipe for each thread.
            self.pipes = []
            slave_pipes = []
            for i in range(self.nchains):
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
                            'delta':delta,
                            'K':K,
                            'nchains':self.nchains,
                            'nGMM_err':nGMM_err,
                            'pi_err':pi_err}
            for i, p in enumerate(self.pipes):
                init_kwargs = init_kwargs0.copy()
                init_kwargs['rng'] = np.random.RandomState(seed+i)
                p.send({'task':'init',
                        'init_args':init_kwargs})
        else:
            self._chains = []
            for i in range(self.nchains):
                self._chains.append(Chain(xArr, yArr, xsigArr, ysigArr, xycovArr, delta, K, nGMM_err, pi_err, self.nchains))
                self._chains[-1].initial_guess()

    def _get_psi(self):
        if self.parallelize:
            print('get_psi')
            for p in self.pipes:
                p.send({'task':'fetch',
                        'key':'chain'})
            chains = [p.recv() for p in self.pipes]
            self.pipes[0].send({'task':'fetch',
                                'key':'ichain'})
            
            ndraw = int(self.pipes[0].recv()/2)
            
        else:
            chains = [c.chain for c in self._chains]
            ndraw = int(self._chains[0].ichain/2)
        psi = np.empty((ndraw, self.nchains, 6), dtype=float)
        psi[:, :, 0] = np.vstack([c['alpha'][0:ndraw] for c in chains]).T
        beta = np.vstack([c['beta'][0:ndraw] for c in chains]).T
        psi[:, :, 1] = beta
        sig_0 = np.vstack([c['sig_0'][0:ndraw] for c in chains]).T # lester
#        psi[:, :, 2] = np.log(sigsqr)
        ximean = np.vstack([np.sum(c['pi'][0:ndraw] * c['mu'][0:ndraw], axis=1)
                            for c in chains]).T
        psi[:, :, 3] = ximean
        xivar = np.vstack([np.sum(c['pi'][0:ndraw] * (c['tausqr'][0:ndraw] + c['mu'][0:ndraw]**2),
                                  axis=1)
                           for c in chains]).T - ximean**2
        psi[:, :, 4] = xivar
        psi[:, :, 5] = np.arctanh(beta * np.sqrt(xivar / (beta**2 * xivar + sig_0))) # lester
        return psi

    def _get_Rhat(self):
        psi = self._get_psi()
        ndraw = psi.shape[0]
        psibarj = np.sum(psi, axis=0)/ndraw
        psibar = np.mean(psibarj, axis=0)
        sjsqr = np.sum((psi-psibarj)**2 / (ndraw-1.0), axis=(0, 1))
        where_are_NaNs = np.isnan(sjsqr) # lester added 
        sjsqr[where_are_NaNs] = 0 # lester added 
        Bvar = ndraw / (self.nchains-1.01) * np.sum((psibarj-psibar)**2, axis=0) # lester changed 1.0 to 1.01
        where_are_NaNs = np.isnan(Bvar) # lester added 
        Bvar[where_are_NaNs] = 0 # lester added 
        Wvar = sjsqr / self.nchains
        where_are_infs = np.isinf(Wvar) # lester added 
        Wvar[where_are_infs] = 1e3 # lester added         
        varplus = (1.0 - 1.0 / ndraw) * Wvar + Bvar / ndraw
        Rhat = np.sqrt(varplus / Wvar)
        return Rhat

    def _initialize_chains(self, miniter, n_burn):
        if self.parallelize:
            for p in self.pipes:
                p.send({'task':'init_chain',
                        'miniter':miniter,'n_burn':n_burn})
        else:
            for c in self._chains:
                c.initialize_chain(miniter,n_burn)

    def _step(self, niter):
        if self.parallelize:
            for p in self.pipes:
                p.send({'task':'step',
                        'niter':niter})
        else:
            for c in self._chains:
                c.step(niter)

    def _extend(self, niter):
        if self.parallelize:
            for p in self.pipes:
                p.send({'task':'extend',
                        'niter':niter})
        else:
            for c in self._chains:
                c.extend(niter)

    def _build_chain(self, ikeep):
        if self.parallelize:
            for p in self.pipes:
                p.send({'task':'fetch',
                        'key':'chain'})
            self.chain = np.hstack([p.recv()[ikeep:] for p in self.pipes])
        else:
            self.chain = np.hstack([c.chain[ikeep:] for c in self._chains])

    def run_mcmc(self, miniter=5000, maxiter=100000, silent=False, n_burn=1000):
        """ Run the Markov Chain Monte Carlo for the LinMix object.

        Bayesian inference is employed, and a Markov chain containing random draws from the
        posterior is developed.  Convergence of the MCMC to the posterior is monitored using the
        potential scale reduction factor (RHAT, Gelman et al. 2004). In general, when RHAT < 1.1
        then approximate convergence is reached.  After convergence is reached, the second halves
        of all chains are concatenated and stored in the `.chain` attribute as a numpy recarray.

        Args:
            miniter(int): The minimum number of iterations to use.
            maxiter(int): The maximum number of iterations to use.
            silent(bool): If true, then suppress updates during sampling.
        """
        checkiter = 100
        # lester
        checkiter = miniter
        
        self._initialize_chains(miniter,n_burn)
        self.n_burn = n_burn
        for i in range(0, miniter, checkiter):
            self._step(checkiter)
          
#            Rhat = self._get_Rhat()
#            
#
#            if not silent:
#                print()
#                print("Iteration: ", i+checkiter)
#                print ("Rhat values for alpha, beta, log(sigma^2)"
#                       ", mean(xi), log(var(xi)), atanh(corr(xi, eta)):")
#                print(Rhat)
#                print()
#
#        i += checkiter
#        while not np.all(Rhat < 1.1) and (i < maxiter):
#            self._extend(checkiter)
#            self._step(checkiter)
#
#            Rhat = self._get_Rhat()
#            if not silent:
#                print()
#                print("Iteration: ", i+checkiter)
#                print ("Rhat values for alpha, beta, log(sigma^2)"
#                       ", mean(xi), log(var(xi)), atanh(corr(xi, eta)):")
#                print(Rhat)
#                print()
#                i += checkiter

        # Throw away first half of each chain
#        self._build_chain(int(i/2))
#        self._build_chain(self.n_burn)

        self._build_chain(0)
        # Clean up threads
        if self.parallelize:
            for p in self.pipes:
                p.send({'task':'kill'})



        


# =============================================================================
# INPUTS
# =============================================================================
#
##sbf = './limix_inputs/' # subfolder
#sbf = './limix_inputs_106_x47_mass5p0/' # subfolder
#
#
#pi_err      = np.load(sbf+'pi_err.npy')         # 3x probability of each posterior gaussian
#GMMx        = np.load(sbf+'GMMx.npy')           # 3x posterior means per mass
#GMMy        = np.load(sbf+'GMMy.npy')           # 3x posterior means per sfr
#GMMxsig     = np.load(sbf+'GMMxsig.npy')        # 3x posterior sigmas per mass
#GMMysig     = np.load(sbf+'GMMysig.npy')        # 3x posterior sigmas per sfr
#GMMxycov    = np.load(sbf+'GMMxycov.npy')       # 3x posterior covar per mass-sfr pair
#
#nK          = np.load(sbf+'nK.npy')             # 3 #gaussians modelling xi
#nGauss      = np.load(sbf+'nGauss.npy')         # 3 #gaussians modelling BEAGLE posterior
#
#nChains     = np.load(sbf+'nChains.npy')        # 2
#minIter     = np.load(sbf+'minIter.npy')        # 3000
#maxIter     = np.load(sbf+'maxIter.npy')        # 3000
#
##print(GMMx)
#
#nChains     = 1
#minIter     = 30
#maxIter     = 30
#
#lm = LinMix(GMMx, GMMy, GMMxsig, GMMysig, xycovArr=GMMxycov, K=nK, nGMM_err=nGauss, pi_err=pi_err, nchains=nChains)
#lm.run_mcmc(miniter=minIter, maxiter=maxIter)
#
## https://github.com/jmeyers314/linmix
#print(' ')
#print("{}, {}".format(lm.chain['alpha'].mean(), lm.chain['alpha'].std()))
#print("{}, {}".format(lm.chain['beta'].mean(), lm.chain['beta'].std()))



#print(lm)

#
#lm_chain_str = str(round(np.min(GMMx), 2)).replace(".","p")
#np.save(lm_chain_str + '_lm_chain', lm.chain)
#
#

# =============================================================================
# PLOTS
# =============================================================================
'''
import matplotlib.pyplot as plt

plt.title('alpha')
plt.plot(lm.chain['alpha'])
plt.show()

plt.title('beta')
plt.plot(lm.chain['beta'])
plt.show()

plt.title('sig0')
plt.plot(lm.chain['sig_0'])
plt.show()

plt.title('k')
plt.plot(lm.chain['k'])
plt.show()


# need to consider plotting similar for xi and eta for individual objects

plt.figure(figsize=(15, 10))
for i in range(len(lm.chain['xi'][0,21:33])):
    plt.plot(lm.chain['xi'][:,i])
plt.show()

'''
