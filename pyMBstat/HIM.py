class HIM: 
    """
    Class for creating the Quantum Ground State of the one dimensional 
    Bosonic Harmonic Interaction Model (HIM)
    see J. Math. Phys. 26(12), 3105, (1985).
    
    Parameters
    -----------
    N - Number of Bosons.
    omega - Frequency of the single particle harmonic confinement potential.
    lambda0 - Two body interaction parameter (= gamma^2 / 2 in the paper).
    """
    def __init__(self, N = 10, omega = 1, lambda0 = 0.1):
        self.N = N
        self.dimensions = 1
        self.omega = omega
        self.lambda0 = lambda0
        self.deltaN = np.sqrt(np.power(omega,2)+N*2*lambda0)
        self.KN = np.power(self.deltaN/np.pi,(self.dimensions/4)*(self.N-1))\
                * np.power(self.omega/np.pi,(self.dimensions/4))
        self.A = (1/(4*self.N))*((self.N-1)*self.deltaN + self.omega)
        self.B = (1/(4*self.N))*(self.omega - self.deltaN)
        # parameters for Rho1
        self.a1 = ((self.N-1.0)*(np.power(self.omega,2) + np.power(self.deltaN,2))
                + 2.0 * (np.power(self.N,2) - self.N + 1)*self.omega*self.deltaN)\
                / ((self.N-1.0)*self.omega+self.deltaN) / (4.0*self.N)
        self.a2 = (self.N-1)*np.power(self.omega-self.deltaN,2)\
                /((self.N-1.0)*self.omega+self.deltaN) / (2.0*self.N)
        # parameters for Rho2
        self.b1 = (1 / (4 * self.N)) * ((self.N - 2) * self.omega**2 + (3 * self.N - 2) * self.deltaN**2
                + 2 * (self.N**2 - 2 * self.N + 2) * self.omega * self.deltaN)\
                / ((self.N - 2) * self.omega + 2 * self.deltaN)
        self.b2 = (1 / (2 * self.N)) * ((self.N - 2) * self.omega**2 - (self.N + 2) * self.deltaN**2
                + 4 * self.omega * self.deltaN)\
                / ((self.N - 2) * self.omega + 2 * self.deltaN)
        self.b3 = (1 / (2 * self.N)) * ((self.N - 2) * (self.omega - self.deltaN)**2)\
                / ((self.N - 2) * self.omega + 2 * self.deltaN)
    def conditonalP(self,x,R,m):
        Pcond = np.sqrt(1.0/(2 * np.pi * self.var(m)))\
              * np.exp(-(x - self.mu(m)*R)**2/(2*self.var(m)))
        return Pcond
    def mu(self,m):
        mu = (self.deltaN - self.omega) / ((m - 1) * self.deltaN + (1 - m + self.N) * self.omega)
        return mu
    def var(self,m):
        var = m * (self.deltaN - self.omega) + self.N * self.omega
        var = var / (2 * self.deltaN * ((m - 1) * self.deltaN + (1 - m + self.N) * self.omega))
        return var
    def sigma(self,m):
        sigma = np.sqrt(self.var(m))    
        return sigma
    def psiN(self,x):
        """ 
        The many-body wavefunction Psi(x1,x2,...,xN) 
        ---------------------------------------------
        Input - N-dimensional postion vector 
        Output - value of wavefunction 
        
        """
        if (self.N != len(x)):
            return('Wrong Input. The length of x must be {}'.format(self.N))
        wavefunction = self.KN * np.exp(- 2 * self.A * np.dot(x,x) 
                                   - 4 * self.B * (np.sum(np.sum(np.outer(x, x) - np.diag(x*x))))/2)
        return wavefunction
    def rho1(self, x, xp = 0.0, diag = True, normN = False): 
        """
        The one-body reduced density matrix
        -----------------------------------
        rho (x, xp)
        
        Input - two position values or mesh grid matrices
        Output - the density matrix over the input
        
        diag - flag to return only the boson density, i.e., x = xp
        normN - flag to normalize the boson density to N. 
        """
        if diag == True: # will return only the diagonal of the density matrix as a vector
            rho = np.exp(-self.a1*(np.power(x,2)+np.power(x,2))+self.a2*x*x)
        else: # return the full first order density matrix
            rho = np.exp(-self.a1*(np.power(x,2)+np.power(xp,2))+self.a2*x*xp)          
        rho = rho * np.power(self.deltaN * self.N * self.omega
                                              / np.pi / ((self.N-1)*self.omega+self.deltaN),self.dimensions/2)           
        if normN == True:
            rho = rho * self.N
        return rho
    def rho2(self, x1, x2, x1p = 0.0, x2p = 0.0, diag = True, normN = False): 
        """
        The two-body reduced density matrix
        -----------------------------------
        rho2 (x1, x2, x1p, x2p)
        
        Input - four position values
        Output - the second order reduced density over the input
        
        diag - flag to return only the diagonal of the two-body density, i.e., x1 = x1p, x2 = x2p
        normN - flag to normalize the second order density to N(N-1)/2. 
        """
        if diag == True: # will return only the diagonal of the density matrix as a vector
            rho = np.exp(-2 * self.b1*(np.power(x1,2) + np.power(x2,2))
                        - 2 * self.b2 * (x1 * x2) + self.b3 * (np.power(x1,2) + np.power(x2,2) +2 * x1 * x2))
        else: # return the full first order density matrix
            rho = 0
        rho = rho * np.power(self.deltaN**2 * self.N * self.omega
                                              / np.pi**2 / ((self.N-2)*self.omega+2*self.deltaN),self.dimensions/2) 
        if normN == True:
            rho = rho * self.N * (self.N - 1) / 2
        return rho
    def rhoN(self, x, xp = 0.0, diag = True): 
        """
        The N-body density matrix
        -------------------------
        rhoN(x1,...,xM,x1p,...,xNp)
        
        Input - two N-dimensional vectors 
        Output - value of the densiy
        
        diag - flag to output the diagonal element only, i.e., x1 = x1p, ..., xN = xNp
        """
        if diag == True:
            rhoN = np.power(self.psiN(x), 2)
        else:
            rhoN = np.power(self.KN, 2) * np.exp(-2 * self.A * (np.power(x,2)+np.power(xp,2))
                                                 - 4 * self.B * ((np.sum(np.sum(np.outer(x, x) - np.diag(x*x))))/2)
                                                 + ((np.sum(np.sum(np.outer(xp, xp) - np.diag(xp*xp))))/2))
        return rhoN 
    def singleshot(self, numshots = 1, xmin = -10.0, xmax = +10.0):
        """
        Single-shot implementation directly from the N-body density
        using rejection sampling. The random samples are taken from 
        the multivariate normal distribution. For more than a few particles
        this doesn't work as the sampling efficiency drops to zero !!! 
        """
        #covmat = np.eye(self.N)/(self.deltaN)
        covmat = np.eye(self.N)
        rv = multivariate_normal(mean=np.zeros(self.N),
                                 cov=covmat, allow_singular=False)
        M = np.power(self.KN,2)
        sample_size =  (max(1,int(1.0/M))*numshots,self.N)
        #coordinate_sample = (xmax - xmin) * np.random.random_sample(sample_size) + xmin
        coordinate_sample = rv.rvs(sample_size[0])
        p = np.random.random_sample(sample_size[0])
        f = np.zeros_like(p)
        for i in range(sample_size[0]):
        #    f[i] = self.rhoN(coordinate_sample[i,:])
            f[i] = self.rhoN(coordinate_sample[i]) / rv.pdf(coordinate_sample[i])
        M = np.power(self.KN,2)*(np.power(2*np.pi,self.N/2)*np.sqrt(np.linalg.det(covmat)))
        idx = (p <= f/M) # rejection sampling
        sshot = coordinate_sample[idx]
        print('Estimated Sampling Efficency={:3.2f}%'.format(100*sshot.shape[0]/coordinate_sample.shape[0]))
        return sshot
    def check_condP(self):
        """
        Chech implementation of the conditional probability by comparing P2(x2|x1) with rho2(x2,x1)/rho1(x1)
        """
        x = np.linspace(-4,4,128,endpoint=True)
        [X1,X2] = np.meshgrid(x,x)
        RHO2 = self.rho2(X1,X2)
        RHO1 = self.rho1(X1)
        PCOND2 = self.conditonalP(X2,X1,2)
        fig, axs = plt.subplots(1,2,figsize = (12,6))
        s = 'Verifying the implementation of the conditional probability\n'
        s = s + 'Maximal difference betweeen the two implementations is {:3.2e}'
        fig.suptitle(s
                  .format(np.max(PCOND2-RHO2/RHO1)),fontsize = 16)
        axs[0].pcolor(X1,X2,RHO2/RHO1);
        axs[1].pcolor(X1,X2,PCOND2);
    def rho1hist(self, numshots = 1, xmin = -5.0, xmax = +5.0, random_dist = 'uniform'):
        """
        
        """
        if random_dist == 'uniform':
            M = np.power(self.KN,2)
            sample_size =  (max(1,int(1.0/M))*numshots)
            coordinate_sample = (xmax - xmin) * np.random.random_sample(sample_size) + xmin
            f = self.rho1(coordinate_sample)
        elif random_dist == 'normal':
            sigma = np.sqrt(1.0/self.deltaN)
            sigma = 5
            M = np.power(self.KN,2)*np.sqrt(2*np.pi)*sigma
            print(M)
            sample_size =  (max(1,int(1.0/M))*numshots)
            coordinate_sample = np.random.normal(0, sigma, sample_size)
            f = self.rho1(coordinate_sample)/self.gaussian(0, sigma, coordinate_sample)
        p = np.random.random_sample(sample_size)
        idx = (p <= f/M) # rejection sampling
        sshot = coordinate_sample[idx]
        print('Estimated Sampling Efficency={:3.2f}%'.format(100*sshot.shape[0]/coordinate_sample.shape[0]))
#        return sshot
        plt.hist(sshot,normed = 1, bins = 32, alpha = 0.5, histtype='stepfilled', linewidth = 3, color = 'b');
        x = np.linspace(-3,3,512,endpoint=True)
        plt.plot(x,  self.rho1(x), linewidth = 3, color = 'k');
    def rho2hist(self, numshots = 1, xmin = -5.0, xmax = +5.0, ymin = -5.0, ymax = +5.0):
        """
        
        """

        #x_sample = (xmax - xmin) * np.random.random_sample(sample_size) + xmin
        #y_sample = (ymax - ymin) * np.random.random_sample(sample_size) + ymin
        #x_sample = np.random.normal(0, 1, sample_size)
        #y_sample = np.random.normal(0, 1, sample_size)
        meanmat=np.zeros(2)
        covmat = np.eye(2)
        M = self.rho2(0,0) * (2*np.pi*np.sqrt(np.linalg.det(covmat)))
        sample_size =  (max(1,int(1.0/M))*numshots,2)
        #x_sample, y_sample = np.random.multivariate_normal(mean, covmat, sample_size).T
        rv = multivariate_normal(mean=meanmat,
                                 cov=covmat, allow_singular=False)
        coordinate_sample = rv.rvs(sample_size[0])
        p = np.random.random_sample(sample_size[0])
        f = np.zeros_like(p)
        for i in range(sample_size[0]):
            f[i] = self.rho2(coordinate_sample[i][0],coordinate_sample[i][1]) / rv.pdf(coordinate_sample[i])
        #f = self.rho2(x_sample,y_sample)/(self.gaussian(0.0, 1.0, x_sample)*self.gaussian(0, 1, y_sample))
        #f = self.rho2(x_sample,y_sample)
        idx = (p <= f/M) # rejection sampling
        #sshot = np.asarray([x_sample[idx],y_sample[idx]])
        print(coordinate_sample.shape)
        sshot = coordinate_sample[idx]
        #print('Estimated Sampling Efficency={:3.2f}%'.format(100*sshot.shape[0]/x_sample.shape[0]))
        print('Estimated Sampling Efficency={:3.2f}%'.format(100*sshot.shape[0]/coordinate_sample.shape[0]))
        return sshot
    def gaussian(self, mu, sigma, x):
        gauss = (1.0/(sigma * np.sqrt(2 * np.pi))) * np.exp( - (x - mu)**2 / (2 * sigma**2)) 
        return gauss
    def singleshot_chainrule(self, numshots = 1):
        """
        Single-shot from the N-body density expanded as a product of conditional probabilities 
        and the one body density. 
        """
        sshot = np.zeros((numshots,self.N))
        Psshot = np.zeros(numshots)
        coordinate_sample = np.zeros(self.N)
        for i in range(numshots):
            mu = 0.0
            sigma = np.sqrt (1 / (2 * (2 * self.a1 - self.a2)))
            coordinate_sample[:] = 0
            coordinate_sample[0] = np.random.normal(mu, sigma)
            Psshot[i] = self.gaussian(mu,sigma,coordinate_sample[0])
            for m in range(2,system.N+1):
                R = np.sum(coordinate_sample)
                mu = self.mu(m)*R
                sigma = self.sigma(m)
                coordinate_sample[m-1] = np.random.normal(mu, sigma)
                Psshot[i] = Psshot[i] * self.gaussian(mu,sigma,coordinate_sample[m-1])
            sshot[i,:] = coordinate_sample
        return sshot, Psshot
