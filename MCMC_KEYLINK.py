''' MCMC FOR KEYLINK (FIRST VERSION) CALIBRATION
    Last write 20/11/2019
# Python version of Markov Chain Monte Carlo (MCMC) Line 
# David Cameron 11/11/2008 dcam@ceh.ac.uk
# modified for KEYLINK model 3/2018 Gaby Deckmyn
'''
import numpy as nm
from numpy import matlib as ml
from numpy import random as ra
import scipy
from scipy import stats
import Keylinkbayesian as bayes

def export_pools(filename, array):
    """Save an array in a readable format, compatible with R"""
    nm.savetxt(filename + '.txt', array)

chainLength    = 10                                 

caldatafile=open('KL_calibData.txt')


titls=caldatafile.readline() # first line are titles
numdata=int(caldatafile.readline()) # second line is number of data to read
data=nm.zeros((numdata,4))

for i in range (0,numdata):
    zap=caldatafile.readline().split()    
    data[i,0]=float(zap[0])
    data[i,1]=float(zap[1])
    data[i,2]=float(zap[2])
    data[i,3]=float(zap[3])
 
Minima = [0, 0,0,0,0,0,0,0,0] ; pMinima=nm.reshape(Minima,[9]); Maxima = [10,10,6,6,6,0.5,0.5,0.5,0.5]; pMaxima=nm.reshape(Maxima,[9]) # max and min for prior

vcovProposal   = ml.diag( (0.005*(pMaxima-pMinima))**2 )  # variance of the proposal                        

Values        = [1.06,1.14,1.04,1.540,1.188,0.038,0.072,0.028,0.14] ; pValues=nm.reshape(Values,[9])  # starting values of the two model parameters (gmax[0] and gmax[1])  
priorChain         = nm.zeros( [chainLength, len(pValues)],'d' ) ; priorChain[0,:] = pValues # start the chain 
posteriorChain         = nm.zeros( [chainLength, len(pValues)],'d' ) ; posteriorChain[0,:] = pValues # start the chain 

logPrior0 = nm.sum( nm.log( stats.uniform.pdf( pValues, pMinima, pMaxima ) ) )                 

# this is the call to the model with the initial parameter values
y=bayes.KeylinkModel(pValues) #returns biomsses of all pools on all days in y


# calculate first log likelihood (log of summed differences betwen measured and simulated and /error)
ii=0
logLi=[]
logLseries=[]
for i in data[:,0]:  #choose the data (value and error)from the right days 
    logLi.append(-0.5*((y[int(i-1),data[ii,1]]-data[ii,2])/data[ii,3])**2 - nm.log(data[ii,3])) #without Pi, because it will be cancelled with the one from LogL1 (to simplify)
    ii=ii+1 
logL0 = nm.sum( logLi )/ numdata                                
psetMAP = pValues ; logMAP = logPrior0 + logL0 #MAP is max fit point, save parameter and likelyhood of best run
logLseries.append(logL0)

sim=nm.zeros([chainLength, numdata], 'd')

# loop over chain length
for c in range(1,chainLength):

    candidatepValues     = ra.multivariate_normal(pValues, vcovProposal )  # candidate values from a multivariate normal
    
    for j in range(5): #Reflection algorithm
        refMin = min(0,candidatepValues[j]-pMinima[j]) #reflection from minimum
        refMax = max(0,candidatepValues[j]-pMaxima[j]) #reflection from maximum
        candidatepValues[j]=candidatepValues[j]-2*refMin-2*refMax
        
    Prior1 = nm.product( stats.uniform.pdf( candidatepValues, pMinima, pMaxima ) ) # uniform prior     
    
    if Prior1 > 0: #if the parameter you want to try is in the range between min and max
       # this is the model  run with the 'candidate parameter' 
       y=bayes.KeylinkModel(candidatepValues)
       # calculate the log likelihood
       ii=0
       logLi=[]
       simulation=[]
       for i in data[:,0]:
           simvalue=y[int(i-1),int(data[ii,1])]
           simulation.append(simvalue)
           logLi.append(-0.5*((simvalue-data[ii,2])/data[ii,3])**2 - nm.log(data[ii,3])) #without Pi, because it will be cancelled with the one from LogL0
           ii=ii+1
       logL1 = (nm.sum( logLi ) / numdata)
       sim[c,:]=simulation
       # form the ratio of this step to previous and accept/reject from this (log a/b = log a-log b)       
       logalpha          = (nm.log(Prior1)+logL1) - (logPrior0+logL0)  
       lograndom= nm.log(ra.random())
       logLseries.append(logL1)
       if lograndom < logalpha: 
          pValues = candidatepValues ; logPrior0 = nm.log(Prior1) ; logL0 = logL1 # if accepted move to this point
          posteriorChain[c,:]=pValues # add step to the chain	  
          if (logPrior0+logL0)>logMAP:
             logMAP = (logPrior0+logL0) ; psetMAP = pValues # update most likely parameter set
    
    priorChain[c,:]=candidatepValues
    export_pools('keylinkBayesianPrior', priorChain)
    export_pools('keylinkBayesianPosterior', posteriorChain)
    export_pools('keylinkBayesianOutputLogli', logLseries)
    
    
mp                = ml.mean(posteriorChain) ; print('mp',mp)
pCovMatrix        = ml.cov(posteriorChain)  ; print('pCovMatrix',pCovMatrix)
sp                = nm.sqrt(nm.transpose(ml.diag(pCovMatrix))) ; print('sp',sp)
pCorrMatrix       = scipy.corrcoef(posteriorChain) ; print('pCorrMatrix', pCorrMatrix)
print('psetMAP', psetMAP)
 
