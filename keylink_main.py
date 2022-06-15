'''
NEW main for calling KEYLINK model as a function, either single call, bayesian optimisation,
 using posterior (parameter set) or using parameter distributions. 
'''


import numpy as nm
from numpy import matlib as ml
from numpy import random as ra
import scipy
from scipy import stats
import keylink_core as core
import pandas as pd

"runmode can be single, bayesian, posterior, testing (or distribution, not yet implimented)"
runmode='bootstrap'

if (runmode=='single'):
#read the gmax values
 param = core.import_pools('KL_FaunalParams')   
 gmax=param[0,:] 
 Cpools, PWt, PVt, pores=core.KeylinkModel(gmax) 
 core.export_pools('keylinkoutput', Cpools)
 print(Cpools)
 core.show_plot(Cpools, PWt, PVt)  

 
# pValues        =  [1.06,1.14,1.04,1.540,1.188,0.038,0.072,0.028,0.14]    

if(runmode=='bayesian'):
 ''' MCMC FOR KEYLINK (FIRST VERSION) CALIBRATION
    Last write 3/6/2021
 # Python version of Markov Chain Monte Carlo (MCMC) Line 
 # David Cameron 11/11/2008 dcam@ceh.ac.uk
 # modified for KEYLINK model 3/2018 Gaby Deckmyn
 '''
 chainLength    = 10000                          
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

 param = core.import_pools('KL_FaunalParams')   
 gmax=param[0,:] #Original values in the parameterset are the starting values
 #Values        = [1.06,1.14,1.04,1.540,1.188,0.038,0.072,0.028,0.14] ; pValues=nm.reshape(Values,[9])  # starting values of the 9 model being optimised (gmax)  
 priorChain         = nm.zeros( [chainLength, len(gmax)],'d' ) ; priorChain[0,:] = gmax # start the chain 
 posteriorChain         = nm.zeros( [chainLength, len(gmax)],'d' ) ; posteriorChain[0,:] = gmax # start the chain 

 logPrior0 = nm.sum( nm.log( stats.uniform.pdf( gmax, pMinima, pMaxima ) ) )                 

 # this is the call to the model with the initial parameter values
 y,PWt, PVt, pores=core.KeylinkModel(gmax) #returns biomsses of all pools on all days in y


 # calculate first log likelihood (log of summed differences betwen measured and simulated and /error)
 ii=0
 logLi=[]
 logLseries=[]
 for i in data[:,0]:  #choose the data (value and error)from the right days 
    logLi.append(-0.5*((y[int(i-1),int(data[ii,1])]-data[ii,2])/data[ii,3])**2 - nm.log(data[ii,3])) #without Pi, because it will be cancelled with the one from LogL1 (to simplify)
    ii=ii+1 
 logL0 = nm.sum( logLi )/ numdata                                
 psetMAP = gmax ; logMAP = logPrior0 + logL0 #MAP is max fit point, save parameter and likelyhood of best run
 logLseries.append(logL0)

 sim=nm.zeros([chainLength, numdata], 'd')

 # loop over chain length
 for c in range(1,chainLength):

    candidategmax     = ra.multivariate_normal(gmax, vcovProposal )  # candidate values from a multivariate normal
    
    for j in range(5): #Reflection algorithm
        refMin = min(0,candidategmax[j]-pMinima[j]) #reflection from minimum
        refMax = max(0,candidategmax[j]-pMaxima[j]) #reflection from maximum
        candidategmax[j]=candidategmax[j]-2*refMin-2*refMax
        
    Prior1 = nm.product( stats.uniform.pdf( candidategmax, pMinima, pMaxima ) ) # uniform prior     
    
    if Prior1 > 0: #if the parameter you want to try is in the range between min and max
       # this is the model  run with the 'candidate parameter' 
       y, PWt, PVt, pores=core.KeylinkModel(candidategmax)
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
          gmax = candidategmax ; logPrior0 = nm.log(Prior1) ; logL0 = logL1 # if accepted move to this point
          posteriorChain[c,:]=gmax # add step to the chain	  
          if (logPrior0+logL0)>logMAP:
             logMAP = (logPrior0+logL0) ; psetMAP = gmax # update most likely parameter set
    
    priorChain[c,:]=candidategmax
    core.export_pools('keylinkBayesianPrior', priorChain)
    core.export_pools('keylinkBayesianPosterior', posteriorChain)
    core.export_pools('keylinkBayesianOutputLogli', logLseries)
    
    
 mp                = ml.mean(posteriorChain) ; print('mp',mp)
 pCovMatrix        = ml.cov(posteriorChain)  ; print('pCovMatrix',pCovMatrix)
 #sp                = nm.sqrt(nm.transpose(ml.diag(pCovMatrix))) ; print('sp',sp)
 pCorrMatrix       = scipy.corrcoef(posteriorChain) ; print('pCorrMatrix', pCorrMatrix)
 print('psetMAP', psetMAP)
 
 #call best fit again and show graphs
 Cpools, PWt, PVt, pores=core.KeylinkModel(psetMAP) 
 core.export_pools('keylinkoutput', Cpools)
 core.show_plot(Cpools, PWt, PVt)  
 
if runmode == 'posterior':   #the model is run for all the given combinations of gmax values, 
#read line by line, then average and STDEV are calculated
    
 pools = 22 #number of pools
 tStop=3653 #selects the amount of days
 Nvec=10 #number of parameter vectors to test from the Bayesian calibration outputn ormal = 100
 statistics_A = nm.zeros([Nvec, 2*pools], 'd') #matrix to record the averages and standard deviations of population biomass and fluxes
 statistics_B = nm.zeros([Nvec, 2*pools], 'd') #matrix to record the minimum and maximum values of population biomass and fluxes
 soilm = nm.zeros([Nvec, 6], 'd') #matrix to record the average volumens of each pore size class and the average soil water content
 nsim = 0 #number of simulations currently done
 raw = nm.zeros([Nvec*int(tStop), pools], 'd') #matrix for all raw data output from daily simulation of pools
 
 bayesGmax = core.import_pools('KL_gmax') #gmax from (Bayesian) calibration 
 for u in range(0,Nvec):
    gmax=bayesGmax[u,] #vector of gmax randomly selected from Bayesian calibration of gmax
    Cpools, PWt, PVt, pores=core.KeylinkModel(gmax)   
    raw[nsim,:] = Cpools[1, :] #needs to be cheked
    soilm[u,5] = (sum(PWt)/tStop) #average value of soil water content within the simulation period
    for r in range(5): #average values of volumes in each pore size class within the simulation period
        soilm[u,r]=(sum(pores[:,r])/tStop)
    statA=[] #vector for statistic values (each pair of columns have mean and sd for each population)
    statB=[] #vector for minimum and maximum values of biomass in each pool
    for q in range(pools):
        mean=(sum(Cpools[:,q])/len(Cpools[:,q])) #mean biomass on each population
        statA.append(mean)
        sqd=[] #squared differences (for standard deviation formula)
        for w in Cpools[:,q]:
            z=(w-mean)**2
            sqd.append(z)
        sd=(sum(sqd)/(core.n-1))**0.5 #standard deviation of biomass on each population
        statA.append(sd)
        statB.append(min(Cpools[:,q]))
        statB.append(max(Cpools[:,q]))
    statistics_A[u,:]=statA
    statistics_B[u,:]=statB
    core.export_pools('KL_statistics_mean_sd', statistics_A)
    core.export_pools('KL_statistics_min_max', statistics_B)
    core.export_pools('KL_soil_matrix', soilm)
    core.export_pools('KL_raw_data', raw)
  
    #Does a specified amount of runs, to be used in conjunction with testtempmodel being on in core.
    #it will then store all the Carbon pool data into one big csv
if runmode ==  'bootstrap':      
    param = core.import_pools('KL_FaunalParams')   #read the gmax values
    gmax=param[0,:] 
    sampleSize = 100 #amount of runs
    groupAmount = 14
    tStop=3653 #selects the amount of days
    seasonLength = int(nm.floor(tStop/(4*tStop/365.25)))
    seasonAmount = int(nm.floor(4*tStop/365.25))
 
    for i in range(0,sampleSize):
        Cpools, PWt, PVt, temps=core.KeylinkModel(gmax) #Mimic the simple approach, grab the data the model generates
        CpoolSeasonalData = [[0] * groupAmount for i in range(seasonAmount)]
        CpoolSeasonalTemps = [[0] * groupAmount for i in range(seasonAmount)]
        
        for s in range(1,seasonAmount+1): #Grab 40 data points from day 91 onwards in 91 day intervals until 10 years have passed and collect it into a new array
            CpoolSeasonalData[s-1] = Cpools[s*seasonLength][0:groupAmount]
            CpoolSeasonalTemps[s-1] = temps[s*seasonLength]            
        currentCpoolData=core.mk_dataframe(CpoolSeasonalData, CpoolSeasonalTemps, i, tStop, groupAmount) #store into a compatible dataframe
        
        if i == 0:
            currentCpoolData.to_csv('./data/testData.csv') #export data to csv
        else:
            fullCpoolData=pd.concat([core.importFromCSV('./data/testData.csv'),currentCpoolData])
            fullCpoolData.to_csv('./data/testData.csv') #export data to csv

        
        
        
#        if i > 0:
#            fullCpoolData=core.append_df(fullCpoolData, CpoolSeasonalData) #append last dataframe to the full dataframe
        
        