import sys
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.style as style
style.use('tableau-colorblind10')
from scipy.stats import beta, gamma, sem

print("#--------------------------------------------------------------------------------------#")
print("Start of program")
print("#--------------------------------------------------------------------------------------#")

preFix = ""
des = ""
inputPath=""+preFix
outputPath=""+preFix
netcdfPath = inputPath+"mostLikelyVals"+des+".nc"

# Read data for different sliding laws
def readData(iInFile,inFileOptions):
    print(inFileOptions[iInFile]+" was selected")
    if inFileOptions[iInFile] in ["weertmanN0","weertmanNpi"]:
        NetCDFfile = inputPath+"weertmanGridSearch"+des+".nc"
    elif inFileOptions[iInFile] == "coulomb":
        NetCDFfile = inputPath+"coulombGridSearch"+des+".nc"
    elif inFileOptions[iInFile] == "budd":
        NetCDFfile = inputPath+"buddGridSearch"+des+".nc"
    elif inFileOptions[iInFile] == "tsaibuddCFC":
        NetCDFfile = inputPath+"tsaiBuddGridSearchCFC"+des+".nc"
    elif inFileOptions[iInFile] == "tsaibuddBFP":
        NetCDFfile = inputPath+"tsaiBuddGridSearchBFP"+des+".nc"
    elif inFileOptions[iInFile] == "schoofCS":
        print("useCS = True: varying the Schoof friction coefficient")
        NetCDFfile = inputPath+"schoofGridSearchCS"+des+".nc"
    elif inFileOptions[iInFile] == "schoofCMAX":
        print("useCS = False: varying beta (local maximum up-slope of the bedrock) and, therefore, the maximum value C_max=tau_b/N=tan(beta)")
        NetCDFfile = inputPath+"schoofGridSearchCMAX"+des+".nc"
    elif inFileOptions[iInFile] == "zoet":
        print("useCFC = True: varying the Coulomb friction coefficient")
        NetCDFfile = inputPath+"zoetGridSearchCFC"+des+".nc"
    elif inFileOptions[iInFile] == "zoetuTnoN":
        print("useCFC = False: varying the transition speed uT (without dependence on N)")
        NetCDFfile = inputPath+"zoetGridSearchuTnoN"+des+".nc"

    print("Read NetCDF file: "+NetCDFfile)
    print("")

    # Open the NetCDF file
    d = Dataset(NetCDFfile, 'r')

    # Read the variables
    porosityList = d.variables['porosityList'][:]
    grainDiameterList = d.variables['grainDiameterList'][:]
    dof = d.variables['dof'][:]
    chiSq = d.variables['chiSq'][:]
    if inFileOptions[iInFile] not in ["weertmanN0","weertmanNpi"]:
        SLVList = d.variables['SLVList'][:]
    else:
        SLVList = [np.nan]
    if inFileOptions[iInFile] in ["weertmanN0"]:
        chiSq = d.variables['chiSqZero'][:] # CAREFUL!!! this now overwrites the regular chiSq
    if inFileOptions[iInFile] in ["schoofCS","schoofCMAX","zoet","zoetuTnoN"]:
        effectivePressureFailsWeight = d.variables['nEffectivePressureFailsWeight'][:]

    # Close the NetCDF file
    d.close()

    # to actually get chiSq, divide by number of samples
    chiSq = chiSq / dof 

    return chiSq, porosityList, grainDiameterList, SLVList

# aaa
# define prior functions
def priorUniform(parameterList):
    n = len(parameterList)
    uniformPrior = np.ones(n)/np.sum(np.ones(n))
    return uniformPrior

def priorTrapezoidal(parameterList, flatStart, flatEnd, lowestProbability):
    n = len(parameterList)
    centralStart = int(n * flatStart) # set to 0 if you don't want left slope
    centralEnd = int(n * flatEnd) # set to 0 if you don't want right slope
    trapezoidalPrior = np.ones(n)

    # Linearly decrease the probability towards the edges to lowestProbability (0 to 1) of the central value
    if centralStart > 0:
        leftSlope = (1 - lowestProbability) / centralStart
        for i in range(centralStart):
            trapezoidalPrior[i] = lowestProbability + leftSlope * i
    if centralEnd > 0:
        rightSlope = (1 - lowestProbability) / centralEnd
        for i in range(centralEnd):
            trapezoidalPrior[-(i+1)] = lowestProbability + rightSlope * i
    trapezoidalPrior /= np.sum(trapezoidalPrior)
    return trapezoidalPrior

def priorGamma(parameterList, shape, scale):
    # shape, scale = 2.0, 0.1  are shape and scale parameters for the gamma distribution
    gammaPrior = gamma.pdf(parameterList, shape, scale=scale)
    gammaPrior /= np.sum(gammaPrior)
    return gammaPrior

def priorGrainDist(parameterList):
    grainDistPrior = np.zeros(len(parameterList))
    # this expects units in mm
    sandSize = [1/16, 2]
    siltSize = [1/256, 1/16]
    claySize = [1/256]
    grainDistPrior[(parameterList<claySize[0])] = 20
    grainDistPrior[(parameterList>=siltSize[0]) & (parameterList<=siltSize[1])] = 50
    grainDistPrior[(parameterList>sandSize[0])] = 30
    grainDistPrior /= np.sum(grainDistPrior)
    return grainDistPrior

# bbb
# set priors for individual parameters
def defPriors(iInFile, inFileOptions, uniformPriors, porosityList, grainDiameterList, SLVList):
    if uniformPriors == True:
        print("uniform priors selected")
        porosityPrior = np.ones(len(porosityList))/np.sum(np.ones(len(porosityList)))
        grainDiameterPrior = np.ones(len(grainDiameterList))/np.sum(np.ones(len(grainDiameterList)))
        if inFileOptions[iInFile] not in ["rob","weertman","weertmanN0","weertmanNpi"]:
            SLVPrior = np.ones(len(SLVList))/np.sum(np.ones(len(SLVList)))
    else:
        print("custom prior selected")
        # porosity
        centralStart = 0
        centralEnd = (np.where(np.round(porosityList,3)==0.45)[0])/len(porosityList)
        lowestProbability = 0.3
        porosityPrior = priorTrapezoidal(porosityList, centralStart, centralEnd, lowestProbability)

        # grain diameter
        grainDiameterPrior = priorGrainDist(grainDiameterList*1e3)

        # all SLVList priors
        if inFileOptions[iInFile] in ["schoofCS"]:
            SLVList *= 1e-6 * 1e-1 # Pa to MPa and adjust to match buddFrictionParmList range
        if inFileOptions[iInFile] in ["zoetuTnoN"]:
            SLVList *= 1e6 * (60*60*24*365)

        if inFileOptions[iInFile] in ["coulomb","tsaibuddCFC","zoet"]:
            shapeGamma = 3.2
            scaleGamma = 0.1
            SLVPrior = np.flip(priorGamma(SLVList, shapeGamma, scaleGamma))
        elif inFileOptions[iInFile] in ["budd","tsaibuddBFP","schoofCS"]:
            SLVPrior = priorUniform(SLVList)
        elif inFileOptions[iInFile] in ["schoofCMAX"]:
            shapeGamma = 3.2
            scaleGamma = 0.1
            SLVPrior = priorGamma(SLVList, shapeGamma, scaleGamma)
        elif inFileOptions[iInFile] in ["zoetuTnoN"]:
            SLVPrior = priorUniform(SLVList)
        elif inFileOptions[iInFile] in ["weertmanN0","weertmanNpi"]:
            pass
        else:
            print("You forgot to define a prior")

    # calculate n-dimensional parameter prior
    if inFileOptions[iInFile] not in ["weertmanN0","weertmanNpi"]:
        parametersPrior = np.zeros((len(porosityPrior),len(grainDiameterPrior),len(SLVPrior)))
        for i in range(len(porosityPrior)):
            for j in range(len(grainDiameterPrior)):
                for k in range(len(SLVPrior)):
                    parametersPrior[i,j,k] = porosityPrior[k]*grainDiameterPrior[j]*SLVPrior[k]
    else:
        parametersPrior = np.zeros((len(porosityPrior),len(grainDiameterPrior)))
        for i in range(len(porosityPrior)):
            for j in range(len(grainDiameterPrior)):
                parametersPrior[i,j] = porosityPrior[i]*grainDiameterPrior[j]
    
    parametersPrior /= np.sum(parametersPrior)
    print(f"sum(parametersPrior): {np.sum(parametersPrior):.2f}")
    return parametersPrior

# ccc
# Plotting the posterior comparison
def plotPosteriorComparison(x, y, xTickLabels, plotPathPostComp):
    print("create posterior comparison plot")

    # Adjust font size
    plt.close()
    plt.rcParams.update({'font.size': 28})

    # colorblind friendly colors
    c = plt.rcParams['axes.prop_cycle'].by_key()['color']*3
    c = [c[1],c[0],c[2]]
    lw = 3
    
    nrows = 1
    ncols = 1
    my_dpi = 96
    fig, axs = plt.subplots(nrows, ncols, figsize=(1920/my_dpi, 1080/my_dpi), dpi=my_dpi, sharex=False, sharey=False, constrained_layout=True)

    axL = axs
    
    pModels = 1/len(x) # probability of each sliding law (model) at the beginning
    axL.hlines(y=pModels, xmin=np.min(x), xmax=np.max(x), color='k', lw=lw, ls='--')
    for xn in range(len(x)):
        ymin = np.min([pModels,y[xn]])
        ymax = np.max([pModels,y[xn]])
        axL.vlines(x=xn, ymin=ymin, ymax=ymax, color=c[1], lw=lw)

        ydif = round((y[xn]-pModels)/pModels*100,1)
        if ydif <= 0:
            va='top'
            ha='center'
            offset = -pModels*0.05
            sign = "-"
        else:
            va='bottom'
            ha='center'
            offset = pModels*0.05
            sign = "+"
        axL.text(xn, pModels+(y[xn]-pModels)+offset, f"{sign}{abs(ydif):.1f}%", color=c[1], ha=ha, va=va, fontsize=24)

    axL.scatter(x, np.ones(len(x))*pModels, s=200, color=c[1],  edgecolor=c[1], zorder=9, label='prior model probability')
    psm = axL.scatter(x, y, s=200, color=c[0],  edgecolor=c[0], zorder=10, label='posterior model probability')
    axL.set_xlabel('Sliding law')
    axL.set_ylabel('Probability')
    axL.set_xticks(x)
    axL.set_xticklabels(xTickLabels, rotation=45)
    axL.grid(True)
    axL.set_axisbelow(True)
    axL.set_ylim(y[0]-pModels*0.15,y[-1]+pModels*0.15)
    axL.legend(loc='lower right')

    # save plot
    print("Save plot at ",plotPathPostComp)
    fig.savefig(plotPathPostComp)

# ddd
# combine all of the above to determine the posterior probabilities and most likely values
inFileOptions = ["weertmanN0","weertmanNpi","coulomb","budd","tsaibuddCFC","tsaibuddBFP","schoofCMAX","zoet","zoetuTnoN"]
uniformPriors = False

priorsMin = np.zeros((len(inFileOptions)))
likelihoodsMin = np.copy(priorsMin)
posteriorsMin = np.copy(priorsMin)
likelihoods = np.copy(priorsMin)
posteriors = np.copy(priorsMin)
mostLikelyVals = np.zeros((len(inFileOptions),3))*np.nan
for s in range(len(inFileOptions)):
    chiSq, porosityList, grainDiameterList, SLVList = readData(s,inFileOptions)
    porosityListOrg = np.copy(porosityList) # keep input list with original units
    grainDiameterListOrg = np.copy(grainDiameterList)
    SLVListOrg = np.copy(SLVList)
    parametersPrior = np.transpose(defPriors(s, inFileOptions, uniformPriors, porosityList, grainDiameterList, SLVList))

    # Get minimum values
    minIndices = np.unravel_index(np.nanargmin(chiSq, axis=None), chiSq.shape)
    priorsMin[s] = parametersPrior[minIndices]
    likelihoodsMin[s] = np.exp(-0.5*chiSq[minIndices]) # 0.5 assumes a gaussian distribution of the errors of the data
    posteriorsMin[s] = likelihoodsMin[s] * priorsMin[s]

    # Calculate posterior probability
    likelihoods[s] = np.sum(chiSq)
    posteriorDist = np.exp(-0.5*chiSq)*parametersPrior # 0.5 assumes a gaussian distribution of the errors of the data
    posteriors[s] = np.sum(posteriorDist) 

    # Get most likely values
    if inFileOptions[s] in ["weertmanN0","weertmanNpi"]:
        mostLikelyInds = np.unravel_index(np.argmax(posteriorDist),posteriorDist.shape)
        mostLikelyVals[s,1:] = [grainDiameterListOrg[mostLikelyInds[0]],porosityListOrg[mostLikelyInds[1]]]
    else:
        mostLikelyInds = np.unravel_index(np.argmax(posteriorDist),posteriorDist.shape)
        mostLikelyVals[s,:] = [SLVListOrg[mostLikelyInds[0]],grainDiameterListOrg[mostLikelyInds[1]],porosityListOrg[mostLikelyInds[2]]]
    print("#--------------------------------------------------------------------------------------#")

print("Most likely values: SLV,grainDiameter,porosity")
print(mostLikelyVals)

print("Save most likely SlV, grain diameter, and porosity values in",netcdfPath)
ncFile = Dataset(netcdfPath, 'w', format='NETCDF4')

nRows = len(inFileOptions)
nCols = mostLikelyVals.shape[1]
maxStrLen = max([len(identifier) for identifier in inFileOptions])

ncFile.createDimension('identifier', nRows)
ncFile.createDimension('cols', nCols)

idVar = ncFile.createVariable('identifier', str, ('identifier',))
dataVar = ncFile.createVariable('mostLikelyVals', 'f4', ('identifier', 'cols'))
dataVar[:,:] = mostLikelyVals
dataProb = ncFile.createVariable('posteriorProbabilities', 'f4', ('identifier'))
dataProb[:] = posteriors/np.nansum(posteriors)
for i, identifier in enumerate(inFileOptions):
    idVar[i] = identifier
ncFile.close()
print("#--------------------------------------------------------------------------------------#")

# re-name the sliding laws
inFileOptions =np.asarray(inFileOptions)
inFileOptions =np.where(inFileOptions=="weertmanN0",r"Weertman ($N = 0~\mathrm{Pa}$)",inFileOptions)
inFileOptions =np.where(inFileOptions=="weertmanNpi",r"Weertman ($N = P_\mathrm{i}$)",inFileOptions)
inFileOptions =np.where(inFileOptions=="coulomb",r"Coulomb ($\mu$)",inFileOptions)
inFileOptions =np.where(inFileOptions=="budd",r"Budd ($C_{\mathrm{B}}$)",inFileOptions)
inFileOptions =np.where(inFileOptions=="tsaibuddCFC",r"Tsai-Budd ($\mu$)",inFileOptions)
inFileOptions =np.where(inFileOptions=="tsaibuddBFP",r"Tsai-Budd ($C_{\mathrm{B}}$)",inFileOptions)
inFileOptions =np.where(inFileOptions=="schoofCS",r"Schoof ($C_{\mathrm{B}}$)",inFileOptions)
inFileOptions =np.where(inFileOptions=="schoofCMAX",r"Schoof ($C_{\mathrm{max}}$)",inFileOptions)
inFileOptions =np.where(inFileOptions=="zoet",r"Zoet-Iverson ($\mu$)",inFileOptions)
inFileOptions =np.where(inFileOptions=="zoetuTnoN",r"Zoet-Iverson ($u_{\mathrm{t,noN}}$)",inFileOptions)

sortedIndices = np.argsort(posteriors)
sInFileOptions = inFileOptions[sortedIndices]
sLikelihoods = likelihoods[sortedIndices]
sPosteriors = posteriors[sortedIndices]/np.nansum(posteriors)
#print(sInFileOptions)
#print(sPosteriors)
#print(np.nansum(sPosteriors))
#print(sLikelihoods)

sortedIndicesMin = np.argsort(posteriorsMin)
sInFileOptionsMin = inFileOptions[sortedIndicesMin]
sPosteriorsMin = posteriorsMin[sortedIndicesMin]/np.nansum(posteriorsMin)

figName = "bayesTheoremPostComp"
figNameMin= "bayesTheoremPostCompMin"
addOn = 'NobCS'
if uniformPriors == True:
    figName = figName + "UniformPriors"
    figNameMin = figNameMin + "UniformPriors"
else:
    figName = figName + "CustomPriors"
    figNameMin = figNameMin + "CustomPriors"
plotPosteriorComparison(np.linspace(0,len(inFileOptions)-1,len(inFileOptions)), sPosteriors, sInFileOptions, outputPath+figName+des+addOn+".pdf")
plotPosteriorComparison(np.linspace(0,len(inFileOptions)-1,len(inFileOptions)), sPosteriorsMin, sInFileOptionsMin, outputPath+figNameMin+des+addOn+".png")

print("#--------------------------------------------------------------------------------------#")
print("End of program")
print("#--------------------------------------------------------------------------------------#")
