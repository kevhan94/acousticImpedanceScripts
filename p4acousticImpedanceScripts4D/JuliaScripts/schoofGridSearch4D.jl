@eval Base println(xs...) = begin res = println(stdout::IO, xs...); flush(stdout); return res end;
@eval Base println(io::IO) =  begin res =  print(io, '\n'); flush(io); return res end;
using Dates
println("Start time: ",now())

using MAT,ImageFiltering,Interpolations,CSV,DataFrames,Proj4,Dierckx,NetCDF,Roots

des = ARGS[1];

# Include functions in other files
include("geog_to_pol_wgs84_71S.jl");
include("acousticImpedance.jl");

println("#--------------------------------------------------------------------------------------#")
println("Start of program")
println("#--------------------------------------------------------------------------------------#")

# Load MATLAB file
println("Load MATLAB file")
data = matread("Inverse_1km_Both_BasinsV1.mat");
Gh = data["Gh"];

# Get basal drag
tauBedSmooth3ptTri = Gh["tauMagBedN"] ./ Gh["mask"];

# Get basal sliding velocity
uBedSmooth3ptTri = Gh["bedSpeedN"] ./ Gh["mask"];
uBedSmooth3ptTri .= uBedSmooth3ptTri ./ (60*60*24*365) # m/a to m/s
uBx = size(uBedSmooth3ptTri, 1);
uBy = size(uBedSmooth3ptTri, 2);

# Calculate overburden pressure
h=Gh["h"];             # get ice thickness in m
h[.!Gh["ice"]] .= NaN; # set to NaN where there is no ice
h = h ./ Gh["mask"];
g = 9.81; # m/s**2
rhoIce = 910; # kg/m**3
overburdenPressure = h .* g .* rhoIce;

# Set up the parameter space
phiList=(-1:0.01:10);
porosityList=(0.3:0.002:0.6);
grainDiameterList = (1e-3) .* 2 .^ (-phiList);
nPorosity=length(porosityList);
nGrainDiam=length(grainDiameterList);

# Select and define variable Schoof sliding law parameters (CS and cMAX)
println("varying the Schoof friction coefficient and the maximum value of C_max=tau_b/N which is bounded by the local maximum up-slope of the bedrock")
NetCDFfile = "schoofGridSearchCSCMAX"*des*".nc";
startExp = 6;
endExp = 13;
elementsPerOrder = 20;
totalElements = (endExp - startExp) * elementsPerOrder;
schoofFrictionCoeffList = 10.0 .^ range(startExp, stop=endExp, length=totalElements);
nCS = length(schoofFrictionCoeffList);
cMAXLow = parse(Float64, ARGS[2]);
cMAXUp = parse(Float64, ARGS[3]);
cMAXcoeffList = (cMAXLow:0.01:cMAXUp);
ncMAX = length(cMAXcoeffList);

# Define Schoof sliding law parameters
m = 1. / 3.;

# Set up output arrays (misfit)
dof=zeros(nPorosity,nGrainDiam,nCS,ncMAX);
chiSq=zeros(nPorosity,nGrainDiam,nCS,ncMAX);

# Define data sites
sites = ["07", "08", "13", "15", "17"];
nSites = length(sites);
x_all=[];
y_all=[];

# Loop over all sites
icount = 0;
ncount = nPorosity * nGrainDiam * nCS * ncMAX * nSites;
nEffectivePressureNaNsites = zeros(nCS,ncMAX);
println("Start loop over site, time: ",now())
for iSite in 1:nSites
    println("Site: ", iSite)
    
    # Read CSV file
    distLonLat = CSV.read("istar$(sites[iSite])a_lonlatinterp.csv", DataFrame,header=false);
    
    # Load MATLAB file
    matfile = matread("zbout_istar$(sites[iSite]).mat");
    zbout = matfile["zbout"];
    
    dist = distLonLat[:, 1];
    lon = distLonLat[:, 2];
    lat = distLonLat[:, 3];
    
    # Call the geog_to_pol_wgs84_71S function
    x, y = geog_to_pol_wgs84_71S(lat, lon);
    
    # Append x and y values to x_all and y_all respectively
    global x_all = vcat(x_all, x[:]);
    global y_all = vcat(y_all, y[:]);

    # Find indices of non-NaN values in the second column of zbout
    fNotNan = findall(x -> !isnan(x), zbout[:, 2]);

    # Subsample indices with a step size of subSamp
    subSamp = 1;
    fNotNanSubsampled = fNotNan[1:subSamp:end];
    
    # Compute the length of the subsampled indices
    thisDoF = length(fNotNanSubsampled);
    
    # Loop over sliding law parameter (CS or cMAX)
    for iCS=1:nCS
    	for icMAX=1:ncMAX

    	    # Define sliding law function
	    function brondexSlidingLaw(N, uB, tauB)
		    return tauB - schoofFrictionCoeffList[iCS] * uB^m / (1 + (schoofFrictionCoeffList[iCS]/(cMAXcoeffList[icMAX] * N))^(1/m) * uB )^m
	    end

	    # Calculate effective pressure
	    effectivePressure=zeros(uBx,uBy);
	    for i in 1:uBx
            	for j in 1:uBy
	    	    # find_zero() returns initial guess if uB or tauB is NaN -> set effective pressure to NaN if this is the case
                    if isnan(uBedSmooth3ptTri[i, j]) || isnan(tauBedSmooth3ptTri[i, j])
                       effectivePressure[i,j] = NaN;
                       continue
                    end
                    try
                        effectivePressure[i,j] = find_zero(N -> brondexSlidingLaw(N, uBedSmooth3ptTri[i, j], tauBedSmooth3ptTri[i, j]), 1.0e3)
                    catch e
                      	effectivePressure[i,j] = NaN;
		      	# Effective pressure calculation is the same for all sites
		      	# Count how often the effective pressure solve fails
                      	if iSite == 1
                           nEffectivePressureNaNsites[iCS,icMAX] = nEffectivePressureNaNsites[iCS,icMAX] + 1;
                      	end
        	    end
                end
            end
	
	    # Set overburden pressure as upper limit for effective pressure
	    pressureExceeds = effectivePressure .> overburdenPressure;
            if count(pressureExceeds) > 0
               println("Number of points where effective pressure exceeds overburden pressure: ", count(pressureExceeds))
               println("Site: ", iSite, ", CS: ", schoofFrictionCoeffList[iCS], ", cMAX: ", cMAXcoeffList[icMAX])
               println("Effective pressure was limited to overburden pressure")
               effectivePressure = min.(effectivePressure, overburdenPressure);
           end

    	   # Set up interpolation
    	   F = Interpolations.interpolate((Gh["xx"][:,1], Gh["yy"][1,:]), effectivePressure, Gridded(Linear()));

	   # Evaluate interpolator for effective pressure
    	   effectivePressureSites = F(x, y);
	   effectivePressureSites = [effectivePressureSites[i, i] for i in 1:minimum(size(effectivePressureSites))];

	   # Check if interpolated effective pressure contains NaN or negative values
           if size(findall(isnan, effectivePressureSites))[1] > 0
              println("ERROR!!! Interpolated effective pressure contains NaNs!!!")
              println("iSite=$iSite, iCS=$iCS, icMAX=$icMAX")
           elseif minimum(effectivePressureSites)<0.
              println("ERROR!!! Interpolated effective pressure is smaller than 0!!!")
              println("iSite=$iSite, iCS=$iCS, icMAX=$icMAX")
           end

	   # Loop over VGS parameter space
    	   for iPorosity=1:nPorosity
    	       for iGrainDiam=1:nGrainDiam

	    	   # Set parameter values
	    	   porosity=porosityList[iPorosity];
            	   grainDiameter=grainDiameterList[iGrainDiam];
            
		   # Calculate acoustic impedance based on VGS theory
		   ZP_sites,ZS_sites = acousticImpedance(effectivePressureSites,porosity,grainDiameter);
            
		   # Calculate acoustic impedance misfit and sum it for all selected data points of this side (fNotNanSubsampled)
		   thisChiSq=sum(((zbout[fNotNanSubsampled,2]-Spline1D(dist,ZP_sites; k=1)(zbout[fNotNanSubsampled,1])) ./ zbout[fNotNanSubsampled,3]).^2);
            
		   # Store misfit
		   dof[iPorosity,iGrainDiam,iCS,icMAX]=dof[iPorosity,iGrainDiam,iCS,icMAX]+thisDoF;
            	   chiSq[iPorosity,iGrainDiam,iCS,icMAX]=chiSq[iPorosity,iGrainDiam,iCS,icMAX]+thisChiSq;
            
		   global icount = icount + 1;
               end
	       # Print out progress
	       println(icount/ncount*100," % done")
	   end
        end
    end
end

# Calculate failed effective pressure calculation statistics
println("#--------------------------------------------------------------------------------------#")
effectivePressureFails = sum(nEffectivePressureNaNsites)/(nCS+ncMAX);
effectivePressureFailsPercent = effectivePressureFails/size(findall(!isnan,uBedSmooth3ptTri))[1]*100;
nEffectivePressureFailsWeight = nEffectivePressureNaNsites ./ size(findall(!isnan,uBedSmooth3ptTri))[1]
println("On average, the effective pressure solve failed $effectivePressureFails times per CS/cMAX value ($effectivePressureFailsPercent % of all available uB-tauB pairs)")
println("#--------------------------------------------------------------------------------------#")
    
# Create NetCDF output file
println("Save NetCDF output at ",NetCDFfile)
isfile(NetCDFfile) && rm(NetCDFfile);

outputVars = ["dof","chiSq"];
nVars = length(outputVars);
outputData = zeros(nPorosity,nGrainDiam,nCS,ncMAX,nVars);
outputData[:,:,:,:,1] = dof;
outputData[:,:,:,:,2] = chiSq;
for iVars=1:nVars
    nccreate(NetCDFfile, outputVars[iVars], "porosity", nPorosity, "grainDiameter", nGrainDiam, "CS", nCS, "cMAX", ncMAX);
    ncwrite(outputData[:,:,:,:,iVars], NetCDFfile, outputVars[iVars]);
end

nccreate(NetCDFfile, "porosityList", "porosity", nPorosity);
nccreate(NetCDFfile, "grainDiameterList", "grainDiameter", nGrainDiam);
nccreate(NetCDFfile, "CSList", "CS", nCS);
nccreate(NetCDFfile, "cMAXList", "cMAX", ncMAX);
nccreate(NetCDFfile, "nEffectivePressureFailsWeight", "nCS", nCS, "ncMAX", ncMAX);

ncwrite(collect(porosityList), NetCDFfile, "porosityList");
ncwrite(collect(grainDiameterList), NetCDFfile, "grainDiameterList");
ncwrite(collect(schoofFrictionCoeffList), NetCDFfile, "CSList");
ncwrite(collect(cMAXcoeffList), NetCDFfile, "cMAXList");
ncwrite(nEffectivePressureFailsWeight, NetCDFfile, "nEffectivePressureFailsWeight");

println("")
ncinfo(NetCDFfile)
println("")

println("#--------------------------------------------------------------------------------------#")
println("End of program")
println("#--------------------------------------------------------------------------------------#")
println("End time: ",now())