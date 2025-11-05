@eval Base println(xs...) = begin res = println(stdout::IO, xs...); flush(stdout); return res end;
@eval Base println(io::IO) =  begin res =  print(io, '\n'); flush(io); return res end;
using Dates
println("Start time: ",now())

using MAT,ImageFiltering,Interpolations,CSV,DataFrames,Proj4,Dierckx,NetCDF,Roots

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

# Define Budd friction parameter
startExp = -1;		# Starting exponent (10^?)
endExp = 6;		# Ending exponent (10^?)
elementsPerOrder = 20;
totalElements = (endExp - startExp) * elementsPerOrder;
buddFrictionParmList = 10.0 .^ range(startExp, stop=endExp, length=totalElements);
nSLV = length(buddFrictionParmList); # SLV = sliding law variable

# Define Budd sliding law parameters
m = 1.0/3.0;
q = 1.0;

# Set up output arrays (misfit)
dof=zeros(nPorosity,nGrainDiam,nSLV);
chiSq=zeros(nPorosity,nGrainDiam,nSLV);

# Define data sites
sites = ["07", "08", "13", "15", "17"];
nSites = length(sites);
x_all=[];
y_all=[];

# Loop over all sites
icount = 0;
ncount = nPorosity * nGrainDiam * nSLV * nSites;
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
    
    # Loop over sliding law parameter (Cb)
    for iSLV=1:nSLV

    	# Calculate effective pressure
	effectivePressure = (tauBedSmooth3ptTri ./ (buddFrictionParmList[iSLV] .* uBedSmooth3ptTri.^m) ).^(1/q);

	# Set overburden pressure as upper limit for effective pressure
	pressureExceeds = effectivePressure .> overburdenPressure;
        if count(pressureExceeds) > 0
           println("Number of points where effective pressure exceeds overburden pressure: ", count(pressureExceeds))
           println("Site: ", iSite, ", SLV: ", buddFrictionParmList[iSLV])
           println("Effective pressure was limited to overburden pressure")
           effectivePressure = min.(effectivePressure, overburdenPressure);
        end

    	# Set up interpolation
    	F = Interpolations.interpolate((Gh["xx"][:,1], Gh["yy"][1,:]), effectivePressure, Gridded(Linear()));

	# Evaluate interpolator for effective pressure
    	effectivePressureSites = F(x, y);
	effectivePressureSites = [effectivePressureSites[i, i] for i in 1:minimum(size(effectivePressureSites))];

	# Check if interpolated effective pressure is negative
	if minimum(effectivePressureSites)<0.
	   println("ERROR!!! Interpolated effective pressure is smaller than 0!!!")
	   println("iSite=",iSite)
	   println("iSLV=",iSLV)
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
		dof[iPorosity,iGrainDiam,iSLV]=dof[iPorosity,iGrainDiam,iSLV]+thisDoF;
            	chiSq[iPorosity,iGrainDiam,iSLV]=chiSq[iPorosity,iGrainDiam,iSLV]+thisChiSq;
            
		global icount = icount + 1;
            end
	    # Print out progress
	    println(icount/ncount*100," % done")
        end
    end
end
    
# Create NetCDF output file
NetCDFfile = "buddGridSearch.nc";
println("Save NetCDF output at ",NetCDFfile)
isfile(NetCDFfile) && rm(NetCDFfile);

outputVars = ["dof","chiSq"];
nVars = length(outputVars);
outputData = zeros(nPorosity,nGrainDiam,nSLV,nVars);
outputData[:,:,:,1] = dof;
outputData[:,:,:,2] = chiSq;
for iVars=1:nVars
    nccreate(NetCDFfile, outputVars[iVars], "porosity", nPorosity, "grainDiameter", nGrainDiam, "SLV", nSLV);
    ncwrite(outputData[:,:,:,iVars], NetCDFfile, outputVars[iVars]);
end

nccreate(NetCDFfile, "porosityList", "porosity", nPorosity);
nccreate(NetCDFfile, "grainDiameterList", "grainDiameter", nGrainDiam);
nccreate(NetCDFfile, "SLVList", "SLV", nSLV);

ncwrite(collect(porosityList), NetCDFfile, "porosityList");
ncwrite(collect(grainDiameterList), NetCDFfile, "grainDiameterList");
ncwrite(collect(buddFrictionParmList), NetCDFfile, "SLVList");

println("")
ncinfo(NetCDFfile)
println("")

println("#--------------------------------------------------------------------------------------#")
println("End of program")
println("#--------------------------------------------------------------------------------------#")
println("End time: ",now())