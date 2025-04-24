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

# Define the convolution kernels
#kernel = [0.25, 0.5, 0.25] * [0.25 0.5 0.25];

# Get basal drag
println("Perform the convolution")
#tauBedSmooth3ptTri = imfilter( Gh["tauMagBedN"], reflect(centered(kernel)), Fill(0)) ./ imfilter( Gh["mask"], reflect(centered(kernel)), Fill(0));
tauBedSmooth3ptTri = Gh["tauMagBedN"] ./ Gh["mask"];

# Get basal sliding velocity
#uBedSmooth3ptTri = imfilter( Gh["bedSpeedN"], reflect(centered(kernel)), Fill(0)) ./ imfilter( Gh["mask"], reflect(centered(kernel)), Fill(0));
uBedSmooth3ptTri = Gh["bedSpeedN"] ./ Gh["mask"];
uBedSmooth3ptTri .= uBedSmooth3ptTri ./ (60*60*24*365) # m/a to m/s
uBx = size(uBedSmooth3ptTri, 1);
uBy = size(uBedSmooth3ptTri, 2);

# Calculate overburden pressure
h=Gh["h"];             # get ice thickness in m
h[.!Gh["ice"]] .= NaN; # set to NaN where there is no ice
#h = imfilter( h, reflect(centered(kernel)), Fill(0)) ./ imfilter( Gh["mask"], reflect(centered(kernel)), Fill(0));
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

# Select and define variable Zoet-Iverson sliding law parameter (mu or uTnoN)
useCFC = true;
if useCFC == true
   println("useCFC = true: varying the Coulomb friction coefficient")
   NetCDFfile = "zoetGridSearchCFC.nc";
   coulombFrictionCoeffList = (0.01:0.01:0.7);
   SLVList = coulombFrictionCoeffList;
   nSLV = length(coulombFrictionCoeffList); # SLV = sliding law variable
else
   println("useCFC = false: varying the transition speed uT (without dependence on N)")
   NetCDFfile = "zoetGridSearchuTnoN.nc";
   startExp = -13;
   endExp = -10;
   elementsPerOrder = 25;
   totalElements = (endExp - startExp) * elementsPerOrder;
   uTcoeffList = 10.0 .^ range(startExp, stop=endExp, length=totalElements);
   SLVList = uTcoeffList;
   nSLV = length(uTcoeffList);
end

# Define Zoet-Iverson sliding law parameters
eta = 3.2e12; # Pa s
k = 0.1; # range 0.1 to 0.45
Nf = 33; # range 26 to 40
a = 0.25;
R = 0.0153; # m
L =  3.0e8; # J m^(-3)
K = 2.55; # W m^(-1) K^(-1)
Cp = 7.4e-8; # K Pa^(-1)
C1 = Cp * K / L;
k0 = 2.0 * pi / (4.0 * R) ;
uT = ( (1/(eta * (R*a)^2 * k0^3) + (4*C1)/((R*a)^2 * k0)) * Nf ) / (2 + Nf*k); # m/(s*Pa)
Phi = 32.0 * pi / 180.0; # first number is in degrees -> then converted to radians
mu = 0.5;
p = 5.0;

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
nEffectivePressureNaNsites = zeros(nSLV);
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
    
    # Loop over sliding law parameter (mu or uTnoN)
    for iSLV=1:nSLV

    	# Define sliding law function
	function zoetSlidingLaw(N, uB, tauB)
		 if useCFC == true
		    return tauB - N * coulombFrictionCoeffList[iSLV] * (uB / (uB + uT*N))^(1/p)
		 else
		    return tauB - N * mu * (uB / (uB + uTcoeffList[iSLV]*N))^(1/p)
		 end
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
			effectivePressure[i,j] = find_zero(N -> zoetSlidingLaw(N, uBedSmooth3ptTri[i, j], tauBedSmooth3ptTri[i, j]), 1.0e3)
		catch e
		      effectivePressure[i,j] = NaN;
		      # Effective pressure calculation is the same for all sites
		      # Count how often the effective pressure solve fails
		      if iSite == 1
		      	 nEffectivePressureNaNsites[iSLV] = nEffectivePressureNaNsites[iSLV] + 1;
		      end
		end
	    end
        end

	# Set overburden pressure as upper limit for effective pressure
	pressureExceeds = effectivePressure .> overburdenPressure;
        if count(pressureExceeds) > 0
           println("Number of points where effective pressure exceeds overburden pressure: ", count(pressureExceeds))
	   println("Site: ", iSite, ", SLV: ", SLVList[iSLV])
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
           println("iSite=$iSite, iSLV=$iSLV")
	elseif minimum(effectivePressureSites)<0.
           println("ERROR!!! Interpolated effective pressure is smaller than 0!!!")
	   println("iSite=$iSite, iSLV=$iSLV")
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

# Calculate failed effective pressure calculation statistics
println("#--------------------------------------------------------------------------------------#")
effectivePressureFails = sum(nEffectivePressureNaNsites)/nSLV;
effectivePressureFailsPercent = effectivePressureFails/size(findall(!isnan,uBedSmooth3ptTri))[1]*100;
nEffectivePressureFailsWeight = nEffectivePressureNaNsites ./ size(findall(!isnan,uBedSmooth3ptTri))[1]
println("On average, the effective pressure solve failed $effectivePressureFails times per SLV value ($effectivePressureFailsPercent % of all available uB-tauB pairs)")
println("#--------------------------------------------------------------------------------------#")
    
# Create NetCDF output file
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
nccreate(NetCDFfile, "nEffectivePressureFailsWeight", "NaNEffectivePressureWeight", nSLV);

ncwrite(collect(porosityList), NetCDFfile, "porosityList");
ncwrite(collect(grainDiameterList), NetCDFfile, "grainDiameterList");
if useCFC == true
   ncwrite(collect(coulombFrictionCoeffList), NetCDFfile, "SLVList");
else
   ncwrite(collect(uTcoeffList), NetCDFfile, "SLVList");
end
ncwrite(nEffectivePressureFailsWeight, NetCDFfile, "nEffectivePressureFailsWeight");

println("")
ncinfo(NetCDFfile)
println("")

println("#--------------------------------------------------------------------------------------#")
println("End of program")
println("#--------------------------------------------------------------------------------------#")
println("End time: ",now())