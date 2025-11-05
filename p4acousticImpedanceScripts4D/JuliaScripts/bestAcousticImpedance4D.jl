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
println("Perform the convolution")
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

# best values
inFileOptions = ["weertmanN0","weertmanNpi","budd","coulomb","tsaibuddCFCBFP","schoofCSCMAX","zoetCFCuTnoN"]
nInFileOptions = length(inFileOptions);

# Define in- and output paths
preFix = "";
des="";
inputPath = ""*preFix;
inputPath4D = ""*preFix;
outputPath = ""*preFix;

# Define data sites
sites = ["07", "08", "13", "15", "17"];
nSites = length(sites);
x_all=[];
y_all=[];
dist_all=[];

# Set up output arrays
nLoc = 60;
aImodels = ones(nSites, nLoc, nInFileOptions) .* NaN;
aIobs = ones(nSites, nLoc, 2) .* NaN;
obsDist = ones(nSites, nLoc) .* NaN;

minVals = ones(nInFileOptions, 4) .* NaN;
Neff = ones(uBx, uBy, nInFileOptions) .* NaN; # effective pressure fields
NeffSites = ones(nSites, nLoc, nInFileOptions) .* NaN;
AIfields = ones(uBx, uBy, nInFileOptions) .* NaN; # acoustic impedance fields

# Loop over all sites
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
    global dist_all = vcat(dist_all, dist[:]);

    # Find indices of non-NaN values in the second column of zbout
    fNotNan = findall(x -> !isnan(x), zbout[:, 2]);

    # Subsample indices with a step size of subSamp
    subSamp = 1;
    fNotNanSubsampled = fNotNan[1:subSamp:end];
    
    # Compute the length of the subsampled indices
    thisDoF = length(fNotNanSubsampled);

    # Loop over sliding laws
    for iModel = 1:nInFileOptions
	if inFileOptions[iModel] == "weertmanN0"
   	   NetCDFfile = inputPath * "weertmanGridSearch"*des*".nc"
	   porosityList = ncread(NetCDFfile, "porosityList");
	   grainDiameterList = ncread(NetCDFfile, "grainDiameterList");
	   chiSq = ncread(NetCDFfile, "chiSqZero");
	   (iBestPorosity, iBestGrainDiam) = argmin(chiSq).I;
   	   grainDiameter = grainDiameterList[iBestGrainDiam];
   	   porosity = porosityList[iBestPorosity];
   	   effectivePressure = zeros(size(tauBedSmooth3ptTri)).*tauBedSmooth3ptTri; # .*tauBedSmooth3ptTri to stay consistent with NaNs
	   minVals[iModel,1:2] = [grainDiameter,porosity];

	elseif inFileOptions[iModel] == "weertmanNpi"
   	   NetCDFfile = inputPath * "weertmanGridSearch"*des*".nc"
	   porosityList = ncread(NetCDFfile, "porosityList");
	   grainDiameterList = ncread(NetCDFfile, "grainDiameterList");
	   chiSq = ncread(NetCDFfile, "chiSq");
	   (iBestPorosity, iBestGrainDiam) = argmin(chiSq).I;
   	   grainDiameter = grainDiameterList[iBestGrainDiam];
   	   porosity = porosityList[iBestPorosity];
   	   effectivePressure = overburdenPressure;
	   minVals[iModel,1:2] = [grainDiameter,porosity];

	elseif inFileOptions[iModel] == "coulomb"
	   NetCDFfile = inputPath * "coulombGridSearch"*des*".nc"
	   porosityList = ncread(NetCDFfile, "porosityList");
	   grainDiameterList = ncread(NetCDFfile, "grainDiameterList");
	   SLVList = ncread(NetCDFfile, "SLVList");
	   chiSq = ncread(NetCDFfile, "chiSq");
	   (iBestPorosity, iBestGrainDiam, iBestSLV) = argmin(chiSq).I;
   	   grainDiameter = grainDiameterList[iBestGrainDiam];
   	   porosity = porosityList[iBestPorosity];
	   SLV = SLVList[iBestSLV];
	   effectivePressure = tauBedSmooth3ptTri / SLV;
	   minVals[iModel,1:3] = [grainDiameter,porosity,SLV];

	elseif inFileOptions[iModel] == "budd"
	   NetCDFfile = inputPath * "buddGridSearch"*des*".nc"
	   porosityList = ncread(NetCDFfile, "porosityList");
	   grainDiameterList = ncread(NetCDFfile, "grainDiameterList");
	   SLVList = ncread(NetCDFfile, "SLVList");
	   chiSq = ncread(NetCDFfile, "chiSq");
	   (iBestPorosity, iBestGrainDiam, iBestSLV) = argmin(chiSq).I;
   	   grainDiameter = grainDiameterList[iBestGrainDiam];
   	   porosity = porosityList[iBestPorosity];
	   SLV = SLVList[iBestSLV];
	   m = 1.0/3.0;
	   q = 1.0;
	   effectivePressure = (tauBedSmooth3ptTri ./ (SLV .* uBedSmooth3ptTri.^m) ).^(1/q);
	   minVals[iModel,1:3] = [grainDiameter,porosity,SLV];

	elseif inFileOptions[iModel] == "tsaibuddCFC"
	   NetCDFfile = inputPath * "tsaiBuddGridSearchCFC"*des*".nc"
	   porosityList = ncread(NetCDFfile, "porosityList");
	   grainDiameterList = ncread(NetCDFfile, "grainDiameterList");
	   SLVList = ncread(NetCDFfile, "SLVList");
	   chiSq = ncread(NetCDFfile, "chiSq");
	   (iBestPorosity, iBestGrainDiam, iBestSLV) = argmin(chiSq).I;
   	   grainDiameter = grainDiameterList[iBestGrainDiam];
   	   porosity = porosityList[iBestPorosity];
	   SLV = SLVList[iBestSLV];
	   Cb = 37.01;
	   m = 1.0/3.0;
	   mu = 0.5;
	   effectivePressure = tauBedSmooth3ptTri ./ min.(SLV,Cb .* uBedSmooth3ptTri.^m); # this assumes q=1
	   minVals[iModel,1:3] = [grainDiameter,porosity,SLV];

	elseif inFileOptions[iModel] == "tsaibuddBFP"
	   NetCDFfile = inputPath * "tsaiBuddGridSearchBFP"*des*".nc"
	   porosityList = ncread(NetCDFfile, "porosityList");
	   grainDiameterList = ncread(NetCDFfile, "grainDiameterList");
	   SLVList = ncread(NetCDFfile, "SLVList");
	   chiSq = ncread(NetCDFfile, "chiSq");
	   (iBestPorosity, iBestGrainDiam, iBestSLV) = argmin(chiSq).I;
   	   grainDiameter = grainDiameterList[iBestGrainDiam];
   	   porosity = porosityList[iBestPorosity];
	   SLV = SLVList[iBestSLV];
	   Cb = 37.01;
	   m = 1.0/3.0;
	   mu = 0.5;
	   effectivePressure = tauBedSmooth3ptTri ./ min.(mu,SLV .* uBedSmooth3ptTri.^m); # this assumes q=1
	   minVals[iModel,1:3] = [grainDiameter,porosity,SLV];

	elseif inFileOptions[iModel] == "tsaibuddCFCBFP"
	   NetCDFfile = inputPath4D * "tsaiBuddGridSearchCFCBFP.4D.nc"
	   porosityList = ncread(NetCDFfile, "porosityList");
	   grainDiameterList = ncread(NetCDFfile, "grainDiameterList");
	   BFPList = ncread(NetCDFfile, "BFPList");
	   CFCList = ncread(NetCDFfile, "CFCList");
	   chiSq = ncread(NetCDFfile, "chiSq");
	   (iBestPorosity, iBestGrainDiam, iBestCFC, iBestBFP) = argmin(chiSq).I;
   	   grainDiameter = grainDiameterList[iBestGrainDiam];
   	   porosity = porosityList[iBestPorosity];
	   Cb = BFPList[iBestBFP];
	   m = 1.0/3.0;
	   mu = CFCList[iBestCFC];
	   effectivePressure = tauBedSmooth3ptTri ./ min.(mu,Cb .* uBedSmooth3ptTri.^m); # this assumes q=1
	   minVals[iModel,:] = [grainDiameter,porosity,mu,Cb];

	elseif inFileOptions[iModel] == "schoofCMAX"
	   NetCDFfile = inputPath * "schoofGridSearchCMAX"*des*".nc"
	   porosityList = ncread(NetCDFfile, "porosityList");
	   grainDiameterList = ncread(NetCDFfile, "grainDiameterList");
	   SLVList = ncread(NetCDFfile, "SLVList");
	   chiSq = ncread(NetCDFfile, "chiSq");
	   (iBestPorosity, iBestGrainDiam, iBestSLV) = argmin(chiSq).I;
   	   grainDiameter = grainDiameterList[iBestGrainDiam];
   	   porosity = porosityList[iBestPorosity];
	   SLV = SLVList[iBestSLV];
	   m = 1. / 3.;
	   cS = 1.0e9;
	   cMAX = 0.2;

	   function schoofSlidingLawCMAX(N, uB, tauB)
	   	    return tauB - cS * uB^m / (1 + (cS/(SLV * N))^(1/m) * uB )^m
	   end

	   effectivePressure=zeros(uBx,uBy);
           for i in 1:uBx
               for j in 1:uBy
		   # find_zero() returns initial guess if uB or tauB is NaN -> set effective pressure to NaN if this is the case
                   if isnan(uBedSmooth3ptTri[i, j]) || isnan(tauBedSmooth3ptTri[i, j])
                      effectivePressure[i,j] = NaN;
                      continue
                   end
                   try
			effectivePressure[i,j] = find_zero(N -> schoofSlidingLawCMAX(N, uBedSmooth3ptTri[i, j], tauBedSmooth3ptTri[i, j]), 1.0e3)
                   catch e
                   	 effectivePressure[i,j] = NaN;
		   end
	       end
	   end
   	   minVals[iModel,1:3] = [grainDiameter,porosity,SLV];

	elseif inFileOptions[iModel] == "schoofCS"
	   NetCDFfile = inputPath * "schoofGridSearchCS"*des*".nc"
	   porosityList = ncread(NetCDFfile, "porosityList");
	   grainDiameterList = ncread(NetCDFfile, "grainDiameterList");
	   SLVList = ncread(NetCDFfile, "SLVList");
	   chiSq = ncread(NetCDFfile, "chiSq");
	   chiSq = replace(chiSq, NaN => Inf) # otherwise the NaN values get picked up as mins
	   (iBestPorosity, iBestGrainDiam, iBestSLV) = argmin(chiSq).I;
   	   grainDiameter = grainDiameterList[iBestGrainDiam];
   	   porosity = porosityList[iBestPorosity];
	   SLV = SLVList[iBestSLV];
	   m = 1. / 3.;
	   cS = 1.0e9;
	   cMAX = 0.2;

	   function schoofSlidingLawCS(N, uB, tauB)
	   	    return tauB - SLV * uB^m / (1 + (SLV/(cMAX * N))^(1/m) * uB )^m
	   end

	   effectivePressure=zeros(uBx,uBy);
           for i in 1:uBx
               for j in 1:uBy
	       	   # find_zero() returns initial guess if uB or tauB is NaN -> set effective pressure to NaN if this is the case
                   if isnan(uBedSmooth3ptTri[i, j]) || isnan(tauBedSmooth3ptTri[i, j])
                      effectivePressure[i,j] = NaN;
                      continue
                   end
                   try
			effectivePressure[i,j] = find_zero(N -> schoofSlidingLawCS(N, uBedSmooth3ptTri[i, j], tauBedSmooth3ptTri[i, j]), 1.0e3)
                   catch e
                   	 effectivePressure[i,j] = NaN;
		   end
	       end
	   end
   	   minVals[iModel,1:3] = [grainDiameter,porosity,SLV];

	elseif inFileOptions[iModel] == "schoofCSCMAX"
	   NetCDFfile = inputPath4D * "schoofGridSearchCSCMAX.4D.nc"
	   porosityList = ncread(NetCDFfile, "porosityList");
	   grainDiameterList = ncread(NetCDFfile, "grainDiameterList");
	   CSList = ncread(NetCDFfile, "CSList");
	   cMAXList = ncread(NetCDFfile, "cMAXList");
	   chiSq = ncread(NetCDFfile, "chiSq");
	   chiSq = replace(chiSq, NaN => Inf) # otherwise the NaN values get picked up as mins
	   (iBestPorosity, iBestGrainDiam, iBestCS, iBestcMAX) = argmin(chiSq).I;
   	   grainDiameter = grainDiameterList[iBestGrainDiam];
   	   porosity = porosityList[iBestPorosity];
	   cS = CSList[iBestCS];
	   m = 1.0/3.0;
	   cMAX = cMAXList[iBestcMAX];

	   function schoofSlidingLawCScMAX(N, uB, tauB)
	   	    return tauB - cS * uB^m / (1 + (cS/(cMAX * N))^(1/m) * uB )^m
	   end

	   effectivePressure=zeros(uBx,uBy);
           for i in 1:uBx
               for j in 1:uBy
	       	   # find_zero() returns initial guess if uB or tauB is NaN -> set effective pressure to NaN if this is the case
                   if isnan(uBedSmooth3ptTri[i, j]) || isnan(tauBedSmooth3ptTri[i, j])
                      effectivePressure[i,j] = NaN;
                      continue
                   end
                   try
			effectivePressure[i,j] = find_zero(N -> schoofSlidingLawCScMAX(N, uBedSmooth3ptTri[i, j], tauBedSmooth3ptTri[i, j]), 1.0e3)
                   catch e
                   	 effectivePressure[i,j] = NaN;
		   end
	       end
	   end

	   minVals[iModel,:] = [grainDiameter,porosity,cMAX,cS];

	elseif inFileOptions[iModel] == "zoet"
	   NetCDFfile = inputPath * "zoetGridSearchCFC"*des*".nc"
	   porosityList = ncread(NetCDFfile, "porosityList");
	   grainDiameterList = ncread(NetCDFfile, "grainDiameterList");
	   SLVList = ncread(NetCDFfile, "SLVList");
	   chiSq = ncread(NetCDFfile, "chiSq");
	   (iBestPorosity, iBestGrainDiam, iBestSLV) = argmin(chiSq).I;
   	   grainDiameter = grainDiameterList[iBestGrainDiam];
   	   porosity = porosityList[iBestPorosity];
	   SLV = SLVList[iBestSLV];
	   eta = 3.2e12; # Pa s
	   k = 0.1; # range 0.1 to 0.45
	   Nf = 33; # range 26 to 40
	   a = 0.25;
	   R = 0.015; # m
	   L =  3.0e8; # J m^(-3)
	   K = 2.55; # W m^(-1) K^(-1)
	   Cp = 7.4e-8; # K Pa^(-1)
	   C1 = Cp * K / L;
	   k0 = 2.0 * pi / (4.0 * R) ;
	   uT = ( (1/(eta * (R*a)^2 * k0^3) + (4*C1)/((R*a)^2 * k0)) * Nf ) / (2 + Nf*k); # m/(s*Pa)
	   mu = 0.5;
	   p = 5.0;

	   function zoetSlidingLaw(N, uB, tauB)
	   	    return tauB - N * SLV * (uB / (uB + uT*N))^(1/p)
	   end

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
		   end
	       end
	   end
	   minVals[iModel,1:3] = [grainDiameter,porosity,SLV];

	elseif inFileOptions[iModel] == "zoetuTnoN"
	   NetCDFfile = inputPath * "zoetGridSearchuTnoN"*des*".nc"
	   porosityList = ncread(NetCDFfile, "porosityList");
	   grainDiameterList = ncread(NetCDFfile, "grainDiameterList");
	   SLVList = ncread(NetCDFfile, "SLVList");
	   chiSq = ncread(NetCDFfile, "chiSq");
	   (iBestPorosity, iBestGrainDiam, iBestSLV) = argmin(chiSq).I;
   	   grainDiameter = grainDiameterList[iBestGrainDiam];
   	   porosity = porosityList[iBestPorosity];
	   SLV = SLVList[iBestSLV];
	   eta = 3.2e12; # Pa s
	   k = 0.1; # range 0.1 to 0.45
	   Nf = 33; # range 26 to 40
	   a = 0.25;
	   R = 0.015; # m
	   L =  3.0e8; # J m^(-3)
	   K = 2.55; # W m^(-1) K^(-1)
	   Cp = 7.4e-8; # K Pa^(-1)
	   C1 = Cp * K / L;
	   k0 = 2.0 * pi / (4.0 * R) ;
	   uT = ( (1/(eta * (R*a)^2 * k0^3) + (4*C1)/((R*a)^2 * k0)) * Nf ) / (2 + Nf*k); # m/(s*Pa)
	   mu = 0.5;
	   p = 5.0;

	   function zoetSlidingLawuTnoN(N, uB, tauB)
	   	    return tauB - N * mu * (uB / (uB + SLV*N))^(1/p)
	   end

	   effectivePressure=zeros(uBx,uBy);
           for i in 1:uBx
               for j in 1:uBy
	       	   # find_zero() returns initial guess if uB or tauB is NaN -> set effective pressure to NaN if this is the case
                   if isnan(uBedSmooth3ptTri[i, j]) || isnan(tauBedSmooth3ptTri[i, j])
                      effectivePressure[i,j] = NaN;
                      continue
                   end
                   try
			effectivePressure[i,j] = find_zero(N -> zoetSlidingLawuTnoN(N, uBedSmooth3ptTri[i, j], tauBedSmooth3ptTri[i, j]), 1.0e3)
                   catch e
                   	 effectivePressure[i,j] = NaN;
		   end
	       end
	   end
	   minVals[iModel,1:3] = [grainDiameter,porosity,SLV];

	elseif inFileOptions[iModel] == "zoetCFCuTnoN"
	   NetCDFfile = inputPath4D * "zoetGridSearchCFCuTnoN.4D.nc"
	   porosityList = ncread(NetCDFfile, "porosityList");
	   grainDiameterList = ncread(NetCDFfile, "grainDiameterList");
	   CFCList = ncread(NetCDFfile, "CFCList");
	   uTList = ncread(NetCDFfile, "uTList");
	   chiSq = ncread(NetCDFfile, "chiSq");
	   chiSq = replace(chiSq, NaN => Inf) # otherwise the NaN values get picked up as mins
	   (iBestPorosity, iBestGrainDiam, iBestCFC, iBestuT) = argmin(chiSq).I;
   	   grainDiameter = grainDiameterList[iBestGrainDiam];
   	   porosity = porosityList[iBestPorosity];
	   mu = CFCList[iBestCFC];
	   uTnoN = uTList[iBestuT];
	   p = 5.0;

	   function zoetSlidingLawCFCuTnoN(N, uB, tauB)
	   	    return tauB - N * mu * (uB / (uB + uTnoN*N))^(1/p)
	   end

	   effectivePressure=zeros(uBx,uBy);
           for i in 1:uBx
               for j in 1:uBy
	       	   # find_zero() returns initial guess if uB or tauB is NaN -> set effective pressure to NaN if this is the case
                   if isnan(uBedSmooth3ptTri[i, j]) || isnan(tauBedSmooth3ptTri[i, j])
                      effectivePressure[i,j] = NaN;
                      continue
                   end
                   try
			effectivePressure[i,j] = find_zero(N -> zoetSlidingLawCFCuTnoN(N, uBedSmooth3ptTri[i, j], tauBedSmooth3ptTri[i, j]), 1.0e3)
                   catch e
                   	 effectivePressure[i,j] = NaN;
		   end
	       end
	   end
	   minVals[iModel,:] = [grainDiameter,porosity,mu,uTnoN];

	end

	# Set overburden pressure as upper limit for effective pressure
	pressureExceeds = effectivePressure .> overburdenPressure;
        if count(pressureExceeds) > 0
           println("Number of points where effective pressure exceeds overburden pressure: ", count(pressureExceeds))
           println("Site: ", iSite, ", Model: ", inFileOptions[iModel])
           println("Effective pressure was limited to overburden pressure")
           effectivePressure = min.(effectivePressure, overburdenPressure);
        end

	# Store effective pressure and acoustic impedance fields (same for all sites)
	if iSite==1
           Neff[:,:,iModel] = effectivePressure;
           ZPfield,ZSfield = acousticImpedance(effectivePressure,porosity,grainDiameter)
           AIfields[:,:,iModel] = ZPfield;
        end

	# Set up interpolation
	F = Interpolations.interpolate((Gh["xx"][:,1], Gh["yy"][1,:]), effectivePressure, Gridded(Linear()));

    	# Calculate effective pressure at site coordinates
    	effectivePressureSites = F(x, y);
    	effectivePressureSites = [effectivePressureSites[i, i] for i in 1:minimum(size(effectivePressureSites))];
	NeffSites[iSite,fNotNanSubsampled,iModel] = Spline1D(dist,effectivePressureSites; k=1)(zbout[fNotNanSubsampled,1]);

    	# Calculate acoustic impedance at site coordinates
    	ZPsite,ZSsite = acousticImpedance(effectivePressureSites,porosity,grainDiameter);

    	# Calculate acoustic impedance at data coordinates
    	ZP = Spline1D(dist,ZPsite; k=1)(zbout[fNotNanSubsampled,1]);
    	aImodels[iSite,fNotNanSubsampled,iModel] = ZP;

    end

    aIobs[iSite,fNotNanSubsampled, :] = zbout[fNotNanSubsampled,2:3];
    obsDist[iSite,fNotNanSubsampled] = zbout[fNotNanSubsampled,1];

end 

# Create NetCDF output file
addOn = ""
NetCDFfile = outputPath*"acousticImpedanceModelObsComp"*des*addOn*".nc";
println("Save NetCDF output at ",NetCDFfile)
isfile(NetCDFfile) && rm(NetCDFfile);

nccreate(NetCDFfile, "aImodels", "sites", nSites, "locations", nLoc, "models", nInFileOptions);
nccreate(NetCDFfile, "aIobs", "sites", nSites, "locations", nLoc, "data", 2);
nccreate(NetCDFfile, "obsDist", "sites", nSites, "locations", nLoc);
nccreate(NetCDFfile, "sitesList", "sites", nSites);
nccreate(NetCDFfile, "locationsList", "locations", nLoc);
#nccreate(NetCDFfile, "modelsList", "models", nInFileOptions);
#nccreate(NetCDFfile, "dataList", "data", 2);
nccreate(NetCDFfile, "minVals", "models", nInFileOptions, "parameters", 4);
nccreate(NetCDFfile, "uB", "xData", uBx, "yData", uBy);
nccreate(NetCDFfile, "tauB", "xData", uBx, "yData", uBy);
nccreate(NetCDFfile, "Neff", "xData", uBx, "yData", uBy, "models", nInFileOptions);
nccreate(NetCDFfile, "NeffSites", "sites", nSites, "locations", nLoc, "models", nInFileOptions);
nccreate(NetCDFfile, "AIfields", "xData", uBx, "yData", uBy, "models", nInFileOptions);
nccreate(NetCDFfile, "xx", "xData", uBx, "yData", uBy);
nccreate(NetCDFfile, "yy", "xData", uBx, "yData", uBy);

ncwrite(aImodels, NetCDFfile, "aImodels");
ncwrite(aIobs, NetCDFfile, "aIobs");
ncwrite(obsDist, NetCDFfile, "obsDist");
ncwrite(parse.(Int, sites), NetCDFfile, "sitesList");
ncwrite(collect(1:1:nLoc), NetCDFfile, "locationsList");
#ncwrite(inFileOptions, NetCDFfile, "modelsList");
#ncwrite(["data", "std"], NetCDFfile, "dataList");
ncwrite(minVals, NetCDFfile, "minVals");
ncwrite(uBedSmooth3ptTri, NetCDFfile, "uB");
ncwrite(tauBedSmooth3ptTri, NetCDFfile, "tauB");
ncwrite(Neff, NetCDFfile, "Neff");
ncwrite(NeffSites, NetCDFfile, "NeffSites");
ncwrite(AIfields, NetCDFfile, "AIfields");
ncwrite(Gh["xx"], NetCDFfile, "xx");
ncwrite(Gh["yy"], NetCDFfile, "yy");

println("")
ncinfo(NetCDFfile)
println("")

println("#--------------------------------------------------------------------------------------#")
println("End of program")
println("#--------------------------------------------------------------------------------------#")
println("End time: ",now())