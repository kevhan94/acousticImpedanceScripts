@eval Base println(xs...) = begin res = println(stdout::IO, xs...); flush(stdout); return res end;
@eval Base println(io::IO) =  begin res =  print(io, '\n'); flush(io); return res end;
using Dates
println("Start time: ",now())

using MAT,ImageFiltering,Interpolations,CSV,DataFrames,Proj4,Dierckx,NetCDF

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

# Calculate overburden pressure
h=Gh["h"];             # get ice thickness in m
h[.!Gh["ice"]] .= NaN; # set to NaN where there is no ice
h = h ./ Gh["mask"];
g = 9.81; # m/s**2
rhoIce = 910; # kg/m**3
overburdenPressure = h .* g .* rhoIce;

# Set up interpolation
println("Set up interpolation")
G = Interpolations.interpolate((Gh["xx"][:,1], Gh["yy"][1,:]), overburdenPressure, Gridded(Linear()));

# Set up the parameter space
phiList=(-1:0.01:10);
porosityList=(0.3:0.002:0.6);
grainDiameterList = (1e-3) .* 2 .^ (-phiList);
nPorosity=length(porosityList);
nGrainDiam=length(grainDiameterList);

# Set up output arrays (misfit)
dof=zeros(nPorosity,nGrainDiam);
chiSq=zeros(nPorosity,nGrainDiam);
chiSqZero=zeros(nPorosity,nGrainDiam);

# Define data sites
sites = ["07", "08", "13", "15", "17"];
nSites = length(sites);
x_all=[];
y_all=[];

# Loop over all sites
icount = 0;
ncount = nPorosity * nGrainDiam * nSites;
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
    
    # Evaluate interpolator for overburden pressure
    overburdenPressureSites = G(x, y);
    overburdenPressureSites = [overburdenPressureSites[i, i] for i in 1:minimum(size(overburdenPressureSites))];

    # Check if interpolated overburden pressure is negative
    if minimum(overburdenPressureSites)<0.
           println("ERROR!!! Interpolated overburden pressure is smaller than 0!!!")
           println("iSite=",iSite)
    end

    # Append x and y values to x_all and y_all respectively
    global x_all = vcat(x_all, x[:]);
    global y_all = vcat(y_all, y[:]);

    # Find indices of non-NaN values in the second column of zbout (measured acoustic impedance)
    fNotNan = findall(x -> !isnan(x), zbout[:, 2]);

    # Subsample indices with a step size of subSamp
    subSamp = 1; # use 1 for all data points, 10 for speed up
    fNotNanSubsampled = fNotNan[1:subSamp:end];
    
    # Compute the length of the subsampled indices
    thisDoF = length(fNotNanSubsampled);
    
    # Loop over VGS parameter space
    for iPorosity=1:nPorosity
    	for iGrainDiam=1:nGrainDiam

	    # Set parameter values
	    porosity=porosityList[iPorosity];
            grainDiameter=grainDiameterList[iGrainDiam];
            
	    # Calculate acoustic impedance based on VGS theory
            ZP_sites,ZS_sites = acousticImpedance(overburdenPressureSites,porosity,grainDiameter);			# N = Pi assumption
            ZP0_sites,ZS0_sites = acousticImpedance(zeros(size(overburdenPressureSites)),porosity,grainDiameter);	# N = 0 Pa assumption
               
	    # Calculate acoustic impedance misfit and sum it for all selected data points of this side (fNotNanSubsampled)
            thisChiSq=sum(((zbout[fNotNanSubsampled,2]-Spline1D(dist,ZP_sites; k=1)(zbout[fNotNanSubsampled,1])) ./ zbout[fNotNanSubsampled,3]).^2);
            thisChiSqZero=sum(((zbout[fNotNanSubsampled,2]-Spline1D(dist,ZP0_sites; k=1)(zbout[fNotNanSubsampled,1])) ./ zbout[fNotNanSubsampled,3]).^2);
            
	    # Store misfit
            dof[iPorosity,iGrainDiam]=dof[iPorosity,iGrainDiam]+thisDoF;
            chiSq[iPorosity,iGrainDiam]=chiSq[iPorosity,iGrainDiam]+thisChiSq;
            chiSqZero[iPorosity,iGrainDiam]=chiSqZero[iPorosity,iGrainDiam]+thisChiSqZero;
            
            global icount = icount + 1;
        end
	# Print out progress
	println(icount/ncount*100," % done")
    end
end
    
# Create NetCDF output file
NetCDFfile = "weertmanGridSearch.nc";
println("Save NetCDF output at ",NetCDFfile)
isfile(NetCDFfile) && rm(NetCDFfile);

outputVars = ["dof","chiSq","chiSqZero"];
nVars = length(outputVars);
outputData = zeros(nPorosity,nGrainDiam,nVars);
outputData[:,:,1] = dof;
outputData[:,:,2] = chiSq;
outputData[:,:,3] = chiSqZero;
for iVars=1:nVars
    nccreate(NetCDFfile, outputVars[iVars], "porosity", nPorosity, "grainDiameter", nGrainDiam);
    ncwrite(outputData[:,:,iVars], NetCDFfile, outputVars[iVars]);
end

nccreate(NetCDFfile, "porosityList", "porosity", nPorosity);
nccreate(NetCDFfile, "grainDiameterList", "grainDiameter", nGrainDiam);

ncwrite(collect(porosityList), NetCDFfile, "porosityList");
ncwrite(collect(grainDiameterList), NetCDFfile, "grainDiameterList");

println("")
ncinfo(NetCDFfile)
println("")

println("#--------------------------------------------------------------------------------------#")
println("End of program")
println("#--------------------------------------------------------------------------------------#")
println("End time: ",now())