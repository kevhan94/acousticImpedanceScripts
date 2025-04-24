function acousticImpedance(effectivePressure,porosity,grainDiameter)
    # Calculate acoustic impedance (ZP) and shear impedance (ZS) based on VGS theory
    gravity=9.81;
    densityGrains=2730;
    densityPores=1005;
    KPores=2.374e9;
    KGrains=3.6e10;
    frequencyHz = 100;

    depth=200;

    gammaP0=3.888e8;
    gammaS0=4.588e7;

    strainHardeningIndex=0.0851;

    refT=1;
    refGrainDiameter=1e-3;
    refDepth=0.3;
    refPorosity=0.377;
    refEffectivePressure=(1-refPorosity)*(densityGrains-densityPores)*gravity*refDepth;

    gammaP=gammaP0*((effectivePressure.*grainDiameter)./(refEffectivePressure*refGrainDiameter)).^(1/3);
    gammaS=gammaS0*((effectivePressure.*grainDiameter)./(refEffectivePressure*refGrainDiameter)).^(2/3);

    Density0 = porosity.*densityPores.+(1 .-porosity).*densityGrains;
    K0 = 1 ./(porosity./KPores+(1 .-porosity)./KGrains);

    c0=sqrt.(K0./Density0);

    tauP=0.012; # s 
    poreFluidViscosity = (1. + 1. / (im*2*pi*frequencyHz*tauP))^(-1+strainHardeningIndex);

    cP=c0./real.((1 .+((gammaP.+(4/3).*gammaS)./(Density0.*c0.^2)).*(im      *2*pi*frequencyHz*refT).^strainHardeningIndex .*poreFluidViscosity).^(-0.5));
    #alphaP=-(2*pi*frequencyHz./c0).*imag((1 .+((gammaP.+(4/3).*gammaS)./(Density0.*c0.^2)).*(im      *2*pi*frequencyHz*refT).^strainHardeningIndex .*poreFluidViscosity).^(-0.5));

    # Buckingham (2014), tauS>>tauP, so lim tauS -> inf, so poreFluidViscosity(tauS)=1, => no effect at all on the shear wave
    cS=sqrt.(gammaS./Density0)*((2*pi*frequencyHz*refT)^(strainHardeningIndex/2))/cos(strainHardeningIndex*pi/4);
    #alphaS=2*pi*frequencyHz*tan(strainHardeningIndex*pi/4)./cS;

    ZP=cP.*Density0;
    ZS=cS.*Density0;

    return ZP,ZS    
end