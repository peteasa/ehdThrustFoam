Pphys_raw = Pamb + rho * p;

if (startupT)
{
    startupLPFrac = startupPHalfwidth;
    startupHPFrac = startupPHalfwidth;
    if (runTime.value() < startupT)
    {
        if (enableDetailedLogs && Pstream::master())
            Info << "Slow Start: startupLPFrac: " << startupLPFrac << " until time: " << startupT << nl;
        pStartup = runTime.timeIndex();
    }
    else
    {
        dimensionedScalar minPphys = gMin(Pphys_raw.internalField());
        dimensionedScalar maxPphys = gMax(Pphys_raw.internalField());
        if ((minPphys < Pamb * (1.0 - startupPHalfwidth - runLPFrac)/2) ||
            (Pamb * (1.5 + startupPHalfwidth + runHPFrac) < maxPphys))
        {
            pStartup = runTime.timeIndex();
        }

        int count = runTime.timeIndex() - pStartup;
        scalar pScale = 1.0;
        if (count < startupPMaxCount)
        {
            pScale = count / startupPMaxCount;
        }

        startupLPFrac += runLPFrac * pScale;
        startupHPFrac += runHPFrac * pScale;

        if ((startupLPFrac < 0.98) && enableDetailedLogs && Pstream::master())
            Info << "startupLPFrac: " << startupLPFrac << nl;
    }

    smallP = Pamb * (1.0 - startupLPFrac);
    highP = Pamb * (1.0 + startupHPFrac);
}

Pphys = max(Pphys_raw, smallP);
Pphys = min(Pphys, highP);
Peff = Pphys * (T0TransRef / T);
delta = ( Pphys / P0TransRef ) * ( T0TransRef / T );
Tc = T - Tzero;
Psat = Psat0 * exp( 17.27 * Tc / (Tc + dimensionedScalar("mtA", dimTemperature, 237.3)) );
PH2O = RH * Psat;
tmp<volScalarField> xH2O_raw = PH2O / Pphys;
xH2O = min(xH2O_raw(), 1.0);

if (enableDetailedLogs)
{
    scalar minPphys = gMin(Pphys_raw.internalField());
    scalar maxPphys = gMax(Pphys_raw.internalField());

    scalar minPeff = gMin(Peff.internalField());
    scalar maxPeff = gMax(Peff.internalField());
    if (Pstream::master()) Info << "min/max raw Pphys: " << minPphys << " " << maxPphys
                                << " min/max Peff: " << minPeff << " " << maxPeff << nl;

    scalar minT = gMin(T.internalField());
    scalar maxT = gMax(T.internalField());
    scalar minxH2O = gMin(xH2O.internalField());
    scalar maxxH2O = gMax(xH2O.internalField());
    if (Pstream::master()) Info << "min/max T: " << minT << " " << maxT
                                << " min/max xH2O: " << minxH2O << " " << maxxH2O << nl;
}

magE = mag(E);
magE_NTd = magE * kB * T0TransRef * 1e21 / Peff;

if (0 && enableDetailedLogs)
{
    label maxCell = findMax(magE);
    if (Pstream::master()) Info<< "cell with max magE = " << maxCell << endl;
}

mu_e = mu_e0 * (1.0 - beta * xH2O) / Peff;
mu_N2p = mu_N2p0 * (1.0 - beta * xH2O) / Peff;
mu_O2m = mu_O2m0 * (1.0 - beta * xH2O) / Peff;

if (enableDetailedLogs)
{
    scalar minmu_e = gMin(mu_e.internalField());
    scalar maxmu_e = gMax(mu_e.internalField());
    scalar minmu_N2p = gMin(mu_N2p.internalField());
    scalar maxmu_N2p = gMax(mu_N2p.internalField());
    scalar minmu_O2m = gMin(mu_O2m.internalField());
    scalar maxmu_O2m = gMax(mu_O2m.internalField());
    if (Pstream::master()) Info << "min/max mu_e: " << minmu_e << " " << maxmu_e
                                << " min/max mu_N2p: " << minmu_N2p << " " << maxmu_N2p
                                << " min/max mu_O2m: " << minmu_O2m << " " << maxmu_O2m
                                << nl;

    scalar minDelta = gMin(delta.internalField());
    scalar maxDelta = gMax(delta.internalField());

    scalar minMagE = gMin(magE.internalField());
    scalar maxMagE = gMax(magE.internalField());

    tmp<volScalarField> delta_magE = delta / max(magE, smallE);
    scalar mindelta_magE = gMin(delta_magE());
    scalar maxdelta_magE = gMax(delta_magE());

    scalar minEN = gMin(magE_NTd.internalField());
    scalar maxEN = gMax(magE_NTd.internalField());

    if (Pstream::master()) Info << "delta min/max: " << minDelta << " " << maxDelta
                                << " magE min/max: " << minMagE << " " << maxMagE
                                << nl << "delta/magE min/max: " << mindelta_magE << " " << maxdelta_magE
                                << " magE/N min/max: " << minEN << " " << maxEN
                                << nl;
}

// ionisation N2 -> N2+ + e-
// calculate source ion term (RionN2p) using the Townsend coefficient alpha
alpha = AN20 * Peff * exp(-BN20 * Peff / max(magE, smallE));

// RionN2p is stable when scaled
RionN2p = ionScale * alpha * mu_e * magE * ne;

if (enableDetailedLogs)
{
    scalar minRionN2p = gMin(RionN2p.internalField());
    scalar maxRionN2p = gMax(RionN2p.internalField());
    if (Pstream::master()) Info << "internal min/max RionN2p: " << minRionN2p << " " << maxRionN2p << nl;
}

// attachment
// e- + O2 + M -> O2- + M
eta3b = (C1 + C2 * xH2O) * sqr(Peff) * exp( -D * Peff / max(magE, smallE) );
eta2b = Ad * exp( - pow( ( magE_NTd - EdTd ) / DeltaETd , 2) );
eta = eta3b + eta2b;
RattO2m = eta * mu_e * magE * ne;

if (enableDetailedLogs)
{
    scalar mineta3b = gMin(eta3b);
    scalar maxeta3b = gMax(eta3b);
    scalar mineta2b = gMin(eta2b);
    scalar maxeta2b = gMax(eta2b);
    scalar mineta = gMin(eta);
    scalar maxeta = gMax(eta);
    if (Pstream::master()) Info << "eta3b min/max: " << mineta3b << " " << maxeta3b
                                << " eta2b min/max: " << mineta2b << " " << maxeta2b
                                << " eta min/max: " << mineta << " " << maxeta
                                << nl;

    scalar minSO2 = gMin(RattO2m.internalField());
    scalar maxSO2 = gMax(RattO2m.internalField());
    if (Pstream::master()) Info << "internal min/max RattO2m: " << minSO2 << " " << maxSO2 << nl;

    //label maxCell = findMax(RattO2m);
    //if (Pstream::master()) Info<< "cell with max RattO2m = " << maxCell << nl;
}

// detachment O2- -> O2 + e-
//RdetO2e = kdetO2e * nO2m;
nuDetO2 = AdetO2O2m * Peff * exp(-Eact_eV * eCharge.value() / (kB * T));
RdetO2e =  nuDetO2 * nO2m;

if (enableDetailedLogs)
{
    scalar minnuDetO2 = gMin(nuDetO2.internalField());
    scalar maxnuDetO2 = gMax(nuDetO2.internalField());

    scalar minRdetO2e = gMin(RdetO2e.internalField());
    scalar maxRdetO2e = gMax(RdetO2e.internalField());
    if (Pstream::master()) Info << "internal min/max RdetO2e: " << minRdetO2e << " " << maxRdetO2e
                                << " min/max nuDetO2: " << minnuDetO2 << " " << maxnuDetO2
                                << nl;

    //label maxCell = findMax(RdetO2e);
    //if (Pstream::master()) Info<< "cell with max RdetO2e = " << maxCell << nl;
}

// recombination
// e- + N2+ -> N2
krecN2pLocal = krecN2p0 * pow(T0ChemRef / TeK, 0.39);
RrecN2p = krecN2pLocal * ne * nN2p;

if (enableDetailedLogs)
{
    scalar minRrecN2p = gMin(RrecN2p.internalField());
    scalar maxRrecN2p = gMax(RrecN2p.internalField());
    if (Pstream::master()) Info << "internal min/max RrecN2p: " << minRrecN2p << " " << maxRrecN2p
                                << nl;

    //label maxCell = findMax(RrecN2p);
    //if (Pstream::master()) Info<< "cell with max RrecN2p = " << maxCell << nl;
}

// N2+ + O2- -> N2 + O2
kiiN2pLocal = kiiN2p0 * pow( T0ChemRef / T, 0.5 );
RiiN2p = kiiN2pLocal * nN2p * nO2m;

if (enableDetailedLogs)
{
    scalar minRiiN2p = gMin(RiiN2p.internalField());
    scalar maxRiiN2p = gMax(RiiN2p.internalField());
    if (Pstream::master()) Info << "internal min/max RiiN2p: " << minRiiN2p << " " << maxRiiN2p
                                << nl;

    //label maxCell = findMax(RiiN2p);
    //if (Pstream::master()) Info<< "cell with max RiiN2p = " << maxCell << nl;
}

chemReactionRate = mag(RionN2p) + mag(RattO2m) + mag(RdetO2e) + mag(RrecN2p) + mag(RiiN2p);
maxReactionRate = gMax(chemReactionRate);
if (Pstream::master() && enableDetailedLogs) Info << "maxReactionRate: " << maxReactionRate << nl;
