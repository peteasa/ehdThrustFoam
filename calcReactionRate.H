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

        if ((startupLPFrac < runLPFrac) && enableDetailedLogs && Pstream::master())
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
nH2O = PH2O / (kB * T);
nO2 = xO2 * ( Pphys - PH2O) / (kB * T);
nN2 = xN2 * ( Pphys - PH2O) / (kB * T);
nM = nN2 + nO2 + nH2O;

if (enableDetailedLogs)
{
    scalar minPphys = gMin(Pphys_raw.internalField());
    scalar maxPphys = gMax(Pphys_raw.internalField());

    scalar minPeff = gMin(Peff.internalField());
    scalar maxPeff = gMax(Peff.internalField());
    if (Pstream::master()) Info << "min/max raw Pphys: " << minPphys << " " << maxPphys
                                << " min/max frac Pamb: " << minPphys / Pamb.value() << " " << maxPphys / Pamb.value()
                                << " min/max Peff: " << minPeff << " " << maxPeff << nl;

    scalar minT = gMin(T.internalField());
    scalar maxT = gMax(T.internalField());
    scalar minxH2O = gMin(xH2O.internalField());
    scalar maxxH2O = gMax(xH2O.internalField());
    if (Pstream::master()) Info << "min/max T: " << minT << " " << maxT
                                << " min/max xH2O: " << minxH2O << " " << maxxH2O << nl;
}

magE_NTdOld = magE_NTd;

magE = mag(E);
magE_NTd = magE * kB * T0TransRef * 1e21 / Peff;

if (0 && enableDetailedLogs)
{
    label maxCell = findMax(magE);
    if (Pstream::master()) Info<< "cell with max magE = " << maxCell << endl;
}

if (1)
{
    scalar minENTd = gMin(magE_NTd);
    scalar maxENTd = gMax(magE_NTd);

    tmp<volScalarField> magDENTdDt = mag(magE_NTd - magE_NTdOld);
    scalar maxDENTdDt = gMax(magDENTdDt());
    maxDENTdDt /= runTime.deltaTValue();
    maxDENTdDtRate = maxDENTdDt / maxENTd;

    if (enableDetailedLogs && Pstream::master()) Info << "min/max magE/N: " << minENTd << " " << maxENTd
                                                << " maxDENTdDt: " << maxDENTdDt
                                                << " maxDENTdDtRate: " << maxDENTdDtRate
                                                << nl;
}

mu_e = mu_e0 * (1.0 - beta * xH2O) / Peff;
mu_N2p = mu_N2p0 * (1.0 - beta * xH2O) / Peff;
mu_O4p = mu_O4p0 * (1.0 - beta * xH2O) / Peff;
mu_O2p = mu_O2p0 * (1.0 - beta * xH2O) / Peff;
mu_O2m = mu_O2m0 * (1.0 - beta * xH2O) / Peff;

if (enableDetailedLogs)
{
    scalar minmu_e = gMin(mu_e.internalField());
    scalar maxmu_e = gMax(mu_e.internalField());
    scalar minmu_N2p = gMin(mu_N2p.internalField());
    scalar maxmu_N2p = gMax(mu_N2p.internalField());
    scalar minmu_O4p = gMin(mu_O4p.internalField());
    scalar maxmu_O4p = gMax(mu_O4p.internalField());
    scalar minmu_O2p = gMin(mu_O2p.internalField());
    scalar maxmu_O2p = gMax(mu_O2p.internalField());
    scalar minmu_O2m = gMin(mu_O2m.internalField());
    scalar maxmu_O2m = gMax(mu_O2m.internalField());
    if (Pstream::master()) Info << "min/max mu_e: " << minmu_e << " " << maxmu_e
                                << " min/max mu_N2p: " << minmu_N2p << " " << maxmu_N2p
                                << " min/max mu_O2m: " << minmu_O2m << " " << maxmu_O2m
                                << nl << "min/max mu_O4p: " << minmu_O4p << " " << maxmu_O4p
                                << " min/max mu_O2p: " << minmu_O2p << " " << maxmu_O2p
                                << nl;

    scalar minDelta = gMin(delta.internalField());
    scalar maxDelta = gMax(delta.internalField());

    scalar minMagE = gMin(magE.internalField());
    scalar maxMagE = gMax(magE.internalField());

    tmp<volScalarField> delta_magE = delta / max(magE, smallE);
    scalar mindelta_magE = gMin(delta_magE());
    scalar maxdelta_magE = gMax(delta_magE());

    if (Pstream::master()) Info << "delta min/max: " << minDelta << " " << maxDelta
                                << " magE min/max: " << minMagE << " " << maxMagE
                                << nl << "delta/magE min/max: " << mindelta_magE << " " << maxdelta_magE
                                << nl;
}

// ionisation N2 -> N2+ + e-
// calculate source ion term (RionN2p) using the Townsend coefficient alpha
alphaN2 = AN20 * Peff * exp(-BN20 * Peff / max(magE, smallE));

// RionN2p is stable when scaled
RionN2p = ionScale * alphaN2 * mu_e * magE * ne;

// ionisation O2 -> O2+ + e-
// calculate source ion term (RionO2p) using the Townsend coefficient alpha
alphaO2 = AO20 * Peff * exp(-BO20 * Peff / max(magE, smallE));

// RionO2p is stable when scaled
RionO2p = ionScale * alphaO2 * mu_e * magE * ne;

if (enableDetailedLogs)
{
    scalar minRionN2p = gMin(RionN2p.internalField());
    scalar maxRionN2p = gMax(RionN2p.internalField());
    scalar minRionO2p = gMin(RionO2p.internalField());
    scalar maxRionO2p = gMax(RionO2p.internalField());
    if (Pstream::master()) Info << "min/max RionN2p: " << minRionN2p << " " << maxRionN2p
                                << " min/max RionO2p: " << minRionO2p << " " << maxRionO2p
                                << nl;
}

// charge transfer O2 + N2+ -> O2+ + N2
kCTLocal = kCT0 * pow(T0ChemRef / T, 0.5);
RCT = kCTLocal * nN2p * nO2;

kcluster04pLocal = kcluster04p0 * pow(T0ChemRef / T, 3.2);
RclustO4p = kcluster04pLocal * nO2p * nO2 * nM;

if (enableDetailedLogs)
{
    scalar minRCT = gMin(RCT.internalField());
    scalar maxRCT = gMax(RCT.internalField());
    scalar minRclustO4p = gMin(RclustO4p.internalField());
    scalar maxRclustO4p = gMax(RclustO4p.internalField());
    if (Pstream::master()) Info << "min/max RCT: " << minRCT << " " << maxRCT
                                << " min/max RclustO4p: " << minRclustO4p << " " << maxRclustO4p
                                << nl;
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
    if (Pstream::master()) Info << "min/max RattO2m: " << minSO2 << " " << maxSO2 << nl;
}

// detachment O2- + O2 -> 2O2 + e-
// detachment O2- + N2 -> N2 + O2 + e-
kdetLocal = kdetA0 * exp( - pow(kdetB0 / (kdetC0 + magE_NTd), 2) );
Rdet = kdetLocal * nO2m * (nO2 + nN2);

if (enableDetailedLogs)
{
    scalar minRdet = gMin(Rdet.internalField());
    scalar maxRdet = gMax(Rdet.internalField());
    if (Pstream::master()) Info << "min/max Rdet: " << minRdet << " " << maxRdet
                                << nl;
}

// recombination
// e- + O4+ -> 2O2
krecO4pLocal = krecO4p0 * pow(T0ChemRef / TeK, 0.48);
RrecO4p = krecO4pLocal * ne * nO4p;

// e- + O2+ -> O2
krecO2pLocal = krecO2p0 * pow(T0ChemRef / TeK, 0.7);
RrecO2p = krecO2pLocal * ne * nO2p;

// e- + N2+ -> N2
krecN2pLocal = krecN2p0 * pow(T0ChemRef / TeK, 0.39);
RrecN2p = krecN2pLocal * ne * nN2p;

if (enableDetailedLogs)
{
    scalar minRrecO4p = gMin(RrecO4p.internalField());
    scalar maxRrecO4p = gMax(RrecO4p.internalField());

    scalar minRrecO2p = gMin(RrecO2p.internalField());
    scalar maxRrecO2p = gMax(RrecO2p.internalField());

    scalar minRrecN2p = gMin(RrecN2p.internalField());
    scalar maxRrecN2p = gMax(RrecN2p.internalField());
    if (Pstream::master()) Info << "min/max RrecN2p: " << minRrecN2p << " " << maxRrecN2p
                                << " min/max RrecO4p: " << minRrecO4p << " " << maxRrecO4p
                                << " min/max RrecO2p: " << minRrecO2p << " " << maxRrecO2p
                                << nl;
}

// O4+ + O2- -> 3O2
kiiO4pLocal = kiiO4p0 * pow(T0ChemRef / T, 0.5);
RiiO4p = kiiO4pLocal * nO4p * nO2m;

// O2+ + O2- -> 2O2
kiiO2pLocal = kiiO2p0 * pow(T0ChemRef / T, 0.5);
RiiO2p = kiiO2pLocal * nO2p * nO2m;

// N2+ + O2- -> N2 + O2
kiiN2pLocal = kiiN2p0 * pow(T0ChemRef / T, 0.5);
RiiN2p = kiiN2pLocal * nN2p * nO2m;

if (enableDetailedLogs)
{
    scalar minRiiO4p = gMin(RiiO4p.internalField());
    scalar maxRiiO4p = gMax(RiiO4p.internalField());

    scalar minRiiO2p = gMin(RiiO2p.internalField());
    scalar maxRiiO2p = gMax(RiiO2p.internalField());

    scalar minRiiN2p = gMin(RiiN2p.internalField());
    scalar maxRiiN2p = gMax(RiiN2p.internalField());
    if (Pstream::master()) Info << "min/max RiiO4p: " << minRiiO4p << " " << maxRiiO4p
                                << " min/max RiiO2p: " << minRiiO2p << " " << maxRiiO2p
                                << " min/max RiiN2p: " << minRiiN2p << " " << maxRiiN2p
                                << nl;

    //label maxCell = findMax(RiiN2p);
    //if (Pstream::master()) Info<< "cell with max RiiN2p = " << maxCell << nl;
}

chemReactionRateControl = mag(RionN2p) + mag(RionO2p) + mag(Rdet) + mag(RattO2m) + mag(RrecO4p);
maxReactionRate = gMax(chemReactionRateControl);
chemReactionRateFull = mag(RionN2p) + mag(RionO2p) +
    mag(RattO2m) + mag(Rdet) +
    mag(RrecO4p) + mag(RrecO2p) + mag(RrecN2p) + mag(RrecDissO) + mag(RrecDissN) +
    mag(RiiO4p) + mag(RiiO2p) + mag(RiiN2p);

scalar maxChemReactionRateFull = gMax(chemReactionRateFull);
if (Pstream::master() && enableDetailedLogs) Info << "maxReactionRate: " << maxReactionRate
                                                  << " maxChemReactionRateFull: " << maxChemReactionRateFull
                                                  << nl;
