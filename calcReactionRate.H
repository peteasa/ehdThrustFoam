
magE = mag(E);

if (0 && enableDetailedLogs)
{
    label maxCell = findMax(magE);
    if (Pstream::master()) Info<< "cell with max magE = " << maxCell << endl;
}

Pphys_raw = Pamb + rho * p;
Pphys = max(Pphys_raw, smallP);
Peff_raw = Pphys_raw * (T0TransRef / T);
Peff = max(Peff_raw, smallP);
delta = ( Pphys / P0TransRef ) * ( T0TransRef / T );
Tc = T - Tzero;
Psat = Psat0 * exp( 17.27 * Tc / (Tc + dimensionedScalar("mtA", dimTemperature, 237.3)) );
PH20 = RH * Psat;
tmp<volScalarField> xH20_raw = PH20 / Pphys;
xH20 = min(xH20_raw(), 1.0);

if (enableDetailedLogs)
{
    scalar minPphys = gMin(Pphys_raw.internalField());
    scalar maxPphys = gMax(Pphys_raw.internalField());

    scalar minPeff = gMin(Peff_raw.internalField());
    scalar maxPeff = gMax(Peff_raw.internalField());
    if (Pstream::master()) Info << "min/max raw Pphys: " << minPphys << " " << maxPphys
                                << " min/max raw Peff: " << minPeff << " " << maxPeff << nl;

    scalar minxH20 = gMin(xH20.internalField());
    scalar maxxH20 = gMax(xH20.internalField());
    if (Pstream::master()) Info << "min/max xH20: " << minxH20 << " " << maxxH20 << nl;
}

mu_e = mu_e0 * (1.0 - beta * xH20) / Peff;
mu_N2p = mu_N2p0 * (1.0 - beta * xH20) / Peff;
mu_O2m = mu_O2m0 * (1.0 - beta * xH20) / Peff;

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

    scalar minEN = T0TransRef.value() * kB.value()/(P0TransRef.value() * maxdelta_magE * 1e-21);
    scalar maxEN = T0TransRef.value() * kB.value()/(P0TransRef.value() * mindelta_magE * 1e-21);

    if (Pstream::master()) Info << "delta min/max: " << minDelta << " " << maxDelta
                                << " magE min/max: " << minMagE << " " << maxMagE << nl
                                << "delta/magE min/max: " << mindelta_magE << " " << maxdelta_magE
                                << " magE/N min/max: " << minEN << " " << maxEN
                                << nl;
}

// ionisation N2 -> N2+ + e-
// calculate source ion term (RionN2p) using the Townsend coefficient alpha
alpha = A0 * Peff * exp(-B0 * Peff / max(magE, smallE));

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
eta3b = (C1 + C2 * xH20) * sqr(Peff) * exp( -D * Peff / max(magE, smallE) );
eta2b = Ad * exp( - pow( ( magE * kB * T0TransRef / Peff - Ed ) / DeltaE , 2) );
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
//RdetO2e = kO2e * nO2m;
nuDet = Adet * Peff * exp(-Eact_eV * eCharge.value() / (kB * T));
RdetO2e =  nuDet * nO2m;

if (enableDetailedLogs)
{
    scalar minnuDet = gMin(nuDet.internalField());
    scalar maxnuDet = gMax(nuDet.internalField());

    scalar minRdetO2e = gMin(RdetO2e.internalField());
    scalar maxRdetO2e = gMax(RdetO2e.internalField());
    if (Pstream::master()) Info << "internal min/max RdetO2e: " << minRdetO2e << " " << maxRdetO2e
                                << " min/max nuDet: " << minnuDet << " " << maxnuDet
                                << nl;

    //label maxCell = findMax(RdetO2e);
    //if (Pstream::master()) Info<< "cell with max RdetO2e = " << maxCell << nl;
}

// recombination
// e- + N2+ -> N2
krecN2pLocal = krecN2p0 * pow(T0ChemRef / T, 0.5) * (Pphys / P0TransRef);
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
kiiN2pLocal = kiiN2p0 * ( T0ChemRef / T ) * ( 1.0 + humidityCorrection * xH20 );
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
