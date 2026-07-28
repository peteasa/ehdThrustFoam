
magE = mag(E);

if (0 && enableDetailedLogs)
{
    label maxCell = findMax(magE);
    if (Pstream::master()) Info<< "cell with max magE = " << maxCell << endl;
}

Pphys = Pamb + rho * p;
Peff = Pphys * (T0TransRef / T);
Peff_safe = max(Peff, smallP);
delta = ( Pphys / P0TransRef ) * ( T0TransRef / T );
Tc = T - Tzero;
Psat = Psat0 * exp( 17.27 * Tc / (Tc + dimensionedScalar("mtA", dimTemperature, 237.3)) );
PH20 = RH * Psat;
xH20 = PH20 / Pphys;

if (enableDetailedLogs)
{
    scalar minPphys = gMin(Pphys.internalField());
    scalar maxPphys = gMax(Pphys.internalField());

    scalar minPeff = gMin(Peff.internalField());
    scalar maxPeff = gMax(Peff.internalField());
    if (Pstream::master()) Info << "min/max Pphys: " << minPphys << " " << maxPphys
                                << " min/max Peff: " << minPeff << " " << maxPeff << nl;

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
// calculate source ion term (Sion) using the Townsend coefficient alpha
alpha = A0 * Peff_safe * exp(-B0 * Peff_safe / max(magE, smallE));

// Sion is stable when scaled
Sion = ionScale * alpha * mu_e * magE * ne;

if (enableDetailedLogs)
{
    scalar minSion = gMin(Sion.internalField());
    scalar maxSion = gMax(Sion.internalField());
    if (Pstream::master()) Info << "internal min/max Sion: " << minSion << " " << maxSion << nl;
}

// attachment
// e- + O2 + M -> O2- + M
eta3b = (C1 + C2 * xH20) * sqr(Peff_safe) * exp( -D * Peff_safe / max(magE, smallE) );
eta2b = Ad * exp( - pow( ( magE * kB * T0TransRef / Peff_safe - Ed ) / DeltaE , 2) );
eta = eta3b + eta2b;
SO2ion = eta * mu_e * magE * ne;

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

    scalar minSO2 = gMin(SO2ion.internalField());
    scalar maxSO2 = gMax(SO2ion.internalField());
    if (Pstream::master()) Info << "internal min/max SO2ion: " << minSO2 << " " << maxSO2 << nl;

    //label maxCell = findMax(SO2ion);
    //if (Pstream::master()) Info<< "cell with max SO2ion = " << maxCell << nl;
}

// detachment O2- -> O2 + e-
//SO2e = kO2e * nO2m;
nuDet = Adet * Peff * exp(-Eact_eV * eCharge.value() / (kB * T));
SO2e =  nuDet * nO2m;

if (enableDetailedLogs)
{
    scalar minnuDet = gMin(nuDet.internalField());
    scalar maxnuDet = gMax(nuDet.internalField());

    scalar minSO2e = gMin(SO2e.internalField());
    scalar maxSO2e = gMax(SO2e.internalField());
    if (Pstream::master()) Info << "internal min/max SO2e: " << minSO2e << " " << maxSO2e
                                << " min/max nuDet: " << minnuDet << " " << maxnuDet
                                << nl;

    //label maxCell = findMax(SO2e);
    //if (Pstream::master()) Info<< "cell with max SO2e = " << maxCell << nl;
}

// recombination
// e- + N2+ -> N2
keiLocal = kei0 * pow(T0ChemRef / T, 0.5) * (Pphys / P0TransRef);
Rei = keiLocal * ne * nN2p;

if (enableDetailedLogs)
{
    scalar minRei = gMin(Rei.internalField());
    scalar maxRei = gMax(Rei.internalField());
    if (Pstream::master()) Info << "internal min/max Rei: " << minRei << " " << maxRei
                                << nl;

    //label maxCell = findMax(Rei);
    //if (Pstream::master()) Info<< "cell with max Rei = " << maxCell << nl;
}

// N2+ + O2- -> N2 + O2
kiiLocal = kii0 * ( T0ChemRef / T ) * ( 1.0 + humidityCorrection * xH20 );
Rii = kiiLocal * nN2p * nO2m;

if (enableDetailedLogs)
{
    scalar minRii = gMin(Rii.internalField());
    scalar maxRii = gMax(Rii.internalField());
    if (Pstream::master()) Info << "internal min/max Rii: " << minRii << " " << maxRii
                                << nl;

    //label maxCell = findMax(Rii);
    //if (Pstream::master()) Info<< "cell with max Rii = " << maxCell << nl;
}

chemReactionRate = mag(Sion) + mag(SO2ion) + mag(SO2e) + mag(Rei) + mag(Rii);
maxReactionRate = gMax(chemReactionRate);
if (Pstream::master() && enableDetailedLogs) Info << "maxReactionRate: " << maxReactionRate << nl;
