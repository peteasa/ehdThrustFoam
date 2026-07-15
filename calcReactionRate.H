
magE = mag(E);

Pphys = Pref + rho * p;
Peff = Pphys * (TC0ref / T);
delta = ( Pphys / Psearef ) * ( TC0ref / T );
Tc = T - Tzero;
Psat = Psat0 * exp( 17.27 * Tc / (Tc + dimensionedScalar("mtA", dimTemperature, 237.3)) );
PH20 = RH * Psat;
xH20 = PH20 / Pphys;

mu_e = mu_e0 * (1.0 - beta * xH20) / Peff;
mu_N2p = mu_N2p0 * (1.0 - beta * xH20) / Peff;

if (enableDetailedLogs)
{
    scalar minDelta = gMin(delta.internalField());
    scalar maxDelta = gMax(delta.internalField());

    scalar minMagE = gMin(magE.internalField());
    scalar maxMagE = gMax(magE.internalField());

    if (Pstream::master()) Info << "delta min/max: " << minDelta << " " << maxDelta
                                << " magE min/max: " << minMagE << " " << maxMagE
                                << nl;
}

// ionisation N2 -> N2+ + e-
// calculate source ion term (Sion) using the Townsend coefficient alpha
alpha = A * exp(-B / max(magE, smallE));

// Sion is stable when scaled
Sion = ionScale * alpha * mu_e * magE * ne;

if (enableDetailedLogs)
{
    scalar minSion = gMin(Sion.internalField());
    scalar maxSion = gMax(Sion.internalField());
    if (Pstream::master()) Info << "internal min/max Sion: " << minSion << " " << maxSion << endl;
}

// recombination
// e- + N2+ -> N2
Rei = kei * ne * nN2p;

if (enableDetailedLogs)
{
    scalar minRei = gMin(Rei.internalField());
    scalar maxRei = gMax(Rei.internalField());
    if (Pstream::master()) Info << "internal min/max Rei: " << minRei << " " << maxRei << endl;

    //label maxCell = findMax(Rei);
    //if (Pstream::master()) Info<< "cell with max Rei = " << maxCell << endl;
}

reactionRate = Sion + Rei;

tmp<volScalarField> magRR = mag(reactionRate);
maxReactionRate = gMax(magRR());
if (Pstream::master() && enableDetailedLogs) Info << "maxReactionRate: " << maxReactionRate << endl;
