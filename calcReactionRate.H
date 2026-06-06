
volScalarField magE
(
    "magE",
    mag(E)
);

if (enableDetailedLogs)
{
    scalar minE = gMin(magE);
    scalar maxE = gMax(magE);
    if (Pstream::master()) Info << "min/max magE: " << minE << " " << maxE << endl;
}

// ionisation N2 -> N2+ + e-
// calculate source ion term (Sion) using the Townsend coefficient alpha
// alpha = A * exp(-B / magE);
// volScalarField Sion = alpha * mu_e * magE * ne;

dimensionedScalar smallE("smallE", magE.dimensions(), 1e3);
dimensionedScalar Ecrit("Ecrit", magE.dimensions(), 2.5e7);
dimensionedScalar Smax("Smax", dimless/dimVolume/dimTime, 1e20);
dimensionedScalar neMax("neMax", ne.dimensions(), 1e5); // reducing to 1e16 makes no difference!
dimensionedScalar ionScale("ionScale", dimless, 1e-4);

// bounded electric field response has a big effect on ne nN2p nO2m but has no effect on phiE and E
volVectorField Eeff
(
    "Eeff",
    E / (1.0 + magE/Ecrit)
);

volScalarField magEeff
(
    "magEeff",
    mag(Eeff)
);

if (enableDetailedLogs)
{
    scalar minEeff = gMin(magEeff);
    scalar maxEeff = gMax(magEeff);
    if (Pstream::master()) Info << "min/max magEeff: " << minEeff << " " << maxEeff << endl;
}

// Townsend coefficient
volScalarField alphaRaw
(
    "alphaRaw",
    A * Foam::exp(-B / max(magEeff, smallE))
);

alphaRaw = min(alphaRaw, dimensionedScalar("alphaMax", alphaRaw.dimensions(), 1e5));
if (enableDetailedLogs)
{
    scalar minAlphaRaw = gMin(alphaRaw);
    scalar maxAlphaRaw = gMax(alphaRaw);
    if (Pstream::master()) Info << "min/max alphaRaw: " << minAlphaRaw << " " << maxAlphaRaw << endl;

    label maxCell = findMax(alphaRaw);
    if (Pstream::master()) Info<< "cell with max alphaRaw = " << maxCell << endl;
}

volScalarField ionSwitch  // ionSwitch may not be needed
(
    "ionSwitch",
    magE / (magE + Ecrit)
);

// base ionization source
volScalarField SionRaw
(
    "SionRaw",
    ionScale * alphaRaw * mu_e * magEeff * ne // magEeff may not be needed
);

volScalarField Sion
(
    "Sion",
    min(SionRaw * ionSwitch, Smax)
);

if (enableDetailedLogs)
{
    scalar minSion = gMin(Sion.internalField());
    scalar maxSion = gMax(Sion.internalField());
    if (Pstream::master()) Info << "internal min/max Sion: " << minSion << " " << maxSion << endl;
}

// recombination
// e- + N2+ -> N2
volScalarField Rei
(
    IOobject("Rei", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    kei * ne * nN2p
);

if (enableDetailedLogs)
{
    scalar minRei = gMin(Rei.internalField());
    scalar maxRei = gMax(Rei.internalField());
    if (Pstream::master()) Info << "internal min/max Rei: " << minRei << " " << maxRei << endl;
}

volScalarField reactionRate
(
    IOobject("reactionRate", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    Sion + Rei
);

volScalarField magRR
(
     "magRR",
     mag(reactionRate)
);

scalar maxReactionRate = gMax(magRR);
if (Pstream::master() && enableDetailedLogs) Info << "maxReactionRate: " << maxReactionRate << endl;
