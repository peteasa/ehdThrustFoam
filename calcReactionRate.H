
volScalarField magE
(
    "magE",
    mag(E)
);

if (0 && enableDetailedLogs)
{
    scalar minE = gMin(magE);
    scalar maxE = gMax(magE);
    if (Pstream::master()) Info << "min/max magE: " << minE << " " << maxE << endl;

    label maxCell = findMax(magE);
    if (Pstream::master()) Info<< "cell with max magE = " << maxCell << endl;
}

// ionisation N2 -> N2+ + e-
// calculate source ion term (Sion) using the Townsend coefficient alpha
dimensionedScalar smallE("smallE", magE.dimensions(), 1e3);
volScalarField alpha
(
    "alpha",
    A * exp(-B / max(magE, smallE))
);

dimensionedScalar ionScale("ionScale", dimless, 1e-1);

// Sion is stable when scaled
volScalarField Sion
(
    "Sion",
    ionScale * alpha * mu_e * magE * ne
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

if (0 && enableDetailedLogs)
{
    scalar minRei = gMin(Rei.internalField());
    scalar maxRei = gMax(Rei.internalField());
    if (Pstream::master()) Info << "internal min/max Rei: " << minRei << " " << maxRei << endl;

    label maxCell = findMax(Rei);
    if (Pstream::master()) Info<< "cell with max Rei = " << maxCell << endl;
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

maxReactionRate = gMax(magRR);
if (Pstream::master() && enableDetailedLogs) Info << "maxReactionRate: " << maxReactionRate << endl;
