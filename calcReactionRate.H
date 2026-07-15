
magE = mag(E);

if (enableDetailedLogs)
{
    scalar minE = gMin(magE);
    scalar maxE = gMax(magE);
    if (Pstream::master()) Info << "min/max magE: " << minE << " " << maxE << endl;

    label maxCell = findMax(magE);
    if (0 && Pstream::master()) Info<< "cell with max magE = " << maxCell << endl;
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

    label maxCell = findMax(Rei);
    if (0 && Pstream::master()) Info<< "cell with max Rei = " << maxCell << endl;
}

reactionRate = Sion + Rei;

tmp<volScalarField> magRR = mag(reactionRate);
maxReactionRate = gMax(magRR());
if (Pstream::master() && enableDetailedLogs) Info << "maxReactionRate: " << maxReactionRate << endl;
