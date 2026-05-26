
volScalarField magE
(
    "magE",
    mag(E)
);
Info << "min/max magE: " << min(magE).value() << " " << max(magE).value() << endl;
Info << "internal min/max magE: " << min(mag(E.internalField())).value() << " " << max(mag(E.internalField())).value() << endl;


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

Info << "min/max magEeff: " << min(magEeff).value() << " " << max(magEeff).value() << endl;

// Townsend coefficient
volScalarField alphaRaw
(
    "alphaRaw",
    A * Foam::exp(-B / max(magEeff, smallE))
);

alphaRaw = min(alphaRaw, dimensionedScalar("alphaMax", alphaRaw.dimensions(), 1e5));
Info << "min/max alphaRaw: " << min(alphaRaw).value() << " " << max(alphaRaw).value() << endl;

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

Info << "internal min/max Sion: " << min(Sion.internalField()).value() << " " << max(Sion.internalField()).value() << endl;

// recombination
// e- + N2+ -> N2
volScalarField Rei
(
    IOobject("Rei", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    kei * ne * nN2p
);

Info << "internal min/max Rei: " << min(Rei.internalField()).value() << " " << max(Rei.internalField()).value() << endl;

volScalarField reactionRate
(
    IOobject("reactionRate", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    Sion + Rei
);

scalar maxReactionRate = gMax(mag(reactionRate.primitiveField()));
Info << "maxReactionRate: " << maxReactionRate << endl;
