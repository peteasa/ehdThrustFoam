if (Pstream::master()) Info<< "Creating plasma species\n" << nl;
volScalarField ne
(
    IOobject("ne", runTime.name(), mesh,
             IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);

volScalarField nN2p
(
    IOobject("nN2p", runTime.name(), mesh,
             IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);

volScalarField nO4p
(
    IOobject("nO4p", runTime.name(), mesh,
             IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);

volScalarField nO2p
(
    IOobject("nO2p", runTime.name(), mesh,
             IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);

volScalarField nO2m
(
    IOobject("nO2m", runTime.name(), mesh,
             IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);

volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    eCharge * ( nN2p + nO4p + nO2p - nO2m - ne )
);

dimensionedScalar beta
(
    physicalProperties.lookup("beta")
);

dimensionSet dimVoltage(1, 2, -3, 0, 0, -1, 0);

// dimArea / dimVoltage / dimTime
dimensionedScalar mu_e0
(
    physicalProperties.lookup("mu_e0")
);

volScalarField mu_e
(
    IOobject ( "mu_e", runTime.name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
    mu_e0 * (1.0 - beta * xH2O) / Peff
);

dimensionedScalar mu_N2p0
(
    physicalProperties.lookup("mu_N2p0")
);

volScalarField mu_N2p
(
    IOobject ( "mu_N2p", runTime.name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
    mu_N2p0 * (1.0 - beta * xH2O) / Peff
);

dimensionedScalar mu_O4p0
(
    physicalProperties.lookup("mu_O4p0")
);

volScalarField mu_O4p
(
    IOobject ( "mu_O4p", runTime.name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
    mu_O4p0 * (1.0 - beta * xH2O) / Peff
);

dimensionedScalar mu_O2p0
(
    physicalProperties.lookup("mu_O2p0")
);

volScalarField mu_O2p
(
    IOobject ( "mu_O2p", runTime.name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
    mu_O2p0 * (1.0 - beta * xH2O) / Peff
);

dimensionedScalar mu_O2m0
(
    physicalProperties.lookup("mu_O2m0")
);

volScalarField mu_O2m
(
    IOobject ( "mu_O2m", runTime.name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
    mu_O2m0 * (1.0 - beta * xH2O) / Peff
);

// 2 eV dimTemperature
dimensionedScalar Te
(
    physicalProperties.lookup("Te")
);

dimensionedScalar TeK = Te * eCharge / kB;

volScalarField De
(
    IOobject ( "De", runTime.name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
    mu_e * Te
);

volScalarField DN2p
(
    IOobject ( "DN2p", runTime.name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
    mu_N2p * Te
);

volScalarField DO4p
(
    IOobject ( "DO4p", runTime.name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
    mu_O4p * Te
);

volScalarField DO2p
(
    IOobject ( "DO2p", runTime.name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
    mu_O2p * Te
);

volScalarField DO2m
(
    IOobject ( "DO2m", runTime.name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
    mu_O2m * Te
);

// calculate particle fluxes
volVectorField GammaE
(
    IOobject ( "GammaE", runTime.name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
    U * ne
  - mu_e * E * ne
  - De * fvc::grad(ne)
);

volVectorField GammaN2p
(
    IOobject ( "GammaN2p", runTime.name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
    U * nN2p
  + mu_N2p * E * nN2p
  - DN2p * fvc::grad(nN2p)
);

volVectorField GammaO4p
(
    IOobject ( "GammaO4p", runTime.name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
    U * nO4p
  + mu_O4p * E * nO4p
  - DO4p * fvc::grad(nO4p)
);

volVectorField GammaO2p
(
    IOobject ( "GammaO2p", runTime.name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
    U * nO2p
  + mu_O2p * E * nO2p
  - DO2p * fvc::grad(nO2p)
);

volVectorField GammaO2m
(
    IOobject ( "GammaO2m", runTime.name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
    U * nO2m
  - mu_O2m * E * nO2m
  - DO2m * fvc::grad(nO2m)
);

volVectorField Je
(
    IOobject ( "Je", runTime.name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
  - eCharge * GammaE
);

volVectorField JN2p
(
    IOobject ( "JN2p", runTime.name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
    eCharge * GammaN2p
);

volVectorField JO4p
(
    IOobject ( "JO4p", runTime.name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
    eCharge * GammaO4p
);

volVectorField JO2p
(
    IOobject ( "JO2p", runTime.name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
    eCharge * GammaO2p
);

volVectorField JO2m
(
    IOobject ( "JO2m", runTime.name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
  - eCharge * GammaO2m
);

volVectorField J
(
    IOobject
    (
        "J",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Je + JN2p + JO4p + JO2p + JO2m
);

dimensionedScalar E0("E0", E.dimensions(), 2.5e6);
volScalarField magE
(
    IOobject("magE", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mag(E)
);

volScalarField magE_NTd
(
    IOobject("magE_NTd", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    magE * kB * T0TransRef * 1e21 / Peff
);

// reference temperature for chemical reactions
dimensionedScalar T0ChemRef
(
    physicalProperties.lookup("T0ChemRef")
);

// ionisation - for Townsend ionisation N2 -> N2+ + e-
dimensionedScalar AN2T0
(
    physicalProperties.lookup("AN2T0")
);

dimensionedScalar AN20 = AN2T0 * dimensionedScalar("ToR_T0", dimTemperature, 298.15) / T0TransRef;

dimensionedScalar BN2T0
(
    physicalProperties.lookup("BN2T0")
);

dimensionedScalar BN20 = BN2T0 * dimensionedScalar("ToR_T0", dimTemperature, 298.15) / T0TransRef;

dimensionedScalar smallE("smallE", magE.dimensions(), 1e3);
volScalarField alphaN2
(
    IOobject("alphaN2", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    AN20 * Peff * exp(-BN20 * Peff / max(magE, smallE))
);

// ionisation - for Townsend ionisation O2 -> O2+ + e-
dimensionedScalar AO2T0
(
    physicalProperties.lookup("AO2T0")
);

dimensionedScalar AO20 = AO2T0 * dimensionedScalar("ToR_T0", dimTemperature, 298.15) / T0TransRef;

dimensionedScalar BO2T0
(
    physicalProperties.lookup("BO2T0")
);

dimensionedScalar BO20 = BO2T0 * dimensionedScalar("ToR_T0", dimTemperature, 298.15) / T0TransRef;

volScalarField alphaO2
(
    IOobject("alphaO2", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    AO20 * Peff * exp(-BO20 * Peff / max(magE, smallE))
);

// charge transfer
dimensionedScalar kCT0
(
    physicalProperties.lookup("kCT0")
);

// O2+ clustering
dimensionedScalar kcluster04p0
(
    physicalProperties.lookup("kcluster04p0")
);

// attachment e- + O2 + M -> O2- + M
dimensionedScalar C1
(
    physicalProperties.lookup("C1")
);

dimensionedScalar C2
(
    physicalProperties.lookup("C2")
);

dimensionedScalar D
(
    physicalProperties.lookup("D")
);

dimensionedScalar Ad
(
    physicalProperties.lookup("Ad")
);

dimensionedScalar EdTd
(
    physicalProperties.lookup("EdTd")
);

dimensionedScalar DeltaETd
(
    physicalProperties.lookup("DeltaETd")
);

// detachment
dimensionedScalar kdetA0
(
    physicalProperties.lookup("kdetA0")
);

dimensionedScalar kdetB0
(
    physicalProperties.lookup("kdetB0")
);

dimensionedScalar kdetC0
(
    physicalProperties.lookup("kdetC0")
);

// recombination e- + O4+ -> 2O2
dimensionedScalar krecO4p0
(
    physicalProperties.lookup("krecO4p0")
);

// recombination e- + O2+ -> O2
dimensionedScalar krecO2p0
(
    physicalProperties.lookup("krecO2p0")
);

// recombination e- + N2+ -> N2
dimensionedScalar krecN2p0
(
    physicalProperties.lookup("krecN2p0")
);

// dissociative recombination e- + O2+ -> 2O
dimensionedScalar krecDissO2p0
(
    physicalProperties.lookup("krecDissO2p0")
);

// dissociative recombination e- + N2+ -> 2N
dimensionedScalar krecDissN2p0
(
    physicalProperties.lookup("krecDissN2p0")
);

// for recombination O4+ + O2- -> 3O2
dimensionedScalar kiiO4p0
(
    physicalProperties.lookup("kiiO4p0")
);

// for recombination O2+ + O2- -> N2 + O2
dimensionedScalar kiiO2p0
(
    physicalProperties.lookup("kiiO2p0")
);

// for recombination N2+ + O2- -> N2 + O2
dimensionedScalar kiiN2p0
(
    physicalProperties.lookup("kiiN2p0")
);

dimensionedScalar ionScale("ionScale", dimless, 1e-1);        // 1e-4
volScalarField RionN2p
(
    IOobject("RionN2p", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    ionScale * alphaN2 * mu_e * magE * ne
);

volScalarField RionO2p
(
    IOobject("RionO2p", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    ionScale * alphaO2 * mu_e * magE * ne
);

volScalarField kCTLocal
(
    IOobject("kCTLocal", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    kCT0 * pow(T0ChemRef / T, 0.5)
);

volScalarField RCT
(
    IOobject("RCT", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    kCTLocal * nN2p * nO2
);

volScalarField kcluster04pLocal
(
    IOobject("kcluster04pLocal", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    kcluster04p0 * pow(T0ChemRef / T, 3.2)
);

volScalarField RclustO4p
(
    IOobject("RclustO4p", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    kcluster04p0 * nO2p * nO2 * nM
);

volScalarField eta3b
(
    IOobject("eta3b", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    (C1 + C2 * xH2O) * sqr(Peff) * exp( -D * Peff / max(magE, smallE))
);

volScalarField eta2b
(
    IOobject("eta2b", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    Ad * exp( - pow( ( magE_NTd - EdTd ) / DeltaETd , 2) )
);

volScalarField eta
(
    IOobject("eta", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    eta3b + eta2b
);

volScalarField RattO2m
(
    IOobject("RattO2m", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    eta * mu_e * magE * ne
);

volScalarField kdetLocal
(
    IOobject("kdetLocal", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    kdetA0 * exp( - pow(kdetB0 / (kdetC0 + magE_NTd), 2) )
);

volScalarField Rdet
(
    IOobject("Rdet", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    kdetLocal * nO2m * (nO2 + nN2)
);

dimensionedScalar krecO4pLocal = krecO4p0 * pow(T0ChemRef / TeK, 0.48);

volScalarField RrecO4p
(
    IOobject("RrecO4p", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    krecO4pLocal * ne * nO4p
);

dimensionedScalar krecO2pLocal = krecO2p0 * pow(T0ChemRef / TeK, 0.7);

volScalarField RrecO2p
(
    IOobject("RrecO2p", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    krecO2pLocal * ne * nO2p
);

dimensionedScalar krecN2pLocal = krecN2p0 * pow(T0ChemRef / TeK, 0.39);

volScalarField RrecN2p
(
    IOobject("RrecN2p", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    krecN2pLocal * ne * nN2p
);

volScalarField krecDissO2pLocal
(
    IOobject("krecDissO2pLocal", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    krecDissO2p0 * pow(T0ChemRef / T, 0.5) * (Pphys / P0TransRef)
);

volScalarField RrecDissO
(
    IOobject("RrecDissO", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    krecDissO2pLocal * ne * nO2p
);

volScalarField krecDissN2pLocal
(
    IOobject("krecDissN2pLocal", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    krecDissN2p0 * pow(T0ChemRef / T, 0.5) * (Pphys / P0TransRef)
);

volScalarField RrecDissN
(
    IOobject("RrecDissN", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    krecDissN2pLocal * ne * nN2p
);

volScalarField kiiO4pLocal
(
    IOobject("kiiO4pLocal", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    kiiO4p0 * pow(T0ChemRef / T, 0.5)
);

volScalarField RiiO4p
(
    IOobject("RiiO4p", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    kiiO4pLocal * nO4p * nO2m
);

volScalarField kiiO2pLocal
(
    IOobject("kiiO2pLocal", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    kiiO2p0 * pow(T0ChemRef / T, 0.5)
);

volScalarField RiiO2p
(
    IOobject("RiiO2p", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    kiiO2pLocal * nO2p * nO2m
);

volScalarField kiiN2pLocal
(
    IOobject("kiiN2pLocal", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    kiiN2p0 * pow(T0ChemRef / T, 0.5)
);

volScalarField RiiN2p
(
    IOobject("RiiN2p", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    kiiN2pLocal * nN2p * nO2m
);

volScalarField chemReactionRateControl
(
    IOobject("chemReactionRateControl", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mag(RionN2p) + mag(RionO2p) + mag(Rdet) + mag(RattO2m) + mag(RrecO4p)
);

volScalarField chemReactionRateFull
(
    IOobject("chemReactionRateFull", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mag(RionN2p) + mag(RionO2p) +
    mag(RattO2m) + mag(Rdet) +
    mag(RrecO4p) + mag(RrecO2p) + mag(RrecN2p) + mag(RrecDissO) + mag(RrecDissN) +
    mag(RiiO4p) + mag(RiiO2p) + mag(RiiN2p)
);

scalar maxReactionRate = 0;
