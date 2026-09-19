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
    eCharge * ( nN2p - nO2m - ne )
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
    Je + JN2p + JO2m
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

volScalarField magE_NTdOld
(
    IOobject("magE_NTdOld", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    magE_NTd
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
volScalarField alpha
(
    IOobject("alpha", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    AN20 * Peff * exp(-BN20 * Peff / max(magE, smallE))
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

// detachment O2- -> O2 + e-
dimensionedScalar AdetO2O2m
(
    physicalProperties.lookup("AdetO2O2m")
);

dimensionedScalar Eact_eV
(
    physicalProperties.lookup("Eact_eV")
);

// recombination e- + N2+ -> N2
dimensionedScalar krecN2p0
(
    physicalProperties.lookup("krecN2p0")
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
    ionScale * alpha * mu_e * magE * ne
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

volScalarField nuDetO2
(
    IOobject("nuDetO2", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    AdetO2O2m * Peff * exp(-Eact_eV * eCharge.value() / (kB * T))
);

volScalarField RdetO2e
(
    IOobject("RdetO2e", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    nuDetO2 * nO2m
);

dimensionedScalar krecN2pLocal = krecN2p0 * pow(T0ChemRef / TeK, 0.39);

volScalarField RrecN2p
(
    IOobject("RrecN2p", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    krecN2pLocal * ne * nN2p
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

volScalarField chemReactionRate
(
    IOobject("chemReactionRate", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mag(RionN2p) + mag(RattO2m) + mag(RdetO2e) + mag(RrecN2p) + mag(RiiN2p)
);

scalar maxReactionRate = 0;
