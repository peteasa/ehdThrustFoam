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
    mu_e0 * (1.0 - beta * xH20) / Peff
);

dimensionedScalar mu_N2p0
(
    physicalProperties.lookup("mu_N2p0")
);

volScalarField mu_N2p
(
    IOobject ( "mu_N2p", runTime.name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
    mu_N2p0 * (1.0 - beta * xH20) / Peff
);

dimensionedScalar mu_O2m0
(
    physicalProperties.lookup("mu_O2m0")
);

volScalarField mu_O2m
(
    IOobject ( "mu_O2m", runTime.name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
    mu_O2m0 * (1.0 - beta * xH20) / Peff
);

// 2 eV dimTemperature
dimensionedScalar Te
(
    physicalProperties.lookup("Te")
);

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

dimensionedScalar kB
(
    "kB",
    dimEnergy/dimTemperature,
    1.380649e-23
);

// reference temperature for chemical reactions
dimensionedScalar T0ChemRef
(
    physicalProperties.lookup("T0ChemRef")
);

// ionisation - for Townsend ionisation N2 -> N2+ + e-
dimensionedScalar A0
(
    physicalProperties.lookup("A0")
);

dimensionedScalar B0
(
    physicalProperties.lookup("B0")
);

dimensionedScalar smallE("smallE", magE.dimensions(), 1e3);
volScalarField alpha
(
    IOobject("alpha", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    A0 * Peff * exp(-B0 * Peff / max(magE, smallE))
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

dimensionedScalar Ed
(
    EdTd * 1e-21
);

dimensionedScalar DeltaETd
(
    physicalProperties.lookup("DeltaETd")
);

dimensionedScalar DeltaE
(
    DeltaETd * 1e-21
);

// detachment O2- -> O2 + e-
dimensionedScalar Adet
(
    physicalProperties.lookup("Adet")
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

dimensionedScalar humidityCorrection
(
    physicalProperties.lookup("humidityCorrection")
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
    (C1 + C2 * xH20) * sqr(Peff) * exp( -D * Peff / max(magE, smallE))
);

volScalarField eta2b
(
    IOobject("eta2b", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    Ad * exp( - pow( ( magE * kB * T0TransRef / Peff - Ed ) / DeltaE , 2) )
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

volScalarField nuDet
(
    IOobject("nuDet", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    Adet * Peff * exp(-Eact_eV * eCharge.value() / (kB * T))
);

volScalarField RdetO2e
(
    IOobject("RdetO2e", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    nuDet * nO2m
);

volScalarField krecN2pLocal
(
    IOobject("krecN2pLocal", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    krecN2p0 * pow(T0ChemRef / T, 0.5) * (Pphys / P0TransRef)
);

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
    kiiN2p0 * ( T0ChemRef / T ) * ( 1.0 + humidityCorrection * xH20 )
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
