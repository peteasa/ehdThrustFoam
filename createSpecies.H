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
    eCharge * ( nN2p - ne )
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

// 2 eV dimTemperature
dimensionedScalar Te
(
    physicalProperties.lookup("Te")
);

// dimensionedScalar De = mu_e * Te;  dimensionedScalar De("De", dimArea/dimTime, 0.01);
volScalarField De
(
    IOobject ( "De", runTime.name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
    mu_e * Te
);

// dimensionedScalar DN2p = mu_N2p * Te; //dimensionedScalar DN2p("DN2p", dimArea/dimTime, 0.0001);
volScalarField DN2p
(
    IOobject ( "DN2p", runTime.name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
    mu_N2p * Te
);

// calculate particle fluxes
volVectorField GammaE
(
    IOobject ( "GammaE", runTime.name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
    U * ne
  - mu_e * E * ne
  - De * fvc::grad(ne)
);
//if (Pstream::master()) Info << "GammaE dimensions: " << GammaE.dimensions() << nl;

volVectorField GammaN2p
(
    IOobject ( "GammaN2p", runTime.name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
    U * nN2p
  + mu_N2p * E * nN2p
  - DN2p * fvc::grad(nN2p)
);
//if (Pstream::master()) Info << "GammaN2p dimensions: " << GammaN2p.dimensions() << nl;

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
    Je + JN2p
);

dimensionedScalar E0("E0", E.dimensions(), 2.5e6);
volScalarField magE
(
    IOobject("magE", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mag(E)
);

// ionisation - for Townsend ionisation N2 -> N2+ + e-
// dimless / dimLength
dimensionedScalar A
(
    physicalProperties.lookup("A")
);
// dimVoltage / dimLength
dimensionedScalar B
(
    physicalProperties.lookup("B")
);

dimensionedScalar smallE("smallE", magE.dimensions(), 1e3);
volScalarField alpha
(
    IOobject("alpha", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    A * exp(-B / max(magE, smallE))
);

// recombination (kei * ne * nN2p) -> dimless / dimVolume / dimTime for recombination e- + N2+ -> N2
dimensionedScalar kei
(
    physicalProperties.lookup("kei")
);

dimensionedScalar ionScale("ionScale", dimless, 1e-1);        // 1e-4
volScalarField Sion
(
    IOobject("Sion", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    ionScale * alpha * mu_e * magE * ne
);

volScalarField Rei
(
    IOobject("Rei", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    kei * ne * nN2p
);

volScalarField reactionRate
(
    IOobject("reactionRate", runTime.name(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    Sion + Rei
);

scalar maxReactionRate = 0;
