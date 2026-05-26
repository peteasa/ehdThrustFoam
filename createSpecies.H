Info<< "Creating plasma species\n" << endl;
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

dimensionSet dimVoltage(1, 2, -3, 0, 0, -1, 0);

// dimArea / dimVoltage / dimTime
dimensionedScalar mu_e
(
    physicalProperties.lookup("mu_e")
);

dimensionedScalar mu_N2p
(
    physicalProperties.lookup("mu_N2p")
);

// 2 eV dimTemperature
dimensionedScalar Te
(
    physicalProperties.lookup("Te")
);

dimensionedScalar De = mu_e * Te;     //dimensionedScalar De("De", dimArea/dimTime, 0.01);
dimensionedScalar DN2p = mu_N2p * Te; //dimensionedScalar DN2p("DN2p", dimArea/dimTime, 0.0001);

Info << "De = " << De << endl
     << "DN2p = " << DN2p << endl;

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

// recombination (kei * ne * nN2p) -> dimless / dimVolume / dimTime for recombination e- + N2+ -> N2
dimensionedScalar kei
(
    physicalProperties.lookup("kei")
);

