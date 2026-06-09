
Info<< "Reading physicalProperties\n" << endl;
IOdictionary physicalProperties
(
    IOobject
    (
        "physicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading fvSolutionDict\n" << endl;
IOdictionary fvSolutionDict
(
    IOobject
    (
        "fvSolution",
        runTime.system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    )
);

// Transport properties

Info<< "Reading field kinematic p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

const dictionary& pDict = mesh.solution().subDict("solvers").subDict("p");
label pRefCell = pDict.lookupOrDefault<label>("referenceCell", 0);
scalar pRefValue = pDict.lookupOrDefault<scalar>("referenceValue", 0.0);

Info << "p settings:  pRefCell = " << pRefCell << "  pRefValue = " << pRefValue << endl;

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

dimensionedScalar rho
(
    physicalProperties.lookup("rho")
);

Info<< "Reading field phi\n" << endl;
//#include "createPhi.H"
surfaceScalarField phi
(
    IOobject(
        "phi",
        runTime.name(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(U)     // fvc::interpolate(U) & mesh.Sf()
);

const dictionary& potentialDict = fvSolutionDict.subDict("potentialFlow");
int nNonOrthogonalPotCorrectors = potentialDict.lookupOrDefault<int>("nNonOrthogonalCorrectors", 50);

Info << "potential init settings:" << endl
     << "  nNonOrthogonalPotCorrectors      = " << nNonOrthogonalPotCorrectors << endl;

const dictionary& pimpleDict = fvSolutionDict.subDict("PIMPLE");
// New to ehdMultiSpeciesFoam
int nPhiECorrectors      = pimpleDict.lookupOrDefault<int>("nPhiECorrectors", 2);

// Existing PIMPLE configuration See $FOAM_SRC/finiteVolume/cfdTools/general/solutionControl/pimpleControl
int nOuterCorrectors = pimpleDict.lookupOrDefault<int>("nOuterCorrectors", 1);
int nCorrectors      = pimpleDict.lookupOrDefault<int>("nCorrectors", 2);
int nNonOrthogonalCorrectors = pimpleDict.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);
bool momentumPredictor = pimpleDict.lookupOrDefault("momentumPredictor", true);

Info << "PIMPLE settings:" << endl
     << "  nPhiECorrectors      = " << nPhiECorrectors << endl
     << "  nOuterCorrectors      = " << nOuterCorrectors << endl
     << "  nCorrectors      = " << nCorrectors << endl
     << "  nNonOrthogonalCorrectors = " << nNonOrthogonalCorrectors << endl
     << "  momentumPredictor = " << (momentumPredictor ? "true" : "false") << endl;

autoPtr<viscosityModel> viscosity(viscosityModel::New(mesh));
autoPtr<incompressible::momentumTransportModel> momentumTransport
(
    incompressible::momentumTransportModel::New(U, phi, viscosity)
);

// Electrostatic properties

dimensionedScalar epsilon0
(
    physicalProperties.lookup("epsilon0")
);

dimensionedScalar eCharge
(
    physicalProperties.lookup("eCharge")
);

Info<< "Reading field phiE\n" << endl;
volScalarField phiE
(
    IOobject
    (
        "phiE",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField E
(
    IOobject
    (
        "E",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    -fvc::grad(phiE)
);

