
if (Pstream::master()) Info << "Reading physicalProperties\n" << nl;
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

if (Pstream::master()) Info << "Reading fvSolutionDict\n" << nl;
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

if (Pstream::master()) Info << "Reading field kinematic p\n" << nl;
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

if (Pstream::master()) Info << "p settings:  pRefCell = " << pRefCell << "  pRefValue = " << pRefValue << nl;

if (Pstream::master()) Info << "Reading field U\n" << nl;
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

if (Pstream::master()) Info << "Reading field phi\n" << nl;
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

// Manage a gradual change in the electrostatic and convection interactions
int interactionStart = 0;

// Stabilise phiE electrostatic potential
const dictionary& potentialDict = fvSolutionDict.subDict("potentialFlow");
int nNonOrthogonalPotCorrectors = potentialDict.lookupOrDefault<int>("nNonOrthogonalCorrectors", 50);

if (Pstream::master()) Info << "potential init settings:" << endl
     << "  nNonOrthogonalPotCorrectors      = " << nNonOrthogonalPotCorrectors << nl;

const dictionary& solversDict = fvSolutionDict.subDict("solvers");
const dictionary& solversPDict = solversDict.subDict("p");
int nCellsInCoarsestLevel = solversPDict.lookupOrDefault<int>("nCellsInCoarsestLevel", 10);
label total_cells = returnReduce(mesh.nCells(), sumOp<label>());
label maxRank = Pstream::nProcs();
{
    // nCellsInCoarsestLevel tells GAMG when to stop agglomerating
    // best to choose a number where GAMGSolver uses say 5 levels (level: 0 -> 4)
    // with between 10 -> 1000 cells in the last level depending on total_cells
    // GAMG agglomeration: grouping of fine cells into coarse "super-cells"
    //     a good agglomeration produces coarse cells that still represent the physics
    //     and the connectivity of the original mesh
    //     too many levels and the calculation time increases
    // Look at the cell reduction per level:
    //     good: factors of ~5 –> 10
    //      bad: tiny reductions followed by a huge collapse
    // Look at the GAMG iteration count:
    //     fewer iterations per time step/nonlinear iteration is better
    //     if increasing levels makes GAMG slower, it is not worth it
    // Check the coarse level:
    //     avoid nCells:1 or a handful of cells unless the case is small
    //     something in the range 50–500 cells is usually reasonable
    // agglomeration quality matters more than the reduction factor:
    //     a forced x5 reduction between each level creates that poor coarse meshes
    //     can be worse than a x2.2 reduction with excellent agglomeration

    if (Pstream::master()) Info << (((int(total_cells / (2.2 * 2.2 * 2.2 * 2.2 * maxRank)) - 1 < nCellsInCoarsestLevel)
                                     || nCellsInCoarsestLevel < 10 ) ? "WARNING " : "")
                                << "nCellsInCoarsestLevel: " << nCellsInCoarsestLevel
                                << " total_cells / maxRank: " << total_cells / maxRank
                                << " total_cells / (23.5 * maxRank) - 1: " << total_cells / (2.2 * 2.2 * 2.2 * 2.2 * maxRank) - 1
                                << nl;
    if (((int(total_cells / (2.2 * 2.2 * 2.2 * 2.2 * maxRank)) - 1 < nCellsInCoarsestLevel)
          || nCellsInCoarsestLevel < 10 ) && Pstream::master())
        Info << "        to tune nCellsInCoarsestLevel enable DebugSwitches/GAMG" << nl;
    if (Pstream::master()) Info << "total_cells: " << total_cells << " nProcs (maxRank): " << maxRank << nl;
}

const dictionary& pimpleDict = fvSolutionDict.subDict("PIMPLE");
// New to ehdMultiSpeciesFoam
int nPhiECorrectors      = pimpleDict.lookupOrDefault<int>("nPhiECorrectors", 2);

// Existing PIMPLE configuration See $FOAM_SRC/finiteVolume/cfdTools/general/solutionControl/pimpleControl
int nOuterCorrectors = pimpleDict.lookupOrDefault<int>("nOuterCorrectors", 1);
int nCorrectors      = pimpleDict.lookupOrDefault<int>("nCorrectors", 2);
int nNonOrthogonalCorrectors = pimpleDict.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);
bool momentumPredictor = pimpleDict.lookupOrDefault("momentumPredictor", true);

if (Pstream::master()) Info << "PIMPLE settings:" << nl
     << "  nPhiECorrectors      = " << nPhiECorrectors << nl
     << "  nOuterCorrectors      = " << nOuterCorrectors << nl
     << "  nCorrectors      = " << nCorrectors << nl
     << "  nNonOrthogonalCorrectors = " << nNonOrthogonalCorrectors << nl
     << "  momentumPredictor = " << (momentumPredictor ? "true" : "false") << nl;

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

if (Pstream::master()) Info << "Reading field phiE\n" << nl;
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

label pPatchID = mesh.boundaryMesh().findIndex("PELEMENT");
if (pPatchID < 0)
{
    FatalErrorInFunction
        << "Patch PELEMENT not found"
        << exit(FatalError);
}

label nPatchID = mesh.boundaryMesh().findIndex("NELEMENT");
if (nPatchID < 0)
{
    FatalErrorInFunction
        << "Patch NELEMENT not found"
        << exit(FatalError);
}
