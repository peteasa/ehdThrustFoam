
surfaceScalarField phi_e
(
    "phi_e",
    ( fvc::interpolate(-mu_e * E) & mesh.Sf() )
);

//Info << "phi_e dimensions: " << phi_e.dimensions() << endl;

surfaceScalarField phi_N2p
(
    "phi_N2p",
    ( fvc::interpolate(mu_N2p * E) & mesh.Sf() )
);

//Info << "De dimensions: " << De.dimensions() << endl;

fvScalarMatrix neEqn
(
      fvm::ddt(ne)
    + fvm::div(phi_e, ne)
    - fvm::laplacian(De, ne)
 ==
    Sion
);
// neEqn.solve();

fvScalarMatrix nN2pEqn
(
      fvm::ddt(nN2p)
    + fvm::div(phi_N2p, nN2p)
    - fvm::laplacian(DN2p, nN2p)
 ==
      Sion
);

// nN2pEqn.solve();

// Solve species
neEqn.relax();
SolverPerformance<scalar> nePerf = neEqn.solve();

// Print electron density residual
if (Pstream::master() && enableDetailedLogs) Info << "Residual ne = " << nePerf.finalResidual() << nl;

nN2pEqn.relax();
SolverPerformance<scalar> nN2pPerf = nN2pEqn.solve();

// Print ion density residual
if (Pstream::master() && enableDetailedLogs) Info << "Residual nN2p = " << nN2pPerf.finalResidual() << nl;

if (0)
{
    // Experimental
    // suppress un-physical negative densities
    // balance the minimum charge: nN2p - ne = 0
    minNe = gMin(ne);
    maxNe = gMax(ne);
    minN2 = gMin(nN2p);
    maxN2 = gMax(nN2p);

    dimensionedScalar minNeLimD("minNeLimD", ne.dimensions(), minNeLim);
    ne = 0.5 * (ne + sqrt(sqr(ne) + sqr(minNeLimD)));

    dimensionedScalar minN2LimD("minN2LimD", ne.dimensions(), minN2Lim);
    nN2p = 0.5 * (nN2p + sqrt(sqr(nN2p) + sqr(minN2LimD)));
}

// update the reference values
ne.correctBoundaryConditions();
nN2p.correctBoundaryConditions();

if (enableDetailedLogs)
{
    scalar minNe = gMin(ne);
    scalar maxNe = gMax(ne);
    if (0 < maxNe && Pstream::master()) Info << "min/max ne: " << minNe << " " << maxNe << " min/max: " << minNe / maxNe << nl;
    scalar minN2 = gMin(nN2p);
    scalar maxN2 = gMax(nN2p);
    if (0 < maxN2 && Pstream::master()) Info << "min/max nN2p: " << minN2 << " " << maxN2 << " min/max: " << minN2 / maxN2 << nl;
}

if (enableDetailedLogs)
{
    label nePatchID = mesh.boundaryMesh().findIndex("NELEMENT");
    if (nePatchID < 0)
    {
        FatalErrorInFunction
            << "Patch NELEMENT not found"
            << exit(FatalError);
    }

    const fvPatchScalarField& nepatch = ne.boundaryField()[nePatchID];
    scalar minNeP = gMin(nepatch);
    scalar maxNeP = gMax(nepatch);
    if (Pstream::master()) Info << "min/max ne at NELEMENT boundary: " << minNeP << " " << maxNeP << nl;

    label pePatchID = mesh.boundaryMesh().findIndex("PELEMENT");
    if (pePatchID < 0)
    {
        FatalErrorInFunction
            << "Patch PELEMENT not found"
            << exit(FatalError);
    }

    const fvPatchScalarField& pepatch = ne.boundaryField()[pePatchID];
    scalar minPeP = gMin(pepatch);
    scalar maxPeP = gMax(pepatch);
    if (Pstream::master()) Info << "min/max ne at PELEMENT boundary: " << minPeP << " " << maxPeP << nl;
}

// balancing of charge causes big error in ne
volScalarField chargeError = nN2p - ne;
if (enableDetailedLogs)
{
    scalar minChErr = gMin(chargeError);
    scalar maxChErr = gMax(chargeError);
    if (Pstream::master()) Info << " chargeError min/max: " << minChErr << " " << maxChErr << nl;
}

volScalarField rhoEOld
(
    "rhoEOld",
    rhoE
);

// rhoE = eCharge * (- ne);
rhoE = eCharge * (nN2p - ne);

if (enableDetailedLogs)
{
    scalar minRhoE = gMin(rhoE);
    scalar maxRhoE = gMax(rhoE);
    if (Pstream::master()) Info << "min/max rhoE: " << minRhoE << " " << maxRhoE << nl;

    volScalarField magChange
    (
        "magChange",
        mag(rhoE - rhoEOld)
    );
    scalar rhoChange = gMax(magChange);
    if (Pstream::master()) Info << "max change in rhoE: " << rhoChange << nl;
}

