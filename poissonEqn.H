
volScalarField phiEOld
(
    "phiEOld",
    phiE
);

fvScalarMatrix phiEqn
(
    fvm::laplacian(phiE)
    ==
  - rhoE / epsilon0
);

// Solve Poisson
phiEqn.relax();
SolverPerformance<scalar> phiEPerf = phiEqn.solve();

// Print Poisson residual
if (Pstream::master() && enableDetailedLogs) Info << "Residual phiE = " << phiEPerf.finalResidual() << nl;

if (enableDetailedLogs)
{
    scalar minPhiE = gMin(phiE);
    scalar maxPhiE = gMax(phiE);
    scalar avPhiE = gAverage(phiE);
    if (Pstream::master()) Info << "min/max phiE: " << minPhiE << " " << maxPhiE << " average phiE: " << avPhiE << " phiE range = " << mag(maxPhiE - minPhiE) << nl;

    volScalarField magChange
    (
        "magChange",
        mag(phiE - phiEOld)
    );
    scalar phiEChange = gMax(magChange);
    if (Pstream::master()) Info << "max change in phiE: " << phiEChange << nl;
}

if (enableDetailedLogs)
{
    // display the the reference value
    label ppPatchID = mesh.boundaryMesh().findIndex("PELEMENT");
    if (ppPatchID < 0)
    {
        FatalErrorInFunction
            << "Patch PELEMENT not found"
            << exit(FatalError);
    }

    const fvPatchScalarField& pphiEpatch = phiE.boundaryField()[ppPatchID];
    scalar minPPhiE = gMin(pphiEpatch);
    scalar maxPPhiE = gMax(pphiEpatch);
    if (Pstream::master()) Info << "precorrection min/max phiE at PELEMENT boundary: " << minPPhiE << " " << maxPPhiE << nl;
}
