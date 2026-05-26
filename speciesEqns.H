volVectorField Ue = -mu_e * E;

//surfaceScalarField phiUe = fvc::interpolate(Ue) & mesh.Sf();
//surfaceScalarField phiEFlux_e = fvc::flux(-mu_e * E);
//surfaceScalarField phiEFlux_e = fvc::flux(Ue);
//surfaceScalarField phiEFlux_N2p = fvc::flux(mu_N2p * E);

surfaceScalarField phid_e
(
    "phid_e",
    ( fvc::interpolate(-mu_e * Eeff) & mesh.Sf() )
);

//Info << "phid_e dimensions: " << phi_e.dimensions() << endl;

surfaceScalarField phid_N2p
(
    "phid_N2p",
    ( fvc::interpolate(mu_N2p * Eeff) & mesh.Sf() )
);

//Info << "De dimensions: " << De.dimensions() << endl;

fvScalarMatrix neEqn
(
      fvm::ddt(ne)
    + fvm::div(phid_e, ne)
    - fvm::laplacian(De, ne)
 ==
    Sion
);
// neEqn.solve();

fvScalarMatrix nN2pEqn
(
      fvm::ddt(nN2p)
    + fvm::div(phid_N2p, nN2p)
    - fvm::laplacian(DN2p, nN2p)
 ==
      Sion
);

// nN2pEqn.solve();

// Solve species
neEqn.relax();
SolverPerformance<scalar> nePerf = neEqn.solve();

// update the reference values
ne.correctBoundaryConditions();

// Print electron density residual
Info << "Residual ne = " << nePerf.finalResidual() << endl;

nN2pEqn.relax();
SolverPerformance<scalar> nN2pPerf = nN2pEqn.solve();

// update the reference values
nN2p.correctBoundaryConditions();

// Print ion density residual
Info << "Residual nN2p = " << nN2pPerf.finalResidual() << endl;

if (1)
{
    Info << "internal min/max ne: " << min(ne.internalField()).value() << " " << max(ne.internalField()).value() << endl;
    Info << "internal min/max nN2p: " << min(nN2p.internalField()).value() << " " << max(nN2p.internalField()).value() << endl;
}

if (1)
{
    label nePatchID = mesh.boundaryMesh().findIndex("NELEMENT");
    if (nePatchID < 0)
    {
        FatalErrorInFunction
            << "Patch NELEMENT not found"
            << exit(FatalError);
    }

    const fvPatchScalarField& nepatch = ne.boundaryField()[nePatchID];
    Info << "min/max ne at NELEMENT boundary: " << gMin(nepatch) << " " << gMax(nepatch) << endl;

    label pePatchID = mesh.boundaryMesh().findIndex("PELEMENT");
    if (pePatchID < 0)
    {
        FatalErrorInFunction
            << "Patch PELEMENT not found"
            << exit(FatalError);
    }

    const fvPatchScalarField& pepatch = ne.boundaryField()[pePatchID];
    Info << "min/max ne at PELEMENT boundary: " << gMin(pepatch) << " " << gMax(pepatch) << endl;
}

// balancing of charge causes big error in ne
volScalarField chargeError = nN2p - ne;
Info << " chargeError min/max: " << min(chargeError.internalField()).value() << " " << max(chargeError.internalField()).value() << endl;

if (0)
{
    ne += 0.1 * chargeError;
    Info << "balanced charge ne min/max: " << min(ne.internalField()).value() << " " << max(ne.internalField()).value() << endl;
}

// rhoE = eCharge * (- ne);
rhoE = eCharge * (nN2p - ne);

// TODO Find a way to remove this clamp
if (0) {
    Info << "DEBUG: clamp rhoE" << endl;
    rhoE = min(max(rhoE,
                   dimensionedScalar("rhoMin", rhoE.dimensions(), -1e-6)),
               dimensionedScalar("rhoMax", rhoE.dimensions(), 1e-6));
}

Info << "min/max rhoE: " << min(rhoE).value() << " " << max(rhoE).value() << endl;

