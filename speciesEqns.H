
De = mu_e * Te;
DN2p = mu_N2p * Te;
DO4p = mu_O4p * Te;
DO2p = mu_O2p * Te;
DO2m = mu_O2m * Te;

if (enableDetailedLogs)
{
    scalar minDe = gMin(De.internalField());
    scalar maxDe = gMax(De.internalField());

    scalar minDN2p = gMin(DN2p.internalField());
    scalar maxDN2p = gMax(DN2p.internalField());

    scalar minDO4p = gMin(DO4p.internalField());
    scalar maxDO4p = gMax(DO4p.internalField());

    scalar minDO2p = gMin(DO2p.internalField());
    scalar maxDO2p = gMax(DO2p.internalField());

    scalar minDO2m = gMin(DO2m.internalField());
    scalar maxDO2m = gMax(DO2m.internalField());

    if (Pstream::master()) Info << "De min/max: " << minDe << " " << maxDe << nl
                            << "DN2p min/max: " << minDN2p << " " << maxDN2p << nl
                            << "DO4p min/max: " << minDO4p << " " << maxDO4p << nl
                            << "DO2p min/max: " << minDO2p << " " << maxDO2p << nl
                            << "DO2m min/max: " << minDO2m << " " << maxDO2m << nl;
}

// convection term generated from the bulk neutral gas velocity
surfaceScalarField phi_gas
(
    "phi_gas",
    fvc::flux(U)
);

if (startupT)
{
    if (1 && (runTime.value() < startupT))
    {
        if (enableDetailedLogs && Pstream::master()) Info << "Slow Start: Separate convection interaction until time: " << startupT << nl;
        phi_gas = dimensionedScalar("zv", phi_gas.dimensions(), 0);
        interactionStart = runTime.timeIndex();
    }
    else
    {
        int count = runTime.timeIndex() - interactionStart;
        if (count < startupMaxCount)
        {
            // gradually increase impact of movement
            dimensionedScalar gScale("gScale", dimless, startupIncrement * count);
            if (enableDetailedLogs && Pstream::master()) Info << "gScale: " << startupIncrement * count << nl;
            phi_gas *= gScale;
        }
    }
}

// drift terms
surfaceScalarField phi_de
(
    "phi_de",
    fvc::flux(-mu_e * E)
);

surfaceScalarField phi_dO4p
(
    "phi_dO4p",
    fvc::flux(mu_O4p * E)
);

surfaceScalarField phi_dO2p
(
    "phi_dO2p",
    fvc::flux(mu_O2p * E)
);

surfaceScalarField phi_dN2p
(
    "phi_dN2p",
    fvc::flux(mu_N2p * E)
);

surfaceScalarField phi_dO2m
(
    "phi_dO2m",
    fvc::flux(-mu_O2m * E)
);

// combine with the diffusion terms and the source terms
fvScalarMatrix neEqn
(
      fvm::ddt(ne)
    + fvm::div(phi_gas, ne)
    + fvm::div(phi_de, ne)
    - fvm::laplacian(De, ne)
 ==
      RionN2p + RionO2p + Rdet - RrecO4p - RrecN2p - RrecO2p - RattO2m
);

fvScalarMatrix nN2pEqn
(
      fvm::ddt(nN2p)
    + fvm::div(phi_gas, nN2p)
    + fvm::div(phi_dN2p, nN2p)
    - fvm::laplacian(DN2p, nN2p)
 ==
      RionN2p - RrecN2p - RiiN2p - RCT
);

fvScalarMatrix nO4pEqn
(
      fvm::ddt(nO4p)
    + fvm::div(phi_gas, nO4p)
    + fvm::div(phi_dO4p, nO4p)
    - fvm::laplacian(DO4p, nO4p)
 ==
      RclustO4p - RrecO4p - RiiO4p
);

fvScalarMatrix nO2pEqn
(
      fvm::ddt(nO2p)
    + fvm::div(phi_gas, nO2p)
    + fvm::div(phi_dO2p, nO2p)
    - fvm::laplacian(DO2p, nO2p)
 ==
      RionO2p + RCT - RrecO2p - RiiO2p - RclustO4p
);

fvScalarMatrix nO2mEqn
(
      fvm::ddt(nO2m)
    + fvm::div(phi_gas, nO2m)
    + fvm::div(phi_dO2m, nO2m)
    - fvm::laplacian(DO2m, nO2m)
 ==
      RattO2m - RiiN2p - RiiO2p - RiiO4p - Rdet
);

// Solve species
neEqn.relax();
SolverPerformance<scalar> nePerf = neEqn.solve();

// Print electron density residual
// if (enableDetailedLogs && Pstream::master()) Info << "Residual ne = " << nePerf.finalResidual() << nl;

nN2pEqn.relax();
SolverPerformance<scalar> nN2pPerf = nN2pEqn.solve();

nO4pEqn.relax();
SolverPerformance<scalar> nO4pPerf = nO4pEqn.solve();

nO2pEqn.relax();
SolverPerformance<scalar> nO2pPerf = nO2pEqn.solve();

// Print ion density residual
// if (enableDetailedLogs && Pstream::master()) Info << "Residual nN2p = " << nN2pPerf.finalResidual() << nl;

nO2mEqn.relax();
SolverPerformance<scalar> nO2mPerf = nO2mEqn.solve();

// Print ion density residual
// if (enableDetailedLogs && Pstream::master()) Info << "Residual nO2m = " << nO2mPerf.finalResidual() << nl;

if (0)
{
    // Experimental
    // suppress un-physical negative densities
    minNe = gMin(ne);
    maxNe = gMax(ne);
    minN2p = gMin(nN2p);
    maxN2p = gMax(nN2p);
    minO4p = gMin(nO4p);
    maxO4p = gMax(nO4p);
    minO2p = gMin(nO2p);
    maxO2p = gMax(nO2p);
    minO2m = gMin(nO2m);
    maxO2m = gMax(nO2m);

    if (1)
    {
        if (0 && enableDetailedLogs && 0 < maxNe && Pstream::master()) Info << "before clamp min/max ne: " << minNe << " " << maxNe
                                                                            << " min/max: " << minNe / maxNe << nl;
        dimensionedScalar minNeLimD("minNeLimD", ne.dimensions(),
                                    (minNe < 0? max(1e-10, min(maxNe * 0.1, minNeLim)) : 0));
        ne = 0.5 * (ne + sqrt(sqr(ne) + sqr(minNeLimD)));

        if (0 && enableDetailedLogs && 0 < maxN2p && Pstream::master()) Info << "before clamp min/max nN2p: " << minN2p << " " << maxN2p
                                                                             << " min/max: " << minN2p / maxN2p << nl;
        dimensionedScalar minN2pLimD("minN2pLimD", ne.dimensions(),
                                     (minN2p < 0? max(1e-10, min(maxN2p * 0.1, minN2pLim)) : 0));
        nN2p = 0.5 * (nN2p + sqrt(sqr(nN2p) + sqr(minN2pLimD)));

        if (0 && enableDetailedLogs && 0 < maxO4p && Pstream::master()) Info << "before clamp min/max nO4p: " << minO4p << " " << maxO4p
                                                                             << " min/max: " << minN2p / maxN2p << nl;
        dimensionedScalar minO4pLimD("minO4pLimD", ne.dimensions(),
                                     (minO4p < 0? max(1e-10, min(maxO4p * 0.1, minO4pLim)) : 0));
        nO4p = 0.5 * (nO4p + sqrt(sqr(nO4p) + sqr(minO4pLimD)));

        if (0 && enableDetailedLogs && 0 < maxO2p && Pstream::master()) Info << "before clamp min/max nO2p: " << minO2p << " " << maxO2p
                                                                             << " min/max: " << minN2p / maxN2p << nl;
        dimensionedScalar minO2pLimD("minO2pLimD", ne.dimensions(),
                                     (minO2p < 0? max(1e-10, min(maxO2p * 0.1, minO2pLim)) : 0));
        nO2p = 0.5 * (nO2p + sqrt(sqr(nO2p) + sqr(minO2pLimD)));

        if (0 && enableDetailedLogs && 0 < maxO2m && Pstream::master()) Info << "before clamp min/max nO2m: " << minO2m << " " << maxO2m
                                                                             << " min/max: " << minO2m / maxO2m << nl;
        dimensionedScalar minO2mLimD("minO2mLimD", ne.dimensions(),
                                     (minO2m < 0? max(1e-10, min(maxO2m * 0.1, minO2mLim)) : 0));
        nO2m = 0.5 * (nO2m + sqrt(sqr(nO2m) + sqr(minO2mLimD)));
    }

    if (0)
    {
        dimensionedScalar minD("minD", ne.dimensions(), 0.0);
        ne = max(minD, ne);
        nN2p = max(minD, nN2p);
        nO4p = max(minD, nO4p);
        nO2p = max(minD, nO2p);
        nO2m = max(minD, nO2m);
    }
}

// update the reference values
ne.correctBoundaryConditions();
nN2p.correctBoundaryConditions();
nO4p.correctBoundaryConditions();
nO2p.correctBoundaryConditions();
nO2m.correctBoundaryConditions();

if (enableDetailedLogs)
{
    scalar minNe = gMin(ne);
    scalar maxNe = gMax(ne);
    if (0 < maxNe && Pstream::master()) Info << "min/max ne: " << minNe << " " << maxNe
                                             << " min/max: " << minNe / maxNe << nl;
    scalar minN2p = gMin(nN2p);
    scalar maxN2p = gMax(nN2p);
    if (0 < maxN2p && Pstream::master()) Info << "min/max nN2p: " << minN2p << " " << maxN2p
                                             << " min/max: " << minN2p / maxN2p << nl;
    scalar minO4p = gMin(nO4p);
    scalar maxO4p = gMax(nO4p);
    if (0 < maxO4p && Pstream::master()) Info << "min/max nO4p: " << minO4p << " " << maxO4p
                                             << " min/max: " << minO4p / maxO4p << nl;
    scalar minO2p = gMin(nO2p);
    scalar maxO2p = gMax(nO2p);
    if (0 < maxO2p && Pstream::master()) Info << "min/max nO2p: " << minO2p << " " << maxO2p
                                             << " min/max: " << minO2p / maxO2p << nl;
    scalar minO2m = gMin(nO2m);
    scalar maxO2m = gMax(nO2m);
    if (0 < maxO2m && Pstream::master()) Info << "min/max nO2m: " << minO2m << " " << maxO2m
                                             << " min/max: " << minO2m / maxO2m << nl;
}

if (enableDetailedLogs)
{
    const fvPatchScalarField& nepatch = ne.boundaryField()[nPatchID];
    scalar minNeP = gMin(nepatch);
    scalar maxNeP = gMax(nepatch);
    if (Pstream::master()) Info << "min/max ne at NELEMENT boundary: " << minNeP << " " << maxNeP << nl;

    const fvPatchScalarField& pepatch = ne.boundaryField()[pPatchID];
    scalar minPeP = gMin(pepatch);
    scalar maxPeP = gMax(pepatch);
    if (Pstream::master()) Info << "min/max ne at PELEMENT boundary: " << minPeP << " " << maxPeP << nl;
}

// balancing of charge causes big error in ne
volScalarField chargeError = nN2p + nO4p + nO2p - nO2m - ne;
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
rhoE = eCharge * (nN2p + nO4p + nO2p - nO2m - ne);

// calculate particle fluxes
GammaE = U * ne - mu_e * E * ne - De * fvc::grad(ne);
GammaN2p = U * nN2p + mu_N2p * E * nN2p - DN2p * fvc::grad(nN2p);
GammaO4p = U * nO4p + mu_O4p * E * nO4p - DO4p * fvc::grad(nO4p);
GammaO2p = U * nO2p + mu_O2p * E * nO2p - DO2p * fvc::grad(nO2p);
GammaO2m = U * nO2m - mu_O2m * E * nO2m - DO2m * fvc::grad(nO2m);
Je = - eCharge * GammaE;
JN2p = eCharge * GammaN2p;
JO4p = eCharge * GammaO4p;
JO2p = eCharge * GammaO2p;
JO2m = - eCharge * GammaO2m;
J = Je + JN2p + JO4p + JO2p + JO2m;

if (enableDetailedLogs)
{
    auto currentOnPatch = []
    (
        const Foam::fvMesh& mesh,
        const volVectorField& Jfield,
        const label patchID
    )
    {
        return gSum
        (
            Jfield.boundaryField()[patchID]
            &
            mesh.Sf().boundaryField()[patchID]
        );
    };

    scalar eAnodeI = currentOnPatch(mesh, Je, pPatchID);
    scalar eCathodeI = currentOnPatch(mesh, Je, nPatchID);

    scalar N2pAnodeI = currentOnPatch(mesh, JN2p, pPatchID);
    scalar N2pCathodeI = currentOnPatch(mesh, JN2p, nPatchID);

    scalar O4pAnodeI = currentOnPatch(mesh, JO4p, pPatchID);
    scalar O4pCathodeI = currentOnPatch(mesh, JO4p, nPatchID);

    scalar O2pAnodeI = currentOnPatch(mesh, JO2p, pPatchID);
    scalar O2pCathodeI = currentOnPatch(mesh, JO2p, nPatchID);

    scalar O2mAnodeI = currentOnPatch(mesh, JO2m, pPatchID);
    scalar O2mCathodeI = currentOnPatch(mesh, JO2m, nPatchID);

    scalar anodeI = eAnodeI + N2pAnodeI + O4pAnodeI + O2pAnodeI + O2mAnodeI;
    scalar cathodeI = eCathodeI + N2pCathodeI + O4pCathodeI + O2pCathodeI + O2mCathodeI;

    if (Pstream::master()) Info << "anode currents (e/N2p/O4p/O2p/O2m): " << eAnodeI << " "
                                << N2pAnodeI << " " << O4pAnodeI << " " << O2pAnodeI << " " << O2mAnodeI << nl;
    if (Pstream::master()) Info << "cathode currents (e/N2p/O4p/O2p/O2m): " << eCathodeI << " "
                                << N2pCathodeI << " " << O4pCathodeI << " " << O2pCathodeI << " " << O2mCathodeI << nl;
    if (Pstream::master()) Info << "anode/cathode currents: " << anodeI << " " << cathodeI << nl;
}

if (1)
{
    scalar minRhoE = gMin(rhoE);
    scalar maxRhoE = gMax(rhoE);
    scalar maxMagRhoE = max(mag(minRhoE), mag(maxRhoE));

    volScalarField magDRhoEDt
    (
        "magDRhoEDt",
        mag(rhoE - rhoEOld)/runTime.deltaTValue()
    );

    scalar maxDRhoEDt = gMax(magDRhoEDt);

    if (0 < maxMagRhoE && maxDRhoEDtRateDec < 1000)
    {
        maxDRhoEDtRate = maxDRhoEDt / maxMagRhoE;
    }
    else
    {
        maxDRhoEDtRate = 0;
    }

    if (enableDetailedLogs && Pstream::master()) Info << "min/max rhoE: " << minRhoE << " " << maxRhoE
                                                << " maxDRhoEDt: " << maxDRhoEDt
                                                << " maxDRhoEDtRate: " << maxDRhoEDtRate
                                                << " maxMagRhoE: " << maxMagRhoE
                                                << nl;
}

