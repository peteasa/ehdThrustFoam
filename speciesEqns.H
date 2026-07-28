
De = mu_e * Te;
DN2p = mu_N2p * Te;
DO2m = mu_O2m * Te;

if (enableDetailedLogs)
{
    scalar minDe = gMin(De.internalField());
    scalar maxDe = gMax(De.internalField());

    scalar minDN2p = gMin(DN2p.internalField());
    scalar maxDN2p = gMax(DN2p.internalField());

    scalar minDO2m = gMin(DO2m.internalField());
    scalar maxDO2m = gMax(DO2m.internalField());

    if (Pstream::master()) Info << "De min/max: " << minDe << " " << maxDe << nl
                            << "DN2p min/max: " << minDN2p << " " << maxDN2p << nl
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
        if (count < 1/startupIncrement)
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
      Sion - Rei - SO2ion + SO2e
);
// neEqn.solve();

fvScalarMatrix nN2pEqn
(
      fvm::ddt(nN2p)
    + fvm::div(phi_gas, nN2p)
    + fvm::div(phi_dN2p, nN2p)
    - fvm::laplacian(DN2p, nN2p)
 ==
      Sion - Rei - Rii
);
// nN2pEqn.solve();

fvScalarMatrix nO2mEqn
(
      fvm::ddt(nO2m)
    + fvm::div(phi_gas, nO2m)
    + fvm::div(phi_dO2m, nO2m)
    - fvm::laplacian(DO2m, nO2m)
 ==
      SO2ion - Rii - SO2e
);
// nO2mEqn.solve();

// Solve species
neEqn.relax();
SolverPerformance<scalar> nePerf = neEqn.solve();

// Print electron density residual
if (enableDetailedLogs && Pstream::master()) Info << "Residual ne = " << nePerf.finalResidual() << nl;

nN2pEqn.relax();
SolverPerformance<scalar> nN2pPerf = nN2pEqn.solve();

// Print ion density residual
if (enableDetailedLogs && Pstream::master()) Info << "Residual nN2p = " << nN2pPerf.finalResidual() << nl;

nO2mEqn.relax();
SolverPerformance<scalar> nO2mPerf = nO2mEqn.solve();

// Print ion density residual
if (enableDetailedLogs && Pstream::master()) Info << "Residual nO2m = " << nO2mPerf.finalResidual() << nl;

if (0)
{
    // Experimental
    // suppress un-physical negative densities
    // balance the minimum charge: nN2p - nO2m - ne = 0
    minNe = gMin(ne);
    maxNe = gMax(ne);
    minN2 = gMin(nN2p);
    maxN2 = gMax(nN2p);
    minO2 = gMin(nO2m);
    maxO2 = gMax(nO2m);

    dimensionedScalar minNeLimD("minNeLimD", ne.dimensions(), minNeLim);
    ne = 0.5 * (ne + sqrt(sqr(ne) + sqr(minNeLimD)));

    dimensionedScalar minN2LimD("minN2LimD", ne.dimensions(), minN2Lim);
    nN2p = 0.5 * (nN2p + sqrt(sqr(nN2p) + sqr(minN2LimD)));

    dimensionedScalar minO2LimD("minO2LimD", ne.dimensions(), minO2Lim);
    nO2m = 0.5 * (nO2m + sqrt(sqr(nO2m) + sqr(minO2LimD)));
}

// update the reference values
ne.correctBoundaryConditions();
nN2p.correctBoundaryConditions();

if (enableDetailedLogs)
{
    scalar minNe = gMin(ne);
    scalar maxNe = gMax(ne);
    if (0 < maxNe && Pstream::master()) Info << "min/max ne: " << minNe << " " << maxNe
                                             << " min/max: " << minNe / maxNe << nl;
    scalar minN2 = gMin(nN2p);
    scalar maxN2 = gMax(nN2p);
    if (0 < maxN2 && Pstream::master()) Info << "min/max nN2p: " << minN2 << " " << maxN2
                                             << " min/max: " << minN2 / maxN2 << nl;
    scalar minO2 = gMin(nO2m);
    scalar maxO2 = gMax(nO2m);
    if (0 < maxO2 && Pstream::master()) Info << "min/max nO2m: " << minO2 << " " << maxO2
                                             << " min/max: " << minO2 / maxO2 << nl;
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
volScalarField chargeError = nN2p - nO2m - ne;
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
rhoE = eCharge * (nN2p - nO2m - ne);

// calculate particle fluxes
GammaE = U * ne - mu_e * E * ne - De * fvc::grad(ne);
GammaN2p = U * nN2p + mu_N2p * E * nN2p - DN2p * fvc::grad(nN2p);
GammaO2m = U * nO2m - mu_O2m * E * nO2m - DO2m * fvc::grad(nO2m);
Je = - eCharge * GammaE;
JN2p = eCharge * GammaN2p;
JO2m = - eCharge * GammaO2m;
J = Je + JN2p + JO2m;

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

    scalar O2mAnodeI = currentOnPatch(mesh, JO2m, pPatchID);
    scalar O2mCathodeI = currentOnPatch(mesh, JO2m, nPatchID);

    scalar anodeI = eAnodeI + N2pAnodeI + O2mAnodeI;
    scalar cathodeI = eCathodeI + N2pCathodeI + O2mCathodeI;

    if (Pstream::master()) Info << "anode currents (e/N2/O2): " << eAnodeI << " "
                                << N2pAnodeI << " " << O2mAnodeI << nl;
    if (Pstream::master()) Info << "cathode currents (e/N2/O2): " << eCathodeI << " "
                                << N2pCathodeI << " " << O2mCathodeI << nl;
    if (Pstream::master()) Info << "anode/cathode currents: " << anodeI << " " << cathodeI << nl;
}

if (1)
{
    scalar minRhoE = gMin(rhoE);
    scalar maxRhoE = gMax(rhoE);
    scalar maxMagRhoE = max(mag(minRhoE ), mag(maxRhoE ));

    volScalarField magDRhoEDt
    (
        "magDRhoEDt",
        mag(rhoE - rhoEOld)/runTime.deltaTValue()
    );

    scalar maxDRhoEDt = gMax(magDRhoEDt);

    if (0 < maxMagRhoE)
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

