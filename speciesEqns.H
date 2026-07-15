
De = mu_e * Te;
DN2p = mu_N2p * Te;

if (enableDetailedLogs)
{
    scalar minDe = gMin(De.internalField());
    scalar maxDe = gMax(De.internalField());

    scalar minDN2p = gMin(DN2p.internalField());
    scalar maxDN2p = gMax(DN2p.internalField());

    if (Pstream::master()) Info << "De min/max: " << minDe << " " << maxDe << nl
                            << "DN2p min/max: " << minDN2p << " " << maxDN2p << nl;
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
        if (enableDetailedLogs && Pstream::master()) Info << "Slow Start: Separate convection interaction @ time: " << runTime.name() << nl;
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

//Info << "phi_de dimensions: " << phi_e.dimensions() << endl;

surfaceScalarField phi_dN2p
(
    "phi_dN2p",
    fvc::flux(mu_N2p * E)
);

//Info << "De dimensions: " << De.dimensions() << endl;

// combine with the diffusion terms and the source terms
fvScalarMatrix neEqn
(
      fvm::ddt(ne)
    + fvm::div(phi_gas, ne)
    + fvm::div(phi_de, ne)
    - fvm::laplacian(De, ne)
 ==
    Sion
);
// neEqn.solve();

fvScalarMatrix nN2pEqn
(
      fvm::ddt(nN2p)
    + fvm::div(phi_gas, nN2p)
    + fvm::div(phi_dN2p, nN2p)
    - fvm::laplacian(DN2p, nN2p)
 ==
      Sion
);

// nN2pEqn.solve();

// Solve species
neEqn.relax();
SolverPerformance<scalar> nePerf = neEqn.solve();

// Print electron density residual
if (enableDetailedLogs && Pstream::master()) Info << "Residual ne = " << nePerf.finalResidual() << nl;

nN2pEqn.relax();
SolverPerformance<scalar> nN2pPerf = nN2pEqn.solve();

// Print ion density residual
if (enableDetailedLogs && Pstream::master()) Info << "Residual nN2p = " << nN2pPerf.finalResidual() << nl;

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

// calculate particle fluxes
GammaE = U * ne - mu_e * E * ne - De * fvc::grad(ne);
GammaN2p = U * nN2p + mu_N2p * E * nN2p - DN2p * fvc::grad(nN2p);
Je = - eCharge * GammaE;
JN2p = eCharge * GammaN2p;
J = Je + JN2p;

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

    scalar anodeI = eAnodeI + N2pAnodeI;
    scalar cathodeI = eCathodeI + N2pCathodeI;

    if (Pstream::master()) Info << "anode currents (e/N2): " << eAnodeI << " "
                                << N2pAnodeI << nl;
    if (Pstream::master()) Info << "cathode currents (e/N2): " << eCathodeI << " "
                                << N2pCathodeI << nl;
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

