/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    ehdThrustFoam

Description
    Solver for ehd thruster with electrostatics and transient turbulent flow of
    incompressible isothermal fluids

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"

#include "volFields.H"
#include "surfaceFields.H"
#include "fvMesh.H"
#include "fvMatrices.H"

#include "fvc.H"
#include "fvcFlux.H"
#include "fvcSnGrad.H"
#include "fvcDdt.H"

#include "fvm.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"

#include "viscosityModel.H"
#include "incompressibleMomentumTransportModel.H"
#include "adjustPhi.H"
#include "constrainHbyA.H"
#include "pimpleControl.H"
#include "functionObjectList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "createFields.H"
    #include "createSpecies.H"

    // support for functionObjectList.H
    Foam::functionObjectList functions(runTime);
    functions.start();

    // support for continuityErrs.H
    scalar cumulativeContErr = 0.0;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting iteration loop\n" << nl;

    // depending on the environment
    scalar factMulti = 8.0;
    int pwr = -10;
    scalar factor = factMulti * pow(10.0, pwr);
    bool factorChange = false;

    // maxEeCo above 0.3: the simulation is likely to terminate early
    scalar ecoHyperThresh = 0.3;
    // maxEeCo above 0.2: reducing factor seems to reduce maxEeCo 0.2 -> 0.1
    scalar ecoUpperThresh = 0.2;
    // ecoLowerThresh lower than the expected normal reduction caused by factor changes
    scalar ecoLowerThresh = 0.08;
    scalar maxEeCo = (ecoUpperThresh + ecoLowerThresh) / 2;
    scalar maxEeCo_old = maxEeCo;
    // eeCoRateLimit is likely to be exceeded in the cycle after the factor is increased
    scalar eeCoRateLimit = 0.5;

    scalar dtUpperLimit = 1e-6;

    scalar minNeLimI = 1.0e-2;
    scalar minN2LimI = 1.0e-2;
    scalar minNeLim = minNeLimI;
    scalar minN2Lim = minN2LimI;

    scalar ratioThr = -0.01;
    scalar recoveryRatio = -0.00001;
    scalar minRatioThr = ratioThr;

    Info<< "currentTime = " << runTime.name() << nl;
    Info<< "endTime     = " << runTime.endTime().value() << nl;
    Info<< "deltaT      = " << runTime.deltaTValue() << nl;
    Info<< "runTime.run() = " << runTime.run() << nl;
    Info<< "runTime.loop() first check = " << runTime.loop() << nl;

    int nInitialIterations = 50;
    scalar PPhiE_old = 0;
    while (runTime.loop())
    {
        int iterPerLogs = 1000;
        int enableDetailedLogs = !(runTime.timeIndex() % iterPerLogs);
        if (nInitialIterations || runTime.timeIndex() < 10)
        {
            enableDetailedLogs = true;
        }

        if (1)
        {
            label ppPatchID = mesh.boundaryMesh().findIndex("PELEMENT");
            if (ppPatchID < 0)
            {
                FatalErrorInFunction
                    << "Patch PELEMENT not found"
                    << exit(FatalError);
            }

            const fvPatchScalarField& pphiEpatch = phiE.boundaryField()[ppPatchID];
            scalar PPhiE = gMax(pphiEpatch);
            if (!nInitialIterations && (1e-10 < mag(PPhiE - PPhiE_old)))
            {
                if (Pstream::master()) Info << "PPhiE: " << PPhiE << " PPhiE - PPhiE_old: " << PPhiE - PPhiE_old << nl;
                nInitialIterations = 50;
                minRatioThr = ratioThr;
                enableDetailedLogs = true;
            }

            PPhiE_old = PPhiE;
        }

        scalar minNe = gMin(ne);
        scalar maxNe = gMax(ne);
        scalar minN2 = gMin(nN2p);
        scalar maxN2 = gMax(nN2p);
        if (minNeLim == minNeLimI && 0 < maxNe && 0 < maxN2)
        {
            // Experimental
            minNeLim = minNeLimI * 1e-10;
            minN2Lim = minN2LimI * 1e-10;
        }

        if (0 < maxNe && 0 < maxN2)
        {
            // Experimental
            scalar old = minRatioThr;
            // hysteresis for density correction recoveryRatio -> minRatioThr
            // This is an experiment to reduce negative densities
            minRatioThr = ( recoveryRatio < (minNe / maxNe)
                               && recoveryRatio < (minN2 / maxN2) )
                ? ratioThr : minRatioThr;

            if (old != minRatioThr && ratioThr == minRatioThr && Pstream::master())
                Info << runTime.timeIndex() << ": min/max ne: " << minNe / maxNe << " nN2p: " << minN2 / maxN2
                     << nl;

            // strong clamping
            minRatioThr = ( (minNe / maxNe) < -0.1
                               || (minN2 / maxN2) < -0.1 )
                ? ratioThr : minRatioThr;
        }

        bool potCorrection = true;
        while (potCorrection)
        {
            potCorrection = false;

            if (nInitialIterations)
            {
                // Experimental
                // suppress un-physical negative densities
                // balance the minimum charge: nN2p - ne = 0
                if (enableDetailedLogs && 0 < maxNe && Pstream::master()) Info << "before clamp min/max ne: " << minNe << " " << maxNe
                                                                               << " min/max: " << minNe / maxNe << nl;
                dimensionedScalar minNeLimD("minNeLimD", ne.dimensions(), minNeLim);
                ne = 0.5 * (ne + sqrt(sqr(ne) + sqr(minNeLimD)));

                if (enableDetailedLogs && 0 < maxN2 && Pstream::master()) Info << "before clamp min/max nN2p: " << minN2 << " " << maxN2
                                                                               << " min/max: " << minN2 / maxN2 << nl;
                dimensionedScalar minN2LimD("minN2LimD", ne.dimensions(), minN2Lim);
                nN2p = 0.5 * (nN2p + sqrt(sqr(nN2p) + sqr(minN2LimD)));

                // implement hysteresis using an arbitarily large negative number
                minRatioThr = -VGREAT;
            }

            if (nInitialIterations)
            {
                scalar minPhiE = gMax(phiE);
                scalar maxPhiE = gMin(phiE);
                if (Pstream::master()) Info << "Initialising phiE solving the Poisson equation at time " << runTime.name()
                                            << ": phiE extent = " << mag(maxPhiE - minPhiE) << nl;

                for (int i = 0; i < nNonOrthogonalPotCorrectors && nInitialIterations; i++)
                {
                    fvScalarMatrix phiInitEqn
                    (
                        fvm::laplacian(phiE)
                        ==
                      - rhoE / epsilon0
                    );

                    SolverPerformance<scalar> pPerf = phiInitEqn.solve();
                    nInitialIterations = pPerf.nIterations();
                }

                // ensure E is well behaved at startup
                phiE.correctBoundaryConditions();
                E = -fvc::grad(phiE);
                E.correctBoundaryConditions();

                minPhiE = gMin(phiE);
                maxPhiE = gMax(phiE);
                if (Pstream::master()) Info << "phiE initialised: min/max = " << minPhiE << " / " << maxPhiE
                                            << ": phiE extent = " << mag(maxPhiE - minPhiE) << nl;
            }

            if (2 < runTime.timeIndex())
            {
                // allow time for simulation to settle before
                // calculating the drift CourantNo
                #include "mpECourantNo.H"
            }

            if (1) {
                #include "calcReactionRate.H"
                scalar eeCoRate = (maxEeCo - maxEeCo_old) / maxEeCo_old;
                scalar factorCh_old = factorChange;
                factorChange = false;
                if (ecoLowerThresh < maxEeCo_old && eeCoRateLimit < eeCoRate)
                {
                    nInitialIterations = 50;
                    minRatioThr = ratioThr;
                    potCorrection = true;
                    pwr -= 2;
                    factor = factMulti * pow(10.0, pwr);
                    factorChange = true;
                }
                else if (factorCh_old && ecoUpperThresh < maxEeCo && eeCoRateLimit < eeCoRate)
                {
                    // previous change has had little effect
                    // reaction rate is changing rapidly still
                    nInitialIterations = 50;
                    factMulti -= 1.0;
                    if (factMulti < 1.0)
                    {
                        pwr -= 1;
                        factMulti = 9.0;
                    }

                    factor = factMulti * pow(10.0, pwr);
                    factorChange = true;
                }
                else if (ecoHyperThresh < maxEeCo)
                {
                    pwr -= 1;
                    factor = factMulti * pow(10.0, pwr);
                    factorChange = true;
                }
                else if (ecoUpperThresh < maxEeCo)
                {
                    factMulti -= 1.0;
                    if (factMulti < 1.0)
                    {
                        pwr -= 1;
                        factMulti = 9.0;
                    }

                    factor = factMulti * pow(10.0, pwr);
                    factorChange = true;
                }
                else if (maxEeCo < ecoLowerThresh)
                {
                    factMulti += 1.0;
                    if (9.0 < factMulti)
                    {
                        pwr += 1;
                        factMulti = 1.0;
                    }

                    factor = factMulti * pow(10.0, pwr);
                    factorChange = true;
                }

                if (factorChange)
                {
                    enableDetailedLogs = true;
                    if (ratioThr == minRatioThr && 0 < maxNe && 0 < maxN2
                        && ( (minNe / maxNe) < minRatioThr
                             || (minN2 / maxN2) < minRatioThr ))
                    {
                        // Experimental
                        nInitialIterations = 50;
                        minRatioThr = ratioThr;
                        potCorrection = true;
                    }
                }

                if (Pstream::master() && enableDetailedLogs)
                    Info << "maxEeCo: " << maxEeCo << " eeCoRate: " << eeCoRate
                         << " factMulti: " << factMulti << " pwr: " << pwr
                         << " prev maxEeCo: " << maxEeCo_old
                         << " minRatioThr: " << minRatioThr << nl;
                maxEeCo_old = maxEeCo;

                // Now update the DeltaT value
                scalar Vmax = gMax(mesh.V());
                if (Pstream::master() && enableDetailedLogs) Info << "Vmax: " << Vmax << nl;
                scalar meshDeltaX = min(Foam::exp(Foam::log(Vmax)/3.0), GREAT);
                if (Pstream::master() && enableDetailedLogs) Info << "meshDeltaX: " << meshDeltaX << nl;
                volScalarField magUe
                (
                    "magUe",
                    mag(mu_e * E)
                );
                scalar maxDriftVelocity = gMax(magUe);
                if (Pstream::master() && enableDetailedLogs) Info << "maxDriftVelocity: " << maxDriftVelocity << nl;
                scalar dt1 = meshDeltaX / (50 * (maxDriftVelocity + SMALL));
                if (Pstream::master() && enableDetailedLogs) Info << "dt1: " << dt1 << nl;

                maxReactionRate = max(maxReactionRate, SMALL);
                scalar dt2 = factor / maxReactionRate;
                if (Pstream::master() && enableDetailedLogs) Info << "dt2: " << dt2 << nl;
                scalar deltaT = min(dt1, dt2);
                deltaT = min(deltaT, dtUpperLimit);
                if (Pstream::master() && enableDetailedLogs) Info << "computed new deltaT: " << deltaT << nl;

                runTime.setDeltaT(deltaT);
            }
        }

        if (enableDetailedLogs)
        {
            #include "mpCourantNo.H"
            // #include "setDeltaT.H" : instead opt for a dynamic DeltaT based on the reaction rate
        }

        if (Pstream::master() && enableDetailedLogs) Info<< "Iteration = " << runTime.name() << " index: " << runTime.timeIndex() << nl << nl;
        while (pimple.loop())
        {
            for (int corr=0; corr<nPhiECorrectors; corr++)
            {
                // update species densities and calculate new rhoE
                #include "calcReactionRate.H"
                #include "speciesEqns.H"

                // update phiE using updated rhoE
                #include "poissonEqn.H"

                // Update E
                E = -fvc::grad(phiE);

                if (0)
                {
                    // TEMP work out how to remove this
                    // clamp E
                    scalar Emax = 1e6;
                    forAll(E, i)
                    {
                        scalar magE = mag(E[i]);
                        if (magE > Emax)
                        {
                            E[i] *= Emax / magE;
                        }
                    }
                }

                E.correctBoundaryConditions();
                if (enableDetailedLogs)
                {
                    label nePatchID = mesh.boundaryMesh().findIndex("NELEMENT");
                    if (nePatchID < 0)
                    {
                        FatalErrorInFunction
                            << "Patch NELEMENT not found"
                            << exit(FatalError);
                    }

                    const fvPatchVectorField& nEpatch = E.boundaryField()[nePatchID];
                    const tmp<vectorField> tnHat = mesh.boundary()[nePatchID].nf();
                    const vectorField& nHat = tnHat();
                    scalarField nEn = nEpatch & nHat;
                    scalar minNEn = gMin(nEn);
                    scalar maxNEn = gMax(nEn);
                    if (Pstream::master()) Info << "min/max En at NELEMENT boundary: " << minNEn << " " << maxNEn << nl;

                    label pePatchID = mesh.boundaryMesh().findIndex("PELEMENT");
                    if (pePatchID < 0)
                    {
                        FatalErrorInFunction
                            << "Patch PELEMENT not found"
                            << exit(FatalError);
                    }

                    const fvPatchVectorField& pEpatch = E.boundaryField()[pePatchID];
                    const tmp<vectorField> tpHat = mesh.boundary()[pePatchID].nf();
                    const vectorField& pHat = tpHat();
                    scalarField pEn = pEpatch & pHat;
                    scalar minPEn = gMin(pEn);
                    scalar maxPEn = gMax(pEn);
                    if (Pstream::master()) Info << "min/max En at PELEMENT boundary: " << minPEn << " " << maxPEn << nl;
                }

                if (enableDetailedLogs)
                {
                    volScalarField magE
                    (
                        "magE",
                        mag(E)
                    );

                    scalar minMagE = gMin(magE);
                    scalar maxMagE = gMax(magE);
                    if (Pstream::master()) Info << "min/max magE: " << minMagE << " " << maxMagE << nl;

                    scalar minPhiE = gMin(phiE);
                    scalar maxPhiE = gMax(phiE);
                    if (Pstream::master()) Info << "phiE min/max = " << minPhiE << " / " << maxPhiE << nl;
                }
            }

            // update U
            #include "momentumEqns.H"

            // update rAU in preparation for solving for pressure
            #include "pressureEqns.H"

            while (pimple.correct())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);

                while (pimple.correctNonOrthogonal())
                {
                    pEqn.solve();

                    if (pimple.finalNonOrthogonalIter())
                    {
                        // calculate phi
                        phi = phiHbyA - pEqn.flux();
                    }
                }

                // recompute HbyA for better stability
                HbyA = rAU * UEqn.H();

                // correct velocity with updated rAU
                U = HbyA - rAU * fvc::grad(p);
                U.correctBoundaryConditions();

                // update the momentum with the updated U
                momentumTransport->correct();

                if (enableDetailedLogs)
                {
                    #include "continuityErrs.H"
                    volScalarField magU = mag(U);
                    scalar maxMagU = gMax(magU);
                    if (Pstream::master()) Info << "Pressure corrector: max(U) = " << maxMagU << nl;
                }
            }

            // recompute phi with updated U
            phi = fvc::flux(U);

            if (0)
            {
                // TEMP work out how to remove this
                // clamp U
                scalar Umax = 2000;
                forAll(U, i)
                {
                    scalar magUi = mag(U[i]);
                    if (magUi > Umax)
                    {
                        U[i] *= Umax / magUi;
                    }
                }

                Info << "after clamp min/max U: " << min(mag(U.internalField())).value() << " " << max(mag(U.internalField())).value() << nl;

                // update the momentum with the updated U
                momentumTransport->correct();
            }
        }

        // support for functionObjectList.H
        functions.execute();

        runTime.write();

        if (Pstream::master() && enableDetailedLogs) Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << nl;
    }

    Info<< nl << "End" << nl << nl;
    return 0;
}


// ************************************************************************* //
