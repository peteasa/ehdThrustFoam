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
    scalar factorMulti = 8.0;
    int pwr = -9;
    scalar factor = factorMulti * pow(10.0, pwr);

    // maxEeCo above 0.2: the simulation is likely to terminate early
    scalar ecoHyperThresh = 0.2;
    // maxEeCo above 0.1: reducing factor seems to reduce maxEeCo 0.1 -> 0.05
    scalar ecoUpperThresh = 0.1;
    // ecoLowerThresh lower than the expected normal reduction caused by factor changes
    scalar ecoLowerThresh = 0.04;
    scalar maxEeCo = (ecoUpperThresh + ecoLowerThresh) / 2;
    scalar maxEeCo_old = maxEeCo;
    // eeCoRateLimit is likely to be exceeded in the cycle after the factor is increased
    scalar eeCoRateLimit = 0.5;

    scalar dtUpperLimit = 1e-6;

    Info<< "currentTime = " << runTime.name() << nl;
    Info<< "endTime     = " << runTime.endTime().value() << nl;
    Info<< "deltaT      = " << runTime.deltaTValue() << nl;
    Info<< "runTime.run() = " << runTime.run() << nl;
    Info<< "runTime.loop() first check = " << runTime.loop() << nl;

    int nInitialIterations = 50;
    while (runTime.loop())
    {
        int logsPerIter = 1000;
        int enableDetailedLogs = !(runTime.timeIndex() % logsPerIter);
        if (runTime.timeIndex() < 10) enableDetailedLogs = 1;
        if (nInitialIterations)
        {
            // Balance the charge at the start of the run: nN2p - ne = 0
            if (enableDetailedLogs)
            {
                scalar minNe = gMin(ne);
                scalar maxNe = gMax(ne);
                if (Pstream::master()) Info << "before clamp min/max ne: " << minNe << " " << maxNe << endl;
            }

            ne = max(ne, dimensionedScalar("minNe", ne.dimensions(), 9.8e-3));
            ne.correctBoundaryConditions();

            if (enableDetailedLogs)
            {
                scalar minN2 = gMin(nN2p);
                scalar maxN2 = gMax(nN2p);
                if (Pstream::master()) Info << "before clamp min/max nN2p: " << minN2 << " " << maxN2 << endl;
            }

            nN2p = max(nN2p, dimensionedScalar("minN2p", nN2p.dimensions(), 1.0e-2));
            nN2p.correctBoundaryConditions();

            rhoE = eCharge * (nN2p - ne);
        }

        if (nInitialIterations)
        {
            if (enableDetailedLogs)
            {
                scalar minPhiE = gMax(phiE);
                scalar maxPhiE = gMin(phiE);
                if (Pstream::master()) Info << "Initialising phiE solving the Laplace equation at time 0: phiE extent = "
                                            << mag(maxPhiE - minPhiE) << endl;
            }

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

            if (enableDetailedLogs)
            {
                scalar minPhiE = gMin(phiE);
                scalar maxPhiE = gMax(phiE);
                if (Pstream::master()) Info << "phiE initialised: min/max = " << minPhiE << " / " << maxPhiE << endl;
            }
        }

        if (1) {
            #include "calcReactionRate.H"
            scalar Vmax = gMax(mesh.V());
            if (Pstream::master() && enableDetailedLogs) Info << "Vmax: " << Vmax << endl;
            scalar meshDeltaX = min(Foam::exp(Foam::log(Vmax)/3.0), GREAT);
            if (Pstream::master() && enableDetailedLogs) Info << "meshDeltaX: " << meshDeltaX << endl;
            volScalarField magUe
            (
                "magUe",
                mag(mu_e * E)
            );
            scalar maxDriftVelocity = gMax(magUe);
            if (Pstream::master() && enableDetailedLogs) Info << "maxDriftVelocity: " << maxDriftVelocity << endl;
            scalar dt1 = meshDeltaX / (50 * (maxDriftVelocity + SMALL));
            if (Pstream::master() && enableDetailedLogs) Info << "dt1: " << dt1 << endl;

            scalar eeCoRate = (maxEeCo - maxEeCo_old) / maxEeCo_old;
            if (ecoLowerThresh < maxEeCo_old && eeCoRateLimit < eeCoRate)
            {
                pwr -= 2;
                factor = factorMulti * pow(10.0, pwr);
                enableDetailedLogs = 1;
            }
            else if (ecoHyperThresh < maxEeCo)
            {
                pwr -= 1;
                factor = factorMulti * pow(10.0, pwr);
                enableDetailedLogs = 1;
            }
            else if (ecoUpperThresh < maxEeCo)
            {
                factorMulti -= 1.0;
                if (factorMulti < 1.0)
                {
                    pwr -= 1;
                    factorMulti = 9.0;
                }

                factor = factorMulti * pow(10.0, pwr);
                enableDetailedLogs = 1;
            }
            else if (maxEeCo < ecoLowerThresh)
            {
                factorMulti += 1.0;
                if (9.0 < factorMulti)
                {
                    pwr += 1;
                    factorMulti = 1.0;
                }

                factor = factorMulti * pow(10.0, pwr);
                enableDetailedLogs = 1;
            }

            if (Pstream::master() && enableDetailedLogs) Info << "maxEeCo: " << maxEeCo << " eeCoRate: " << eeCoRate << " factorMulti: " << factorMulti << " pwr: " << pwr << " previous maxEeCo: " << maxEeCo_old << endl;
            maxEeCo_old = maxEeCo;

            maxReactionRate = max(maxReactionRate, SMALL);
            scalar dt2 = factor / maxReactionRate;
            if (Pstream::master() && enableDetailedLogs) Info << "dt2: " << dt2 << endl;
            scalar deltaT = min(dt1, dt2);
            deltaT = min(deltaT, dtUpperLimit);
            if (Pstream::master() && enableDetailedLogs) Info << "computed new deltaT: " << deltaT << endl;

            runTime.setDeltaT(deltaT);
        }

        if (0)
        {
            // Note to use this header re-write for parallel operation
            #include "CourantNo.H"
            // #include "setDeltaT.H" : instead opt for a dynamic DeltaT based on the reaction rate
        }

        if (Pstream::master() && enableDetailedLogs) Info<< "Iteration = " << runTime.name() << " index: " << runTime.timeIndex() << nl << endl;
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
                    if (Pstream::master()) Info << "min/max En at NELEMENT boundary: " << minNEn << " " << maxNEn << endl;

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
                    if (Pstream::master()) Info << "min/max En at PELEMENT boundary: " << minPEn << " " << maxPEn << endl;
                }


                {
                    surfaceScalarField phiEFlux_e = fvc::flux(-mu_e * E);
                    surfaceScalarField phiEFlux_N2p = fvc::flux(mu_N2p * E);

                    maxEeCo = 0.0;
                    scalar maxEpCo = 0.0;
                    const scalar dt = runTime.deltaTValue();

                    forAll(phiEFlux_e, faceI)
                    {
                        scalar flux = mag(phiEFlux_e[faceI]);

                        label own = mesh.owner()[faceI];
                        label nei = mesh.neighbour()[faceI];

                        scalar coOwn = flux * dt / mesh.V()[own];
                        scalar coNei = flux * dt / mesh.V()[nei];

                        maxEeCo = max(maxEeCo, max(coOwn, coNei));
                    }

                    maxEeCo = returnReduce(maxEeCo, maxOp<scalar>());

                    forAll(phiEFlux_N2p, faceI)
                    {
                        scalar flux = mag(phiEFlux_N2p[faceI]);

                        label own = mesh.owner()[faceI];
                        label nei = mesh.neighbour()[faceI];

                        scalar coOwn = flux * dt / mesh.V()[own];
                        scalar coNei = flux * dt / mesh.V()[nei];

                        maxEpCo = max(maxEpCo, max(coOwn, coNei));
                    }

                    maxEpCo = returnReduce(maxEpCo, maxOp<scalar>());

                    scalar totalFlux =
                        gSum(mag(phiEFlux_e)())
                        + gSum(mag(phiEFlux_N2p)());

                    scalar meanECo = (totalFlux * dt) / gSum(mesh.V());

                    if (enableDetailedLogs)
                    {
                        volScalarField magE
                        (
                             "magE",
                             mag(E)
                        );

                        scalar minMagE = gMin(magE);
                        scalar maxMagE = gMax(magE);
                        if (Pstream::master()) Info << "E min/max: " << minMagE << " / " << maxMagE << endl;

                        volScalarField magUe
                        (
                             "magUe",
                             mag(Ue)
                        );

                        scalar minMagUe = gMin(magUe);
                        scalar maxMagUe = gMax(magUe);
                        if (Pstream::master()) Info << "Ue min/max: " << minMagUe << " / " << maxMagUe << endl;

                        scalar minPhiEFlux_e = gMin(phiEFlux_e);
                        scalar maxPhiEFlux_e = gMax(phiEFlux_e);
                        if (Pstream::master()) Info << "phiEFlux_e min/max: "
                                                    << minPhiEFlux_e << " / " << maxPhiEFlux_e << endl;

                        if (Pstream::master()) Info << "Mean drift Courant Number = " << meanECo << endl;
                        if (Pstream::master()) Info << "Max drift Courant Numbers: electron: " << maxEeCo
                                                    << " positive ion: " << maxEpCo << endl;
                    }
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
                    if (Pstream::master()) Info << "Pressure corrector: max(U) = " << maxMagU << endl;
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

                Info << "after clamp min/max U: " << min(mag(U.internalField())).value() << " " << max(mag(U.internalField())).value() << endl;

                // update the momentum with the updated U
                momentumTransport->correct();
            }
        }

        // support for functionObjectList.H
        functions.execute();

        runTime.write();

        if (Pstream::master() && enableDetailedLogs) Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< nl << "End" << nl << endl;
    return 0;
}


// ************************************************************************* //
