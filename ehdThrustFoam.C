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

    scalar dtUpperLimit = 1e-6;
    List<scalar> dtLowerLimits(11, 1.25892541179417e-9);
    List<scalar> reactionRateLimits(11, 7.94328234724283e+0);
    dtLowerLimits[9] = 1.58489319246111e-9;
    reactionRateLimits[9] = 6.30957344480194e+0;
    dtLowerLimits[8] = 1.99526231496888e-9;
    reactionRateLimits[8] = 5.01187233627273e+0;
    dtLowerLimits[7] = 2.51188643150958e-9;
    reactionRateLimits[7] = 3.98107170553498e+0;
    dtLowerLimits[6] = 3.16227766016838e-9;
    reactionRateLimits[6] = 3.16227766016838e+0;
    dtLowerLimits[5] = 3.98107170553497e-9;
    reactionRateLimits[5] = 2.51188643150958e+0;
    dtLowerLimits[4] = 5.01187233627272e-9;
    reactionRateLimits[4] = 1.99526231496888e+0;
    dtLowerLimits[3] = 6.30957344480193e-9;
    reactionRateLimits[3] = 1.58489319246111e+0;
    dtLowerLimits[2] = 5.56029764306997e-9;
    reactionRateLimits[2] = 1.25892541179417e+0;
    dtLowerLimits[1] = 5e-9;
    reactionRateLimits[1] = 1e+0;
    dtLowerLimits[0] = 3.7767762353825e-9;
    reactionRateLimits[0] = 7.94328234724282e-1;
    scalar reactionRateScale = 1.258925411794;

    Info<< "currentTime = " << runTime.name() << nl;
    Info<< "endTime     = " << runTime.endTime().value() << nl;
    Info<< "deltaT      = " << runTime.deltaTValue() << nl;
    Info<< "runTime.run() = " << runTime.run() << nl;
    Info<< "runTime.loop() first check = " << runTime.loop() << nl;

    int nInitialIterations = 50;
    while (runTime.loop())
    {
        if (nInitialIterations)
        {
            // Balance the charge at the start of the run: nN2p - ne = 0
            Info << "internal before clamp min/max ne: " << min(ne.internalField()).value() << " " << max(ne.internalField()).value() << endl;
            ne = max(ne, dimensionedScalar("minNe", ne.dimensions(), 9.8e-3));
            ne.correctBoundaryConditions();
            Info << "internal before clamp min/max nN2p: " << min(nN2p.internalField()).value() << " " << max(nN2p.internalField()).value() << endl;
            nN2p = max(nN2p, dimensionedScalar("minN2p", nN2p.dimensions(), 1.0e-2));
            nN2p.correctBoundaryConditions();

            rhoE = eCharge * (nN2p - ne);
        }

        if (nInitialIterations)
        {
            Info << "Initialising phiE solving the Laplace equation at time 0: phiE extent = "
                 << mag(gMax(phiE) - gMin(phiE)) << endl;

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

            Info << "phiE initialised: min/max = "
                 << min(phiE).value() << " / "
                 << max(phiE).value() << endl;
        }

        #include "CourantNo.H"
        // #include "setDeltaT.H" : instead opt for a dynamic DeltaT based on the reaction rate
        if (1) {
            #include "calcReactionRate.H"
            scalar Vmax = gMax(mesh.V().primitiveField());
            Info << "Vmax: " << Vmax << endl;
            scalar meshDeltaX = min(Foam::exp(Foam::log(Vmax)/3.0), GREAT);
            Info << "meshDeltaX: " << meshDeltaX << endl;
            volScalarField magUe
            (
                "magUe",
                mag(mu_e * E)
            );
            scalar maxDriftVelocity = max(magUe).value();
            Info << "maxDriftVelocity: " << maxDriftVelocity << endl;
            scalar dt1 = meshDeltaX / (50 * (maxDriftVelocity + SMALL));
            Info << "dt1: " << dt1 << endl;

            scalar reactionRateMax = reactionRateLimits[0];
            scalar dtLowerLimit = dtLowerLimits[0];
            for (int i = 1; i < dtLowerLimits.size(); i++)
            {
                reactionRateMax = reactionRateLimits[i] / reactionRateScale < maxReactionRate?
                    reactionRateLimits[i] : reactionRateMax;
                dtLowerLimit = reactionRateMax == reactionRateLimits[i]?
                    dtLowerLimits[i] : dtLowerLimit;
            }

            Info << "dtLowerLimit: " << dtLowerLimit
                 << " reactionRateMax: " << reactionRateMax << endl;

            maxReactionRate = max(maxReactionRate, SMALL);
            maxReactionRate = min(maxReactionRate, reactionRateMax);
            scalar dt2 = dtLowerLimit * reactionRateMax / maxReactionRate;
            Info << "dt2: " << dt2 << endl;
            scalar deltaT = max(min(dt1, dt2), dtLowerLimit);
            deltaT = min(deltaT, dtUpperLimit);
            Info << "computed new deltaT: " << deltaT << endl;

            runTime.setDeltaT(deltaT);
        }

        Info<< "Iteration = " << runTime.name() << nl << endl;
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
                if (1)
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
                    Info << "min/max En at NELEMENT boundary: " << gMin(nEn) << " " << gMax(nEn) << endl;

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
                    Info << "min/max En at PELEMENT boundary: " << gMin(pEn) << " " << gMax(pEn) << endl;
                }


                {
                    surfaceScalarField phiEFlux_e = fvc::flux(-mu_e * E);
                    surfaceScalarField phiEFlux_N2p = fvc::flux(mu_N2p * E);

                    scalar maxEeCo = 0.0;
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

                    forAll(phiEFlux_N2p, faceI)
                    {
                        scalar flux = mag(phiEFlux_N2p[faceI]);

                        label own = mesh.owner()[faceI];
                        label nei = mesh.neighbour()[faceI];

                        scalar coOwn = flux * dt / mesh.V()[own];
                        scalar coNei = flux * dt / mesh.V()[nei];

                        maxEpCo = max(maxEpCo, max(coOwn, coNei));
                    }

                    scalar totalFlux =
                        gSum(mag(phiEFlux_e)())
                        + gSum(mag(phiEFlux_N2p)());

                    scalar meanECo = (totalFlux * dt) / gSum(mesh.V());

                    Info << "E min/max: " << min(mag(E)).value()
                         << " / " << max(mag(E)).value() << endl;

                    Info << "Ue min/max: " << min(mag(Ue)).value()
                         << " / " << max(mag(Ue)).value() << endl;

                    Info << "phiEFlux_e min/max: "
                         << min(phiEFlux_e).value()
                         << " / " << max(phiEFlux_e).value() << endl;

                    Info << "Mean drift Courant Number = " << meanECo << endl;
                    Info << "Max drift Courant Numbers: electron: " << maxEeCo
                         << " positive ion: " << maxEpCo << endl;
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

                #include "continuityErrs.H"
                Info << "Pressure corrector: max(U) = " << max(mag(U)) << endl;
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

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< nl << "End" << nl << endl;
    return 0;
}


// ************************************************************************* //
