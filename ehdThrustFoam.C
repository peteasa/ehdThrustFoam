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

    scalar factorMulti = 8.0;
    int pwr = -9;
    scalar factor = factorMulti * pow(10.0, pwr);
    scalar ecoUpperThresh = 0.2;
    scalar ecoLowerThresh = 0.08;
    scalar maxEeCo = (ecoUpperThresh + ecoLowerThresh) / 2;

    scalar dtUpperLimit = 1e-6;

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

            if (ecoUpperThresh < maxEeCo)
            {
                factorMulti -= 1.0;
                if (factorMulti < 1.0)
                {
                    pwr -= 1;
                    factorMulti = 9.0;
                }

                factor = factorMulti * pow(10.0, pwr);
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
            }

            Info << "maxEeCo " << maxEeCo << " factorMulti: " << factorMulti << " pwr: " << pwr << endl;

            maxReactionRate = max(maxReactionRate, SMALL);
            scalar dt2 = factor / maxReactionRate;
            Info << "dt2: " << dt2 << endl;
            scalar deltaT = min(dt1, dt2);
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
