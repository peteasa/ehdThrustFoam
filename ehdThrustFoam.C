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
    int pwr = -10;
    scalar factor = factorMulti * pow(10.0, pwr);

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

    scalar minNeLim = 1.0e-2;
    scalar minN2Lim = 1.0e-2;

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

        if (9 < runTime.timeIndex())
        {
            // after startup allow zero densities
            minNeLim = 0;
            minN2Lim = 0;
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
                enableDetailedLogs = true;
            }

            PPhiE_old = PPhiE;
        }

        scalar minNe = gMin(ne);
        scalar maxNe = gMax(ne);
        scalar minN2 = gMin(nN2p);
        scalar maxN2 = gMax(nN2p);
        scalar minThresh = -0.01;
        if (1 && 0 < maxNe && 0 < maxN2
            && ( minNe / maxNe < minThresh
                 || minN2 / maxN2 < minThresh ))
        {
            // Experimental
            nInitialIterations = 50;
        }

        if (nInitialIterations)
        {
            if (1 && 0 < maxNe && 0 < maxN2)
            {
                // Experimental:
                // reduce the unphysical negative density
                scalar minRelax = 0.9;  // 1.0
                minNeLim = 0 < maxNe && minNe / maxNe < minThresh ? minRelax * minNe : minNe;
                minN2Lim = 0 < maxN2 && minN2 / maxN2 < minThresh ? minRelax * minN2 : minN2;
            }

            // First iteration (nNe, nN2p = 0) use arbitary non-zero startup values
            // The thrust is generated by positive ions moving to the cathode at the rear
            // Balance the charge at the start of the run: nN2p - ne = 0
            if (0 < maxNe && Pstream::master()) Info << runTime.timeIndex() << ": before clamp min/max ne: " << minNe << " " << maxNe << " min/max: " << minNe / maxNe << nl;

            ne = max(ne, dimensionedScalar("minNe", ne.dimensions(), minNeLim));
            ne.correctBoundaryConditions();

            if (0 < maxN2 && Pstream::master()) Info << runTime.timeIndex() << ": before clamp min/max nN2p: " << minN2 << " " << maxN2 << " min/max: " << minN2 / maxN2 << nl;

            nN2p = max(nN2p, dimensionedScalar("minN2p", nN2p.dimensions(), minN2Lim));
            nN2p.correctBoundaryConditions();

            rhoE = eCharge * (nN2p - ne);
        }

        bool potCorrection = true;
        while (potCorrection)
        {
            potCorrection = false;
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
                if (ecoLowerThresh < maxEeCo_old && eeCoRateLimit < eeCoRate)
                {
                    nInitialIterations = 50;
                    potCorrection = true;
                    pwr -= 2;
                    factor = factorMulti * pow(10.0, pwr);
                    enableDetailedLogs = true;
                }
                else if (ecoHyperThresh < maxEeCo)
                {
                    pwr -= 1;
                    factor = factorMulti * pow(10.0, pwr);
                    enableDetailedLogs = true;
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
                    enableDetailedLogs = true;
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
                    enableDetailedLogs = true;
                }

                if (Pstream::master() && enableDetailedLogs) Info << "maxEeCo: " << maxEeCo << " eeCoRate: " << eeCoRate << " factorMulti: " << factorMulti << " pwr: " << pwr << " previous maxEeCo: " << maxEeCo_old << nl;
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
