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
#include "patchDistWave.H"

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

    if (1 && runTime.value() < runTime.deltaTValue() * 4)
    {
        // create initial conditions only if starting from time zero!

        if (Pstream::master()) Info << "create initial charge density at " << runTime.name() << nl;
        // alternative for this is to use createZones and setFields
        // however setFields only expects geometry defined zones so is not
        // as flexible
        label ppPatchID = mesh.boundaryMesh().findIndex("PELEMENT");

        if (ppPatchID == -1)
        {
            FatalErrorInFunction
                << "Patch PELEMENT not found"
                << exit(FatalError);
        }

        labelHashSet patchIDs;
        patchIDs.insert(ppPatchID);

        scalarField cellDist(mesh.nCells(), GREAT);

        Foam::patchDistWave::calculate
        (
            mesh,
            patchIDs,
            cellDist
        );

        volScalarField distanceToAnode
        (
            IOobject
            (
                "distanceToAnode",
                runTime.name(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimLength, 0.0)
         );

        // copy internal field
        forAll(cellDist, cellI)
        {
            distanceToAnode[cellI] = cellDist[cellI];
        }

        distanceToAnode.correctBoundaryConditions();

        scalar sheathThickness = 200e-6;
        forAll(mesh.C(), cellI)
        {
            if (distanceToAnode[cellI] < sheathThickness)
            {
                nN2p[cellI]  = 1e15;
                ne[cellI]    = 1e15;
            }
        }
    }

    scalar minNeLimI = 1.0e-2;
    scalar minN2LimI = 1.0e-2;
    scalar minNeLim = minNeLimI;
    scalar minN2Lim = minN2LimI;

    // depending on the environment
    scalar factMulti = 5.0;
    int pwr = 8;
    scalar factor = factMulti * pow(10.0, pwr);
    bool factorChange = false;

    /*************************************************************************/
    /*                       Auto adjust parameters                          */
    // manage changes in maxEeCo:
    // divide factor by 10 if (ecoHyperThresh < maxEeCo) to avoid early termination
    // suggest: between 0.3 -> 1.0
    scalar ecoHyperThresh = 1.0;
    // subtract one from the first none-zero digit of factor if (ecoUpperThresh < maxEeCo)
    // suggest: ecoHyperThresh / 2
    scalar ecoUpperThresh = ecoHyperThresh / 2;
    // add one to factor if (maxEeCo < ecoLowerThresh) :- WARNING: avoid maxEeCo bounce by
    // choosing a threshold lower than ecoUpperThresh: preferably lower then the expected
    // reduction in maxEeCo caused by applying subraction of one as a result of maxEeCo
    // exceeding ecoUpperThresh
    // suggest ecoUpperThresh / 10
    scalar ecoLowerThresh = ecoUpperThresh / 10;
    scalar maxEeCo = (ecoUpperThresh + ecoLowerThresh) / 2;
    scalar maxEeCo_old = maxEeCo;
    // monitor the rate of change of maxEeCo:
    // eeCoRateLimit is likely to be exceeded in the cycle after the factor is increased
    // suggest: between 0.5 -> 1.0
    scalar eeCoRateLimit = 1.0;

    // manage changes in density min/max ratios
    // if ne min/max ratio exceeds -0.9 then the simulation is very likely to terminate early
    // suggest: between -0.9 -> -0.1
    scalar ratioHThr = -0.1;
    // clamp the density min value: suggest -0.01
    scalar ratioThr = -0.01;
    // hysteresis value allows the system to recover from a clamping event
    // suggest: ratioThr / 100
    scalar recoveryRatio = ratioThr / 100;
    scalar minRatioThr = ratioThr;
    // when minRCyDec is 0 any threshold violation will be acted on
    int minRCyDec = 0;

    // track the change in rho
    scalar maxDRhoEDtRate = 0;
    scalar maxDRhoEDtRateThr = 2e8;
    // when maxDRhoEDtRateDec is 0 any threshold violation will be acted on
    int maxDRhoEDtRateDec = 0;
    /*************************************************************************/

    Info<< "currentTime = " << runTime.name() << nl;
    Info<< "endTime     = " << runTime.endTime().value() << nl;
    Info<< "deltaT      = " << runTime.deltaTValue() << nl;
    Info<< "runTime.run() = " << runTime.run() << nl;
    Info<< "runTime.loop() first check = " << runTime.loop() << nl;

    // ensure that Poisson Equation is run at appropriate times
    int intervalCount = 0;
    int maxInterval = 12;
    scalar PPhiE_old = 0;
    while (runTime.loop())
    {
        int iterPerLogs = 1000;
        int enableDetailedLogs = !(runTime.timeIndex() % iterPerLogs);
        if (!(intervalCount % maxInterval) || runTime.timeIndex() < 10 || runTime.deltaTValue() < 1e-40)
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
            if ((intervalCount % maxInterval) && (1e-10 < mag(PPhiE - PPhiE_old)))
            {
                if (Pstream::master()) Info << "PPhiE: " << PPhiE << " PPhiE - PPhiE_old: " << PPhiE - PPhiE_old << nl;
                intervalCount = 0;
                minRatioThr = ratioThr;
                enableDetailedLogs = true;
            }

            PPhiE_old = PPhiE;
        }

        scalar minNe = gMin(ne);
        scalar maxNe = gMax(ne);
        scalar minN2 = gMin(nN2p);
        scalar maxN2 = gMax(nN2p);
        if (0 < maxNe && 0 < maxN2)
        {
            // Experimental
            // allow larger negative densities as the max density increases
            minNeLim = maxNe * (-ratioThr) * 1e-2 * minNeLimI;
            minN2Lim = maxN2 * (-ratioThr) * 1e-2 * minN2LimI;
        }

        if (0 < maxNe && 0 < maxN2)
        {
            // Experimental
            scalar minRatioThrOld = minRatioThr;
            // hysteresis for density correction recoveryRatio -> minRatioThr
            // This is an experiment to reduce negative densities
            minRatioThr = ( recoveryRatio < (minNe / maxNe)
                            && recoveryRatio < (minN2 / maxN2) ) ? ratioThr : minRatioThr;
            if ( ratioThr == minRatioThr )
            {
                minRCyDec = ( (minNe / maxNe) < minRatioThr
                             || (minN2 / maxN2) < minRatioThr ) ? minRCyDec - 1 : 0;
            } else {
                minRCyDec = ( (minNe / maxNe) < recoveryRatio
                              || (minN2 / maxN2) < recoveryRatio ) ? minRCyDec - 1 : 0;
            }

            if (minRatioThrOld != minRatioThr && ratioThr == minRatioThr && Pstream::master())
                Info << runTime.timeIndex() << ": re-arm trigger min/max ne: " << minNe / maxNe << " nN2p: " << minN2 / maxN2
                     << nl;

            // strong clamping
            bool clampnow = ( (minNe / maxNe) < ratioHThr
                || (minN2 / maxN2) < ratioHThr );
            if (clampnow)
            {
                minRatioThr = ratioThr;
            }
        }

        if (maxDRhoEDtRateThr < maxDRhoEDtRate)
        {
            if ((maxDRhoEDtRateDec < 990) && 10 * maxDRhoEDtRateThr < maxDRhoEDtRate)
            {
                if (Pstream::master()) Info << "WARNING: correct for very high density changes" << nl;
                maxDRhoEDtRateDec = 0;

                intervalCount = 0;
                minRatioThr = ratioThr;
                enableDetailedLogs = true;
            }

            maxDRhoEDtRateDec--;
        }
        else
        {
            // provide hysteresis to avoid rapid repeat
            maxDRhoEDtRateDec = ((maxDRhoEDtRateDec < 500) && maxDRhoEDtRate < maxDRhoEDtRateThr * 0.9) ? 0 : maxDRhoEDtRateDec - 1;
        }

        bool potCorrection = true;
        while (potCorrection)
        {
            potCorrection = false;

            if (!(intervalCount % maxInterval))
            {
                // Experimental
                // suppress un-physical negative densities
                // balance the minimum charge: nN2p - ne = 0
                if (enableDetailedLogs && 0 < maxNe && Pstream::master()) Info << "before clamp min/max ne: " << minNe << " " << maxNe
                                                                               << " min/max: " << minNe / maxNe
                                                                               << " minNeLim: " << minNeLim << nl;
                dimensionedScalar minNeLimD("minNeLimD", ne.dimensions(), minNeLim);
                ne = 0.5 * (ne + sqrt(sqr(ne) + sqr(minNeLimD)));

                if (enableDetailedLogs && 0 < maxN2 && Pstream::master()) Info << "before clamp min/max nN2p: " << minN2 << " " << maxN2
                                                                               << " min/max: " << minN2 / maxN2
                                                                               << " minN2Lim: " << minN2Lim << nl;
                dimensionedScalar minN2LimD("minN2LimD", ne.dimensions(), minN2Lim);
                nN2p = 0.5 * (nN2p + sqrt(sqr(nN2p) + sqr(minN2LimD)));

                // implement hysteresis using an arbitarily large negative number
                minRatioThr = -VGREAT;

                if (0 < maxNe && 0 < maxN2)
                {
                    // restart density hysteresis count
                    minRCyDec = ( (minNe / maxNe) < recoveryRatio
                                  || (minN2 / maxN2) < recoveryRatio ) ? 3000 : 0;
                }

                // restart rapid change count
                maxDRhoEDtRateDec = 1000;

                // densities likely to have been updated so update rhoE
                rhoE = eCharge * (nN2p - ne);
            }

            if (!(intervalCount % maxInterval))
            {
                intervalCount = 1;
                factorChange = true;

                scalar minPhiE = gMax(phiE);
                scalar maxPhiE = gMin(phiE);
                if (Pstream::master()) Info << "Initialising phiE solving the Poisson equation at time " << runTime.name()
                                            << ": phiE extent = " << mag(maxPhiE - minPhiE) << nl;

                int nIterations = 50;
                for (int i = 0; i < nNonOrthogonalPotCorrectors && nIterations; i++)
                {
                    fvScalarMatrix phiInitEqn
                    (
                        fvm::laplacian(phiE)
                        ==
                      - rhoE / epsilon0
                    );

                    SolverPerformance<scalar> pPerf = phiInitEqn.solve();
                    nIterations = pPerf.nIterations();
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

                if (factorCh_old && ecoUpperThresh < maxEeCo && eeCoRateLimit < eeCoRate)
                {
                    // previous change has had little effect
                    // reaction rate is changing rapidly still
                    if (Pstream::master()) Info << "WARNING: eeCoRate is high" << nl;
                    factMulti -= 1.0;
                    if (factMulti < 1.0)
                    {
                        pwr -= 1;
                        factMulti = 9.0;
                    }

                    factor = factMulti * pow(10.0, pwr);
                    factorChange = true;

                    intervalCount = 0;
                }
                else if (ecoHyperThresh < maxEeCo)
                {
                    if (Pstream::master()) Info << "WARNING: maxEeCo is high" << nl;
                    pwr -= 1;
                    factor = factMulti * pow(10.0, pwr);
                    factorChange = true;
                }
                else if ( minRCyDec && minRCyDec < 2000 && ( (minNe / maxNe) < ratioHThr
                                                || (minN2 / maxN2) < ratioHThr ) )
                {
                    factMulti -= 1.0;
                    if (factMulti < 1.0)
                    {
                        pwr -= 1;
                        factMulti = 9.0;
                    }

                    factor = factMulti * pow(10.0, pwr);
                    factorChange = true;

                    if ((intervalCount % maxInterval) && (minNe / maxNe) < ratioHThr) intervalCount++;
                    if ((intervalCount % maxInterval) && (minN2 / maxN2) < ratioHThr) intervalCount++;

                    scalar minRatioThrOld = minRatioThr;
                    minRatioThr = ratioThr;
                    if (minRatioThrOld != minRatioThr && Pstream::master())
                        Info << runTime.timeIndex() << ": force re-arm trigger min/max ne: " << minNe / maxNe
                             << " nN2p: " << minN2 / maxN2 << nl;
                }
                else if ( minRCyDec < 0 )
                {
                    if (Pstream::master()) Info << "WARNING: correct for negative density" << nl;
                    factMulti -= 1.0;
                    if (factMulti < 1.0)
                    {
                        pwr -= 1;
                        factMulti = 9.0;
                    }

                    factor = factMulti * pow(10.0, pwr);
                    factorChange = true;

                    intervalCount = 0;
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

                    if (intervalCount % maxInterval) intervalCount++;
                }
                else if ((maxDRhoEDtRate < maxDRhoEDtRateThr)
                         && minRCyDec == 0
                         && maxDRhoEDtRateDec == 0
                         && maxEeCo < ecoLowerThresh)
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
                else if ( (maxDRhoEDtRateDec < 990)
                          && maxDRhoEDtRateThr < maxDRhoEDtRate
                          && ( (minNe / maxNe) < ratioHThr
                               || (minN2 / maxN2) < ratioHThr ) )
                {
                    if (Pstream::master()) Info << "WARNING: correct for multiple threshold violation!" << nl;
                    pwr -= 1;
                    factor = factMulti * pow(10.0, pwr);
                    factorChange = true;

                    intervalCount = 0;
                } else if (ecoLowerThresh < maxEeCo_old && eeCoRateLimit < eeCoRate)
                {
                    if (Pstream::master()) Info << "WARNING: eeCoRate is high" << nl;
                    pwr -= 2;
                    factor = factMulti * pow(10.0, pwr);
                    factorChange = true;

                    intervalCount = 0;
                }
                else if (maxDRhoEDtRateDec < 0)
                {
                    if (Pstream::master()) Info << "WARNING: correct for rapid density changes" << nl;
                    enableDetailedLogs = true;

                    intervalCount = 0;
                    minRatioThr = ratioThr;
                    potCorrection = true;
                }

                if (factorChange)
                {
                    enableDetailedLogs = true;
                    if (ratioThr == minRatioThr && 0 < maxNe && 0 < maxN2
                        && ( (minNe / maxNe) < minRatioThr
                             || (minN2 / maxN2) < minRatioThr ) )
                    {
                        // Experimental
                        intervalCount = 0;
                    }

                    if (!(intervalCount % maxInterval))
                    {
                        minRatioThr = ratioThr;
                        potCorrection = true;
                    }
                }

                if (1 && enableDetailedLogs && (ecoHyperThresh < maxEeCo) && Pstream::master())
                    Info << runTime.timeIndex() << ": THRESH extreamly high maxEeCo: " << maxEeCo << nl;
                else if (1 && enableDetailedLogs && (ecoUpperThresh < maxEeCo) && Pstream::master())
                    Info << runTime.timeIndex() << ": THRESH very high maxEeCo: " << maxEeCo << nl;
                if (enableDetailedLogs && 0 < maxNe && (minNe / maxNe) < ratioHThr && Pstream::master())
                    Info << runTime.timeIndex() << ": THRESH high min/max ne: " << minNe / maxNe << nl;
                if (enableDetailedLogs && 0 < maxN2 && (minN2 / maxN2) < ratioHThr && Pstream::master())
                    Info << runTime.timeIndex() << ": THRESH high min/max nN2p: " << minN2 / maxN2 << nl;
                if (enableDetailedLogs && (maxDRhoEDtRateThr < maxDRhoEDtRate) && Pstream::master())
                    Info << runTime.timeIndex() << ": THRESH rapid density changes maxDRhoEDtRate: " << maxDRhoEDtRate << nl;

                if (enableDetailedLogs && Pstream::master())
                    Info << "maxEeCo: " << maxEeCo << " eeCoRate: " << eeCoRate
                         << " factMulti: " << factMulti << " pwr: " << pwr
                         << " intervalCount: " << intervalCount
                         << " mnRCyDec: " << minRCyDec
                         << " mxDRDtRate: " << maxDRhoEDtRate
                         << " mxDRDtRateDec: " << maxDRhoEDtRateDec
                         << " mnRatioThr: " << minRatioThr << nl;
                maxEeCo_old = maxEeCo;

                // Now update the DeltaT value
                scalar Vmax = gMax(mesh.V());
                if (enableDetailedLogs && Pstream::master()) Info << "Vmax: " << Vmax << nl;
                scalar meshDeltaX = min(Foam::exp(Foam::log(Vmax)/3.0), GREAT);
                if (enableDetailedLogs && Pstream::master()) Info << "meshDeltaX: " << meshDeltaX << nl;
                volScalarField magUe
                (
                    "magUe",
                    mag(mu_e * E)
                );
                scalar maxDriftVelocity = gMax(magUe);
                if (enableDetailedLogs && Pstream::master()) Info << "maxDriftVelocity: " << maxDriftVelocity << nl;
                scalar dt1 = meshDeltaX / (50 * (maxDriftVelocity + SMALL));
                if (enableDetailedLogs && Pstream::master()) Info << "dt1: " << dt1 << nl;

                maxReactionRate = max(maxReactionRate, SMALL);
                scalar dt2 = factor / maxReactionRate;
                if (enableDetailedLogs && Pstream::master()) Info << "dt2: " << dt2 << nl;
                scalar deltaT = min(dt1, dt2);
                deltaT = min(deltaT, dtUpperLimit);
                if (runTime.deltaTValue() < 1e-40) enableDetailedLogs = true;
                if (enableDetailedLogs && Pstream::master()) Info << ((runTime.deltaTValue() < 1e-40) ? "WARNING: " : "")
                                                                  << "computed new deltaT: " << deltaT
                                                                  << " maxDRhoEDtRateThr: " << maxDRhoEDtRateThr << nl;

                runTime.setDeltaT(deltaT);
            }
        }

        if (enableDetailedLogs)
        {
            #include "mpCourantNo.H"
            // #include "setDeltaT.H" : instead opt for a dynamic DeltaT based on the reaction rate
        }

        if (enableDetailedLogs && Pstream::master()) Info<< "Iteration = " << runTime.name() << " index: " << runTime.timeIndex() << nl << nl;
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

        if (enableDetailedLogs && Pstream::master()) Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << nl;
    }

    Info<< nl << "End" << nl << nl;
    return 0;
}


// ************************************************************************* //
