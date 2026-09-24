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
#include "constrainPressure.H"
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

    if (Pstream::master()) Info << "\nStarting iteration loop\n" << nl;

    scalar dtUpperLimit = 1e-6;

    if (1 && runTime.value() < runTime.deltaTValue() * 4)
    {
        // create initial conditions only if starting from time zero!
        scalar sheathThickness = 200e-6;
        scalar initNe = 1e4;
        scalar initNN2p = 1e11;
        scalar initNO2m = 1e4;
        #include "initDensity.H"
    }

    // arbitary non-zero startup values
    // if ne collapses during the simulation then make this as small
    // as possible because at least this charge will be assigned to
    // every cell in the model whilst correcting for negative densities
    scalar minNeLimI = 1e-6;
    scalar minN2pLimI = 1e-6;
    scalar minO2mLimI = 1e-6;
    scalar minNeLim = minNeLimI;
    scalar minN2pLim = minN2pLimI;
    scalar minO2mLim = minO2mLimI;

    // depending on the environment
    scalar factMulti = 5.0;
    int pwr = -2;
    scalar factor = factMulti * pow(10.0, pwr);
    bool factorChange = false;
    bool factorCh_old = factorChange;

    /*************************************************************************/
    // startup time - separate electrostatic and convection interaction
    scalar startupT = 1.5e-6;
    scalar startupMaxCount = 4e3;
    scalar startupIncrement = 1 / startupMaxCount;

    scalar startupPHalfwidth = 0.01;
    scalar runHPFrac = 10;
    scalar startupPMaxCount = 1e3;
    scalar startupPCount = 1;

    scalar runLPFrac = 1.0 - 2.0 * startupPHalfwidth;
    if (1.0 <= startupPHalfwidth && Pstream::master()) Info << "ERROR: startupPHalfwidth must be less than 1.0" << nl;
    if (1.0 <= (startupPHalfwidth + runLPFrac) && Pstream::master())
        Info << "ERROR: runLPFrac must be less than " << 1.0 - startupPHalfwidth << nl;
    if (runHPFrac + startupPHalfwidth < 1.0 && Pstream::master())
        Info << "ERROR: runHPFrac should be greater than " << 1.0 - startupPHalfwidth << nl;

    /*************************************************************************/
    /*                       Auto adjust parameters                          */
    // manage changes in maxEeCo:
    // divide factor by 10 if (ecoHyperThresh < maxEeCo) to avoid early termination
    // suggest: between 0.3 -> 1.0
    scalar ecoHyperInitial = 2.0;
    scalar ecoHyperRunning = 0.5;
    scalar ecoHyperThresh = ecoHyperInitial;
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
    scalar eeCoRateInitial = 1000.0;
    scalar eeCoRateRunning = 0.9;
    scalar eeCoRateLimit = eeCoRateInitial;
    scalar eeCoRunTime = 5.05e-12;

    // manage changes in density min/max ratios
    // hysteresis value allows the system to recover from a clamping event
    // to disable use high negative value
    const scalar recoveryRatio = -5e-6;
    // limit the dynamic recoveryRatio to negative numbers
    const scalar limRecoveryRatio = recoveryRatio / 50;
    scalar NeTargetMin = max(minNeLimI, minNeLim);
    scalar N2pTargetMin = max(minN2pLimI, minN2pLim);
    scalar O2mTargetMin = max(minO2mLimI, minO2mLim);
    // initial the dynamic recoveryRatio values
    scalar neRecoveryRatio = recoveryRatio;
    scalar N2pRecoveryRatio = recoveryRatio;
    scalar O2mRecoveryRatio = recoveryRatio;
    // clamp the density min value
    const scalar densityMulti = 100;
    const scalar ratioThr = recoveryRatio * densityMulti;
    scalar denRatioThr = ratioThr;
    // initial the dynamic ratioThr values
    scalar neRatioThr = ratioThr;
    scalar N2pRatioThr = ratioThr;
    scalar O2mRatioThr = ratioThr;

    // if ne min/max ratio exceeds say -1.0 then the simulation is very likely to terminate early
    // initial the dynamic high negative ratio threshold values
    const scalar densityMultiHigh = 10000;
    scalar neRatioHThr = recoveryRatio * densityMultiHigh;
    scalar N2pRatioHThr = recoveryRatio * densityMultiHigh;
    scalar O2mRatioHThr = recoveryRatio * densityMultiHigh;

    // when minRCyDec is 0 any threshold violation will be acted on
    int minRCyDec = 0;

    // track the change in rho
    scalar maxDRhoEDtRate = 0;
    // to disable use high value
    scalar maxDRhoEDtRateThrInitial = 2e80;
    scalar maxDRhoEDtRateThrRunning = 5e8;
    scalar maxDRhoEDtRateRunTime = 1.4e-6;
    scalar maxDRhoEDtRateThr = maxDRhoEDtRateThrInitial;
    // when maxDRhoEDtRateDec is 0 any threshold violation will be acted on
    int maxDRhoEDtRateDec = 0;

    // track the change in E/N
    scalar maxDENTdDtRateThr = 1e9;
    /*************************************************************************/

    if (Pstream::master()) Info << "currentTime   = " << runTime.name() << nl;
    if (Pstream::master()) Info << "endTime       = " << runTime.endTime().value() << nl;
    if (Pstream::master()) Info << "deltaT        = " << runTime.deltaTValue() << nl;
    if (Pstream::master()) Info << "negative density recoveryRatio: " << recoveryRatio
                                << " ratioThr: " << ratioThr
                                << " RatioHThr: " << neRatioHThr << nl;

    // ensure that Poisson Equation is not run too frequently
    int maxInterval = 12;
    int intervalCount = maxInterval - 1;
    scalar PPhiE_old = 0;
    while (runTime.loop())
    {
        int iterPerLogs = 1000;
        int enableDetailedLogs = !(runTime.timeIndex() % iterPerLogs);
        if (!(intervalCount % maxInterval) || runTime.timeIndex() < 10 || runTime.deltaTValue() < 1e-40
            || (startupPHalfwidth < startupLPFrac && startupLPFrac < startupPHalfwidth + 0.005)
            || (runLPFrac - 0.005 < startupLPFrac && startupLPFrac < runLPFrac ))
        {
            enableDetailedLogs = true;
        }

        if (runTime.value() < min(maxDRhoEDtRateRunTime, startupT))
        {
            if (enableDetailedLogs && Pstream::master()) Info << "Slow Start: using initial maxDRhoEDtRate threshold until time: " << maxDRhoEDtRateRunTime << nl;
            maxDRhoEDtRateThr = maxDRhoEDtRateThrInitial;
        }
        else
        {
            maxDRhoEDtRateThr = maxDRhoEDtRateThrRunning;
        }

        if (runTime.value() < min(eeCoRunTime, startupT))
        {
            if (enableDetailedLogs && Pstream::master()) Info << "Slow Start: using initial eeCoRate until time: " << eeCoRunTime << nl;
            ecoHyperThresh = ecoHyperInitial;
            ecoUpperThresh = ecoHyperThresh / 2;
            ecoLowerThresh = ecoUpperThresh / 10;
            eeCoRateLimit = eeCoRateInitial;
        }
        else
        {
            ecoHyperThresh = ecoHyperRunning;
            ecoUpperThresh = ecoHyperThresh / 2;
            ecoLowerThresh = ecoUpperThresh / 10;
            eeCoRateLimit = eeCoRateRunning;
        }

        if (1)
        {
            const fvPatchScalarField& pphiEpatch = phiE.boundaryField()[pPatchID];
            scalar PPhiE = gMax(pphiEpatch);
            if ((intervalCount % maxInterval) && (1e-6 < mag((PPhiE - PPhiE_old)/PPhiE)))
            {
                if (Pstream::master()) Info << runTime.timeIndex() << ": PPhiE: " << PPhiE << " PPhiE - PPhiE_old: " << PPhiE - PPhiE_old << nl;

                intervalCount = 0;
                denRatioThr = ratioThr;
                neRatioThr = neRecoveryRatio * densityMulti;
                N2pRatioThr = N2pRecoveryRatio * densityMulti;
                O2mRatioThr = O2mRecoveryRatio * densityMulti;
            }

            PPhiE_old = PPhiE;
        }

        // calcReactionRate.H may be included multiple times
        // maxDENTdDtRate is calculated once
        scalar maxDENTdDtRate = 0;

        // sample the densities
        scalar minNe = gMin(ne);
        scalar maxNe = gMax(ne);
        scalar minN2p = gMin(nN2p);
        scalar maxN2p = gMax(nN2p);
        scalar minO2m = gMin(nO2m);
        scalar maxO2m = gMax(nO2m);
        if (0 && 0 < maxNe && 0 < maxN2p && 0 < maxO2m)
        {
            // Experimental
            // allow larger negative densities as the max density increases
            minNeLim = maxNe * (-ratioThr) * 1e-2 * minNeLimI;
            minN2pLim = maxN2p * (-ratioThr) * 1e-2 * minN2pLimI;
            minO2mLim = maxO2m * (-ratioThr) * 1e-2 * minO2mLimI;
        }

        if (0 < maxNe && 0 < maxN2p && 0 < maxO2m)
        {
            // Experimental
            scalar denRatioThrOld = denRatioThr;

            NeTargetMin = min(max(minNeLimI, minNeLim), maxNe);
            N2pTargetMin = min(max(minN2pLimI, minN2pLim), maxN2p);
            O2mTargetMin = min(max(minO2mLimI, minO2mLim), maxO2m);
            neRecoveryRatio = min(limRecoveryRatio,
                                  NeTargetMin / maxNe - (NeTargetMin / maxNe - ratioThr) / densityMulti);
            N2pRecoveryRatio = min(limRecoveryRatio,
                                  N2pTargetMin / maxN2p - (N2pTargetMin / maxN2p - ratioThr) / densityMulti);
            O2mRecoveryRatio = min(limRecoveryRatio,
                                  O2mTargetMin / maxO2m - (O2mTargetMin / maxO2m - ratioThr) / densityMulti);

            neRatioHThr = neRecoveryRatio * densityMultiHigh;
            N2pRatioHThr = N2pRecoveryRatio * densityMultiHigh;
            O2mRatioHThr = O2mRecoveryRatio * densityMultiHigh;

            // hysteresis for density correction recoveryRatio -> denRatioThr
            // This is an experiment to reduce negative densities
            denRatioThr = ( !factorCh_old
                            && neRecoveryRatio < (minNe / maxNe)
                            && N2pRecoveryRatio < (minN2p / maxN2p)
                            && O2mRecoveryRatio < (minO2m / maxO2m) ) ? ratioThr : denRatioThr;
            if ( ratioThr == denRatioThr )
            {
                neRatioThr = neRecoveryRatio * densityMulti;
                N2pRatioThr = N2pRecoveryRatio * densityMulti;
                O2mRatioThr = O2mRecoveryRatio * densityMulti;

                minRCyDec = ( (minNe / maxNe) < neRatioThr
                             || (minN2p / maxN2p) < N2pRatioThr
                             || (minO2m / maxO2m) < O2mRatioThr ) ? minRCyDec - 1 : 0;

                if (Pstream::master() && (neRatioThr < neRatioHThr
                                          || N2pRatioThr < N2pRatioHThr
                                          || O2mRatioThr < O2mRatioHThr)) Info << "ERROR: neRatioHThr: " << neRatioHThr
                                            << " neRatioThr: " << neRatioThr
                                            << " N2pRatioHThr: " << N2pRatioHThr
                                            << " N2pRatioThr: " << N2pRatioThr
                                            << " O2mRatioHThr: " << O2mRatioHThr
                                            << " O2mRatioThr: " << O2mRatioThr
                                            << nl;
            }
            else
            {
                minRCyDec = ( (minNe / maxNe) < neRecoveryRatio
                              || (minN2p / maxN2p) < N2pRecoveryRatio
                              || (minO2m / maxO2m) <  O2mRecoveryRatio ) ? minRCyDec - 1 : 0;

                neRatioThr = denRatioThr;
                N2pRatioThr = denRatioThr;
                O2mRatioThr = denRatioThr;
            }

            if (denRatioThrOld != denRatioThr && ratioThr == denRatioThr && Pstream::master())
                Info << runTime.timeIndex() << ": re-arm trigger min/max ne: " << minNe / maxNe
                     << " nN2p: " << minN2p / maxN2p
                     << " nO2m: " << minO2m / maxO2m << nl;

            // strong clamping
            bool clampnow = ( (minNe / maxNe) < neRatioHThr
                || (minN2p / maxN2p) < N2pRatioHThr
                || (minO2m / maxO2m) < O2mRatioHThr );
            if (clampnow)
            {
                denRatioThr = ratioThr;
                neRatioThr = neRecoveryRatio * densityMulti;
                N2pRatioThr = N2pRecoveryRatio * densityMulti;
                O2mRatioThr = O2mRecoveryRatio * densityMulti;
            }
        }

        if (maxDRhoEDtRateThr < maxDRhoEDtRate)
        {
            if ((maxDRhoEDtRateDec < 990) && 10 * maxDRhoEDtRateThr < maxDRhoEDtRate)
            {
                if (Pstream::master()) Info << "WARNING: correct for very high density changes" << nl;
                maxDRhoEDtRateDec = 0;

                intervalCount = 0;
                denRatioThr = ratioThr;
                neRatioThr = neRecoveryRatio * densityMulti;
                N2pRatioThr = N2pRecoveryRatio * densityMulti;
                O2mRatioThr = O2mRecoveryRatio * densityMulti;
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
        for (int potLoopCount = 0; potCorrection && potLoopCount < 2; potLoopCount++)
        {
            potCorrection = false;

            if (!(intervalCount % maxInterval))
            {
                // Experimental
                // suppress un-physical negative densities
                // balance the minimum charge: nN2p - nO2m - ne = 0
                if (enableDetailedLogs && 0 < maxNe && Pstream::master()) Info << "before clamp min/max ne: " << minNe << " " << maxNe
                                                                               << " min/max: " << minNe / maxNe
                                                                               << " neRecRatio: " << neRecoveryRatio << nl;
                if ((minNe / maxNe) < neRatioHThr)
                {
                    dimensionedScalar minNeLimD("minNeLimD", ne.dimensions(),
                                                (minNe < 0? max(1e-10, min(maxNe * 0.1, minNeLim)) : 0));
                    ne = 0.5 * (ne + sqrt(sqr(ne) + sqr(minNeLimD)));
                }

                if (enableDetailedLogs && 0 < maxN2p && Pstream::master()) Info << "before clamp min/max nN2p: " << minN2p << " " << maxN2p
                                                                               << " min/max: " << minN2p / maxN2p
                                                                               << " N2pRecRatio: " << N2pRecoveryRatio << nl;
                if ((minN2p / maxN2p) < N2pRatioHThr)
                {
                    dimensionedScalar minN2pLimD("minN2pLimD", ne.dimensions(),
                                                 (minN2p < 0? max(1e-10, min(maxN2p * 0.1, minN2pLim)) : 0));
                    nN2p = 0.5 * (nN2p + sqrt(sqr(nN2p) + sqr(minN2pLimD)));
                }

                if (enableDetailedLogs && 0 < maxO2m && Pstream::master()) Info << "before clamp min/max nO2m: " << minO2m << " " << maxO2m
                                                                               << " min/max: " << minO2m / maxO2m
                                                                               << " O2mRecRatio: " << O2mRecoveryRatio << nl;
                if ((minO2m / maxO2m) < O2mRatioHThr)
                {
                    dimensionedScalar minO2mLimD("minO2mLimD", ne.dimensions(),
                                                 (minO2m < 0? max(1e-10, min(maxO2m * 0.1, minO2mLim)) : 0));
                    nO2m = 0.5 * (nO2m + sqrt(sqr(nO2m) + sqr(minO2mLimD)));
                }

                ne.correctBoundaryConditions();
                nN2p.correctBoundaryConditions();
                nO2m.correctBoundaryConditions();

                // implement hysteresis using an arbitarily large negative number
                denRatioThr = -VGREAT;
                neRatioThr = denRatioThr;
                N2pRatioThr = denRatioThr;
                O2mRatioThr = denRatioThr;

                if (0 < maxNe && 0 < maxN2p && 0 < maxO2m)
                {
                    // restart density hysteresis count
                    minRCyDec = ( (minNe / maxNe) < neRecoveryRatio
                                  || (minN2p / maxN2p) < N2pRecoveryRatio
                                  || (minO2m / maxO2m) < O2mRecoveryRatio ) ? 3000 : 0;
                }

                // restart rapid change count
                maxDRhoEDtRateDec = 1000;

                // densities likely to have been updated so update rhoE
                rhoE = eCharge * (nN2p - nO2m - ne);
            }

            if (!(intervalCount % maxInterval))
            {
                intervalCount = 1;
                factorChange = true;

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
                factorCh_old = factorChange;
                factorChange = false;

                if (factorCh_old && ecoUpperThresh < maxEeCo && eeCoRateLimit < eeCoRate)
                {
                    // previous change has had little effect
                    // reaction rate is changing rapidly still
                    if (Pstream::master()) Info << "WARNING: eeCoRate is high" << nl;
                    enableDetailedLogs = true;
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
                    enableDetailedLogs = true;
                    pwr -= 1;
                    factor = factMulti * pow(10.0, pwr);
                    factorChange = true;

                    if (intervalCount % maxInterval) intervalCount++;
                }
                else if ( minRCyDec && minRCyDec < 2000 && ( (minNe / maxNe) < neRatioHThr
                                                || (minN2p / maxN2p) < N2pRatioHThr
                                                || (minO2m / maxO2m) < O2mRatioHThr ) )
                {
                    factMulti -= 1.0;
                    if (factMulti < 1.0)
                    {
                        pwr -= 1;
                        factMulti = 9.0;
                    }

                    factor = factMulti * pow(10.0, pwr);
                    factorChange = true;

                    if ((intervalCount % maxInterval) && (minNe / maxNe) < neRatioHThr) intervalCount++;
                    if ((intervalCount % maxInterval) && (minN2p / maxN2p) < N2pRatioHThr) intervalCount++;
                    if ((intervalCount % maxInterval) && (minO2m / maxO2m) < O2mRatioHThr) intervalCount++;

                    scalar denRatioThrOld = denRatioThr;
                    denRatioThr = ratioThr;
                    neRatioThr = neRecoveryRatio * densityMulti;
                    N2pRatioThr = N2pRecoveryRatio * densityMulti;
                    O2mRatioThr = O2mRecoveryRatio * densityMulti;
                    if (denRatioThrOld != denRatioThr && Pstream::master())
                        Info << runTime.timeIndex() << ": force re-arm trigger min/max ne: " << minNe / maxNe
                             << " nN2p: " << minN2p / maxN2p
                             << " nO2m: " << minO2m / maxO2m << nl;
                }
                else if ( minRCyDec < 0 )
                {
                    if (Pstream::master()) Info << "WARNING: correct for negative density" << nl;
                    enableDetailedLogs = true;
                    factMulti -= 1.0;
                    if (factMulti < 1.0)
                    {
                        pwr -= 1;
                        factMulti = 9.0;
                    }

                    factor = factMulti * pow(10.0, pwr);
                    factorChange = true;

                    if (intervalCount % maxInterval) intervalCount++;
                    if (intervalCount % maxInterval) intervalCount++;
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
                else if ( (maxDRhoEDtRateDec < 990)
                          && maxDRhoEDtRateThr < maxDRhoEDtRate
                          && ( (minNe / maxNe) < neRatioHThr
                               || (minN2p / maxN2p) < N2pRatioHThr
                               || (minO2m / maxO2m) < O2mRatioHThr ) )
                {
                    if (Pstream::master()) Info << "WARNING: correct for multiple threshold violation!" << nl;
                    enableDetailedLogs = true;
                    pwr -= 1;
                    factor = factMulti * pow(10.0, pwr);
                    factorChange = true;

                    intervalCount = 0;
                }
                else if (maxDRhoEDtRateDec < 0)
                {
                    if (Pstream::master()) Info << "WARNING: correct for rapid density changes" << nl;
                    enableDetailedLogs = true;

                    intervalCount = 0;
                    denRatioThr = ratioThr;
                    neRatioThr = neRecoveryRatio * densityMulti;
                    N2pRatioThr = N2pRecoveryRatio * densityMulti;
                    O2mRatioThr = O2mRecoveryRatio * densityMulti;
                    potCorrection = true;
                }
                else if (ecoLowerThresh < maxEeCo_old && eeCoRateLimit < eeCoRate)
                {
                    if (Pstream::master()) Info << "WARNING: eeCoRate is high" << nl;
                    enableDetailedLogs = true;
                    pwr -= 2;
                    factor = factMulti * pow(10.0, pwr);
                    factorChange = true;

                    intervalCount = 0;
                }
                else if (!(intervalCount+1 % maxInterval)
                         && maxDENTdDtRateThr < maxDENTdDtRate
                         && ecoLowerThresh < maxEeCo)
                {
                    if (Pstream::master()) Info << "WARNING: maxDENTdDtRate is high" << nl;
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
                else if (!factorCh_old
                         && ( neRatioThr < (minNe / maxNe)
                         && N2pRatioThr < (minN2p / maxN2p)
                         && O2mRatioThr < (minO2m / maxO2m) )
                         && (maxDRhoEDtRate < maxDRhoEDtRateThr)
                         && (maxDENTdDtRate < maxDENTdDtRateThr)
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

                    if ( neRecoveryRatio < (minNe / maxNe)
                         && N2pRecoveryRatio < (minN2p / maxN2p)
                         && O2mRecoveryRatio < (minO2m / maxO2m) ) intervalCount = maxInterval - 1;
                }
                else if (intervalCount < maxInterval - 1) intervalCount++;

                if (factorChange)
                {
                    // for debug enableDetailedLogs = true;
                    if (ratioThr == denRatioThr && 0 < maxNe && 0 < maxN2p && 0 < maxO2m
                        && ( (minNe / maxNe) < neRatioThr
                             || (minN2p / maxN2p) < N2pRatioThr
                             || (minO2m / maxO2m) < O2mRatioThr ) )
                    {
                        // Experimental
                        intervalCount = 0;
                    }

                    if (!(intervalCount % maxInterval))
                    {
                        enableDetailedLogs = true;

                        denRatioThr = ratioThr;
                        neRatioThr = neRecoveryRatio * densityMulti;
                        N2pRatioThr = N2pRecoveryRatio * densityMulti;
                        O2mRatioThr = O2mRecoveryRatio * densityMulti;
                        potCorrection = true;
                    }
                }

                if (1 && enableDetailedLogs && (ecoHyperThresh < maxEeCo) && Pstream::master())
                    Info << runTime.timeIndex() << ": THRESH extreamly high maxEeCo: " << maxEeCo << nl;
                else if (1 && enableDetailedLogs && (ecoUpperThresh < maxEeCo) && Pstream::master())
                    Info << runTime.timeIndex() << ": THRESH very high maxEeCo: " << maxEeCo << nl;
                else if (1 && enableDetailedLogs && (ecoLowerThresh < maxEeCo) && Pstream::master())
                    Info << runTime.timeIndex() << ": THRESH high maxEeCo: " << maxEeCo << nl;
                if (enableDetailedLogs && 0 < maxNe && (minNe / maxNe) < neRatioHThr && Pstream::master())
                    Info << runTime.timeIndex() << ": THRESH high min/max ne: " << minNe / maxNe << nl;
                if (enableDetailedLogs && 0 < maxN2p && (minN2p / maxN2p) < N2pRatioHThr && Pstream::master())
                    Info << runTime.timeIndex() << ": THRESH high min/max nN2p: " << minN2p / maxN2p << nl;
                if (enableDetailedLogs && 0 < maxO2m && (minO2m / maxO2m) < O2mRatioHThr && Pstream::master())
                    Info << runTime.timeIndex() << ": THRESH high min/max nO2m: " << minO2m / maxO2m << nl;
                if (enableDetailedLogs && (maxDRhoEDtRateThr < maxDRhoEDtRate) && Pstream::master())
                    Info << runTime.timeIndex() << ": THRESH rapid density changes maxDRhoEDtRate: " << maxDRhoEDtRate << nl;
                if (enableDetailedLogs && (maxDENTdDtRateThr < maxDENTdDtRate) && Pstream::master())
                    Info << runTime.timeIndex() << ": THRESH rapid E/N changes maxDENTdDtRate: " << maxDENTdDtRate << nl;

                if ((enableDetailedLogs || factorChange) && Pstream::master())
                    Info << runTime.timeIndex() << ": maxEeCo: " << maxEeCo << " eeCoRate: " << eeCoRate
                         << " factMulti: " << factMulti << " pwr: " << pwr
                         << " intervalCount: " << intervalCount
                         << " mnRCyDec: " << minRCyDec
                         << " mxDRDtRate: " << maxDRhoEDtRate
                         << " mxDRDtRateDec: " << maxDRhoEDtRateDec
                         << " mxDRDtRateThr: " << maxDRhoEDtRateThr
                         << " mxDENDtRate: " << maxDENTdDtRate
                         << " denRatioThr: " << denRatioThr
                         << nl;
                maxEeCo_old = maxEeCo;

                // Now update the DeltaT value
                scalar Vmax = gMax(mesh.V());
                scalar meshDeltaX = min(Foam::exp(Foam::log(Vmax)/3.0), GREAT);
                volScalarField magUe
                (
                    "magUe",
                    mag(mu_e * E)
                );
                scalar maxDriftVelocity = gMax(magUe);
                scalar dt1 = meshDeltaX / (50 * (maxDriftVelocity + SMALL));

                maxReactionRate = max(maxReactionRate, SMALL);
                scalar dt2 = factor / maxReactionRate;
                scalar deltaT = min(dt1, dt2);
                deltaT = min(deltaT, dtUpperLimit);
                if (enableDetailedLogs && Pstream::master())
                    Info << "Vmax: " << Vmax << " maxDriftVelocity: " << maxDriftVelocity
                         << " meshDeltaX: " << meshDeltaX
                         << " maxDriftVelocity: " << maxDriftVelocity
                         << " dt1: " << dt1
                         << " dt2: " << dt2
                         << nl;
                if (runTime.deltaTValue() < 1e-40) enableDetailedLogs = true;
                if (enableDetailedLogs && (runTime.deltaTValue() < 1e-40) && Pstream::master())
                    Info << "WARNING: small deltaT: " << deltaT << nl;

                runTime.setDeltaT(deltaT);
            }
        }

        if (enableDetailedLogs || factorChange)
        {
            #include "mpCourantNo.H"
            // #include "setDeltaT.H" : instead opt for a dynamic DeltaT based on the reaction rate
        }

        if (enableDetailedLogs && Pstream::master()) Info << "Iteration: " << runTime.name()
                                                          << " index: " << runTime.timeIndex()
                                                          << " deltaT: " << runTime.deltaTValue()
                                                          << nl << nl;
        while (pimple.loop())
        {
            for (int corr=0; corr < nPhiECorrectors; corr++)
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
                    const fvPatchVectorField& nEpatch = E.boundaryField()[nPatchID];
                    const tmp<vectorField> tnHat = mesh.boundary()[nPatchID].nf();
                    const vectorField& nHat = tnHat();
                    scalarField nEn = nEpatch & nHat;
                    scalar minNEn = gMin(nEn);
                    scalar maxNEn = gMax(nEn);
                    if (Pstream::master()) Info << "min/max En at NELEMENT boundary: " << minNEn << " " << maxNEn << nl;

                    const fvPatchVectorField& pEpatch = E.boundaryField()[pPatchID];
                    const tmp<vectorField> tpHat = mesh.boundary()[pPatchID].nf();
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
                    if (Pstream::master()) Info << "phiE min/max: " << minPhiE << " " << maxPhiE << nl;
                }
            }

            // update U
            #include "momentumEqns.H"

            // update rAU in preparation for solving for pressure
            #include "pressureEqns.H"

            while (pimple.correct())
            {
                constrainPressure(p, U, phiHbyA, rAU);

                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);

                while (pimple.correctNonOrthogonal())
                {
                    SolverPerformance<scalar> pPerf = pEqn.solve();
                    // if (enableDetailedLogs && Pstream::master()) Info << "Residual p = " << pPerf.finalResidual() << nl;

                    if (pimple.finalNonOrthogonalIter())
                    {
                        // calculate phi
                        phi = phiHbyA - pEqn.flux();
                    }
                }

                // recompute HbyA for better stability
                // HbyA = rAU * UEqn.H();
                // HbyA = constrainHbyA(rAU * UEqn.H(), U, p);
                // adjustPhi(phiHbyA, U, p);

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
                    if (Pstream::master()) Info << "Pressure corrector: max(U): " << maxMagU << nl;
                }
            }

            #include "thermalEqns.H"

            if (enableDetailedLogs)
            {
                dimensionedVector thrust = fvc::domainIntegrate(FEHD);

                surfaceVectorField phiU(phi * fvc::interpolate(U));
                vector num = gSum(phiU.primitiveField());
                scalar den = gSum(phi.primitiveField());
                vector Uavg = num / (den + SMALL);

                dimensionedScalar power = fvc::domainIntegrate(FEHD & U);

                // energy based velocity
                // ie velocity required that would produce the same power
                scalar v_eff = 0 < mag(thrust.value()) ? power.value() / mag(thrust.value()) : 0;
                if (Pstream::master()) Info << "thrust: " << thrust << nl
                                            << " Uavg: " << Uavg << nl
                                            << " power: " << power << " v_eff: " << v_eff << nl;
            }

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

        if (enableDetailedLogs && Pstream::master()) Info<< "ExecutionTime: " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime: " << runTime.elapsedClockTime() << " s"
            << nl << nl;
    }

    Info<< nl << "End" << nl << nl;
    return 0;
}


// ************************************************************************* //
