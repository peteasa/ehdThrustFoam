volVectorField FEHD = rhoE * E;
if (enableDetailedLogs)
{
    volScalarField magF
    (
         "magF",
         mag(FEHD)
    );

    scalar minmagF = gMin(magF);
    scalar maxmagF = gMax(magF);
    if (Pstream::master()) Info << "min/max FEHD: " << minmagF << " " << maxmagF << nl;
}

if (startupT)
{
    if (1 && (runTime.value() < startupT))
    {
        if (enableDetailedLogs && Pstream::master()) Info << "Slow Start: Separate electrostatic interaction until time: " << startupT << nl;
        FEHD = dimensionedVector("zv", FEHD.dimensions(), vector::zero);
        interactionStart = runTime.timeIndex();
    }
    else
    {
        int count = runTime.timeIndex() - interactionStart;
        if (count < 1/startupIncrement)
        {
            // gradually increase FEHD
            dimensionedScalar fScale("fScale", dimless, startupIncrement * count);
            if (enableDetailedLogs && Pstream::master()) Info << "fScale: " << startupIncrement * count << nl;
            FEHD *= fScale;
        }
    }
}

fvVectorMatrix UEqn
(
      fvm::ddt(U)
    + fvm::div(phi, U)
    - momentumTransport->divDevTau(U)
      ==
      FEHD / rho
);

UEqn.relax();
if (momentumPredictor)
{
    UEqn.solve();
}

//call correct() after the pressure loop has updated U
//momentumTransport->correct();
