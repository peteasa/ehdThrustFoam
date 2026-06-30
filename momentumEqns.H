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
