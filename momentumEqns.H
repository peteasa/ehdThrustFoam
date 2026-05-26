volVectorField FEHD = rhoE * E;
Info << "min/max FEHD: " << min(mag(FEHD)).value() << " " << max(mag(FEHD)).value() << endl;

if (0)
{
    if (1 && (runTime.value() < 3e-5))
    {
        Info << "DEBUG: Separate electrostatic interaction" << endl;
        FEHD = dimensionedVector("zv", FEHD.dimensions(), vector::zero);
    } else {
        dimensionedScalar fScale("fScale", dimless, 1e-14);
        if (0) FEHD *= fScale;
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
