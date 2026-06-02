
volScalarField rAU
(
    "rAU",
    1.0 / ( UEqn.A().ref() + dimensionedScalar("deltaA", UEqn.A().ref().dimensions(), 1e-12 ) )
);

if (enableDetailedLogs)
{
    scalar minRAU = gMin(rAU);
    scalar maxRAU = gMax(rAU);
    if (Pstream::master()) Info << "rAU min/max: " << minRAU << " / " << maxRAU << endl;

    scalar minA = gMin(UEqn.A().ref());
    scalar maxA = gMax(UEqn.A().ref());
    if (Pstream::master()) Info << "UEqn.A() min/max: " << minA << " / " << maxA << endl;

    vector minH = gMin(UEqn.H().ref());
    vector maxH = gMax(UEqn.H().ref());
    if (Pstream::master()) Info << "UEqn.H() min/max: " << minH << " / " << maxH << endl;
}

rAU.correctBoundaryConditions();

volVectorField HbyA = constrainHbyA(rAU*UEqn.H(), U, p);
surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA)
);

adjustPhi(phiHbyA, U, p);
