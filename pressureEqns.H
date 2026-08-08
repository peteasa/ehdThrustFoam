
rAU = 1.0 / ( UEqn.A().ref() + dimensionedScalar("deltaA", UEqn.A().ref().dimensions(), 1e-12 ) );
rAU.correctBoundaryConditions();

volVectorField HbyA = constrainHbyA(rAU * UEqn.H(), U, p);
surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA)
);

if (enableDetailedLogs)
{
    scalar minRAU = gMin(rAU);
    scalar avgRAU = gSum(mesh.V() * rAU.primitiveField()) / gSum(mesh.V());
    scalar maxRAU = gMax(rAU);
    tmp<volScalarField> magHbyA = mag(HbyA);
    scalar minHbyA = gMin(magHbyA());
    scalar maxHbyA = gMax(magHbyA());
    scalar minPhiHbyA = gMin(phiHbyA);
    scalar maxPhiHbyA = gMax(phiHbyA);
    if (Pstream::master()) Info << "rAU min/mean/max: " << minRAU << " " << avgRAU << " " << maxRAU
                                << " HbyA min/max: " << minHbyA << " " << maxHbyA
                                << " phiHbyA min/max: " << minPhiHbyA << " " << maxPhiHbyA << nl;
}

adjustPhi(phiHbyA, U, p);

