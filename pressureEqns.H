
volScalarField rAU
(
    "rAU",
    1.0 / ( UEqn.A().ref() + dimensionedScalar("deltaA", UEqn.A().ref().dimensions(), 1e-12 ) )
);

Info << "rAU min/max: " << gMin(rAU) << " / " << gMax(rAU) << endl;
Info << "UEqn.A() min/max: " << gMin(UEqn.A().ref()) << " / " << gMax(UEqn.A().ref()) << endl;
Info << "UEqn.H() min/max: " << gMin(UEqn.H().ref()) << " / " << gMax(UEqn.H().ref()) << endl;

rAU.correctBoundaryConditions();

volVectorField HbyA = constrainHbyA(rAU*UEqn.H(), U, p);
surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA)
);

adjustPhi(phiHbyA, U, p);
