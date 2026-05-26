
/*dimensionedScalar lambda
(
    "lambda",
    dimless / sqr(dimLength),
    1e-10
);
*/

fvScalarMatrix phiEqn
(
    fvm::laplacian(phiE)
  // - fvm::Sp(lambda, phiE)  // has no effect! introduced to reduce the growth of phiE
    ==
  - rhoE / epsilon0
    //- rhoEff / epsilon0 introduced to reduce growth of E and phiE but has no effect!
);

// Solve Poisson
phiEqn.relax();
SolverPerformance<scalar> phiEPerf = phiEqn.solve();

// Print Poisson residual
Info << "Residual phiE = " << phiEPerf.finalResidual() << endl;

Info << "min/max phiE: " << min(phiE).value() << " " << max(phiE).value() << " phiE range = " << mag(gMax(phiE) - gMin(phiE)) << endl;
