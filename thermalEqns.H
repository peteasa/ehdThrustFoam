volScalarField alphaEff
(
    "alphaEff",
    momentumTransport->nu() / Pr + momentumTransport->nut() / Prt
);
//if (Pstream::master()) Info << "alphaEff dimensions: " << alphaEff.dimensions() << nl;

Cp = Cp0 + CpSlope * (T - CpRef);

//if (Pstream::master()) Info << "Cp dimensions: " << Cp.dimensions() << nl;

fvScalarMatrix TEqn
(
    fvm::ddt(T)
  + fvm::div(phi, T)
  - fvm::laplacian(alphaEff, T)
 ==
    (J & E) / (rho * Cp)
);

SolverPerformance<scalar> TPerf = TEqn.solve();

if (enableDetailedLogs && Pstream::master()) Info << "Residual T = " << TPerf.finalResidual() << nl;

