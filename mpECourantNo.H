{
    // now calculate the Courant Numbers for e and ions
    surfaceScalarField phiEFlux_e = fvc::flux(-mu_e * E);
    surfaceScalarField phiEFlux_N2p = fvc::flux(mu_N2p * E);

    maxEeCo = 0.0;
    scalar maxEpCo = 0.0;
    const scalar dt = runTime.deltaTValue();

    forAll(phiEFlux_e, faceI)
    {
        scalar flux = mag(phiEFlux_e[faceI]);

        label own = mesh.owner()[faceI];
        label nei = mesh.neighbour()[faceI];

        scalar coOwn = flux * dt / mesh.V()[own];
        scalar coNei = flux * dt / mesh.V()[nei];

        maxEeCo = max(maxEeCo, max(coOwn, coNei));
    }

    maxEeCo = returnReduce(maxEeCo, maxOp<scalar>());

    forAll(phiEFlux_N2p, faceI)
    {
        scalar flux = mag(phiEFlux_N2p[faceI]);

        label own = mesh.owner()[faceI];
        label nei = mesh.neighbour()[faceI];

        scalar coOwn = flux * dt / mesh.V()[own];
        scalar coNei = flux * dt / mesh.V()[nei];

        maxEpCo = max(maxEpCo, max(coOwn, coNei));
    }

    maxEpCo = returnReduce(maxEpCo, maxOp<scalar>());

    scalar totalFlux =
        gSum(mag(phiEFlux_e)())
        + gSum(mag(phiEFlux_N2p)());

    scalar meanECo = (totalFlux * dt) / gSum(mesh.V());

    if (enableDetailedLogs)
    {
        volVectorField Ue
        (
            "Ue",
            -mu_e * E
        );

        volScalarField magUe
        (
            "magUe",
             mag(Ue)
        );

        scalar minMagUe = gMin(magUe);
        scalar maxMagUe = gMax(magUe);
        if (Pstream::master()) Info << "Ue min/max: " << minMagUe << " / " << maxMagUe << nl;

        scalar minPhiEFlux_e = gMin(phiEFlux_e);
        scalar maxPhiEFlux_e = gMax(phiEFlux_e);
        if (Pstream::master()) Info << "phiEFlux_e min/max: "
                                    << minPhiEFlux_e << " / " << maxPhiEFlux_e << nl;

        if (Pstream::master()) Info << "Mean drift Courant Number = " << meanECo << nl;
        if (Pstream::master()) Info << "Max drift Courant Numbers: electron: " << maxEeCo
                                    << " positive ion: " << maxEpCo << nl;
    }
}
