{
    // now calculate the Courant Numbers for e and ions
    surfaceScalarField phiEFlux_e = fvc::flux(-mu_e * E);
    surfaceScalarField phiEFlux_N2p = fvc::flux(mu_N2p * E);
    surfaceScalarField phiEFlux_O4p = fvc::flux(mu_O4p * E);
    surfaceScalarField phiEFlux_O2p = fvc::flux(mu_O2p * E);
    surfaceScalarField phiEFlux_O2m = fvc::flux(-mu_O2m * E);

    maxEeCo = 0.0;
    scalar maxEN2pCo = 0.0;
    scalar maxEO4pCo = 0.0;
    scalar maxEO2pCo = 0.0;
    scalar maxEpCo = 0.0;
    scalar maxEmCo = 0.0;
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

        maxEN2pCo = max(maxEN2pCo, max(coOwn, coNei));
    }

    maxEN2pCo = returnReduce(maxEN2pCo, maxOp<scalar>());

    forAll(phiEFlux_O4p, faceI)
    {
        scalar flux = mag(phiEFlux_O4p[faceI]);

        label own = mesh.owner()[faceI];
        label nei = mesh.neighbour()[faceI];

        scalar coOwn = flux * dt / mesh.V()[own];
        scalar coNei = flux * dt / mesh.V()[nei];

        maxEO4pCo = max(maxEO4pCo, max(coOwn, coNei));
    }

    maxEO4pCo = returnReduce(maxEO4pCo, maxOp<scalar>());

    forAll(phiEFlux_O2p, faceI)
    {
        scalar flux = mag(phiEFlux_O2p[faceI]);

        label own = mesh.owner()[faceI];
        label nei = mesh.neighbour()[faceI];

        scalar coOwn = flux * dt / mesh.V()[own];
        scalar coNei = flux * dt / mesh.V()[nei];

        maxEO2pCo = max(maxEO2pCo, max(coOwn, coNei));
    }

    maxEO2pCo = returnReduce(maxEO2pCo, maxOp<scalar>());

    maxEpCo = max(maxEN2pCo, max(maxEO4pCo, maxEO2pCo));

    forAll(phiEFlux_O2m, faceI)
    {
        scalar flux = mag(phiEFlux_O2m[faceI]);

        label own = mesh.owner()[faceI];
        label nei = mesh.neighbour()[faceI];

        scalar coOwn = flux * dt / mesh.V()[own];
        scalar coNei = flux * dt / mesh.V()[nei];

        maxEmCo = max(maxEmCo, max(coOwn, coNei));
    }

    maxEmCo = returnReduce(maxEmCo, maxOp<scalar>());

    scalar totalFlux =
        gSum(mag(phiEFlux_e)())
        + gSum(mag(phiEFlux_N2p)())
        + gSum(mag(phiEFlux_O2p)())
        + gSum(mag(phiEFlux_O4p)())
        + gSum(mag(phiEFlux_O2m)());

    scalar meanECo = (totalFlux * dt) / gSum(mesh.V());

    if (enableDetailedLogs)
    {
        tmp<volVectorField> Ue = -mu_e * E;
        tmp<volScalarField> magUe = mag(Ue());

        scalar minMagUe = gMin(magUe());
        scalar maxMagUe = gMax(magUe());
        if (Pstream::master()) Info << "Ue min/max: " << minMagUe << " / " << maxMagUe << nl;

        scalar minPhiEFlux_e = gMin(phiEFlux_e);
        scalar maxPhiEFlux_e = gMax(phiEFlux_e);
        if (Pstream::master()) Info << "phiEFlux_e min/max: "
                                    << minPhiEFlux_e << " / " << maxPhiEFlux_e << nl;

        if (Pstream::master()) Info << "Mean drift Courant Number = " << meanECo << nl;
        if (Pstream::master()) Info << "Max drift Courant Numbers: electron: " << maxEeCo
                                    << " positive ion: " << maxEpCo
                                    << " negative ion: " << maxEmCo << nl;
    }
}
