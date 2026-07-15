
{
    if (Pstream::master()) Info << "create initial charge density at " << runTime.name() << nl;
    // alternative for this is to use createZones and setFields
    // however setFields only expects geometry defined zones so is not
    // as flexible
    labelHashSet patchIDs;
    patchIDs.insert(pPatchID);

    scalarField cellDist(mesh.nCells(), GREAT);

    Foam::patchDistWave::calculate
    (
        mesh,
        patchIDs,
        cellDist
    );

    volScalarField distanceToAnode
    (
        IOobject
        (
            "distanceToAnode",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimLength, 0.0)
    );

    // copy internal field
    forAll(cellDist, cellI)
    {
        distanceToAnode[cellI] = cellDist[cellI];
    }

    distanceToAnode.correctBoundaryConditions();

    forAll(mesh.C(), cellI)
    {
        if (distanceToAnode[cellI] < sheathThickness)
        {
            ne[cellI]    = initNe;
            nN2p[cellI]  = initNN2p;
        }
    }
}
