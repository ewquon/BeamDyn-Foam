// this file should be included by solver

Info << "\nReading coupling properties" << endl;
IOdictionary couplingProperties
(
    IOobject
    (
        "couplingProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Switch fluidSolve;
Switch beamSolve;

BD::readInputs(mesh,
               couplingProperties,
               fluidSolve,
               beamSolve);

Info<< "  fluid/beam solvers: " 
    << fluidSolve << " " << beamSolve << endl;
