// Header for globally accessible routines

#ifndef beamDynInterface_H
#define beamDynInterface_H

#include "vectorList.H"
#include "turbulenceModel.H"
#include "dynamicFvMesh.H"

namespace BD
{
#ifndef namespaceFoam
#define namespaceFoam
    using namespace Foam;
#endif

    // access functions

    //-used by solver
    label interfacePatchID();

    //-used by boundary patch
    vectorList& linDisp();
    vectorList& angDisp();
    double* h();
    label N();
    Switch enforce2D();
    label bladeDirection();

    // member functions

    void readInputs( const dynamicFvMesh& mesh, 
                     const IOdictionary& dict,
                     Switch fluidSolve,
                     Switch beamSolve );

    void start( double t0, double dt );
    void calculateShapeFunctions( const pointField& pf );

    void updateSectionLoads( const double* F,  // set constant distributed load, for testing
                             const double* M );
    void updateSectionLoads( const dynamicFvMesh& mesh, 
                             const volScalarField& p, 
                             const incompressible::turbulenceModel& turbulence );

    void update( double t, double dt );
    void updatePrescribedDeflection( double t );

    void updateNodePositions();

    void write( bool writeNow, std::string timeName );

    void stop();

}

#endif
