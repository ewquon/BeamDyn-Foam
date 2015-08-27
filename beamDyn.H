// Header file that declares the BD namespace and defines the member access functions
// This file includes beamDynInterface.H which declares the globally accessible functions

#ifndef beamDyn_H
#define beamDyn_H

#include "className.H"
#include "constants.H"

// need to include these headers so that compilation of 
// beamDynInterfacePointPatchField doesn't crap out!
//#include "fvMesh.H"
//#include "fvMatrices.H"
//#include "vectorList.H"
//#include "turbulenceModel.H"
//#include "dynamicFvMesh.H"

#include <fstream>

#include "beamDynInterface.H" // same as what's included by the other codes!

namespace BD
{
#ifndef namespaceFoam
#define namespaceFoam
    using namespace Foam;
#endif

    NamespaceName("BeamDyn");

    const scalar pi(Foam::constant::mathematical::pi);
    const scalar eps(1.0e-8);

    // inputs from couplingProperties

    word    patchName;
    label   patchID;

    label   bladeDir;
    scalar  bladeR0;
    scalar  bladeR;
    vector  origin;

    Switch  twoD;

    scalar  loadMultiplier; // scaling applied to aero loads, for development
 
    vector  prescribed_max_deflection; // only works for y deflections for now
    vector  prescribed_max_rotation;   // only works for x rotations for now, overwrites bending rotation if specified

    // coupling variables

    scalar  rhoRef;
    scalar  pRef;

    bool first=true; // first iteration?

    std::fstream loadFile;
    std::fstream dispFile;

    int nnodes=0;
    
    double  currentTime = -1;
    double  currentDeltaT = -1;

    vectorList  *pos0_ptr,*rot0_ptr;    // starting position and orientation of beam nodes
    vectorList  *pos_ptr, *rot_ptr;     // current position and orientation of beam nodes
    vectorList  *disp_ptr;              // linear displacement at beam nodes
    vectorList  *adisp_ptr;             // angular displacement at beam nodes
    scalarList  *r_ptr;                 // spanwise coordinates of beam nodes
    double      *h_ptr;                 // shape functions [ # CFD nodes * # BD nodes ]

    // access functions

    //-used by solver
    label interfacePatchID() { return patchID; };

    //-used by boundary patch
    vectorList& linDisp() { return *disp_ptr; };  // linear displacement
    vectorList& angDisp() { return *adisp_ptr; }; // angular displacement
    //scalarList& r() { return *r_ptr; };         // beam coordinate, NOT USED?
    double* h() { return h_ptr; };                // shape functions
    label N() { return nnodes; };                 // number of beamdyn nodes
    Switch enforce2D() { return twoD; };          // force 2D displacements
    label bladeDirection() { return bladeDir; };

}

// only need this is we're not compiling this source as part of a library(?)
//#ifdef NoRepository
//#   include "beamDyn.C"
//#endif

#endif
