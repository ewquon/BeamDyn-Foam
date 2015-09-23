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
#include "IFstream.H"
#include "OFstream.H"

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
    const scalar radToDeg(180.0/pi);

    // inputs from couplingProperties

    word    patchName;
    label   patchID;

    labelList  OFtoIEC(3); // simple coordinate mappings
    labelList  IECtoOF(3);

    label   bladeDir;
    scalar  bladeR0;
    scalar  bladeR;
    vector  origin;

    Switch  twoD;

    // debug inputs
//    scalar  loadMultiplier; // scaling applied to aero loads, for development
// 
//    vector  prescribed_max_deflection; // only works for y deflections for now
//    vector  prescribed_max_rotation;   // only works for x rotations for now, overwrites bending rotation if specified

    // coupling variables

    scalar  rhoRef;
    scalar  pRef;

    bool first=true; // first iteration?
    bool restarted=false;

    std::fstream loadFile;
    //std::fstream dispFile;
    std::fstream posFile;
    std::fstream trackFile;

    int nnodes=0;
    
    double  currentTime = -1;
    double  currentDeltaT = -1;

    // calculated/used once at t=0
    double      *h_ptr;                 // shape functions [ # CFD nodes * # BD nodes ]
    vectorList  *p_ptr;                 // surface offset vectors

    vectorList  *pos0_ptr,*crv0_ptr;    // starting position and orientation of beam nodes (IEC coords)
                                        // - pos0_ptr used to evaluate r to calculate shape functions
                                        // - crv0_ptr is NOT USED

    labelList   *trackedPts_ptr;        // tracked points from pointSet "trackedPoints" (generated with topoSet)

    // updated from BeamDyn in updateNodePositions() after calling beamDynStep()
    vectorList  *pos_ptr, *crv_ptr;     // current position and orientation of beam nodes (IEC coords)
    vectorList  *disp_ptr;              // linear displacement at beam nodes (IEC coords)
    vectorList  *adisp_ptr;             // angular displacement at beam nodes (IEC coords)
    
    // calculated at all t
//    vectorList  *x1_ptr;                // position vector along the beam (x1) axis
//                                        // - initial value calculated by calculateInitialDisplacementVectors()
//                                        // - updated every time step by beamDynInterfacePointPatchField::updateCoeffs()
    scalarList  *r_ptr;                 // spanwise coordinates of beam nodes used to calculate shape functions
                                        // - calculated in updateNodePositions()


    //
    // access functions
    //

    //-used by solver
    label interfacePatchID() { return patchID; };

    //-used by boundary patch
    vectorList& pos() { return *pos_ptr; };  // linear displacement
    vectorList& crv() { return *crv_ptr; }; // angular displacement
    vectorList& crv0() { return *crv0_ptr; }; // angular displacement
    vectorList& linDisp() { return *disp_ptr; };  // linear displacement
    vectorList& angDisp() { return *adisp_ptr; }; // angular displacement
    //scalarList& r() { return *r_ptr; };           // beam coordinate at beamdyn nodes --NOT USED
    double* h() { return h_ptr; };                // shape functions
    label N() { return nnodes; };                 // number of beamdyn nodes
    Switch enforce2D() { return twoD; };          // force 2D displacements
    label bladeDirection() { return bladeDir; };  // in OpenFOAM coordinates
    //label openfoamDir(label dir) { return IECtoOF[dir]; };
    vectorList& p() { return *p_ptr; };           // surface node offset
    //vectorList& x1() { return *x1_ptr; };         // beam axis position
    labelList& trackedPoints() { return *trackedPts_ptr; };

    //
    // utility functions
    //
    inline double crvToRad(double crv) { return 4.0 * Foam::atan( crv/4.0 ); };
    std::string createFname( std::string fname );
}

// only need this is we're not compiling this source as part of a library(?)
//#ifdef NoRepository
//#   include "beamDyn.C"
//#endif

#endif
