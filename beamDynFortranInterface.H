#ifndef BEAMDYNINTERFACE_H
#define BEAMDYNINTERFACE_H

// These Fortran functions are defined in BeamDyn/Library_BeamDyn_AM2.f90
//
// NOTES:
// - Fortran uses call-by-address convention; only pointers should be passed to Fortran.
// - quadruple precision variables (long double) are not exchanged properly

extern "C" 
{
    //
    // solution control
    //
    void beamDynStart( double *t0_foam, double *dt_foam );
    void beamDynEnd();
    void beamDynStep( double *dt_foam );

    void beamDynWriteState( const char *fname );
    void beamDynReadState( const char *fname );

    //
    // set functions
    //
    //void beamDynSetDistributedLoadAtNode( int *inode0, double *F, double *M );
    void beamDynSetDistributedLoad( int *ig, double *F, double *M );
    void beamDynSetDisplacement( int *inode0, double *u );

    void beamDynSetXRotationMatrix( int *inode0, double *ang );    // set the rotation matrix for rotation about x-axis; angle is in radians
    void beamDynSetZRotationMatrix( int *inode0, double *ang );    // set the rotation matrix for rotation about z-axis; angle is in radians

    //
    // access functions
    //
    void beamDynGetNnodes( int *nnodes );
    void beamDynGetDispAtNode( int *ispan, double *lin_disp, double *ang_disp );        // deprecated
    void beamDynGetNode0Displacement( int *ispan, double *lin_disp, double *ang_disp ); // get linear/angular DISPLACEMENT at node (0-indexed)
    void beamDynGetNode0Position( int *ispan, double *pos, double *rot );               // get position/rotation at node (0-indexed)
    void beamDynGetInitNode0Position( int *ispan, double *pos0, double *rot0 );         // get initial position/rotation at node (0-indexed)
    void beamDynGetNode0RotationMatrix( int *ispan, double *R );                        // get rotation matrix[9], stored row-wise

    void beamDynGetShapeFunctions( double *s, double *hx );


    // deprecated:
    void beamDynGetDisp( int *N, double *s, double *u, double *v, double *w, double *kx, double *ky, double *kz );

    // shouldn't be needed:
    void beamDynGetNgp( int *ngp );
    void beamDynGetGLLPts( double *gllp, double *gllw );
    void beamDynGetGaussPts( double *gp, double *gw );

}

#endif
