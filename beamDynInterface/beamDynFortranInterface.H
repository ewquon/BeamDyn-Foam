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
    void beamDynStop();
    void beamDynStep( double *dt_foam );

    void beamDynWriteState( const char *fname );
    void beamDynReadState( const char *fname, int *istat );

    //
    // set functions
    //

    void beamDynSetDistributedLoad( int *ig, double *F, double *M ); // set distributed load at 0-index gauss point
    //void beamDynSetDisplacement( int *inode0, double *u );

    //void beamDynSetXRotationMatrix( int *inode0, double *ang );    // set the rotation matrix for rotation about x-axis; angle is in radians
    //void beamDynSetZRotationMatrix( int *inode0, double *ang );    // set the rotation matrix for rotation about z-axis; angle is in radians

    //
    // access functions
    //

    void beamDynGetNNodes( int *nnodes );
    void beamDynGetNode0InitialPosition( int *ispan, double *pos0, double *crv0 );      // get initial position/rotation at node (0-indexed)
    void beamDynGetNode0Position( int *ispan, double *pos, double *crv );               // get position/rotation at node (0-indexed)
    void beamDynGetNode0Displacement( int *ispan, double *lin_disp, double *ang_disp ); // get linear/angular DISPLACEMENT at node (0-indexed)

    void beamDynGetShapeFunctions( double *s, double *hx );

    // shouldn't be needed:
    //void beamDynGetDispAtNode( int *ispan, double *lin_disp, double *ang_disp );
    //void beamDynGetDisp( int *N, double *s, double *u, double *v, double *w, double *kx, double *ky, double *kz );
    //void beamDynGetNgp( int *ngp );
    //void beamDynGetGLLPts( double *gllp, double *gllw );
    //void beamDynGetGaussPts( double *gp, double *gw );
    //void beamDynGetNode0RotationMatrix( int *ispan, double *R );                        // get rotation matrix[9], stored row-wise

    void beamDynGetRotationMatrix( double *Rmat, double *crv );

}

#endif
