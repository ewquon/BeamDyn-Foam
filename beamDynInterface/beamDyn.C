#include "beamDyn.H"
#include "beamDynFortranInterface.H" // define Fortran functions

#include "fvCFD.H"

#include <limits>

namespace BD
{
    ///////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Routines that directly interface with BeamDyn
    //
    ///////////////////////////////////////////////////////////////////////////////////////////////

    void start( double t0, double dt )
    {
        currentTime = t0;
        currentDeltaT = dt;
        restarted = 0;

        if (Pstream::master())
        {
            Info<< "\n================================" << endl;
            Info<<   "Starting BeamDyn" << endl;
            Info<<   "================================" << endl;
            //Info<< "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n" << endl;
            beamDynStart( &t0, &dt );
            beamDynGetNNodes( &nnodes ); // total number of nodes in beam model
            //Info<< "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
            Info<< "...done\n" << endl;
        }

        // initialize arrays for storing configuration
        Pstream::scatter(nnodes);
        r_ptr    = new scalarList(nnodes, 0.0);
        pos0_ptr = new vectorList(nnodes, vector::zero);
        crv0_ptr = new vectorList(nnodes, vector::zero);
        pos_ptr  = new vectorList(nnodes, vector::zero);
        crv_ptr  = new vectorList(nnodes, vector::zero);
        disp_ptr = new vectorList(nnodes, vector::zero);
        adisp_ptr= new vectorList(nnodes, vector::zero);

        // perform restart read of saved state data if necessary
        // open disp.out and load.out for writing later
        if(Pstream::master())
        {
            if (t0 > 0)
            {
                Info<< "Attempting restart from previous time = " << t0 << endl;

                //std::string rstFile("BeamDynState_" + Foam::Time::timeName(t0) + ".dat");
                std::string rstFile(Foam::Time::timeName(t0) + "/BeamDynState.dat");
                if (FILE *file = fopen(rstFile.c_str(), "r"))
                {
                    fclose(file);
                    beamDynReadState( rstFile.c_str() );
                    restarted = 1;
                }

                if( !restarted ) // try different restart file
                {
                    std::string rstFile("BeamDynState_last.dat");
                    if (FILE *file = fopen(rstFile.c_str(), "r"))
                    {
                        fclose(file);
                        beamDynReadState( rstFile.c_str() );
                        restarted = 1;
                    }
                }

                if( !restarted ) 
                {
                    Info<< "Problem opening restart file " << rstFile << endl;
                }

                //***these outputs are in IEC coordinates***
                // angles are in degrees
                loadFile.open("load.out", std::ios::in | std::ios::out | std::ios::app);
                dispFile.open("disp.out", std::ios::in | std::ios::out | std::ios::app);
            }
            else
            {
                //***these outputs are in IEC coordinates***
                // angles are in degrees
                loadFile.open("load.out", std::ios::out);
                dispFile.open("disp.out", std::ios::out);
            }
            if (!loadFile.is_open()) Info<< "Problem opening load.out???" << endl;
            if (!dispFile.is_open()) Info<< "Problem opening disp.out???" << endl;

            //Info<< "Setting precision to " << std::numeric_limits<double>::digits10 << endl;
            //loadFile.precision(std::numeric_limits<double>::digits10);
            //dispFile.precision(std::numeric_limits<double>::digits10);
            loadFile.precision(8);
            dispFile.precision(8);

            // get initial configuration (IN BEAMDYN COORDINATES)
            double posi[3], crvi[3];
            Info<< "Initial linear/angular position:" << endl;
            for( int inode=0; inode < nnodes; ++inode )
            {
                beamDynGetNode0InitialPosition( &inode, posi, crvi );
                for(int i=0; i < 3; ++i) {
                    (*pos0_ptr)[inode][i] = posi[i];
                    (*crv0_ptr)[inode][i] = crvi[i];
                }
                Info<< " " << posi[0] 
                    << " " << posi[1] 
                    << " " << posi[2]
                    << "  "
                    << " " << crvToRad(crvi[0])*radToDeg
                    << " " << crvToRad(crvi[1])*radToDeg
                    << " " << crvToRad(crvi[2])*radToDeg
                    << endl;
            }

        } //if Pstream master

        updateNodePositions(); // this should write out either the initial configuration (0's)
                               // or the restart configuration

//        // save initial configuration, used by updateNodePositions() in subsequent iterations
//        for( int inode=0; inode < nnodes; ++inode )
//        {
//            (*pos0_ptr)[inode] = (*pos_ptr)[inode];
//            (*crv0_ptr)[inode] = (*crv_ptr)[inode];
//        }

        Info<< "BeamDyn initialization complete.\n\n";

    }

    //*********************************************************************************************

    void readInputs( const dynamicFvMesh& mesh, 
                     const IOdictionary& couplingProperties,
                     Switch& fluidSolve,
                     Switch& beamSolve )
    {
        // Read interface patch info
        // TODO: only one interface patch for now

        beamSolve = 1;
        couplingProperties.lookup("interfacePatch") >> patchName;
        patchID = mesh.boundaryMesh().findPatchID(patchName);
        Info<< "  Found interfacePatchID for " << patchName 
            << " : " << patchID << endl;

        if (patchID < 0)
        {
            //FatalErrorIn(args.executable())
            //    << "Problem with finding interface patch"
            //    << abort(FatalError);
            Info<< "Problem with finding interface patch" << endl;
            //Foam::error::abort();
            beamSolve = 0;
        }

        // Read reference values

        //dimensionedScalar 
        //rhoRef( 
        //        "rhoRef", dimDensity, 
        //        couplingProperties.lookupOrDefault<scalar>("rhoRef",1.0) 
        //);
        //scalar rhoRef( couplingProperties.lookupOrDefault<scalar>("rhoRef",1.0) );
        rhoRef = couplingProperties.lookupOrDefault<scalar>("rhoRef",1.0);
        pRef = readScalar(mesh.solutionDict().subDict("PIMPLE").lookup("pRefValue"));
        Info<< "  Reference density  : " << rhoRef << endl;
        Info<< "  Reference pressure : " << pRef << endl;

        //label bladeDir( couplingProperties.lookupOrDefault<label>("bladeDir",0) );
        bladeDir = couplingProperties.lookupOrDefault<label>("bladeDir",0);
        Info<< "  Blade axis (beamdyn coords): " << bladeDir << endl;

        vector coordMap = couplingProperties.lookupOrDefault<vector>("coordinateMapping",vector(0,1,2));
        for( int i=0; i<3; ++i )
        {
            OFtoIEC[i] = label(coordMap[i]);
            IECtoOF[OFtoIEC[i]] = i;
        }
        Info<< "  Coordinate mapping from OpenFOAM to IEC : " << OFtoIEC << endl;
        Info<< "  Coordinate mapping from IEC to OpenFoam : " << IECtoOF << endl;

        twoD = couplingProperties.lookupOrDefault<Switch>("twoD",0);
        Info<< "  Constrain rotation to about the blade axis (2D) : " << twoD << endl;

        //scalar bladeR( couplingProperties.lookupOrDefault("R",-1) );
        //scalar bladeR0( couplingProperties.lookupOrDefault("R0",-1) );
        bladeR0 = couplingProperties.lookupOrDefault<scalar>("R0",-1); // default < 0 : use minimum input r value
        bladeR  = couplingProperties.lookupOrDefault<scalar>("R",-1);

        //dimensionedVector
        //origin( 
        //        "origin", dimLength, 
        //        couplingProperties.lookupOrDefault<vector>("origin",vector::zero)
        //);
        //vector origin( couplingProperties.lookupOrDefault<vector>("origin",vector::zero) );
        origin = couplingProperties.lookupOrDefault<vector>("origin",vector::zero);
        Info<< "  Origin (rotation/moment reference) : " << origin << endl;

        //loadMultiplier = couplingProperties.lookupOrDefault<scalar>("loadMultiplier",1.0);
        //Info<< "  Load multiplier (FOR DEVELOPMENT) : " << loadMultiplier << endl;

        fluidSolve = couplingProperties.lookupOrDefault<Switch>("fluidSolve",1);
        if(beamSolve) beamSolve = couplingProperties.lookupOrDefault<Switch>("beamSolve",1);
        //if(!fluidSolve) Info<< "  SKIPPING fluid solution" << endl;
        //if(!beamSolve) Info<< "  SKIPPING structural dynamics solution" << endl;

        //prescribed_max_deflection = couplingProperties.lookupOrDefault<vector>
        //(
        //    "prescribed_max_deflection",
        //    vector::zero
        //);
        //prescribed_max_rotation = couplingProperties.lookupOrDefault<vector>
        //(
        //    "prescribed_max_rotation",
        //    vector::zero
        //);
        //Info<< "  Prescribed max deflection (when beamSolve==0): " 
        //    << prescribed_max_deflection << endl;
        //Info<< "  Prescribed max rotation   (when beamSolve==0): " 
        //    << prescribed_max_rotation << endl;
        //prescribed_max_rotation *= Foam::constant::mathematical::pi/180;

    }

    //*********************************************************************************************

    void stop()
    {
        Info<< "================================" << endl;
        Info<< "Stopping BeamDyn" << endl;
        Info<< "================================\n" << endl;

        if(Pstream::master()) beamDynWriteState("BeamDynState_last.dat");

        //Info<< "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" << endl;
        if(Pstream::master()) beamDynStop();
        //Info<< "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n" << endl;

        delete pos0_ptr;
        delete crv0_ptr;
        delete pos_ptr;
        delete crv_ptr;
        delete disp_ptr;
        delete adisp_ptr;
        delete r_ptr;
        delete [] h_ptr;

        if (Pstream::master())
        {
            loadFile.close();
            dispFile.close();
        }
    }

    //*********************************************************************************************

    void update( double t, double dt )
    {
//        Info<< "================================" << endl;
//        Info<< "| Calling BeamDyn update" << endl;
//        Info<< "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" << endl;
        currentTime = t; // do we really need to save this?
        currentDeltaT = dt;
        if(Pstream::master()) 
        {
            beamDynStep( &dt ); // NOTE: dt isn't actually used!
        }
//        Info<< "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << nl << endl;
        updateNodePositions();
    }

    //*********************************************************************************************

//    void updatePrescribedDeflection( double t )/*{{{*/
//    {
//        if(Pstream::master()) 
//        {
//            // prescribed deflection for testing
//            // y = ymax * (x/L)^2
//            //
//            // - at t=1.0, deflection is 100%
//            // - calculate shortening of the blade, normalized by blade length:
//            // b = (3*a**4 + sqrt(9*a**8+2*a**6))**(1./3.)
//            // xmax = b/2**(2./3.)/a**2 - 1/2**(1./3.)/b
//
//            double L = bladeR - bladeR0;
//            double a = 2*prescribed_max_deflection[1]*t; // prescribed deflection only in y-dir for now
//            double ang = prescribed_max_rotation[0]*currentDeltaT; // rotation increment
//            Info<< "current delta t : " << currentDeltaT << endl;
//            if (a != 0)
//            {
//                double b = Foam::pow( 3*Foam::pow(a,4) 
//                                    + Foam::sqrt( 9*Foam::pow(a,8) + 2*Foam::pow(a,6) )
//                                 ,1.0/3.0);
//                double xmax = L * ( b/Foam::pow(2,2.0/3.0)/pow(a,2) - 1.0/(b*pow(2,1.0/3.0)) );
//                double ymax = L * prescribed_max_deflection[1]*t;
//
//                double x[3], tmp[3], u[3];//, p[3];
//                for (int i=0; i<3; ++i)
//                {
//                    u[i] = 0;
//                }
//
//                Info<< "Prescribed tip position at time " << t 
//                    << " : x/L,y/L = " << xmax/L << " " << ymax/L 
//                    << " with twist increment = " << ang
//                    << endl;
//
//                for (int inode=0; inode<nnodes; ++inode)
//                {
//                    //beamDynGetInitNode0Position( &inode, x0, tmp );
//                    beamDynGetNode0Position( &inode, x, tmp );
//                    //double xi = x[bladeDir]/L;
//                    double s = (*pos0_ptr)[inode][bladeDir] / L; //parametric coordinate from 0 to 1
//
//                    // NOTE: prescribed deflection is only in y-dir for now
//                    //       bladeDir is implied to be x-dir
//                    u[bladeDir] = (xmax-L)*s; //xi * xmax - (*pos0_ptr)[inode][bladeDir];
//                    u[1] = ymax*s*s - (*pos0_ptr)[inode][1];
//                    u[2] = 0.0;
//
//    // DEBUG
//    //                Info<< "  initial position is " << (*pos0_ptr)[inode] 
//    //                    << "  (s=" << s << ")"
//    //                    << endl;
//    //
//    //                Info<< "  setting disp for node " << inode+1
//    //                    << " at " << x[0] << " " << x[1] << " " << x[2]
//    //                    << " to " << u[0] << " " << u[1] << " " << u[2]
//    //                    << endl;
//
//                    beamDynSetDisplacement( &inode, u );
//
//                    // Since we're prescribing motion and not actually running the BeamDyn solver,
//                    // the rotation matrix is never actually updated. We manually set it here.
//                    //double dydx = a*xi;
//                    //double ang = Foam::atan(a*xi);
//                    ang = Foam::atan(a*s);
//                    beamDynSetZRotationMatrix( &inode, &ang );
//
//                }
//            }
//
//            // For now, just overwrite the rotation matrix if we're doing a test with x-axis rotation
//            if (ang != 0)
//            {
//                Info<< "Specifying twist increment = " << ang << endl;
//                for (int inode=0; inode<nnodes; ++inode)
//                {
////                    ang = prescribed_max_rotation[0]*t + (*rot0_ptr)[inode][0]; // this is a rotation angle, not the absolute orientation!
////                    ang = prescribed_max_rotation[0]*currentDeltaT + (*rot0_ptr)[inode][0];
////                    Info<< "Prescribing twist angle " << ang*180/pi << endl;
//                    beamDynSetXRotationMatrix( &inode, &ang );
//                }
//            }
//        }
//
//        updateNodePositions();
//    }/*}}}*/

    //*********************************************************************************************

    void write( std::string timeName, bool writeNow=true )
    {
        if (!writeNow || !Pstream::master()) return;

        //std::string fname("BeamDynState_" + timeName + ".dat");
        std::string fname(timeName + "/BeamDynState.dat");
        beamDynWriteState( fname.c_str() );
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Interface calculations
    //
    ///////////////////////////////////////////////////////////////////////////////////////////////

    // Retrieve disp array from the BeamDyn library
    // TODO: clean this up?
    // - updates pos/rot, used to calculate disp
    // - updates disp, accessed through BD::disp() in beamDynInterfacePointPatch::updateCoeffs()
    // - updates r, used by updateSectionLoads()
    // - writes displacement at the starting time step to disp.out
    void updateNodePositions()
    {
        scalarList &r    = *r_ptr;
        vectorList &pos0 = *pos0_ptr;
        //vectorList &rot0 = *rot0_ptr; // not used
        vectorList &pos  = *pos_ptr;
        vectorList &crv  = *crv_ptr;
        vectorList &disp = *disp_ptr;
        vectorList &adisp= *adisp_ptr;

//        Info<< "Retrieving node positions for the next iteration" << endl;
        if(Pstream::master())
        {
            if (first) Info<< "Initial displacements from BeamDyn lib [m,deg]: " << endl;
            else dispFile << currentTime;

            // --loop over nodes in the BeamDyn blade model (assumed single element)
            //   TODO: handle multiple elements
            double lin_disp[3], ang_disp[3];
            //double ang;
//            double R[9]; //rotation matrix from BeamDyn
            for( int inode=0; inode<nnodes; ++inode ) 
            {

                // get node linear/angular displacement [m, crv]
                beamDynGetNode0Displacement( &inode, lin_disp, ang_disp );

// ROTATION SHOULD BE APPLIED FOR THE POINTS ON THE INTERFACE PATCH!!!
//                beamDynGetNode0RotationMatrix( &inode, R );
// DEBUG
// these are identity before the first iteration and approximately the identity matrix afterwards
//                Info << "DEBUG R matrix " << inode << " : \n";
//                Info << " R=[" << R[0] << " " << R[1] << " " << R[2] << "; ...\n";
//                Info << "    " << R[3] << " " << R[4] << " " << R[5] << "; ...\n";
//                Info << "    " << R[6] << " " << R[7] << " " << R[8] << "]\n";
//
//                disp[inode] = vector::zero;
//                for( int j=0; j<3; ++j ) {
//                    for( int i=0; i<3; ++i )
//                    {
//                        disp[inode].component(j) += R[3*j+i] * lin_disp[i];
//                    }
//
//                    //disp[inode].component(j) = lin_disp[j];
//                }
//                if(twoD) disp[inode].component(bladeDir) = 0.0;
                for( int i=0; i<3; ++i )
                {
                    disp[inode].component(i) = lin_disp[i];
                    adisp[inode].component(i)= ang_disp[i];
                }

                // used for sectional loads calculation
                // note: pos0 and disp are in IEC coords
                r[inode] = pos0[inode][bladeDir] + disp[inode].component(bladeDir);

                if (first) // print out initial displaced config, either 0's or (hopefully) repeated on restart
                {
                    //NOTE: for prescribed motion, since the beamdyn solve routine isn't called, the rotation parameters
                    //  are not properly updated so ang_disp will be inaccurate; just use the rotation matrix R instead.
                    Info<< " " << lin_disp[0] 
                        << " " << lin_disp[1] 
                        << " " << lin_disp[2]
                        << " " << crvToRad(ang_disp[0])*radToDeg
                        << " " << crvToRad(ang_disp[1])*radToDeg
                        << " " << crvToRad(ang_disp[2])*radToDeg
                        << endl;
                }
                else // write subsequent displacements to file
                {
                    // note: these are in the beamdyn frame
                    dispFile << " " << lin_disp[0] 
                             << " " << lin_disp[1] 
                             << " " << lin_disp[2];
                    dispFile << " " << crvToRad(ang_disp[0])*radToDeg 
                             << " " << crvToRad(ang_disp[1])*radToDeg 
                             << " " << crvToRad(ang_disp[2])*radToDeg;
                }

            }// loop over beam nodes

            if (first) Info<< endl;
            else dispFile << std::endl;

        }// if Pstream::master

        Pstream::scatter(r);
        Pstream::scatter(pos); // verified that this works
        Pstream::scatter(crv);
        Pstream::scatter(disp);
        Pstream::scatter(adisp);
    }

    //*********************************************************************************************

    void calculateShapeFunctions( const pointField& pf )
    {
        int nSurfNodes = pf.size(); // number of local points on interface patch
        h_ptr = new double[nSurfNodes*nnodes];
        //double &h = *h_ptr;

        // read/write shape function information
        std::string fname("shape.dat");
        if( Pstream::parRun() ) 
        {
            std::stringstream ss;
            ss << "processor" << Pstream::myProcNo() << "/" << fname;
            fname = ss.str();
        }

        if( restarted )
        {
            // Make sure we only calculate this at the beginning using the undeformed surface nodes
            Info<< "Skipping shape function calculation for restarted simulation" << endl;

            // read saved shape function
            std::ifstream hFile(fname, std::ios::in | std::ios::binary);
            //hFile.read((char *) &h_ptr, sizeof(h_ptr)); // C-style
            hFile.read( reinterpret_cast<char*>(&h_ptr), sizeof(h_ptr) );

            return;
        }

        std::ofstream hFile(fname, std::ios::out | std::ios::binary);

        scalarList &r = *r_ptr;
        if( bladeR0 < 0.0 ) bladeR0 = r[0];
        if( bladeR  < 0.0 ) bladeR  = r[nnodes-1];
        Info<< "r: " << r << endl;
        //Pout<< "Blade span : " << bladeR0 << " " << bladeR << endl;

        if( nSurfNodes > 0 )
        {
            Pout << "calculating shape functions for " << nSurfNodes << " surface nodes" << endl;

            // DEBUG
            scalar hmin( 9e9);
            scalar hmax(-9e9);
            scalar smin( 9e9);
            scalar smax(-9e9);

            double s;
            double L_2 = (bladeR-bladeR0)/2.0;
            double hi[nnodes];

            double num, den;
            double GLL[nnodes];
            Info<< "GLL pts : "; // DEBUG
            for( int i=0; i<nnodes; ++i )
            {
                GLL[i] = 2.0*(r[i]-r[0])/(r[nnodes-1]-r[0]) - 1.0;
                Info<< " " << GLL[i];
            }
            Info<< endl;

            forAll( pf, ptI )
            {
                //s = ( pf[ptI].component(bladeDir) - bladeR0 ) / L_2 - 1.0;
                s = ( pf[ptI].component(IECtoOF[bladeDir]) - bladeR0 ) / L_2 - 1.0;
                smin = min(smin,s);
                smax = max(smax,s);
                //beamDynGetShapeFunctions( &s, hi ); // this only works on the master node...
                //vvvvvvvvvv Code snippet from BeamDyn diffmtc subroutine vvvvvvvvvv
                for( int j=0; j<nnodes; ++j )
                {
                    hi[j] = 0.0;
                    num = 1.0;
                    den = 1.0;
                    if( fabs(s-GLL[j]) <= eps )
                    {
                        hi[j] = 1.0;
                    }
                    else
                    {
                        for( int i=0; i<nnodes; ++i )
                        {
                            if( i != j )
                            {
                                den *= (GLL[j] - GLL[i]);
                                num *= (s - GLL[i]);
                            }
                        }
                        hi[j] = num/den;
                    }
                }
                //^^^^^^^^^^^^^^^^^^^ End of code snippet ^^^^^^^^^^^^^^^^^^^^^^^^^^

                scalar hsum(0);
                for( int inode=0; inode < nnodes; ++inode )
                {
                    h_ptr[ptI*nnodes + inode] = hi[inode];
                    hsum += hi[inode];

                    hFile << hi[inode];
                }
                hmin = min(hmin,hsum);
                hmax = max(hmax,hsum);

            }// loop over surface nodes
            Pout<< "-- sum(h) : [ " << hmin << ", " << hmax << " ]" 
                << " for s : [ " << smin << ", " << smax << " ]"
                << endl;

        }//if nSurfNodes > 0

        //hFile.write((char *) &h_ptr, sizeof(h_ptr)); // C-style
        hFile.write( reinterpret_cast<char*>(&h_ptr), sizeof(h_ptr) );

    }

    //*********************************************************************************************

    // for testing...
    void updateSectionLoads( const double* F,
                             const double* M )
    {
        if (!Pstream::master()) return;

        // for safety, define new variables to pass
        double Ftot[3], Mtot[3];
        for( int i=0; i<3; ++i )
        {
            Ftot[i] = F[i];
            Mtot[i] = M[i];
        }

        //for( int ig=0; ig<nnodes-1; ++ig ) 
        //{
        //    beamDynSetDistributedLoad(&ig, Ftot, Mtot);
        //}

        // Only update nodes at x>=0 
        for( int ig=nnodes/2; ig<nnodes-1; ++ig )  // only set interior gauss pts
        {
            //Info<< "gauss pt " << ig << " btwn "
            //    << " " << (*pos0_ptr)[ig][0]
            //    << " " << (*pos0_ptr)[ig][1]
            //    << " " << (*pos0_ptr)[ig][2]
            //    << " and "
            //    << " " << (*pos0_ptr)[ig+1][0]
            //    << " " << (*pos0_ptr)[ig+1][1]
            //    << " " << (*pos0_ptr)[ig+1][2]
            //    << endl;
            beamDynSetDistributedLoad(&ig, Ftot, Mtot);
        }

        first = false;
    }


    // actual interface routine
    void updateSectionLoads( const dynamicFvMesh& mesh, 
                             const volScalarField& p, 
                             const incompressible::turbulenceModel& turbulence )
    {
        Info<< "Calculating section loads for BeamDyn" << endl;

        scalarList &r = *r_ptr;
        scalar p0( pRef / rhoRef );

        // setup arrays, pointers
        double r0, r1;
        const polyPatch& bladePatch = mesh.boundaryMesh()[patchID];
        const vectorField& bladePatchNormals = mesh.Sf().boundaryField()[patchID];

        // calculate shear stress
        //   note: devReff returns the effective stress tensor including the laminar stress
        //   note: face normals point _outside_ the computational domain

//        Info<< "Calculating surface shear stresses" << endl;
//        const volSymmTensorField Reff(turbulence.devReff());
//        vectorField bladePatchShearStress = 
//            (
//                -mesh.Sf().boundaryField()[interfacePatchID]
//                /mesh.magSf().boundaryField()[interfacePatchID]
//            ) & Reff.boundaryField()[interfacePatchID];

        // Face normals point into solid surface, i.e., outward from fluid volume, 
        // i.e. the direction the fluid is pushing on the wall.
        // This matches the 'forces' function object implementation 
        // in Foam::forces::calcForcesMoment() at
        //   ~/OpenFOAM/OpenFOAM-2.3.1/src/postProcessing/functionObjects/forces/forces/forces.C
        // - also, no need to normalize by magSf since we multiply by mag(Sf) later
        const volSymmTensorField Reff(turbulence.devReff());
        //THIS DOESN'T COMPILE WITH CLANG:
        //vectorField bladePatchShearStress = 
        //    mesh.Sf().boundaryField()[interfacePatchID]
        //    & Reff.boundaryField()[interfacePatchID];
        vectorField bladePatchShearStress( bladePatchNormals & Reff.boundaryField()[patchID] );

        //
        // --loop over nodes in the BeamDyn blade model, assumed single element
        //   i.e., nnodes = nodes_elem = order_elem+1 = ngp+1
        //
        //loadFile << currentTime; // at this point, still equal to t at beginning of time step
        loadFile << currentTime + currentDeltaT;
//        Info<< "Integrating sectional loads" << endl;
        if(first) Info<< "Initial info:" << endl;
        for( int ig=0; ig<nnodes-1; ++ig ) 
        {
            vector Fp(vector::zero);
            vector Fv(vector::zero);
            vector Mp(vector::zero);
            vector Mv(vector::zero);

            r0 = r[ig];
            r1 = r[ig+1];
            scalar dr = r1 - r0; // note: this is the width of the integration segment 
                                 //       corresponding to the Gauss pt

            // 
            // --loop over faces on interface patch
            //
            int nFacesFound = 0;
            forAll( bladePatch, faceI )
            {
                vector rc( bladePatch.faceCentres()[faceI] );
                //if( rc[bladeDir] >= r0 && rc[bladeDir] < r1 )
                if( rc[IECtoOF[bladeDir]] >= r0 && rc[IECtoOF[bladeDir]] < r1 )
                {
                    vector Sf( bladePatchNormals[faceI] ); // surface normal
                    vector dm( rc - origin );

                    vector dFp = rhoRef * Sf * (p.boundaryField()[patchID][faceI] - p0) / dr;
                    Fp += dFp;
                    Mp += dm ^ dFp;

                    vector dFv = mag(Sf) * bladePatchShearStress[faceI] / dr;
                    Fv += dFv;
                    Mv += dm ^ dFv;

                    nFacesFound += 1;
                }

            }// end of face loop

            Pstream::gather(Fp, sumOp<vector>()); // These are the DISTRIBUTED loads!
            Pstream::gather(Fv, sumOp<vector>());
            Pstream::gather(Mp, sumOp<vector>());
            Pstream::gather(Mv, sumOp<vector>());

            double Ftot[3], Mtot[3]; // these should be in IEC coordinates
            if(Pstream::master())
            {
                for (int i=0; i<3; ++i) 
                {
                    //Ftot[i] = loadMultiplier*(Fp[i] + Fv[i]);
                    //Mtot[i] = loadMultiplier*(Mp[i] + Mv[i]);
                    //Ftot[i] = Fp[i] + Fv[i];
                    //Mtot[i] = Mp[i] + Mv[i];
                    Ftot[OFtoIEC[i]] = Fp[i] + Fv[i];
                    Mtot[OFtoIEC[i]] = Mp[i] + Mv[i];
                }

                // Update F_foam and M_foam in BeamDyn library
                //beamDynSetDistributedLoadAtNode(&inode, Ftot, Mtot);
                beamDynSetDistributedLoad(&ig, Ftot, Mtot);

                // Output loads are in the beamdyn frame
                loadFile << " " << Ftot[0] 
                         << " " << Ftot[1] 
                         << " " << Ftot[2];
                loadFile << " " << Mtot[0] 
                         << " " << Mtot[1] 
                         << " " << Mtot[2];
            }

            if (first)
            {
                Pstream::gather(nFacesFound, sumOp<int>());
                //Info<< "  seg " << ig
                //    << " with " << nFacesFound << " faces btwn " << r0 << " " << r1 << ":"
                //    << " (" << Ftot[0] << " " << Ftot[1] << " " << Ftot[2] << ") " 
                //    << " (" << Mtot[0] << " " << Mtot[1] << " " << Mtot[2] << ") " 
                //    << endl;
                Info<< "segment " << ig
                    << " with " << nFacesFound << " faces"
                    << " between " << r0 << " " << r1 << endl;
            }

        } // end loop over beamdyn nodes

        loadFile << std::endl;

        first = false;

    }

    // disp is the interpolated displacement vector that will be rotated
    // adisp is the interpolated angular displacement
    void rotateDisplacementVector(double* disp, double* adisp)
    {
        double disp0[3];
        for( int i=0; i<3; ++i) { disp0[i] = disp[i]; }

        double Rmat[9];
        beamDynGetRotationMatrix( adisp, Rmat ); // fortran subroutine will convert from angular disp to CRV params

        for( int i=0; i<3; ++i )
        {
            disp[i] = 0;
            for( int j=0; j<3; ++j )
            {
                disp[i] += Rmat[3*i+j]*disp0[j];
            }
        }

        Info<< "adisp=" 
            << adisp[0] << ","
            << adisp[1] << ","
            << adisp[2] << " : "
            << " disp : " 
            << disp0[0] << ","
            << disp0[1] << ","
            << disp0[2] << " --> "
            << disp[0] << ","
            << disp[1] << ","
            << disp[2] << endl;
    }

} // end of BD namespace
