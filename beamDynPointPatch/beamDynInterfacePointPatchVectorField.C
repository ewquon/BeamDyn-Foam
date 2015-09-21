/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "beamDynInterfacePointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "polyMesh.H"
#include "displacementMotionSolver.H"
#include "vectorList.H"

#include "beamDynInterface.H" // BD namespace

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

beamDynInterfacePointPatchVectorField::
beamDynInterfacePointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF)
{
    Info<< "Created instance of beamDynInterfacePointPatchVectorField (2)" << endl;
}


beamDynInterfacePointPatchVectorField::
beamDynInterfacePointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict)
{
    Info<< "Created instance of beamDynInterfacePointPatchVectorField (3)" << endl;
    if (!dict.found("value"))
    {
        updateCoeffs();
    }
}


beamDynInterfacePointPatchVectorField::
beamDynInterfacePointPatchVectorField
(
    const beamDynInterfacePointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper)
{
    Info<< "Created instance of beamDynInterfacePointPatchVectorField (4)" << endl;
}


beamDynInterfacePointPatchVectorField::
beamDynInterfacePointPatchVectorField
(
    const beamDynInterfacePointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF)
{
    Info<< "Created instance of beamDynInterfacePointPatchVectorField (2,ptf)" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void beamDynInterfacePointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

//    double minTwist=9e9, maxTwist=-9e9;

    vector pos(vector::zero); // interpolated linear displacement
    vector crv(vector::zero); // interpolated angular displacement
    vector p(vector::zero); // position vector from beam curve to a surface node

    //- patch coordinates, if needed
    //const labelList& meshPoints = patch().meshPoints(); // returns polyPatch_.meshPoints(), i.e. node IDs
    const pointField& localPoints = patch().localPoints(); // returns polyPatch_.localPoints(), i.e. DISPLACED node coords
    //Info<< "- TEST: node coords, pt disp vec : " << localPoints[0] << " " << this->operator[](0) << endl;

    //- current time, if needed
    //const polyMesh& mesh = this->dimensionedInternalField().mesh()();
    //const Time& t = mesh.time();

    vectorList& posList = BD::pos(); // position of BD nodes in OpenFOAM coords
    vectorList& crvList = BD::crv(); // orientation of BD nodes in OpenFOAM coords

    vectorList& pList = BD::p();        // surface offset vector in OpenFOAM coords
    //vectorList& x1List = BD::x1();      // position vector along beam axis in OpenFOAM coords

    //
    // --loop over all surface nodes
    //
    forAll(*this, ptI)
    {
        // interpolate displacement from pre-calculated shape function
        pos = vector::zero; // in OpenFoam coords
        crv = vector::zero; // in OpenFoam coords
        for( int inode=0; inode<BD::N(); ++inode )
        {
            for( int i=0; i<3; ++i )
            {
                pos.component(i) += BD::h()[ptI*BD::N()+inode] * posList[inode].component(i);
                crv.component(i) += BD::h()[ptI*BD::N()+inode] * crvList[inode].component(i);
            }
        }

        // calculate vector components
        //vector x1 = x1List[ptI] + u;    // in OpenFoam coords
        p = pList[ptI];                 // in OpenFoam coords

/////////////////////////////////////////////////////////////////////
// TODO: general rotations, retrieve rotation matrix -- need VABS...
/////////////////////////////////////////////////////////////////////
        BD::rotateVector( p, crv );     // rotates p, in OpenFoam coords

// for now, simple rotation about x0-axis (in BD frame) with cross sections remaining planar
//        ang = a.component(BD::bladeDirection());
//        minTwist = min(minTwist,ang);
//        maxTwist = max(maxTwist,ang);
//        v[BD::openfoamDir(0)] = 0;
//        v[BD::openfoamDir(1)] = localPoints[ptI].component(BD::openfoamDir(1)) * Foam::cos(ang) 
//                              - localPoints[ptI].component(BD::openfoamDir(2)) * Foam::sin(ang)
//                              - localPoints[ptI].component(BD::openfoamDir(1));
//        v[BD::openfoamDir(2)] = localPoints[ptI].component(BD::openfoamDir(1)) * Foam::sin(ang) 
//                              + localPoints[ptI].component(BD::openfoamDir(2)) * Foam::cos(ang)
//                              - localPoints[ptI].component(BD::openfoamDir(2));
//

        // Apply displacement
//        //this->operator[](ptI) = u;
//        this->operator[](ptI) = u + v;

//        this->operator[](ptI) = x1 + p - localPoints[ptI];
//        Info<< " x: " << localPoints[ptI] << endl;
//        Info<< " u: " << u << endl;
//        Info<< "x1: " << x1List[ptI] << endl;
//        Info<< " p: " << pList[ptI] << endl << endl;

        // new position is at x = pos + p
        // pointDisplacement values are relative to the current deformed configuration
        this->operator[](ptI) = pos + p - localPoints[ptI];

        // for testing
        if(BD::enforce2D()) this->operator[](ptI).component(BD::bladeDirection()) = 0.0;

        // update displacement vector lists
//        pList[ptI] = p;
//        x1List[ptI] = x1;

    } //loop over all (surface) nodes on patch

//    Info<< "- min/max twist : " 
//        << minTwist*180.0/Foam::constant::mathematical::pi << " " 
//        << maxTwist*180.0/Foam::constant::mathematical::pi << " deg" << endl;

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void beamDynInterfacePointPatchVectorField::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    beamDynInterfacePointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
