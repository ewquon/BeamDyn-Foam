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

Class
    Foam::beamDynInterfacePointPatchVectorField

Description
    Foam::beamDynInterfacePointPatchVectorField

    Created from oscillatingDisplacementPointPatchVectorField
    Originally in src/fvMotionSolver/pointPatchFields/derived/beamDynInterface

SourceFiles
    beamDynInterfacePointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef beamDynInterfacePointPatchVectorField_H
#define beamDynInterfacePointPatchVectorField_H

#include "fixedValuePointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class beamDynInterfacePointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class beamDynInterfacePointPatchVectorField
:
    public fixedValuePointPatchField<vector>
{
    //const scalar pi(Foam::constant::mathematical::pi);

    // Private data

public:

    //- Runtime type information
    TypeName("beamDynInterface");


    // Constructors

        //- Construct from patch and internal field
        beamDynInterfacePointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        beamDynInterfacePointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        beamDynInterfacePointPatchVectorField
        (
            const beamDynInterfacePointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector> > clone() const
        {
            return autoPtr<pointPatchField<vector> >
            (
                new beamDynInterfacePointPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        beamDynInterfacePointPatchVectorField
        (
            const beamDynInterfacePointPatchVectorField&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<vector> > clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<vector> >
            (
                new beamDynInterfacePointPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
