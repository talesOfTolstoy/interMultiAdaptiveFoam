/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of interMultiAdaptiveFoam.

    interMultiAdaptiveFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    interMultiAdaptiveFoam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with interMultiAdaptiveFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#ifndef constantSurfaceTension_H
#define constantSurfaceTension_H

#include "surfaceTensionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace surfaceTensionModels
{

/*---------------------------------------------------------------------------*\
                           Class constant Declaration
\*---------------------------------------------------------------------------*/

class constant
:
    public surfaceTensionModel
{
    // Private data

        //- Surface tension coefficient
        dimensionedScalar sigma_;


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from dictionary and mesh
        constant
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~constant();


    // Member Functions

        //- Surface tension coefficient
        virtual tmp<volScalarField> sigma() const;

        //- Update surface tension coefficient from given dictionary
        virtual bool readDict(const dictionary& dict);

        //- Write in dictionary format
        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceTensionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
