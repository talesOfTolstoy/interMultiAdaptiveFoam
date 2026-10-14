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
#ifndef surfaceTensionModel_H
#define surfaceTensionModel_H

#include "regIOobject.H"
#include "dimensionedTypes.H"
#include "volFieldsFwd.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class surfaceTensionModel Declaration
\*---------------------------------------------------------------------------*/

class surfaceTensionModel
:
    public regIOobject
{
protected:

    // Protected member data

        //- Reference to mesh
        const fvMesh& mesh_;


    // Protected member functions

        static const dictionary& sigmaDict(const dictionary& dict)
        {
            return dict.subDict("sigma");
        }


public:

    //- Runtime type information
    TypeName("surfaceTensionModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            surfaceTensionModel,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (dict, mesh)
        );


    // Static data members

        //- Surface tension coefficient dimensions
        static const dimensionSet dimSigma;


    // Constructors

        // Construct from mesh
        surfaceTensionModel(const fvMesh& mesh);


    //- Destructor
    virtual ~surfaceTensionModel();


    // Selectors

        static autoPtr<surfaceTensionModel> New
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        //- Surface tension coefficient
        virtual tmp<volScalarField> sigma() const = 0;

        //- Update surface tension coefficient from given dictionary
        virtual bool readDict(const dictionary& dict) = 0;

        //- Write in dictionary format
        virtual bool writeData(Ostream& os) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
