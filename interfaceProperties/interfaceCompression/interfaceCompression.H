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

#ifndef interfaceCompression_H
#define interfaceCompression_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class interfaceCompressionWeight Declaration
\*---------------------------------------------------------------------------*/

class interfaceCompressionLimiter
{

public:

    interfaceCompressionLimiter(Istream&)
    {}

    scalar limiter
    (
        const scalar cdWeight,
        const scalar faceFlux,
        const scalar phiP,
        const scalar phiN,
        const vector&,
        const scalar
    ) const
    {
        // Quadratic compression scheme
        //return min(max(4*min(phiP*(1 - phiP), phiN*(1 - phiN)), 0), 1);

        // Quartic compression scheme
        return
            min(max(
            1 - max(sqr(1 - 4*phiP*(1 - phiP)), sqr(1 - 4*phiN*(1 - phiN))),
            0), 1);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
