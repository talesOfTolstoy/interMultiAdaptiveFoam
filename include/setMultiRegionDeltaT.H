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

if (adjustTimeStep)
{
    if (CoNum == -GREAT)
    {
        CoNum = SMALL;
    }

    if (DiNum == -GREAT)
    {
        DiNum = SMALL;
    }

    scalar maxDeltaTFluid = maxCo/(CoNum + SMALL);
    scalar maxDeltaTSolid = maxDi/(DiNum + SMALL);

    scalar deltaTFluid =
        min
        (
            min(maxDeltaTFluid, 1.0 + 0.1*maxDeltaTFluid),
            1.2
        );

    runTime.setDeltaT
    (
        min
        (
            min(deltaTFluid, maxDeltaTSolid)*runTime.deltaT().value(),
            maxDeltaT
        )
    );

    Info<< "deltaT = " <<  runTime.deltaT().value() << endl;
}

// ************************************************************************* //
