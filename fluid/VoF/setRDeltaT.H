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
{
    volScalarField& rDeltaT = trDeltaT.ref();

    const dictionary& pisoDict = piso.dict();

    scalar maxCo
    (
        pisoDict.getOrDefault<scalar>("maxCo", 0.9)
    );

    scalar maxAlphaCo
    (
        pisoDict.getOrDefault<scalar>("maxAlphaCo", 0.2)
    );

    scalar rDeltaTSmoothingCoeff
    (
        pisoDict.getOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.1)
    );

    label nAlphaSpreadIter
    (
        pisoDict.getOrDefault<label>("nAlphaSpreadIter", 1)
    );

    scalar alphaSpreadDiff
    (
        pisoDict.getOrDefault<scalar>("alphaSpreadDiff", 0.2)
    );

    scalar alphaSpreadMax
    (
        pisoDict.getOrDefault<scalar>("alphaSpreadMax", 0.99)
    );

    scalar alphaSpreadMin
    (
        pisoDict.getOrDefault<scalar>("alphaSpreadMin", 0.01)
    );

    label nAlphaSweepIter
    (
        pisoDict.getOrDefault<label>("nAlphaSweepIter", 5)
    );

    scalar rDeltaTDampingCoeff
    (
        pisoDict.getOrDefault<scalar>("rDeltaTDampingCoeff", 1.0)
    );

    scalar maxDeltaT
    (
        pisoDict.getOrDefault<scalar>("maxDeltaT", GREAT)
    );

    volScalarField rDeltaT0("rDeltaT0", rDeltaT);

    // Set the reciprocal time-step from the local Courant number
    rDeltaT.ref() = max
    (
        1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
        fvc::surfaceSum(mag(rhoPhi))()()
       /((2*maxCo)*mesh.V()*rho())
    );

    if (maxAlphaCo < maxCo)
    {
        // Further limit the reciprocal time-step
        // in the vicinity of the interface

        volScalarField alpha1Bar(fvc::average(alpha1));

        rDeltaT.ref() = max
        (
            rDeltaT(),
            pos0(alpha1Bar() - alphaSpreadMin)
           *pos0(alphaSpreadMax - alpha1Bar())
           *fvc::surfaceSum(mag(phi))()()
           /((2*maxAlphaCo)*mesh.V())
        );
    }

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    Info<< "Flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;

    if (rDeltaTSmoothingCoeff < 1.0)
    {
        fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
    }

    if (nAlphaSpreadIter > 0)
    {
        fvc::spread
        (
            rDeltaT,
            alpha1,
            nAlphaSpreadIter,
            alphaSpreadDiff,
            alphaSpreadMax,
            alphaSpreadMin
        );
    }

    if (nAlphaSweepIter > 0)
    {
        fvc::sweep(rDeltaT, alpha1, nAlphaSweepIter, alphaSpreadDiff);
    }

    Info<< "Smoothed flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;

    // Limit rate of change of time scale
    // - reduce as much as required
    // - only increase at a fraction of old time scale
    if
    (
        rDeltaTDampingCoeff < 1.0
     && runTime.timeIndex() > runTime.startTimeIndex() + 1
    )
    {
        rDeltaT = max
        (
            rDeltaT,
            (scalar(1) - rDeltaTDampingCoeff)*rDeltaT0
        );

        Info<< "Damped flow time scale min/max = "
            << gMin(1/rDeltaT.primitiveField())
            << ", " << gMax(1/rDeltaT.primitiveField()) << endl;
    }
}
