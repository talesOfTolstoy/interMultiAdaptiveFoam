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
// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
 fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
 + MRF.DDt(rho, U)
 + turbulence.divDevRhoReff(rho,U)
 ==
 fvOptions(rho, U)
 );
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

Info << "Min(U) = " << min(U).value() << " Max(U) = " << max(U).value() << endl;

