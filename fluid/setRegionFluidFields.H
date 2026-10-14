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
dynamicFvMesh& mesh = fluidRegions[i];

immiscibleIncompressibleTwoPhaseMixture& mixture(mixtureFluid[i]);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();

volScalarField& rho = rhoFluid[i];
volVectorField& U = UFluid[i];
surfaceScalarField& phi = phiFluid[i];
surfaceScalarField& alphaPhi10 = alphaPhi10Fluid[i];

incompressible::turbulenceModel& turbulence = turbulenceFluid[i];
volScalarField& p_rgh = p_rghFluid[i];
volScalarField& p = pFluid[i];

volScalarField& mDot = mDotFluid[i];
volScalarField& qDot = qDotFluid[i];

volScalarField& T = TFluid[i];
volScalarField& amrIndicator = amrIndicatorFluid[i];

volScalarField& kappa = kappaFluid[i];

const bool& alphaRestart = alphaRestartFluid[i];

scalar& cumulativeContErr = cumulativeContErrFluid[i];

volScalarField& gh = ghFluid[i];
surfaceScalarField& ghf = ghfFluid[i];
const dimensionedScalar& ghRef = ghRefFluid[i];

surfaceScalarField& brackFS = brackFSFluid[i];

fv::options& fvOptions = fluidFvOptions[i];

// bool frozenFlow = frozenFlowFluid[i];

const label pRefCell = pRefCellFluid[i];
const scalar pRefValue = pRefValueFluid[i];

dictionary phase1 = mixture.subDict(mixture.phase1Name());
dictionary phase2 = mixture.subDict(mixture.phase2Name());
dictionary pcProperties = mixture.subDict("PhaseChangeProperties");

dimensionedScalar k1
(
    "kappa",
    dimensionSet(1,1,-3,-1,0,0,0),
    readScalar(phase1.lookup("kappa"))
);
dimensionedScalar cp1
(
    "Cp",
    dimensionSet(0,2,-2,-1,0,0,0),
    readScalar(phase1.lookup("Cp"))
);

dimensionedScalar k2
(
    "kappa",
    dimensionSet(1,1,-3,-1,0,0,0),
    readScalar(phase2.lookup("kappa"))
);
dimensionedScalar cp2
(
    "Cp",
    dimensionSet(0,2,-2,-1,0,0,0),
    readScalar(phase2.lookup("Cp"))
);

dimensionedScalar Tsat
(
    "Tsat",
    dimTemperature,
    readScalar(pcProperties.lookup("Tsat"))
);

dimensionedScalar hlv
(
    "hlv",
    dimensionSet(0,2,-2,0,0,0,0),
    readScalar(pcProperties.lookup("hlv"))
);

dimensionedScalar diff
(
 "diff",
 dimensionSet(0,2,0,0,0,0,0),
 readScalar(pcProperties.lookup("diff"))
 );

dimensionedScalar maxQdot
(
 "maxQdot",
 dimensionSet(1,-1,-3,0,0,0,0),
 readScalar(pcProperties.lookup("maxQdot"))
 );

dimensionedScalar R_g
(
 "R_g",
 dimensionSet(0,2,-2,-1,0,0,0),
 readScalar(pcProperties.lookup("R_g"))
 );

// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

// volScalarField rhoCp
// (
//     IOobject
//     (
//         "rhoCp",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT
//     ),
//     alpha1*rho1*cp1 + alpha2*rho2*cp2,
//     alpha1.boundaryField().types()
// );
// rhoCp.oldTime();


//Energy flux

// surfaceScalarField rhoPhiCp
// (
//     IOobject
//     (
//         "rhoPhiCp",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//     fvc::interpolate(rhoCp)*phi
// );

surfaceScalarField alphaPhiUn
(
    IOobject
    (
        "alphaPhiUn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(phi.dimensions(), Zero)
);

dimensionedScalar sigmaEvap
(
 "sigmaEvap",
 dimless,
 readScalar(pcProperties.lookup("sigmaEvap"))
 );

pisoControl& piso(pisoFluid[i]);

bool correctPhi
(
    piso.dict().getOrDefault("correctPhi", mesh.dynamic())
);

// bool checkMeshCourantNo
// (
//     piso.dict().getOrDefault("checkMeshCourantNo", false)
// );

// bool moveMeshOuterCorrectors
// (
//     piso.dict().getOrDefault("moveMeshOuterCorrectors", false)
// );

IOMRFZoneList& MRF = MRFfluid[i];

// surfaceScalarField& alphaPhi0(alphaPhi10Fluid[i]);

// MULES Correction
tmp<surfaceScalarField> talphaPhi1Corr0;

autoPtr<surfaceVectorField> Uf;

if (fluidRegions[i].dynamic())
{
    Info<< "Constructing face velocity Uf\n" << endl;

    Uf.reset
    (
        new surfaceVectorField
        (
            IOobject
            (
                "Uf",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(UFluid[i])
        )
    );
}

volScalarField rAU
(
    IOobject
    (
        "rAU",
        runTime.timeName(),
        fluidRegions[i],
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fluidRegions[i],
    dimensionedScalar("rAU", dimTime/dimDensity, 1)
);
