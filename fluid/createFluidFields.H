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
// Initialise fluid field pointer lists
PtrList<immiscibleIncompressibleTwoPhaseMixture> mixtureFluid
(
    fluidRegions.size()
);

PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
PtrList<surfaceScalarField> alphaPhi10Fluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> hRefFluid(fluidRegions.size());
PtrList<volScalarField> ghFluid(fluidRegions.size());
PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
PtrList<incompressible::turbulenceModel> turbulenceFluid(fluidRegions.size());
PtrList<volScalarField> p_rghFluid(fluidRegions.size());
PtrList<volScalarField> pFluid(fluidRegions.size());
PtrList<volScalarField> rhoFluid(fluidRegions.size());
PtrList<volScalarField> kappaFluid(fluidRegions.size());
PtrList<volScalarField> TFluid(fluidRegions.size());
PtrList<volScalarField> amrIndicatorFluid(fluidRegions.size());
PtrList<volVectorField> qwFluid(fluidRegions.size());

PtrList<volScalarField> mDotFluid(fluidRegions.size());
PtrList<volScalarField> qDotFluid(fluidRegions.size());

PtrList<IOMRFZoneList> MRFfluid(fluidRegions.size());
PtrList<fv::options> fluidFvOptions(fluidRegions.size());

List<label> pRefCellFluid(fluidRegions.size());
List<scalar> pRefValueFluid(fluidRegions.size());
List<scalar> cumulativeContErrFluid(fluidRegions.size());
PtrList<dimensionedScalar> ghRefFluid(fluidRegions.size());
// List<bool> frozenFlowFluid(fluidRegions.size(),false);
List<bool> alphaRestartFluid(fluidRegions.size(),false);

PtrList<pisoControl> pisoFluid(fluidRegions.size());

PtrList <surfaceScalarField> brackFSFluid(fluidRegions.size());

const uniformDimensionedVectorField& g = meshObjects::gravity::New(runTime);

scalarField maxTFluid(fluidRegions.size());

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{
    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;

    Info<< "    Adding to p_rghFluid\n" << endl;
    p_rghFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p_rgh",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );


    Info<< "    Adding to UFluid\n" << endl;
    UFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to phiFluid\n" << endl;
    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(UFluid[i])
                & fluidRegions[i].Sf()
        )
    );

    Info<< "    Adding to hRefFluid\n" << endl;
    hRefFluid.set
    (
        i,
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "hRef",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar("hRef", dimLength, Zero) // uses name
        )
    );

    Info << "    Adding to mDot\n" << endl;
    mDotFluid.set
      (
       i,
       new volScalarField
       (
        IOobject
        (
         "mDot",
         runTime.timeName(),
         fluidRegions[i],
         IOobject::MUST_READ,
         IOobject::AUTO_WRITE
         ),
        fluidRegions[i]
        )
       );

    Info << "    Adding to qDot\n" << endl;
    qDotFluid.set
      (i,
       new volScalarField
       (
        IOobject
        (
         "qDot",
         runTime.timeName(),
         fluidRegions[i],
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
         ),
        fluidRegions[i],
        dimensionedScalar
        (
         "qDot",
         dimensionSet(1,-1,-3,0,0,0,0),
         0.0
         )
        )
       );

    dimensionedScalar ghRef
    (
        mag(g.value()) > SMALL
      ? g & (cmptMag(g.value())/mag(g.value()))*hRefFluid[i]
      : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
    );

    ghRefFluid.set(i, new dimensionedScalar(ghRef));

    Info<< "    Adding to ghFluid\n" << endl;
    ghFluid.set
    (
        i,
        new volScalarField
        (
            "gh",
            (g & fluidRegions[i].C()) - ghRef
        )
    );

    Info<< "    Adding to ghfFluid\n" << endl;
    ghfFluid.set
    (
        i,
        new surfaceScalarField
        (
            "ghf",
            (g & fluidRegions[i].Cf()) - ghRef
        )
    );

    mixtureFluid.set
    (
        i,
        new immiscibleIncompressibleTwoPhaseMixture
        (
            UFluid[i],
            phiFluid[i]
        )
    );


    Info<< "    Adding to turbulenceFluid\n" << endl;
    turbulenceFluid.set
    (
        i,
        incompressible::turbulenceModel::New
        (
            UFluid[i],
            phiFluid[i],
            mixtureFluid[i]
        ).ptr()
    );


    immiscibleIncompressibleTwoPhaseMixture& mixture(mixtureFluid[i]);
    volScalarField& alpha1(mixture.alpha1());
    volScalarField& alpha2(mixture.alpha2());


    //- This also builds the surface reconstruction scheme

    const dimensionedScalar& rho1 = mixture.rho1();
    const dimensionedScalar& rho2 = mixture.rho2();

    dictionary phase1 = mixture.subDict(mixture.phase1Name());
    dictionary phase2 = mixture.subDict(mixture.phase2Name());

    dimensionedScalar k1
    (
        "kappa",
        dimensionSet(1,1,-3,-1,0,0,0),
        readScalar(phase1.lookup("kappa"))
    );

    dimensionedScalar k2
    (
        "kappa",
        dimensionSet(1,1,-3,-1,0,0,0),
        readScalar(phase2.lookup("kappa"))
    );



    Info<< "    Adding to kappaFluid\n" << endl;
    kappaFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "kappa",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            k1*alpha1 + k2*alpha2
        )
    );

    Info<< "    Adding to TFluid\n" << endl;
    TFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to amrIndicatorFluid\n" << endl;
    amrIndicatorFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "amrIndicator",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar("amrIndicator", dimless, 0.0)
        )
    );

    maxTFluid[i] = gMax(TFluid[i]);

    rhoFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            alpha1*rho1 + alpha2*rho2
        )
    );

    pisoFluid.set
    (
        i,
        new pisoControl(fluidRegions[i])
    );

    // pimpleFluid[i].dict().readIfPresent("frozenFlow", frozenFlowFluid[i]);

    pFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            p_rghFluid[i] + rhoFluid[i]*ghFluid[i]
        )
    );

    pRefCellFluid[i] = -1;
    pRefValueFluid[i] = 0.0;
    setRefCell
    (
        pFluid[i],
        p_rghFluid[i],
        pisoFluid[i].dict(),
        pRefCellFluid[i],
        pRefValueFluid[i]
    );



    if (p_rghFluid[i].needReference())
    {
        pFluid[i] += dimensionedScalar
        (
            "p",
            pFluid[i].dimensions(),
            pRefValueFluid[i] - getRefCellValue(pFluid[i], pRefCellFluid[i])
        );
        p_rghFluid[i] = pFluid[i] - rhoFluid[i]*ghFluid[i];
    }

    fluidRegions[i].setFluxRequired(p_rghFluid[i].name());
    fluidRegions[i].setFluxRequired(alpha1.name());


    Info<< "    Adding fvOptions\n" << endl;
    fluidFvOptions.set
    (
        i,
        new fv::options(fluidRegions[i])
    );

    turbulenceFluid[i].validate();

    Info<< "    Adding MRF\n" << endl;
    MRFfluid.set
    (
        i,
        new IOMRFZoneList(fluidRegions[i])
    );

    IOobject alphaPhi10Header
    (
        IOobject::groupName("alphaPhi0", alpha1.group()),
        runTime.timeName(),
        fluidRegions[i],
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    );

    alphaRestartFluid[i] =
        alphaPhi10Header.typeHeaderOk<surfaceScalarField>(true);

    if (alphaRestartFluid[i])
    {
        Info << "Restarting alpha" << endl;
    }

    alphaPhi10Fluid.set
    (
        i,
        new surfaceScalarField
        (
            alphaPhi10Header,
            phiFluid[i]*fvc::interpolate(alpha1)
        )
    );

    brackFSFluid.set
      (
       i,
       new surfaceScalarField
       (
        IOobject
        (
         "brackFS",
         runTime.timeName(),
         fluidRegions[i],
         IOobject::NO_READ,
         IOobject::NO_WRITE
        ),
        fvc::interpolate(mixture.sigmaK())*fvc::snGrad(alpha1)
        )
       );

    uniformDimensionedScalarField cumulativeContErrIO
    (
        IOobject
        (
            "cumulativeContErr",
            runTime.timeName(),
            "uniform",
            fluidRegions[i],
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        dimensionedScalar(dimless, Zero)
    );


    cumulativeContErrFluid[i] = cumulativeContErrIO.value();


    if(mixture.found("writeQwall"))
    {
        qwFluid.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "qw",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                -kappaFluid[i]*fvc::grad(TFluid[i])
            )
        );
    }

}
