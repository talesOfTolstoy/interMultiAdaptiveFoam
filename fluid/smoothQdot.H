{
    mDot = qDot / hlv;

    Info << "Max mDot = " << max(mDot).value() << endl;

    fvScalarMatrix mDotEqn
        (
         fvm::laplacian(diff, mDot)
         - fvm::Sp(scalar(1), mDot)
         ==
         - qDot / hlv
        );

    mDotEqn.solve();

    volScalarField epsilonV = 0.001 - alpha1;
    volScalarField vapourCutoff = pos(epsilonV);

    volScalarField epsilonL = alpha1 - 0.999;
    volScalarField liquidCutoff = pos(epsilonL);

    dimensionedScalar NNum = fvc::domainIntegrate(mDot);

    volScalarField vapourVolGen = mDot * epsilonV * vapourCutoff;
    dimensionedScalar NVDen = fvc::domainIntegrate(vapourVolGen);

    volScalarField liquidVolGen = mDot * epsilonL * liquidCutoff;
    dimensionedScalar NLDen = fvc::domainIntegrate(liquidVolGen);

    dimensionedScalar NL("Nl", dimensionSet(0,0,0,0,0,0,0), 2.0);
    dimensionedScalar NV("Nv", dimensionSet(0,0,0,0,0,0,0), 2.0);

    if(NLDen.value() != 0) { NL = NNum / NLDen; }
    if(NVDen.value() != 0) { NV = NNum / NVDen; }

    Info << "NV = " << NV.value() << ", NL = " << NL.value() << endl;

    vapourVolGen = NV * vapourCutoff * mDot * epsilonV;
    liquidVolGen = - NL * liquidCutoff * mDot * epsilonL;

    mDot = vapourVolGen + liquidVolGen;

    const scalar mEvap = gSum( qDot/hlv * mesh.V());


    Info<< "Timeev= " << runTime.timeName() <<" " <<
        "mEvap= " << mEvap  <<" " <<
        nl << endl;
}

