

dimensionedScalar nu_ave=sum(U.mesh().V()*fluid.nu())/sum(U.mesh().V());

Info << "Average of nu: " << nu_ave << endl;



//scalar ReNum = 0.0;
//scalar meanReNum = 0.0;
//
//{
//    scalarField sumPhi
//    (
//        fvc::surfaceSum(mag(phi))().primitiveField()
//    );
//
//    CoNum = 0.5*gMax(mag(U)/fluid.nu());
//
//    meanCoNum =
//        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
//}
//
//Info<< "Courant Number mean: " << meanCoNum
//    << " max: " << CoNum << endl;


