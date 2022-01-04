
Info<< "Reading inhibition characteristics\n" << endl;

//scalar pH_UL = 7.5;


dimensionedScalar pH_UL
(
    transportProperties.lookup("pH_UL")
);

dimensionedScalar pH_LL
(
    transportProperties.lookup("pH_LL")
);

dimensionedScalar pH_UL2
(
    transportProperties.lookup("pH_UL2")
);

dimensionedScalar pH_LL2
(
    transportProperties.lookup("pH_LL2")
);


dimensionedScalar kI_amm
(
    transportProperties.lookup("kI_amm")
);

dimensionedScalar kI_H2
(
    transportProperties.lookup("kI_H2")
);

dimensionedScalar kI_IN
(
    transportProperties.lookup("kI_IN")
);




Info<< "Reading pH characteristics\n" << endl;
dimensionedScalar k_a1
(
    transportProperties.lookup("k_a1")
);

dimensionedScalar k_a2
(
    transportProperties.lookup("k_a2")
);

dimensionedScalar k_a3
(
    transportProperties.lookup("k_a3")
);

dimensionedScalar k_a4
(
    transportProperties.lookup("k_a4")
);

dimensionedScalar k_a5
(
    transportProperties.lookup("k_a5")
);

dimensionedScalar k_a6
(
    transportProperties.lookup("k_a6")
);

dimensionedScalar k_w
(
    transportProperties.lookup("k_w")
);


dimensionedScalar AminusCplus
(
    transportProperties.lookup("AminusCplus")
);


dimensionedScalar Unit
(
    transportProperties.lookup("Unit")
);

dimensionedScalar timeUnit
(
    transportProperties.lookup("timeUnit")
);


dimensionedScalar muA
(
    transportProperties.lookup("muA")
);

dimensionedScalar muAP
(
    transportProperties.lookup("muAP")
);

dimensionedScalar muAB
(
    transportProperties.lookup("muAB")
);

dimensionedScalar muM
(
    transportProperties.lookup("muM")
);



Info<< "Reading pH field\n" << endl;

			
pH = 3 - (log10(Hplus/Unit));

Info<< "Reading inhibition fields\n" << endl;
volScalarField I_pH
(
(1+2*pow(10,(pH_LL-pH_UL)/pH))/(1+pow(10,2-pH_UL)+pow(10,pH_LL-2))
);






volScalarField I_pH2
(
exp(-3*pow(mag(pH-pH_UL2)/mag(pH_UL2-pH_LL2),2))
);

//pos(pH-pH_UL2)+limit(pH,pH_UL2)


forAll(pH, celli)
{
    if (pH[celli] > pH_UL2.value())
    {
    I_pH2[celli] = 1.;
    }
}
//printf("%f", pH[i]);
//if (mesh.V()[celli] > pH_UL2) {
//I_pH2[celli] = 1;
//}
//}


volScalarField I_amm
(
1/(1+amm/kI_amm)
);

volScalarField I_H2
(
1/(1+H2/kI_H2)
);

volScalarField I_IN
(
1/(1+kI_IN/psiIN)
);


Info<< "Reading test\n" << endl;






