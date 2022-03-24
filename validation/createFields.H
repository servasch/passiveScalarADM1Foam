Info<< "Reading field T\n" << endl;

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



volScalarField Hplus
(
    IOobject
    (
        "Hplus",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField pH
(
    IOobject
    (
        "pH",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField amm
(
    IOobject
    (
        "amm",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField CO2
(
    IOobject
    (
        "CO2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField H2
(
    IOobject
    (
        "H2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



volScalarField psi4
(
    IOobject
    (
        "psi4",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField psi5
(
    IOobject
    (
        "psi5",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField psi6
(
    IOobject
    (
        "psi6",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField psi7
(
    IOobject
    (
        "psi7",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField psi8
(
    IOobject
    (
        "psi8",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField CH4
(
    IOobject
    (
        "CH4",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField psi77
(
    IOobject
    (
        "psi77",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField psi777
(
    IOobject
    (
        "psi777",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField psiIN
(
    IOobject
    (
        "psiIN",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField Xc
(
    IOobject
    (
        "Xc",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField X1
(
    IOobject
    (
        "X1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField X2
(
    IOobject
    (
        "X2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField X3
(
    IOobject
    (
        "X3",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



volScalarField X4
(
    IOobject
    (
        "X4",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField X5
(
    IOobject
    (
        "X5",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField X6
(
    IOobject
    (
        "X6",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField X7
(
    IOobject
    (
        "X7",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField X8
(
    IOobject
    (
        "X8",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField X77
(
    IOobject
    (
        "X77",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField X777
(
    IOobject
    (
        "X777",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField XH2
(
    IOobject
    (
        "XH2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



Info<< "Reading field U\n" << endl;

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


Info<< "Reading diffusivity DT\n" << endl;
dimensionedScalar DT
(
    transportProperties.lookup("DT")
);


Info<< "Reading diffusivity DT4\n" << endl;
dimensionedScalar DT4
(
    transportProperties.lookup("DT4")
);

Info<< "Reading diffusivity DT5\n" << endl;
dimensionedScalar DT5
(
    transportProperties.lookup("DT5")
);

Info<< "Reading diffusivity DT6\n" << endl;
dimensionedScalar DT6
(
    transportProperties.lookup("DT6")
);

Info<< "Reading diffusivity DT7\n" << endl;
dimensionedScalar DT7
(
    transportProperties.lookup("DT7")
);

Info<< "Reading diffusivity DTX\n" << endl;
dimensionedScalar DTX
(
    transportProperties.lookup("DTX")
);

Info<< "Reading diffusivity DT8\n" << endl;
dimensionedScalar DT8
(
    transportProperties.lookup("DT8")
);




Info<< "Reading rate k_decXSu\n" << endl;
dimensionedScalar k_decXSu
(
    transportProperties.lookup("k_decXSu")
);

Info<< "Reading rate k_decXAa\n" << endl;
dimensionedScalar k_decXAa
(
    transportProperties.lookup("k_decXAa")
);

Info<< "Reading rate k_decXFa\n" << endl;
dimensionedScalar k_decXFa
(
    transportProperties.lookup("k_decXFa")
);

Info<< "Reading rate k_decXC4\n" << endl;
dimensionedScalar k_decXC4
(
    transportProperties.lookup("k_decXC4")
);

Info<< "Reading rate k_decXPro\n" << endl;
dimensionedScalar k_decXPro
(
    transportProperties.lookup("k_decXPro")
);

Info<< "Reading rate k_decXAc\n" << endl;
dimensionedScalar k_decXAc
(
    transportProperties.lookup("k_decXAc")
);

Info<< "Reading rate k_decXH2\n" << endl;
dimensionedScalar k_decXH2
(
    transportProperties.lookup("k_decXH2")
);



Info<< "Reading rate f_ChXc\n" << endl;
dimensionedScalar f_ChXc
(
    transportProperties.lookup("f_ChXc")
);

Info<< "Reading rate f_PrXc\n" << endl;
dimensionedScalar f_PrXc
(
    transportProperties.lookup("f_PrXc")
);

Info<< "Reading rate f_LiXc\n" << endl;
dimensionedScalar f_LiXc
(
    transportProperties.lookup("f_LiXc")
);


Info<< "Reading rate f_FaLi\n" << endl;
dimensionedScalar f_FaLi
(
    transportProperties.lookup("f_FaLi")
);


Info<< "Reading rate f_BuSu\n" << endl;
dimensionedScalar f_BuSu
(
    transportProperties.lookup("f_BuSu")
);


Info<< "Reading rate f_ProSu\n" << endl;
dimensionedScalar f_ProSu
(
    transportProperties.lookup("f_ProSu")
);

Info<< "Reading rate f_ProAa\n" << endl;
dimensionedScalar f_ProAa
(
    transportProperties.lookup("f_ProAa")
);

Info<< "Reading rate f_ProVa\n" << endl;
dimensionedScalar f_ProVa
(
    transportProperties.lookup("f_ProVa")
);


Info<< "Reading rate f_AcSu\n" << endl;
dimensionedScalar f_AcSu
(
    transportProperties.lookup("f_AcSu")
);

Info<< "Reading rate f_AcAa\n" << endl;
dimensionedScalar f_AcAa
(
    transportProperties.lookup("f_AcAa")
);

Info<< "Reading rate f_AcFa\n" << endl;
dimensionedScalar f_AcFa
(
    transportProperties.lookup("f_AcFa")
);

Info<< "Reading rate f_AcVa\n" << endl;
dimensionedScalar f_AcVa
(
    transportProperties.lookup("f_AcVa")
);

Info<< "Reading rate f_AcBu\n" << endl;
dimensionedScalar f_AcBu
(
    transportProperties.lookup("f_AcBu")
);

Info<< "Reading rate f_AcPro\n" << endl;
dimensionedScalar f_AcPro
(
    transportProperties.lookup("f_AcPro")
);


Info<< "Reading rate f_VaAa\n" << endl;
dimensionedScalar f_VaAa
(
    transportProperties.lookup("f_VaAa")
);

Info<< "Reading rate f_BuAa\n" << endl;
dimensionedScalar f_BuAa
(
    transportProperties.lookup("f_BuAa")
);

Info<< "Reading rate f_H2Su\n" << endl;
dimensionedScalar f_H2Su
(
    transportProperties.lookup("f_H2Su")
);


Info<< "Reading rate f_H2Aa\n" << endl;
dimensionedScalar f_H2Aa
(
    transportProperties.lookup("f_H2Aa")
);


Info<< "Reading rate f_H2Fa\n" << endl;
dimensionedScalar f_H2Fa
(
    transportProperties.lookup("f_H2Fa")
);

Info<< "Reading rate f_H2Va\n" << endl;
dimensionedScalar f_H2Va
(
    transportProperties.lookup("f_H2Va")
);

Info<< "Reading rate f_H2Bu\n" << endl;
dimensionedScalar f_H2Bu
(
    transportProperties.lookup("f_H2Bu")
);

Info<< "Reading rate f_H2Pro\n" << endl;
dimensionedScalar f_H2Pro
(
    transportProperties.lookup("f_H2Pro")
);




Info<< "Reading rate Y_Su\n" << endl;
dimensionedScalar Y_Su
(
    transportProperties.lookup("Y_Su")
);

Info<< "Reading rate Y_Aa\n" << endl;
dimensionedScalar Y_Aa
(
    transportProperties.lookup("Y_Aa")
);

Info<< "Reading rate Y_Fa\n" << endl;
dimensionedScalar Y_Fa
(
    transportProperties.lookup("Y_Fa")
);

Info<< "Reading rate Y_C4\n" << endl;
dimensionedScalar Y_C4
(
    transportProperties.lookup("Y_C4")
);

Info<< "Reading rate Y_Pro\n" << endl;
dimensionedScalar Y_Pro
(
    transportProperties.lookup("Y_Pro")
);

Info<< "Reading rate Y_Ac\n" << endl;
dimensionedScalar Y_Ac
(
    transportProperties.lookup("Y_Ac")
);

Info<< "Reading rate Y_H2\n" << endl;
dimensionedScalar Y_H2
(
    transportProperties.lookup("Y_H2")
);

Info<< "Reading rate k_HydCh\n" << endl;
dimensionedScalar k_HydCh
(
    transportProperties.lookup("k_HydCh")
);

Info<< "Reading rate k_HydPr\n" << endl;
dimensionedScalar k_HydPr
(
    transportProperties.lookup("k_HydPr")
);

Info<< "Reading rate k_HydLi\n" << endl;
dimensionedScalar k_HydLi
(
    transportProperties.lookup("k_HydLi")
);

Info<< "Reading rate k_mSu\n" << endl;
dimensionedScalar k_mSu
(
    transportProperties.lookup("k_mSu")
);

Info<< "Reading rate k_mAa\n" << endl;
dimensionedScalar k_mAa
(
    transportProperties.lookup("k_mAa")
);

Info<< "Reading rate k_mFa\n" << endl;
dimensionedScalar k_mFa
(
    transportProperties.lookup("k_mFa")
);

Info<< "Reading rate k_mC4\n" << endl;
dimensionedScalar k_mC4
(
    transportProperties.lookup("k_mC4")
);

Info<< "Reading rate k_mPro\n" << endl;
dimensionedScalar k_mPro
(
    transportProperties.lookup("k_mPro")
);

Info<< "Reading rate k_mAc\n" << endl;
dimensionedScalar k_mAc
(
    transportProperties.lookup("k_mAc")
);

Info<< "Reading rate k_mH2\n" << endl;
dimensionedScalar k_mH2
(
    transportProperties.lookup("k_mH2")
);


Info<< "Reading rate k_dis\n" << endl;
dimensionedScalar k_dis
(
    transportProperties.lookup("k_dis")
);


Info<< "Reading rate kS\n" << endl;
dimensionedScalar kS
(
    transportProperties.lookup("kS")
);


Info<< "Reading rate K_SSu\n" << endl;
dimensionedScalar K_SSu
(
    transportProperties.lookup("K_SSu")
);

Info<< "Reading rate K_SAa\n" << endl;
dimensionedScalar K_SAa
(
    transportProperties.lookup("K_SAa")
);

Info<< "Reading rate K_SFa\n" << endl;
dimensionedScalar K_SFa
(
    transportProperties.lookup("K_SFa")
);

Info<< "Reading rate K_SC4\n" << endl;
dimensionedScalar K_SC4
(
    transportProperties.lookup("K_SC4")
);

Info<< "Reading rate K_SPro\n" << endl;
dimensionedScalar K_SPro
(
    transportProperties.lookup("K_SPro")
);

Info<< "Reading rate K_SAc\n" << endl;
dimensionedScalar K_SAc
(
    transportProperties.lookup("K_SAc")
);

Info<< "Reading rate K_SH2\n" << endl;
dimensionedScalar K_SH2
(
    transportProperties.lookup("K_SH2")
);







#include "createPhi.H"

#include "createFvOptions.H"
