/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    scalarTransportFoam

Description
    Solves the steady or transient transport equation for a passive scalar.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"
    #include "inhibitions.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {



            fvScalarMatrix TEqn
            (
                fvm::ddt(T)
              + fvm::div(phi, T)
              - fvm::laplacian(DT, T)
             ==
                fvOptions(T)
            );



            fvScalarMatrix HplusEqn
            (
                fvm::ddt(Hplus)
              + fvm::div(phi, Hplus)
              - fvm::laplacian(DT, Hplus)
             ==
                fvOptions(Hplus)
/*              - (amm/17)/(1+k_a6/Hplus) * timeUnit
              + k_w/Hplus * timeUnit
              + (CO2/44)/(1+Hplus/k_a1+k_a2/Hplus) * timeUnit
              + 2*(CO2/44)/(1+Hplus/k_a2+pow(Hplus,2)/(k_a1*k_a2)) * timeUnit
              + (psi8/17)/(1+Hplus/k_a3) * timeUnit
              + (psi777/74)/(1+Hplus/k_a4) * timeUnit
              + (psi77/88)/(1+Hplus/k_a5) * timeUnit
              - AminusCplus * timeUnit
*/           );


            fvScalarMatrix ammEqn
            (
                fvm::ddt(amm)
              + fvm::div(phi, amm)
              - fvm::laplacian(DT, amm)
             ==
                fvOptions(amm)
            );


            fvScalarMatrix CO2Eqn
            (
                fvm::ddt(CO2)
              + fvm::div(phi, CO2)
              - fvm::laplacian(DT, CO2)
             ==
                fvOptions(CO2)
              + 2.413 * muA * (X4 + X5)
              + 1.01 * muAP * X777
              + 3.303 * muAB * X77
              + 16.726 * muM * X8
             );


            fvScalarMatrix H2Eqn
            (
                fvm::ddt(H2)
              + fvm::div(phi, H2)
              - fvm::laplacian(DT, H2)
             ==
                fvOptions(H2)
              + (1-Y_Su) * f_H2Su * k_mSu * (psi4/(psi4+K_SSu)) * X4   * I_pH * I_IN
              + (1-Y_Aa) * f_H2Aa * k_mAa * (psi5/(psi5+K_SAa)) * X5   * I_pH * I_IN
              + (1-Y_Fa) * f_H2Fa * k_mFa * (psi6/(psi6+K_SFa)) * X6  * I_pH * I_IN * I_H2
              + (1-Y_C4) * f_H2Va * k_mC4 * (psi7/(psi7+K_SC4)) * X7  * (psi7/(psi7+psi77))  * I_pH * I_IN * I_H2
              + (1-Y_C4) * f_H2Bu * k_mC4 * (psi77/(psi77+K_SC4)) * X77 * (psi77/(psi7+psi77))  * I_pH * I_IN * I_H2
              + (1-Y_Pro) * f_H2Pro * k_mPro * (psi777/(psi777+K_SPro)) * X777  * I_pH * I_IN * I_H2
              - k_mH2 * (H2/(H2+K_SH2)) * XH2 * I_pH * I_IN

            );



            fvScalarMatrix psi4Eqn
            (
                fvm::ddt(psi4)
              + fvm::div(phi, psi4)
              - fvm::laplacian(DT4, psi4)
             ==
                fvOptions(psi4)
              + k_HydCh * X1
              + (1-f_FaLi) * k_HydLi * X3
              - k_mSu * (psi4/(psi4+K_SSu)) * X4  * I_pH * I_IN

            );

            fvScalarMatrix psi5Eqn
            (
                fvm::ddt(psi5)
              + fvm::div(phi, psi5)
              - fvm::laplacian(DT5, psi5)
             ==
                fvOptions(psi5)
              + k_HydPr * X2
              - k_mAa * (psi5/(psi5+K_SAa)) * X5 * I_pH * I_IN
            );

            fvScalarMatrix psi6Eqn
            (
                fvm::ddt(psi6)
              + fvm::div(phi, psi6)
              - fvm::laplacian(DT6, psi6)
             ==
                fvOptions(psi6)
              + f_FaLi * k_HydLi * X3
              - k_mFa * (psi6/(psi6+K_SFa)) * X6 * I_pH * I_IN * I_H2

            );

            fvScalarMatrix psi7Eqn
            (
                fvm::ddt(psi7)
              + fvm::div(phi, psi7)
              - fvm::laplacian(DT7, psi7)
             ==
                fvOptions(psi7)
              + (1-Y_Aa) * f_VaAa * k_mAa * (psi5/(psi5+K_SAa)) * X5  * I_pH * I_IN
              - k_mC4 * (psi7/(psi7+K_SC4)) * X7 * (psi7/(psi7+psi77))  * I_pH * I_IN * I_H2

            );


            fvScalarMatrix psi77Eqn
            (
                fvm::ddt(psi77)
              + fvm::div(phi, psi77)
              - fvm::laplacian(DT7, psi77)
             ==
                fvOptions(psi77)
              + (1-Y_Su) * f_BuSu * k_mSu * (psi4/(psi4+K_SSu)) * X4  * I_pH * I_IN
              + (1-Y_Aa) * f_BuAa * k_mAa * (psi5/(psi5+K_SAa)) * X5  * I_pH * I_IN
              - k_mC4 * (psi77/(psi77+K_SC4)) * X77 * (psi77/(psi7+psi77))  * I_pH * I_IN * I_H2

            );

            fvScalarMatrix psi777Eqn
            (
                fvm::ddt(psi777)
              + fvm::div(phi, psi777)
              - fvm::laplacian(DT7, psi777)
             ==
                fvOptions(psi777)
              + (1-Y_Su) * f_ProSu * k_mSu * (psi4/(psi4+K_SSu)) * X4  * I_pH * I_IN
              + (1-Y_Aa) * f_ProAa * k_mAa * (psi5/(psi5+K_SAa)) * X5  * I_pH * I_IN
              + (1-Y_C4) * f_ProVa * k_mC4 * (psi7/(psi7+K_SC4)) * X7 * (psi7/(psi7+psi77))  * I_pH * I_IN * I_H2
              - k_mPro * (psi777/(psi777+K_SPro)) * X777  * I_pH * I_IN * I_H2

            );


            fvScalarMatrix psi8Eqn
            (
                fvm::ddt(psi8)
              + fvm::div(phi, psi8)
              - fvm::laplacian(DT8, psi8)
             ==
                fvOptions(psi8)
              + (1-Y_Su) * f_AcSu * k_mSu * (psi4/(psi4+K_SSu)) * X4  * I_pH * I_IN
              + (1-Y_Aa) * f_AcAa * k_mAa * (psi5/(psi5+K_SAa)) * X5  * I_pH * I_IN
              + (1-Y_Fa) * f_AcFa * k_mFa * (psi6/(psi6+K_SFa)) * X6  * I_pH * I_IN * I_H2
              + (1-Y_C4) * f_AcVa * k_mC4 * (psi7/(psi7+K_SC4)) * X7  * (psi7/(psi7+psi77)) * I_pH * I_IN * I_H2
              + (1-Y_C4) * f_AcBu * k_mC4 * (psi77/(psi77+K_SC4)) * X77 * (psi77/(psi7+psi77))  * I_pH * I_IN * I_H2
              + (1-Y_Pro) * f_AcPro * k_mPro * (psi777/(psi777+K_SPro)) * X777  * I_pH * I_IN * I_H2
              - k_mAc * (psi8/(psi8+K_SAc)) * X8 * I_pH2 * I_IN * I_amm
            );
			
			
            fvScalarMatrix CH4Eqn
            (
                fvm::ddt(CH4)
              + fvm::div(phi, CH4)
              - fvm::laplacian(DT, CH4)
             ==
                fvOptions(CH4)
              + (1-Y_Ac) * k_mAc * (psi8/(psi8+K_SAc)) * X8 * I_pH2 * I_IN * I_amm
              + (1-Y_H2) * k_mH2 * (H2/(H2+K_SH2)) * XH2 * I_pH * I_IN
            );


            fvScalarMatrix psiINEqn
            (
                fvm::ddt(psiIN)
              + fvm::div(phi, psiIN)
              - fvm::laplacian(DTX, psiIN)
             ==
                fvOptions(psiIN)
			  - Y_Su * 0.00625 * k_mSu * (psi4/(psi4+K_SSu)) * X4  * I_pH * I_IN
			  + (0.007-Y_Aa*0.00625) * k_mAa * (psi5/(psi5+K_SAa)) * X5  * I_pH * I_IN
			  - Y_Fa * 0.00625 * k_mFa * (psi6/(psi6+K_SFa)) * X6  * I_pH * I_IN * I_H2
			  - Y_C4 * 0.00625 * k_mC4 * (psi7/(psi7+K_SC4)) * X7  * (psi7/(psi7+psi77)) * I_pH * I_IN * I_H2
			  - Y_C4 * 0.00625 * k_mC4 * (psi77/(psi77+K_SC4)) * X77 * (psi77/(psi7+psi77))  * I_pH * I_IN * I_H2
			  - Y_Pro * 0.00625 * k_mPro * (psi777/(psi777+K_SPro)) * X777  * I_pH * I_IN * I_H2
			  - Y_Ac * 0.00625 * k_mAc * (psi8/(psi8+K_SAc)) * X8 * I_pH2 * I_IN * I_amm
			  - Y_H2 * 0.00625 * k_mH2 * (H2/(H2+K_SH2)) * XH2 * I_pH * I_IN
            );


            fvScalarMatrix XcEqn 
            (
                fvm::ddt(Xc)
              + fvm::div(phi, Xc)
              - fvm::laplacian(DTX, Xc)
             ==
                fvOptions(Xc)
              - k_dis * Xc
              + k_decXSu * X4
              + k_decXAa * X5
              + k_decXFa * X6
              + k_decXC4 * X7
              + k_decXC4 * X77
              + k_decXPro * X777
              + k_decXAc * X8
              + k_decXAc * XH2
            );


            fvScalarMatrix X1Eqn
            (
                fvm::ddt(X1)
              + fvm::div(phi, X1)
              - fvm::laplacian(DTX, X1)
             ==
                fvOptions(X1)
              - k_HydCh * X1
              + f_ChXc * k_dis * Xc
            );

            fvScalarMatrix X2Eqn
            (
                fvm::ddt(X2)
              + fvm::div(phi, X2)
              - fvm::laplacian(DTX, X2)
             ==
                fvOptions(X2)
              - k_HydPr * X2
              + f_PrXc * k_dis * Xc
            );

            fvScalarMatrix X3Eqn
            (
                fvm::ddt(X3)
              + fvm::div(phi, X3)
              - fvm::laplacian(DTX, X3)
             ==
                fvOptions(X3)
              - k_HydLi * X3
              + f_LiXc * k_dis * Xc
            );


            fvScalarMatrix X4Eqn
            (
                fvm::ddt(X4)
              + fvm::div(phi, X4)
              - fvm::laplacian(DTX, X4)
             ==
                fvOptions(X4)
              + Y_Su * k_mSu * (psi4/(psi4+K_SSu)) * X4  * I_pH * I_IN
              - k_decXSu * X4
            );

            fvScalarMatrix X5Eqn
            (
                fvm::ddt(X5)
              + fvm::div(phi, X5)
              - fvm::laplacian(DTX, X5)
             ==
                fvOptions(X5)
              + Y_Aa * k_mAa * (psi5/(psi5+K_SAa)) * X5  * I_pH * I_IN
              - k_decXAa * X5
            );

            fvScalarMatrix X6Eqn
            (
                fvm::ddt(X6)
              + fvm::div(phi, X6)
              - fvm::laplacian(DTX, X6)
             ==
                fvOptions(X6)
              + Y_Fa * k_mFa * (psi6/(psi6+K_SFa)) * X6  * I_pH * I_IN
              - k_decXFa * X6
            );

            fvScalarMatrix X7Eqn
            (
                fvm::ddt(X7)
              + fvm::div(phi, X7)
              - fvm::laplacian(DTX, X7)
             ==
                fvOptions(X7)
              + Y_C4 * k_mC4 * (psi7/(psi7+K_SC4)) * X7 * (psi7/(psi7+psi77))  * I_pH * I_IN * I_H2
              - k_decXC4 * X7
            );

            fvScalarMatrix X8Eqn
            (
                fvm::ddt(X8)
              + fvm::div(phi, X8)
              - fvm::laplacian(DTX, X8)
             ==
                fvOptions(X8)
              + Y_Ac * k_mAc * (psi8/(psi8+K_SAc)) * X8 * I_pH2 * I_IN * I_amm
              - k_decXAc * X8

            );

            fvScalarMatrix X77Eqn
            (
                fvm::ddt(X77)
              + fvm::div(phi, X77)
              - fvm::laplacian(DTX, X77)
             ==
                fvOptions(X77)
              + Y_C4 * k_mC4 * (psi77/(psi77+K_SC4)) * X77 * (psi77/(psi7+psi77))  * I_pH * I_IN * I_H2
              - k_decXC4 * X77

            );

            fvScalarMatrix X777Eqn
            (
                fvm::ddt(X777)
              + fvm::div(phi, X777)
              - fvm::laplacian(DTX, X777)
             ==
                fvOptions(X777)
              + Y_Pro * k_mPro * (psi777/(psi777+K_SPro)) * X777 * I_pH * I_IN * I_H2
              - k_decXPro * X777
            );

            fvScalarMatrix XH2Eqn
            (
                fvm::ddt(XH2)
              + fvm::div(phi, XH2)
              - fvm::laplacian(DTX, XH2)
             ==
                fvOptions(XH2)
              + Y_H2 * k_mH2 * (H2/(H2+K_SH2)) * XH2 * I_pH * I_IN
              - k_decXPro * XH2
            );





            TEqn.relax();
            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(T);


            HplusEqn.relax();
            fvOptions.constrain(HplusEqn);
            HplusEqn.solve();
            fvOptions.correct(Hplus);


            ammEqn.relax();
            fvOptions.constrain(ammEqn);
            ammEqn.solve();
            fvOptions.correct(amm);

            CO2Eqn.relax();
            fvOptions.constrain(CO2Eqn);
            CO2Eqn.solve();
            fvOptions.correct(CO2);


            H2Eqn.relax();
            fvOptions.constrain(H2Eqn);
            H2Eqn.solve();
            fvOptions.correct(H2);

            psi4Eqn.relax();
            fvOptions.constrain(psi4Eqn);
            psi4Eqn.solve();
            fvOptions.correct(psi4);

            psi5Eqn.relax();
            fvOptions.constrain(psi5Eqn);
            psi5Eqn.solve();
            fvOptions.correct(psi5);

            psi6Eqn.relax();
            fvOptions.constrain(psi6Eqn);
            psi6Eqn.solve();
            fvOptions.correct(psi6);

            psi7Eqn.relax();
            fvOptions.constrain(psi7Eqn);
            psi7Eqn.solve();
            fvOptions.correct(psi7);

            psi8Eqn.relax();
            fvOptions.constrain(psi8Eqn);
            psi8Eqn.solve();
            fvOptions.correct(psi8);
			
            CH4Eqn.relax();
            fvOptions.constrain(CH4Eqn);
            CH4Eqn.solve();
            fvOptions.correct(CH4);

            psi77Eqn.relax();
            fvOptions.constrain(psi77Eqn);
            psi77Eqn.solve();
            fvOptions.correct(psi77);

            psi777Eqn.relax();
            fvOptions.constrain(psi777Eqn);
            psi777Eqn.solve();
            fvOptions.correct(psi777);

            psiINEqn.relax();
            fvOptions.constrain(psiINEqn);
            psiINEqn.solve();
            fvOptions.correct(psiIN);

            XcEqn.relax();
            fvOptions.constrain(XcEqn);
            XcEqn.solve();
            fvOptions.correct(Xc);

            X1Eqn.relax();
            fvOptions.constrain(X1Eqn);
            X1Eqn.solve();
            fvOptions.correct(X1);

            X2Eqn.relax();
            fvOptions.constrain(X2Eqn);
            X2Eqn.solve();
            fvOptions.correct(X2);

            X3Eqn.relax();
            fvOptions.constrain(X3Eqn);
            X3Eqn.solve();
            fvOptions.correct(X3);


            X4Eqn.relax();
            fvOptions.constrain(X4Eqn);
            X4Eqn.solve();
            fvOptions.correct(X4);

            X5Eqn.relax();
            fvOptions.constrain(X5Eqn);
            X5Eqn.solve();
            fvOptions.correct(X5);

            X6Eqn.relax();
            fvOptions.constrain(X6Eqn);
            X6Eqn.solve();
            fvOptions.correct(X6);

            X7Eqn.relax();
            fvOptions.constrain(X7Eqn);
            X7Eqn.solve();
            fvOptions.correct(X7);

            X8Eqn.relax();
            fvOptions.constrain(X8Eqn);
            X8Eqn.solve();
            fvOptions.correct(X8);

            X77Eqn.relax();
            fvOptions.constrain(X77Eqn);
            X77Eqn.solve();
            fvOptions.correct(X77);

            X777Eqn.relax();
            fvOptions.constrain(X777Eqn);
            X777Eqn.solve();
            fvOptions.correct(X777);

            XH2Eqn.relax();
            fvOptions.constrain(XH2Eqn);
            XH2Eqn.solve();
            fvOptions.correct(XH2);


        }

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //