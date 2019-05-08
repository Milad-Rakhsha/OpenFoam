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

\*---------------------------------------------------------------------------*/

#include "Jop.H"

#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModelsGranular
{
    defineTypeNameAndDebug(Jop, 0);
    addToRunTimeSelectionTable
    (
    	viscosityModelGranular,
        Jop,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModelsGranular::Jop::calcNu() const
{


	volScalarField  P_P ( max ( p_ , dimensionedScalar ("pSmall", dimLength*dimLength/dimTime/dimTime, 1e-10)) );
//	volScalarField  I ( strainRate() * diam_ / pow( P_P / rho_s_, 0.5 ));
//	volScalarField  mu ( mu_s_+(mu_2_-mu_s_)*I / (I0_+ I) );

//	volScalarField mu (
//						IOobject
//			            (
//			            	"nu_fluid",
//			                U_.time().timeName(),
//			                U_.db(),
//			                IOobject::NO_READ,
//			                IOobject::NO_WRITE
//			            ),
//			            U_.mesh(),
//						dimensionedScalar ("nu_small", dimless, 1e-2)
//
//	);

		//    return  tmp<volScalarField> (mu *P_P/sr);

//	volScalarField nu_max (
//						IOobject
//			            (
//			            	"nu_max",
//			                U_.time().timeName(),
//			                U_.db(),
//			                IOobject::NO_READ,
//			                IOobject::NO_WRITE
//			            ),
//			            U_.mesh(),
//						nu0_
//	);

//	dimensionedScalar sr_small("vSmall", dimless/dimTime, 1e-6);
	volScalarField sr (max(strainRate(),reg_));
//	dimensionedScalar sr_small("vSmall", dimViscosity, 1e-6);

//	volScalarField k (max(nu_max- mu_s_*P_P /sr_small,nu_max*0));


//	return tmp<volScalarField> ((mu* P_P) / sr);

//	return tmp<volScalarField> (
//									min
//									(
//										nu_max,  //For very small shear rate we need a high nu value
//										(mu_s_* P_P) / sr
//
//									)
//								);


    dimensionedScalar tone("tone", dimTime, 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);
//	dimensionedScalar mu_plastic("mp", dimViscosity, mu_s_.value());
//	return tmp<volScalarField> ( mu_s_* P_P / rtone  + mu_s_* P_P / sr *(1-exp(-sr*tone*K)));
	return tmp<volScalarField> ( mu_s_* P_P / sr +  mu_s_* P_P / sr * (1-exp(-sr*tone*K)) );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModelsGranular::Jop::Jop
(
        const word& name,
        const dictionary& viscosityProperties,
        const volVectorField& U,
        const surfaceScalarField& phi,
		const volScalarField& p
):
	viscosityModelGranular(name, viscosityProperties, U, phi, p),
    JopCoeffs_
    (
        viscosityProperties.optionalSubDict(typeName + "Coeffs")
    ),
	diam_("par_diameter", dimLength, JopCoeffs_),
	rho_s_("par_density", dimless, JopCoeffs_),
	I0_("I0", dimless, JopCoeffs_),
	mu_s_("par_friction", dimless, JopCoeffs_),
	mu_2_("par_mu2", dimless, JopCoeffs_),
	reg_("regularization", dimless/dimTime, JopCoeffs_),
    nu0_("nu_max", dimViscosity, JopCoeffs_),
    K("K", dimless, JopCoeffs_),
    nu_
    (
        IOobject
        (
        	"nu_fluid",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_.mesh(),
        nu0_
    )

{

	Info << "created Constitutive equation based on Jop et al." << endl;
	Info << "particle diameter: " << diam_ << endl;
	Info << "particle density: " << rho_s_ << endl;
	Info << "particle friction: " << mu_s_ << endl;
	Info << "particle mu2: " << mu_2_ << endl;
	Info << "particle regularization: " << reg_ << endl;
	Info << "I0: " << I0_ << endl;
	Info << "nu: " << nu0_ << endl;
	Info << "K: " << K << endl;

}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModelsGranular::Jop::read
(
    const dictionary& viscosityProperties
)
{
	viscosityModelGranular::read(viscosityProperties);

    JopCoeffs_ =
        viscosityProperties.optionalSubDict(typeName + "Coeffs");

    Info <<  "reading from " << typeName << "Coeffs" << endl;
    JopCoeffs_.lookup("par_diameter") >> diam_;
    JopCoeffs_.lookup("par_density") >> rho_s_;
    JopCoeffs_.lookup("I0") >> I0_;
    JopCoeffs_.lookup("par_friction") >> mu_s_;
    JopCoeffs_.lookup("par_mu2") >> mu_2_;
    JopCoeffs_.lookup("regularization") >> reg_;
    JopCoeffs_.lookup("nu_max") >> nu0_;
    JopCoeffs_.lookup("K") >> K;


    return true;
}


// ************************************************************************* //
