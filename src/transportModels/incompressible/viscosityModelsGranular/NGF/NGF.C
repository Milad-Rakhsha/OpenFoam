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

#include "addToRunTimeSelectionTable.H"
#include "NGF.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModelsGranular
{
    defineTypeNameAndDebug(NGF, 0);
    addToRunTimeSelectionTable
    (
    	viscosityModelGranular,
		NGF,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModelsGranular::NGF::calc_g()
{

	volScalarField  P_P ( max ( p_ , dimensionedScalar ("pSmall", dimLength*dimLength/dimTime/dimTime, 1e-25)) );
	volScalarField  I ( max ( strainRate() * diam_ / Foam::sqrt( P_P / rho_s_ ), 1e-20));
//	volScalarField  mu ( mu_s_+(mu_2_-mu_s_)*I / (I0_+ I) );
	volScalarField  mu ( mu_s_ + b_ * I  );

	volScalarField  zeta (
							A_ / Foam::sqrt(Foam::mag(mu-mu_s_)) * diam_
						  );

	volScalarField g_loc ( strainRate()/mu );

    tmp<fvScalarMatrix> gEqn
					(
							zeta*zeta*fvm::laplacian(g_)  == fvm::Sp(1.0 , g_)-g_loc
					);

    //	const Foam::fvMesh& mesh = U_.mesh();
    //    fv::options& fvOptions(fv::options::New(mesh));
    //	 fvOptions.constrain(gEqn.ref());
    //	 gEqn.ref().relax();
    solve(gEqn);
    //	 fvOptions.correct(g_);

    return tmp<volScalarField> (g_);


}

Foam::tmp<Foam::volScalarField>
Foam::viscosityModelsGranular::NGF::calcNu()
{
	volScalarField  P_P ( max ( p_ , dimensionedScalar ("pSmall", dimLength*dimLength/dimTime/dimTime, 1e-25)) );
    return tmp<volScalarField> ( 2*P_P/NGF::calc_g()() );
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModelsGranular::NGF::NGF
(
        const word& name,
        const dictionary& viscosityProperties,
        const volVectorField& U,
        const surfaceScalarField& phi,
		const volScalarField& p

):
	viscosityModelGranular(name, viscosityProperties, U, phi, p),
	NGFCoeffs_
    (
        viscosityProperties.optionalSubDict(typeName + "Coeffs")
    ),
	diam_("par_diameter", dimLength, NGFCoeffs_),
	rho_s_("par_density", dimless, NGFCoeffs_),
	I0_("I0", dimless, NGFCoeffs_),
	mu_s_("par_friction", dimless, NGFCoeffs_),
	mu_2_("par_mu2", dimless, NGFCoeffs_),
	b_("b", dimless, NGFCoeffs_),
	A_("A", dimless, NGFCoeffs_),
    nu0_("nu_0", dimViscosity, NGFCoeffs_),

    nu_
    (
        IOobject
        (
        	"nu_fluid",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        nu0_
    ),

	g_
	(
		IOobject
		(
			"g",
            U_.time().timeName(),
            U_.db(),
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
        U_.mesh()
	)


{

	Info << "created Constitutive equation based on Non-local Granular Fluidity" << endl;
	Info << "particle diameter: " << diam_ << endl;
	Info << "particle density: " << rho_s_ << endl;
	Info << "I0: " << I0_ << endl;
	Info << "particle friction: " << mu_s_ << endl;
	Info << "particle mu2: " << mu_2_ << endl;
	Info << "nu: " << nu0_ << endl;
	Info << "b:" << b_ << endl;
	Info << "A:" << A_ << endl;


}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModelsGranular::NGF::read
(
    const dictionary& viscosityProperties
)
{
	viscosityModelGranular::read(viscosityProperties);

	NGFCoeffs_ =
        viscosityProperties.optionalSubDict(typeName + "Coeffs");

    Info <<  "reading from " << typeName << "Coeffs" << endl;
    NGFCoeffs_.lookup("par_diameter") >> diam_;
    NGFCoeffs_.lookup("par_density") >> rho_s_;
    NGFCoeffs_.lookup("I0") >> I0_;
    NGFCoeffs_.lookup("par_friction") >> mu_s_;
    NGFCoeffs_.lookup("par_mu2") >> mu_2_;
    NGFCoeffs_.lookup("b") >> b_;
    NGFCoeffs_.lookup("A") >> A_;
    NGFCoeffs_.lookup("nu_0") >> nu0_;



    return true;
}


// ************************************************************************* //
