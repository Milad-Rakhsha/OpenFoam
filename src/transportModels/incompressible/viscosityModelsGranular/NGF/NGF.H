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

Class
    Foam::viscosityModelsGranular::NGF

Description
    An incompressible NGF non-Newtonian viscosity model.

SourceFiles
    NGF.C

\*---------------------------------------------------------------------------*/

#ifndef NGF_H
#define NGF_H

#include "viscosityModelsGranular/viscosityModelGranular/viscosityModelGranular.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "fvCFD.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModelsGranular
{

/*---------------------------------------------------------------------------*\
                           Class NGF Declaration
\*---------------------------------------------------------------------------*/

class NGF
:
    public viscosityModelGranular
{
    // Private data

        dictionary NGFCoeffs_;

        dimensionedScalar diam_;
        dimensionedScalar rho_s_;
        dimensionedScalar I0_;
        dimensionedScalar mu_s_;
        dimensionedScalar mu_2_;
        dimensionedScalar b_;

        dimensionedScalar A_;
        dimensionedScalar nu0_;

        volScalarField nu_;
        volScalarField g_;


    // Private Member Functions

        //- Calculate and return the laminar viscosity
        tmp<volScalarField> calcNu() ;
        tmp<volScalarField> calc_g() ;


public:

    //- Runtime type information
    TypeName("NGF");


    // Constructors

        //- Construct from components
        NGF
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const volScalarField& p

        );


    //- Destructor
    virtual ~NGF()
    {}


    // Member Functions

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Correct the laminar viscosity
        virtual void correct()
        {
            nu_ = calcNu();

        }

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModelsGranular
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
