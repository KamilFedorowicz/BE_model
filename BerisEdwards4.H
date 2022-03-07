

#ifndef BerisEdwards_H
#define BerisEdwards_H

#include "constitutiveEq2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs
{ 

/*---------------------------------------------------------------------------*\
                           Class BerisEdwards Declaration
\*---------------------------------------------------------------------------*/

class BerisEdwards
:
    public constitutiveEq2
{
    // Private data


        //- Transported nematic stress
        volTensorField tau_;
	volTensorField Q_;
	volTensorField sigma2_;
        //- Transported nematic stress
        volTensorField tauTotal_;

        // Model constants
            //- Density
            dimensionedScalar rho_;
            //- Newtnoian viscosity 
            dimensionedScalar etaS_;

            dimensionedScalar a_;

            dimensionedScalar b_;

            dimensionedScalar c_;

            dimensionedScalar L_;

            dimensionedScalar xi_;

            dimensionedScalar gamma_;
 
    // Private Member Functions

        //- Disallow default bitwise copy construct
        BerisEdwards(const BerisEdwards&);

        //- Disallow default bitwise assignment
        void operator=(const BerisEdwards&);
        
protected:

       //- Return the solvent viscosity
       virtual const dimensionedScalar etaS() const
       {
          return etaS_;
       }
      


public:

    //- Runtime type information
    TypeName("BerisEdwards");

    // Constructors

        //- Construct from components
        BerisEdwards
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~BerisEdwards()
        {}


    // Member Functions

       //- Return the viscoelastic stress tensor
        virtual tmp<volTensorField> tau() const
        {
            return tau_;
        }
        
        //- Return the density
        virtual const dimensionedScalar rho() const
        {
            return rho_;
        }
        
        //- Return true if GNF (non-elastic)
        virtual bool isGNF() const
        {
          return false;
        };

        //- Correct the viscoelastic stress
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace constitutiveEqs 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
