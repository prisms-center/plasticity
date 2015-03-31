// created: 2015-3-5 16:54:09
// version: 0.1.0
// url: git@github.com:prisms-center/IntegrationTools.git
// commit: 947b873eb0296ab1631408af59ecf7768231ff63

#ifndef von_mises_HH
#define von_mises_HH

#include <cmath>
#include <cstdlib>
#include "../../../../utils/IntegrationTools/PFunction.hh"

namespace PRISMS
{
    template< class VarContainer>
    class von_mises_f : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  8.1649658092772603e-01*var[4]+pow( pow( 6.6666666666666663e-01*var[0]+-3.3333333333333331e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+6.6666666666666663e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+-3.3333333333333331e-01*var[1]+6.6666666666666663e-01*var[2],2.0000000000000000e+00),(1.0/2.0))+-8.1649658092772603e-01*var[3];
        };

    public:

        von_mises_f()
        {
            this->_name = "von_mises_f";
        };

        von_mises_f* clone() const
        {
            return new von_mises_f(*this);
        };
    };

    template< class VarContainer>
    class von_mises_grad_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return (1.0/2.0)*( 1.3333333333333333e+00*var[0]+-6.6666666666666663e-01*var[1]+-6.6666666666666663e-01*var[2])*pow( pow( 6.6666666666666663e-01*var[0]+-3.3333333333333331e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+6.6666666666666663e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+-3.3333333333333331e-01*var[1]+6.6666666666666663e-01*var[2],2.0000000000000000e+00),-(1.0/2.0));
        };

    public:

        von_mises_grad_0()
        {
            this->_name = "von_mises_grad_0";
        };

        von_mises_grad_0* clone() const
        {
            return new von_mises_grad_0(*this);
        };
    };

    template< class VarContainer>
    class von_mises_grad_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return (1.0/2.0)*( -6.6666666666666663e-01*var[0]+1.3333333333333333e+00*var[1]+-6.6666666666666663e-01*var[2])*pow( pow( 6.6666666666666663e-01*var[0]+-3.3333333333333331e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+6.6666666666666663e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+-3.3333333333333331e-01*var[1]+6.6666666666666663e-01*var[2],2.0000000000000000e+00),-(1.0/2.0));
        };

    public:

        von_mises_grad_1()
        {
            this->_name = "von_mises_grad_1";
        };

        von_mises_grad_1* clone() const
        {
            return new von_mises_grad_1(*this);
        };
    };

    template< class VarContainer>
    class von_mises_grad_2 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return (1.0/2.0)*( -6.6666666666666663e-01*var[0]+-6.6666666666666663e-01*var[1]+1.3333333333333333e+00*var[2])*pow( pow( 6.6666666666666663e-01*var[0]+-3.3333333333333331e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+6.6666666666666663e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+-3.3333333333333331e-01*var[1]+6.6666666666666663e-01*var[2],2.0000000000000000e+00),-(1.0/2.0));
        };

    public:

        von_mises_grad_2()
        {
            this->_name = "von_mises_grad_2";
        };

        von_mises_grad_2* clone() const
        {
            return new von_mises_grad_2(*this);
        };
    };

    template< class VarContainer>
    class von_mises_grad_3 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return -8.1649658092772603e-01;
        };

    public:

        von_mises_grad_3()
        {
            this->_name = "von_mises_grad_3";
        };

        von_mises_grad_3* clone() const
        {
            return new von_mises_grad_3(*this);
        };
    };

    template< class VarContainer>
    class von_mises_grad_4 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 8.1649658092772603e-01;
        };

    public:

        von_mises_grad_4()
        {
            this->_name = "von_mises_grad_4";
        };

        von_mises_grad_4* clone() const
        {
            return new von_mises_grad_4(*this);
        };
    };

    template< class VarContainer>
    class von_mises_hess_0_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  6.6666666666666663e-01*pow( pow( 6.6666666666666663e-01*var[0]+-3.3333333333333331e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+6.6666666666666663e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+-3.3333333333333331e-01*var[1]+6.6666666666666663e-01*var[2],2.0000000000000000e+00),-(1.0/2.0))-pow( 1.3333333333333333e+00*var[0]+-6.6666666666666663e-01*var[1]+-6.6666666666666663e-01*var[2],2.0)*pow( pow( 6.6666666666666663e-01*var[0]+-3.3333333333333331e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+6.6666666666666663e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+-3.3333333333333331e-01*var[1]+6.6666666666666663e-01*var[2],2.0000000000000000e+00),-(3.0/2.0))/4.0;
        };

    public:

        von_mises_hess_0_0()
        {
            this->_name = "von_mises_hess_0_0";
        };

        von_mises_hess_0_0* clone() const
        {
            return new von_mises_hess_0_0(*this);
        };
    };

    template< class VarContainer>
    class von_mises_hess_0_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  -3.3333333333333331e-01*pow( pow( 6.6666666666666663e-01*var[0]+-3.3333333333333331e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+6.6666666666666663e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+-3.3333333333333331e-01*var[1]+6.6666666666666663e-01*var[2],2.0000000000000000e+00),-(1.0/2.0))-( -6.6666666666666663e-01*var[0]+1.3333333333333333e+00*var[1]+-6.6666666666666663e-01*var[2])*( 1.3333333333333333e+00*var[0]+-6.6666666666666663e-01*var[1]+-6.6666666666666663e-01*var[2])*pow( pow( 6.6666666666666663e-01*var[0]+-3.3333333333333331e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+6.6666666666666663e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+-3.3333333333333331e-01*var[1]+6.6666666666666663e-01*var[2],2.0000000000000000e+00),-(3.0/2.0))/4.0;
        };

    public:

        von_mises_hess_0_1()
        {
            this->_name = "von_mises_hess_0_1";
        };

        von_mises_hess_0_1* clone() const
        {
            return new von_mises_hess_0_1(*this);
        };
    };

    template< class VarContainer>
    class von_mises_hess_0_2 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  -3.3333333333333331e-01*pow( pow( 6.6666666666666663e-01*var[0]+-3.3333333333333331e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+6.6666666666666663e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+-3.3333333333333331e-01*var[1]+6.6666666666666663e-01*var[2],2.0000000000000000e+00),-(1.0/2.0))-( -6.6666666666666663e-01*var[0]+-6.6666666666666663e-01*var[1]+1.3333333333333333e+00*var[2])*( 1.3333333333333333e+00*var[0]+-6.6666666666666663e-01*var[1]+-6.6666666666666663e-01*var[2])*pow( pow( 6.6666666666666663e-01*var[0]+-3.3333333333333331e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+6.6666666666666663e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+-3.3333333333333331e-01*var[1]+6.6666666666666663e-01*var[2],2.0000000000000000e+00),-(3.0/2.0))/4.0;
        };

    public:

        von_mises_hess_0_2()
        {
            this->_name = "von_mises_hess_0_2";
        };

        von_mises_hess_0_2* clone() const
        {
            return new von_mises_hess_0_2(*this);
        };
    };

    template< class VarContainer>
    class von_mises_hess_0_3 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        von_mises_hess_0_3()
        {
            this->_name = "von_mises_hess_0_3";
        };

        von_mises_hess_0_3* clone() const
        {
            return new von_mises_hess_0_3(*this);
        };
    };

    template< class VarContainer>
    class von_mises_hess_0_4 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        von_mises_hess_0_4()
        {
            this->_name = "von_mises_hess_0_4";
        };

        von_mises_hess_0_4* clone() const
        {
            return new von_mises_hess_0_4(*this);
        };
    };

    template< class VarContainer>
    class von_mises_hess_1_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  -3.3333333333333331e-01*pow( pow( 6.6666666666666663e-01*var[0]+-3.3333333333333331e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+6.6666666666666663e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+-3.3333333333333331e-01*var[1]+6.6666666666666663e-01*var[2],2.0000000000000000e+00),-(1.0/2.0))-( -6.6666666666666663e-01*var[0]+1.3333333333333333e+00*var[1]+-6.6666666666666663e-01*var[2])*( 1.3333333333333333e+00*var[0]+-6.6666666666666663e-01*var[1]+-6.6666666666666663e-01*var[2])*pow( pow( 6.6666666666666663e-01*var[0]+-3.3333333333333331e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+6.6666666666666663e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+-3.3333333333333331e-01*var[1]+6.6666666666666663e-01*var[2],2.0000000000000000e+00),-(3.0/2.0))/4.0;
        };

    public:

        von_mises_hess_1_0()
        {
            this->_name = "von_mises_hess_1_0";
        };

        von_mises_hess_1_0* clone() const
        {
            return new von_mises_hess_1_0(*this);
        };
    };

    template< class VarContainer>
    class von_mises_hess_1_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  6.6666666666666663e-01*pow( pow( 6.6666666666666663e-01*var[0]+-3.3333333333333331e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+6.6666666666666663e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+-3.3333333333333331e-01*var[1]+6.6666666666666663e-01*var[2],2.0000000000000000e+00),-(1.0/2.0))-pow( -6.6666666666666663e-01*var[0]+1.3333333333333333e+00*var[1]+-6.6666666666666663e-01*var[2],2.0)*pow( pow( 6.6666666666666663e-01*var[0]+-3.3333333333333331e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+6.6666666666666663e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+-3.3333333333333331e-01*var[1]+6.6666666666666663e-01*var[2],2.0000000000000000e+00),-(3.0/2.0))/4.0;
        };

    public:

        von_mises_hess_1_1()
        {
            this->_name = "von_mises_hess_1_1";
        };

        von_mises_hess_1_1* clone() const
        {
            return new von_mises_hess_1_1(*this);
        };
    };

    template< class VarContainer>
    class von_mises_hess_1_2 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  -3.3333333333333331e-01*pow( pow( 6.6666666666666663e-01*var[0]+-3.3333333333333331e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+6.6666666666666663e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+-3.3333333333333331e-01*var[1]+6.6666666666666663e-01*var[2],2.0000000000000000e+00),-(1.0/2.0))-( -6.6666666666666663e-01*var[0]+1.3333333333333333e+00*var[1]+-6.6666666666666663e-01*var[2])*( -6.6666666666666663e-01*var[0]+-6.6666666666666663e-01*var[1]+1.3333333333333333e+00*var[2])*pow( pow( 6.6666666666666663e-01*var[0]+-3.3333333333333331e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+6.6666666666666663e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+-3.3333333333333331e-01*var[1]+6.6666666666666663e-01*var[2],2.0000000000000000e+00),-(3.0/2.0))/4.0;
        };

    public:

        von_mises_hess_1_2()
        {
            this->_name = "von_mises_hess_1_2";
        };

        von_mises_hess_1_2* clone() const
        {
            return new von_mises_hess_1_2(*this);
        };
    };

    template< class VarContainer>
    class von_mises_hess_1_3 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        von_mises_hess_1_3()
        {
            this->_name = "von_mises_hess_1_3";
        };

        von_mises_hess_1_3* clone() const
        {
            return new von_mises_hess_1_3(*this);
        };
    };

    template< class VarContainer>
    class von_mises_hess_1_4 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        von_mises_hess_1_4()
        {
            this->_name = "von_mises_hess_1_4";
        };

        von_mises_hess_1_4* clone() const
        {
            return new von_mises_hess_1_4(*this);
        };
    };

    template< class VarContainer>
    class von_mises_hess_2_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  -3.3333333333333331e-01*pow( pow( 6.6666666666666663e-01*var[0]+-3.3333333333333331e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+6.6666666666666663e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+-3.3333333333333331e-01*var[1]+6.6666666666666663e-01*var[2],2.0000000000000000e+00),-(1.0/2.0))-( -6.6666666666666663e-01*var[0]+-6.6666666666666663e-01*var[1]+1.3333333333333333e+00*var[2])*( 1.3333333333333333e+00*var[0]+-6.6666666666666663e-01*var[1]+-6.6666666666666663e-01*var[2])*pow( pow( 6.6666666666666663e-01*var[0]+-3.3333333333333331e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+6.6666666666666663e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+-3.3333333333333331e-01*var[1]+6.6666666666666663e-01*var[2],2.0000000000000000e+00),-(3.0/2.0))/4.0;
        };

    public:

        von_mises_hess_2_0()
        {
            this->_name = "von_mises_hess_2_0";
        };

        von_mises_hess_2_0* clone() const
        {
            return new von_mises_hess_2_0(*this);
        };
    };

    template< class VarContainer>
    class von_mises_hess_2_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  -3.3333333333333331e-01*pow( pow( 6.6666666666666663e-01*var[0]+-3.3333333333333331e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+6.6666666666666663e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+-3.3333333333333331e-01*var[1]+6.6666666666666663e-01*var[2],2.0000000000000000e+00),-(1.0/2.0))-( -6.6666666666666663e-01*var[0]+1.3333333333333333e+00*var[1]+-6.6666666666666663e-01*var[2])*( -6.6666666666666663e-01*var[0]+-6.6666666666666663e-01*var[1]+1.3333333333333333e+00*var[2])*pow( pow( 6.6666666666666663e-01*var[0]+-3.3333333333333331e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+6.6666666666666663e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+-3.3333333333333331e-01*var[1]+6.6666666666666663e-01*var[2],2.0000000000000000e+00),-(3.0/2.0))/4.0;
        };

    public:

        von_mises_hess_2_1()
        {
            this->_name = "von_mises_hess_2_1";
        };

        von_mises_hess_2_1* clone() const
        {
            return new von_mises_hess_2_1(*this);
        };
    };

    template< class VarContainer>
    class von_mises_hess_2_2 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  6.6666666666666663e-01*pow( pow( 6.6666666666666663e-01*var[0]+-3.3333333333333331e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+6.6666666666666663e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+-3.3333333333333331e-01*var[1]+6.6666666666666663e-01*var[2],2.0000000000000000e+00),-(1.0/2.0))-pow( -6.6666666666666663e-01*var[0]+-6.6666666666666663e-01*var[1]+1.3333333333333333e+00*var[2],2.0)*pow( pow( 6.6666666666666663e-01*var[0]+-3.3333333333333331e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+6.6666666666666663e-01*var[1]+-3.3333333333333331e-01*var[2],2.0000000000000000e+00)+pow( -3.3333333333333331e-01*var[0]+-3.3333333333333331e-01*var[1]+6.6666666666666663e-01*var[2],2.0000000000000000e+00),-(3.0/2.0))/4.0;
        };

    public:

        von_mises_hess_2_2()
        {
            this->_name = "von_mises_hess_2_2";
        };

        von_mises_hess_2_2* clone() const
        {
            return new von_mises_hess_2_2(*this);
        };
    };

    template< class VarContainer>
    class von_mises_hess_2_3 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        von_mises_hess_2_3()
        {
            this->_name = "von_mises_hess_2_3";
        };

        von_mises_hess_2_3* clone() const
        {
            return new von_mises_hess_2_3(*this);
        };
    };

    template< class VarContainer>
    class von_mises_hess_2_4 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        von_mises_hess_2_4()
        {
            this->_name = "von_mises_hess_2_4";
        };

        von_mises_hess_2_4* clone() const
        {
            return new von_mises_hess_2_4(*this);
        };
    };

    template< class VarContainer>
    class von_mises_hess_3_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        von_mises_hess_3_0()
        {
            this->_name = "von_mises_hess_3_0";
        };

        von_mises_hess_3_0* clone() const
        {
            return new von_mises_hess_3_0(*this);
        };
    };

    template< class VarContainer>
    class von_mises_hess_3_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        von_mises_hess_3_1()
        {
            this->_name = "von_mises_hess_3_1";
        };

        von_mises_hess_3_1* clone() const
        {
            return new von_mises_hess_3_1(*this);
        };
    };

    template< class VarContainer>
    class von_mises_hess_3_2 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        von_mises_hess_3_2()
        {
            this->_name = "von_mises_hess_3_2";
        };

        von_mises_hess_3_2* clone() const
        {
            return new von_mises_hess_3_2(*this);
        };
    };

    template< class VarContainer>
    class von_mises_hess_3_3 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        von_mises_hess_3_3()
        {
            this->_name = "von_mises_hess_3_3";
        };

        von_mises_hess_3_3* clone() const
        {
            return new von_mises_hess_3_3(*this);
        };
    };

    template< class VarContainer>
    class von_mises_hess_3_4 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        von_mises_hess_3_4()
        {
            this->_name = "von_mises_hess_3_4";
        };

        von_mises_hess_3_4* clone() const
        {
            return new von_mises_hess_3_4(*this);
        };
    };

    template< class VarContainer>
    class von_mises_hess_4_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        von_mises_hess_4_0()
        {
            this->_name = "von_mises_hess_4_0";
        };

        von_mises_hess_4_0* clone() const
        {
            return new von_mises_hess_4_0(*this);
        };
    };

    template< class VarContainer>
    class von_mises_hess_4_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        von_mises_hess_4_1()
        {
            this->_name = "von_mises_hess_4_1";
        };

        von_mises_hess_4_1* clone() const
        {
            return new von_mises_hess_4_1(*this);
        };
    };

    template< class VarContainer>
    class von_mises_hess_4_2 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        von_mises_hess_4_2()
        {
            this->_name = "von_mises_hess_4_2";
        };

        von_mises_hess_4_2* clone() const
        {
            return new von_mises_hess_4_2(*this);
        };
    };

    template< class VarContainer>
    class von_mises_hess_4_3 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        von_mises_hess_4_3()
        {
            this->_name = "von_mises_hess_4_3";
        };

        von_mises_hess_4_3* clone() const
        {
            return new von_mises_hess_4_3(*this);
        };
    };

    template< class VarContainer>
    class von_mises_hess_4_4 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        von_mises_hess_4_4()
        {
            this->_name = "von_mises_hess_4_4";
        };

        von_mises_hess_4_4* clone() const
        {
            return new von_mises_hess_4_4(*this);
        };
    };

    template<class VarContainer>
    class von_mises : public PFuncBase< VarContainer, double>
    {
    public:
        using PFuncBase< VarContainer, double>::_name;
        using PFuncBase< VarContainer, double>::_var_name;
        using PFuncBase< VarContainer, double>::_var_description;
        
        PSimpleBase< VarContainer, double> *_val;
        PSimpleBase< VarContainer, double> **_grad_val;
        PSimpleBase< VarContainer, double> ***_hess_val;
        
        von_mises()
        {
            construct();
        }

        von_mises(const von_mises &RHS )
        {
            construct();
        }

        von_mises& operator=(const von_mises &RHS )
        {
            _val = RHS._val;
            
            _grad_val[0] = RHS._grad_val[0];
            _grad_val[1] = RHS._grad_val[1];
            _grad_val[2] = RHS._grad_val[2];
            _grad_val[3] = RHS._grad_val[3];
            _grad_val[4] = RHS._grad_val[4];
            _hess_val[0][0] = RHS._hess_val[0][0];
            _hess_val[0][1] = RHS._hess_val[0][1];
            _hess_val[0][2] = RHS._hess_val[0][2];
            _hess_val[0][3] = RHS._hess_val[0][3];
            _hess_val[0][4] = RHS._hess_val[0][4];
            _hess_val[1][0] = RHS._hess_val[1][0];
            _hess_val[1][1] = RHS._hess_val[1][1];
            _hess_val[1][2] = RHS._hess_val[1][2];
            _hess_val[1][3] = RHS._hess_val[1][3];
            _hess_val[1][4] = RHS._hess_val[1][4];
            _hess_val[2][0] = RHS._hess_val[2][0];
            _hess_val[2][1] = RHS._hess_val[2][1];
            _hess_val[2][2] = RHS._hess_val[2][2];
            _hess_val[2][3] = RHS._hess_val[2][3];
            _hess_val[2][4] = RHS._hess_val[2][4];
            _hess_val[3][0] = RHS._hess_val[3][0];
            _hess_val[3][1] = RHS._hess_val[3][1];
            _hess_val[3][2] = RHS._hess_val[3][2];
            _hess_val[3][3] = RHS._hess_val[3][3];
            _hess_val[3][4] = RHS._hess_val[3][4];
            _hess_val[4][0] = RHS._hess_val[4][0];
            _hess_val[4][1] = RHS._hess_val[4][1];
            _hess_val[4][2] = RHS._hess_val[4][2];
            _hess_val[4][3] = RHS._hess_val[4][3];
            _hess_val[4][4] = RHS._hess_val[4][4];
        }

        ~von_mises()
        {
            delete _val;

            delete _grad_val[0];
            delete _grad_val[1];
            delete _grad_val[2];
            delete _grad_val[3];
            delete _grad_val[4];
            delete [] _grad_val;

            delete _hess_val[0][0];
            delete _hess_val[0][1];
            delete _hess_val[0][2];
            delete _hess_val[0][3];
            delete _hess_val[0][4];
            delete _hess_val[1][0];
            delete _hess_val[1][1];
            delete _hess_val[1][2];
            delete _hess_val[1][3];
            delete _hess_val[1][4];
            delete _hess_val[2][0];
            delete _hess_val[2][1];
            delete _hess_val[2][2];
            delete _hess_val[2][3];
            delete _hess_val[2][4];
            delete _hess_val[3][0];
            delete _hess_val[3][1];
            delete _hess_val[3][2];
            delete _hess_val[3][3];
            delete _hess_val[3][4];
            delete _hess_val[4][0];
            delete _hess_val[4][1];
            delete _hess_val[4][2];
            delete _hess_val[4][3];
            delete _hess_val[4][4];
            delete [] _hess_val[0];
            delete [] _hess_val[1];
            delete [] _hess_val[2];
            delete [] _hess_val[3];
            delete [] _hess_val[4];
            delete [] _hess_val;
        };

        von_mises<VarContainer>* clone() const
        {
            return new von_mises<VarContainer>(*this);
        };

        PSimpleFunction< VarContainer, double> simplefunction() const
        {
            return PSimpleFunction< VarContainer, double>( *_val );
        };

        PSimpleFunction< VarContainer, double> grad_simplefunction(int di) const
        {
            return PSimpleFunction< VarContainer, double>( *_grad_val[di] );
        };

        PSimpleFunction< VarContainer, double> hess_simplefunction(int di, int dj) const
        {
            return PSimpleFunction< VarContainer, double>( *_hess_val[di][dj] );
        };

        double operator()(const VarContainer &var)
        {
            return (*_val)(var);
        };

        double grad(const VarContainer &var, int di)
        {
            return (*_grad_val[di])(var);
        };

        double hess(const VarContainer &var, int di, int dj)
        {
            return (*_hess_val[di][dj])(var);
        };

        void eval(const VarContainer &var)
        {
            (*_val)(var);
        };

        void eval_grad(const VarContainer &var)
        {
            (*_grad_val[0])(var);
            (*_grad_val[1])(var);
            (*_grad_val[2])(var);
            (*_grad_val[3])(var);
            (*_grad_val[4])(var);
        };

        void eval_hess(const VarContainer &var)
        {
            (*_hess_val[0][0])(var);
            (*_hess_val[0][1])(var);
            (*_hess_val[0][2])(var);
            (*_hess_val[0][3])(var);
            (*_hess_val[0][4])(var);
            (*_hess_val[1][0])(var);
            (*_hess_val[1][1])(var);
            (*_hess_val[1][2])(var);
            (*_hess_val[1][3])(var);
            (*_hess_val[1][4])(var);
            (*_hess_val[2][0])(var);
            (*_hess_val[2][1])(var);
            (*_hess_val[2][2])(var);
            (*_hess_val[2][3])(var);
            (*_hess_val[2][4])(var);
            (*_hess_val[3][0])(var);
            (*_hess_val[3][1])(var);
            (*_hess_val[3][2])(var);
            (*_hess_val[3][3])(var);
            (*_hess_val[3][4])(var);
            (*_hess_val[4][0])(var);
            (*_hess_val[4][1])(var);
            (*_hess_val[4][2])(var);
            (*_hess_val[4][3])(var);
            (*_hess_val[4][4])(var);
        };

        double operator()() const
        {
            return (*_val)();
        };

        double grad(int di) const
        {
            return (*_grad_val[di])();
        };

        double hess(int di, int dj) const
        {
            return (*_hess_val[di][dj])();
        };

    private:
        void construct()
        {
            _name = "von_mises";
            _var_name.clear();
            _var_name.push_back("beta1");
            _var_name.push_back("beta2");
            _var_name.push_back("beta3");
            _var_name.push_back("tau_y");
            _var_name.push_back("q");
            _var_description.clear();
            _var_description.push_back("First principle stress");
            _var_description.push_back("Second principle stress");
            _var_description.push_back("Third principle stress");
            _var_description.push_back("Yield stress");
            _var_description.push_back("Conjugate stress-like quantity");
            
            _val = new von_mises_f<VarContainer>();
            
            _grad_val = new PSimpleBase< VarContainer, double>*[5];
            _grad_val[0] = new von_mises_grad_0<VarContainer>();
            _grad_val[1] = new von_mises_grad_1<VarContainer>();
            _grad_val[2] = new von_mises_grad_2<VarContainer>();
            _grad_val[3] = new von_mises_grad_3<VarContainer>();
            _grad_val[4] = new von_mises_grad_4<VarContainer>();
            
            _hess_val = new PSimpleBase< VarContainer, double>**[5];
            _hess_val[0] = new PSimpleBase< VarContainer, double>*[5];
            _hess_val[1] = new PSimpleBase< VarContainer, double>*[5];
            _hess_val[2] = new PSimpleBase< VarContainer, double>*[5];
            _hess_val[3] = new PSimpleBase< VarContainer, double>*[5];
            _hess_val[4] = new PSimpleBase< VarContainer, double>*[5];
            _hess_val[0][0] = new von_mises_hess_0_0<VarContainer>();
            _hess_val[0][1] = new von_mises_hess_0_1<VarContainer>();
            _hess_val[0][2] = new von_mises_hess_0_2<VarContainer>();
            _hess_val[0][3] = new von_mises_hess_0_3<VarContainer>();
            _hess_val[0][4] = new von_mises_hess_0_4<VarContainer>();
            _hess_val[1][0] = new von_mises_hess_1_0<VarContainer>();
            _hess_val[1][1] = new von_mises_hess_1_1<VarContainer>();
            _hess_val[1][2] = new von_mises_hess_1_2<VarContainer>();
            _hess_val[1][3] = new von_mises_hess_1_3<VarContainer>();
            _hess_val[1][4] = new von_mises_hess_1_4<VarContainer>();
            _hess_val[2][0] = new von_mises_hess_2_0<VarContainer>();
            _hess_val[2][1] = new von_mises_hess_2_1<VarContainer>();
            _hess_val[2][2] = new von_mises_hess_2_2<VarContainer>();
            _hess_val[2][3] = new von_mises_hess_2_3<VarContainer>();
            _hess_val[2][4] = new von_mises_hess_2_4<VarContainer>();
            _hess_val[3][0] = new von_mises_hess_3_0<VarContainer>();
            _hess_val[3][1] = new von_mises_hess_3_1<VarContainer>();
            _hess_val[3][2] = new von_mises_hess_3_2<VarContainer>();
            _hess_val[3][3] = new von_mises_hess_3_3<VarContainer>();
            _hess_val[3][4] = new von_mises_hess_3_4<VarContainer>();
            _hess_val[4][0] = new von_mises_hess_4_0<VarContainer>();
            _hess_val[4][1] = new von_mises_hess_4_1<VarContainer>();
            _hess_val[4][2] = new von_mises_hess_4_2<VarContainer>();
            _hess_val[4][3] = new von_mises_hess_4_3<VarContainer>();
            _hess_val[4][4] = new von_mises_hess_4_4<VarContainer>();
        };

    };
}
#endif
