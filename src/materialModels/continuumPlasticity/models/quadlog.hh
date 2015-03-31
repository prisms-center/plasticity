// created: 2015-3-5 16:54:09
// version: 0.1.0
// url: git@github.com:prisms-center/IntegrationTools.git
// commit: 947b873eb0296ab1631408af59ecf7768231ff63

#ifndef quadlog_HH
#define quadlog_HH

#include <cmath>
#include <cstdlib>
#include "../../../../utils/IntegrationTools/PFunction.hh"

namespace PRISMS
{
    template< class VarContainer>
    class quadlog_f : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  ( pow(log(var[2]),2.0000000000000000e+00)+pow(log(var[3]),2.0000000000000000e+00)+pow(log(var[4]),2.0000000000000000e+00))*var[1]+5.0000000000000000e-01*pow( log(var[2])+log(var[4])+log(var[3]),2.0000000000000000e+00)*var[0]+5.0000000000000000e-01*var[5]*pow(var[6],2.0000000000000000e+00);
        };

    public:

        quadlog_f()
        {
            this->_name = "quadlog_f";
        };

        quadlog_f* clone() const
        {
            return new quadlog_f(*this);
        };
    };

    template< class VarContainer>
    class quadlog_grad_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 5.0000000000000000e-01*pow( log(var[2])+log(var[4])+log(var[3]),2.0000000000000000e+00);
        };

    public:

        quadlog_grad_0()
        {
            this->_name = "quadlog_grad_0";
        };

        quadlog_grad_0* clone() const
        {
            return new quadlog_grad_0(*this);
        };
    };

    template< class VarContainer>
    class quadlog_grad_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  pow(log(var[2]),2.0000000000000000e+00)+pow(log(var[3]),2.0000000000000000e+00)+pow(log(var[4]),2.0000000000000000e+00);
        };

    public:

        quadlog_grad_1()
        {
            this->_name = "quadlog_grad_1";
        };

        quadlog_grad_1* clone() const
        {
            return new quadlog_grad_1(*this);
        };
    };

    template< class VarContainer>
    class quadlog_grad_2 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  2.0000000000000000e+00*1.0/(var[2])*log(var[2])*var[1]+( log(var[2])+log(var[4])+log(var[3]))/var[2]*var[0];
        };

    public:

        quadlog_grad_2()
        {
            this->_name = "quadlog_grad_2";
        };

        quadlog_grad_2* clone() const
        {
            return new quadlog_grad_2(*this);
        };
    };

    template< class VarContainer>
    class quadlog_grad_3 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  ( log(var[2])+log(var[4])+log(var[3]))/var[3]*var[0]+2.0000000000000000e+00*1.0/(var[3])*log(var[3])*var[1];
        };

    public:

        quadlog_grad_3()
        {
            this->_name = "quadlog_grad_3";
        };

        quadlog_grad_3* clone() const
        {
            return new quadlog_grad_3(*this);
        };
    };

    template< class VarContainer>
    class quadlog_grad_4 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  ( log(var[2])+log(var[4])+log(var[3]))*var[0]/var[4]+2.0000000000000000e+00*log(var[4])*var[1]*1.0/(var[4]);
        };

    public:

        quadlog_grad_4()
        {
            this->_name = "quadlog_grad_4";
        };

        quadlog_grad_4* clone() const
        {
            return new quadlog_grad_4(*this);
        };
    };

    template< class VarContainer>
    class quadlog_grad_5 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 5.0000000000000000e-01*pow(var[6],2.0000000000000000e+00);
        };

    public:

        quadlog_grad_5()
        {
            this->_name = "quadlog_grad_5";
        };

        quadlog_grad_5* clone() const
        {
            return new quadlog_grad_5(*this);
        };
    };

    template< class VarContainer>
    class quadlog_grad_6 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return var[5]*var[6];
        };

    public:

        quadlog_grad_6()
        {
            this->_name = "quadlog_grad_6";
        };

        quadlog_grad_6* clone() const
        {
            return new quadlog_grad_6(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_0_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        quadlog_hess_0_0()
        {
            this->_name = "quadlog_hess_0_0";
        };

        quadlog_hess_0_0* clone() const
        {
            return new quadlog_hess_0_0(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_0_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        quadlog_hess_0_1()
        {
            this->_name = "quadlog_hess_0_1";
        };

        quadlog_hess_0_1* clone() const
        {
            return new quadlog_hess_0_1(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_0_2 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return ( log(var[2])+log(var[4])+log(var[3]))/var[2];
        };

    public:

        quadlog_hess_0_2()
        {
            this->_name = "quadlog_hess_0_2";
        };

        quadlog_hess_0_2* clone() const
        {
            return new quadlog_hess_0_2(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_0_3 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return ( log(var[2])+log(var[4])+log(var[3]))/var[3];
        };

    public:

        quadlog_hess_0_3()
        {
            this->_name = "quadlog_hess_0_3";
        };

        quadlog_hess_0_3* clone() const
        {
            return new quadlog_hess_0_3(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_0_4 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return ( log(var[2])+log(var[4])+log(var[3]))/var[4];
        };

    public:

        quadlog_hess_0_4()
        {
            this->_name = "quadlog_hess_0_4";
        };

        quadlog_hess_0_4* clone() const
        {
            return new quadlog_hess_0_4(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_0_5 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        quadlog_hess_0_5()
        {
            this->_name = "quadlog_hess_0_5";
        };

        quadlog_hess_0_5* clone() const
        {
            return new quadlog_hess_0_5(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_0_6 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        quadlog_hess_0_6()
        {
            this->_name = "quadlog_hess_0_6";
        };

        quadlog_hess_0_6* clone() const
        {
            return new quadlog_hess_0_6(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_1_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        quadlog_hess_1_0()
        {
            this->_name = "quadlog_hess_1_0";
        };

        quadlog_hess_1_0* clone() const
        {
            return new quadlog_hess_1_0(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_1_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        quadlog_hess_1_1()
        {
            this->_name = "quadlog_hess_1_1";
        };

        quadlog_hess_1_1* clone() const
        {
            return new quadlog_hess_1_1(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_1_2 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 2.0000000000000000e+00*1.0/(var[2])*log(var[2]);
        };

    public:

        quadlog_hess_1_2()
        {
            this->_name = "quadlog_hess_1_2";
        };

        quadlog_hess_1_2* clone() const
        {
            return new quadlog_hess_1_2(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_1_3 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 2.0000000000000000e+00*1.0/(var[3])*log(var[3]);
        };

    public:

        quadlog_hess_1_3()
        {
            this->_name = "quadlog_hess_1_3";
        };

        quadlog_hess_1_3* clone() const
        {
            return new quadlog_hess_1_3(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_1_4 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 2.0000000000000000e+00*log(var[4])*1.0/(var[4]);
        };

    public:

        quadlog_hess_1_4()
        {
            this->_name = "quadlog_hess_1_4";
        };

        quadlog_hess_1_4* clone() const
        {
            return new quadlog_hess_1_4(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_1_5 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        quadlog_hess_1_5()
        {
            this->_name = "quadlog_hess_1_5";
        };

        quadlog_hess_1_5* clone() const
        {
            return new quadlog_hess_1_5(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_1_6 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        quadlog_hess_1_6()
        {
            this->_name = "quadlog_hess_1_6";
        };

        quadlog_hess_1_6* clone() const
        {
            return new quadlog_hess_1_6(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_2_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return ( log(var[2])+log(var[4])+log(var[3]))/var[2];
        };

    public:

        quadlog_hess_2_0()
        {
            this->_name = "quadlog_hess_2_0";
        };

        quadlog_hess_2_0* clone() const
        {
            return new quadlog_hess_2_0(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_2_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 2.0000000000000000e+00*1.0/(var[2])*log(var[2]);
        };

    public:

        quadlog_hess_2_1()
        {
            this->_name = "quadlog_hess_2_1";
        };

        quadlog_hess_2_1* clone() const
        {
            return new quadlog_hess_2_1(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_2_2 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  -2.0000000000000000e+00*1.0/(var[2]*var[2])*log(var[2])*var[1]+1.0/(var[2]*var[2])*var[0]-( log(var[2])+log(var[4])+log(var[3]))/(var[2]*var[2])*var[0]+2.0000000000000000e+00*1.0/(var[2])/var[2]*var[1];
        };

    public:

        quadlog_hess_2_2()
        {
            this->_name = "quadlog_hess_2_2";
        };

        quadlog_hess_2_2* clone() const
        {
            return new quadlog_hess_2_2(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_2_3 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 1.0/var[2]/var[3]*var[0];
        };

    public:

        quadlog_hess_2_3()
        {
            this->_name = "quadlog_hess_2_3";
        };

        quadlog_hess_2_3* clone() const
        {
            return new quadlog_hess_2_3(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_2_4 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 1.0/var[2]*var[0]/var[4];
        };

    public:

        quadlog_hess_2_4()
        {
            this->_name = "quadlog_hess_2_4";
        };

        quadlog_hess_2_4* clone() const
        {
            return new quadlog_hess_2_4(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_2_5 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        quadlog_hess_2_5()
        {
            this->_name = "quadlog_hess_2_5";
        };

        quadlog_hess_2_5* clone() const
        {
            return new quadlog_hess_2_5(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_2_6 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        quadlog_hess_2_6()
        {
            this->_name = "quadlog_hess_2_6";
        };

        quadlog_hess_2_6* clone() const
        {
            return new quadlog_hess_2_6(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_3_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return ( log(var[2])+log(var[4])+log(var[3]))/var[3];
        };

    public:

        quadlog_hess_3_0()
        {
            this->_name = "quadlog_hess_3_0";
        };

        quadlog_hess_3_0* clone() const
        {
            return new quadlog_hess_3_0(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_3_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 2.0000000000000000e+00*1.0/(var[3])*log(var[3]);
        };

    public:

        quadlog_hess_3_1()
        {
            this->_name = "quadlog_hess_3_1";
        };

        quadlog_hess_3_1* clone() const
        {
            return new quadlog_hess_3_1(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_3_2 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 1.0/var[2]/var[3]*var[0];
        };

    public:

        quadlog_hess_3_2()
        {
            this->_name = "quadlog_hess_3_2";
        };

        quadlog_hess_3_2* clone() const
        {
            return new quadlog_hess_3_2(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_3_3 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  -2.0000000000000000e+00*1.0/(var[3]*var[3])*log(var[3])*var[1]+1.0/(var[3]*var[3])*var[0]+2.0000000000000000e+00*1.0/(var[3])/var[3]*var[1]-( log(var[2])+log(var[4])+log(var[3]))/(var[3]*var[3])*var[0];
        };

    public:

        quadlog_hess_3_3()
        {
            this->_name = "quadlog_hess_3_3";
        };

        quadlog_hess_3_3* clone() const
        {
            return new quadlog_hess_3_3(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_3_4 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 1.0/var[3]*var[0]/var[4];
        };

    public:

        quadlog_hess_3_4()
        {
            this->_name = "quadlog_hess_3_4";
        };

        quadlog_hess_3_4* clone() const
        {
            return new quadlog_hess_3_4(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_3_5 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        quadlog_hess_3_5()
        {
            this->_name = "quadlog_hess_3_5";
        };

        quadlog_hess_3_5* clone() const
        {
            return new quadlog_hess_3_5(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_3_6 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        quadlog_hess_3_6()
        {
            this->_name = "quadlog_hess_3_6";
        };

        quadlog_hess_3_6* clone() const
        {
            return new quadlog_hess_3_6(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_4_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return ( log(var[2])+log(var[4])+log(var[3]))/var[4];
        };

    public:

        quadlog_hess_4_0()
        {
            this->_name = "quadlog_hess_4_0";
        };

        quadlog_hess_4_0* clone() const
        {
            return new quadlog_hess_4_0(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_4_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 2.0000000000000000e+00*log(var[4])*1.0/(var[4]);
        };

    public:

        quadlog_hess_4_1()
        {
            this->_name = "quadlog_hess_4_1";
        };

        quadlog_hess_4_1* clone() const
        {
            return new quadlog_hess_4_1(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_4_2 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 1.0/var[2]*var[0]/var[4];
        };

    public:

        quadlog_hess_4_2()
        {
            this->_name = "quadlog_hess_4_2";
        };

        quadlog_hess_4_2* clone() const
        {
            return new quadlog_hess_4_2(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_4_3 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 1.0/var[3]*var[0]/var[4];
        };

    public:

        quadlog_hess_4_3()
        {
            this->_name = "quadlog_hess_4_3";
        };

        quadlog_hess_4_3* clone() const
        {
            return new quadlog_hess_4_3(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_4_4 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  2.0000000000000000e+00*1.0/var[4]*var[1]*1.0/(var[4])+-2.0000000000000000e+00*log(var[4])/(var[4]*var[4])*var[1]-( log(var[2])+log(var[4])+log(var[3]))*var[0]/(var[4]*var[4])+var[0]/(var[4]*var[4]);
        };

    public:

        quadlog_hess_4_4()
        {
            this->_name = "quadlog_hess_4_4";
        };

        quadlog_hess_4_4* clone() const
        {
            return new quadlog_hess_4_4(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_4_5 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        quadlog_hess_4_5()
        {
            this->_name = "quadlog_hess_4_5";
        };

        quadlog_hess_4_5* clone() const
        {
            return new quadlog_hess_4_5(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_4_6 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        quadlog_hess_4_6()
        {
            this->_name = "quadlog_hess_4_6";
        };

        quadlog_hess_4_6* clone() const
        {
            return new quadlog_hess_4_6(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_5_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        quadlog_hess_5_0()
        {
            this->_name = "quadlog_hess_5_0";
        };

        quadlog_hess_5_0* clone() const
        {
            return new quadlog_hess_5_0(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_5_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        quadlog_hess_5_1()
        {
            this->_name = "quadlog_hess_5_1";
        };

        quadlog_hess_5_1* clone() const
        {
            return new quadlog_hess_5_1(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_5_2 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        quadlog_hess_5_2()
        {
            this->_name = "quadlog_hess_5_2";
        };

        quadlog_hess_5_2* clone() const
        {
            return new quadlog_hess_5_2(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_5_3 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        quadlog_hess_5_3()
        {
            this->_name = "quadlog_hess_5_3";
        };

        quadlog_hess_5_3* clone() const
        {
            return new quadlog_hess_5_3(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_5_4 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        quadlog_hess_5_4()
        {
            this->_name = "quadlog_hess_5_4";
        };

        quadlog_hess_5_4* clone() const
        {
            return new quadlog_hess_5_4(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_5_5 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        quadlog_hess_5_5()
        {
            this->_name = "quadlog_hess_5_5";
        };

        quadlog_hess_5_5* clone() const
        {
            return new quadlog_hess_5_5(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_5_6 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return var[6];
        };

    public:

        quadlog_hess_5_6()
        {
            this->_name = "quadlog_hess_5_6";
        };

        quadlog_hess_5_6* clone() const
        {
            return new quadlog_hess_5_6(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_6_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        quadlog_hess_6_0()
        {
            this->_name = "quadlog_hess_6_0";
        };

        quadlog_hess_6_0* clone() const
        {
            return new quadlog_hess_6_0(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_6_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        quadlog_hess_6_1()
        {
            this->_name = "quadlog_hess_6_1";
        };

        quadlog_hess_6_1* clone() const
        {
            return new quadlog_hess_6_1(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_6_2 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        quadlog_hess_6_2()
        {
            this->_name = "quadlog_hess_6_2";
        };

        quadlog_hess_6_2* clone() const
        {
            return new quadlog_hess_6_2(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_6_3 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        quadlog_hess_6_3()
        {
            this->_name = "quadlog_hess_6_3";
        };

        quadlog_hess_6_3* clone() const
        {
            return new quadlog_hess_6_3(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_6_4 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        quadlog_hess_6_4()
        {
            this->_name = "quadlog_hess_6_4";
        };

        quadlog_hess_6_4* clone() const
        {
            return new quadlog_hess_6_4(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_6_5 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return var[6];
        };

    public:

        quadlog_hess_6_5()
        {
            this->_name = "quadlog_hess_6_5";
        };

        quadlog_hess_6_5* clone() const
        {
            return new quadlog_hess_6_5(*this);
        };
    };

    template< class VarContainer>
    class quadlog_hess_6_6 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return var[5];
        };

    public:

        quadlog_hess_6_6()
        {
            this->_name = "quadlog_hess_6_6";
        };

        quadlog_hess_6_6* clone() const
        {
            return new quadlog_hess_6_6(*this);
        };
    };

    template<class VarContainer>
    class quadlog : public PFuncBase< VarContainer, double>
    {
    public:
        using PFuncBase< VarContainer, double>::_name;
        using PFuncBase< VarContainer, double>::_var_name;
        using PFuncBase< VarContainer, double>::_var_description;
        
        PSimpleBase< VarContainer, double> *_val;
        PSimpleBase< VarContainer, double> **_grad_val;
        PSimpleBase< VarContainer, double> ***_hess_val;
        
        quadlog()
        {
            construct();
        }

        quadlog(const quadlog &RHS )
        {
            construct();
        }

        quadlog& operator=(const quadlog &RHS )
        {
            _val = RHS._val;
            
            _grad_val[0] = RHS._grad_val[0];
            _grad_val[1] = RHS._grad_val[1];
            _grad_val[2] = RHS._grad_val[2];
            _grad_val[3] = RHS._grad_val[3];
            _grad_val[4] = RHS._grad_val[4];
            _grad_val[5] = RHS._grad_val[5];
            _grad_val[6] = RHS._grad_val[6];
            _hess_val[0][0] = RHS._hess_val[0][0];
            _hess_val[0][1] = RHS._hess_val[0][1];
            _hess_val[0][2] = RHS._hess_val[0][2];
            _hess_val[0][3] = RHS._hess_val[0][3];
            _hess_val[0][4] = RHS._hess_val[0][4];
            _hess_val[0][5] = RHS._hess_val[0][5];
            _hess_val[0][6] = RHS._hess_val[0][6];
            _hess_val[1][0] = RHS._hess_val[1][0];
            _hess_val[1][1] = RHS._hess_val[1][1];
            _hess_val[1][2] = RHS._hess_val[1][2];
            _hess_val[1][3] = RHS._hess_val[1][3];
            _hess_val[1][4] = RHS._hess_val[1][4];
            _hess_val[1][5] = RHS._hess_val[1][5];
            _hess_val[1][6] = RHS._hess_val[1][6];
            _hess_val[2][0] = RHS._hess_val[2][0];
            _hess_val[2][1] = RHS._hess_val[2][1];
            _hess_val[2][2] = RHS._hess_val[2][2];
            _hess_val[2][3] = RHS._hess_val[2][3];
            _hess_val[2][4] = RHS._hess_val[2][4];
            _hess_val[2][5] = RHS._hess_val[2][5];
            _hess_val[2][6] = RHS._hess_val[2][6];
            _hess_val[3][0] = RHS._hess_val[3][0];
            _hess_val[3][1] = RHS._hess_val[3][1];
            _hess_val[3][2] = RHS._hess_val[3][2];
            _hess_val[3][3] = RHS._hess_val[3][3];
            _hess_val[3][4] = RHS._hess_val[3][4];
            _hess_val[3][5] = RHS._hess_val[3][5];
            _hess_val[3][6] = RHS._hess_val[3][6];
            _hess_val[4][0] = RHS._hess_val[4][0];
            _hess_val[4][1] = RHS._hess_val[4][1];
            _hess_val[4][2] = RHS._hess_val[4][2];
            _hess_val[4][3] = RHS._hess_val[4][3];
            _hess_val[4][4] = RHS._hess_val[4][4];
            _hess_val[4][5] = RHS._hess_val[4][5];
            _hess_val[4][6] = RHS._hess_val[4][6];
            _hess_val[5][0] = RHS._hess_val[5][0];
            _hess_val[5][1] = RHS._hess_val[5][1];
            _hess_val[5][2] = RHS._hess_val[5][2];
            _hess_val[5][3] = RHS._hess_val[5][3];
            _hess_val[5][4] = RHS._hess_val[5][4];
            _hess_val[5][5] = RHS._hess_val[5][5];
            _hess_val[5][6] = RHS._hess_val[5][6];
            _hess_val[6][0] = RHS._hess_val[6][0];
            _hess_val[6][1] = RHS._hess_val[6][1];
            _hess_val[6][2] = RHS._hess_val[6][2];
            _hess_val[6][3] = RHS._hess_val[6][3];
            _hess_val[6][4] = RHS._hess_val[6][4];
            _hess_val[6][5] = RHS._hess_val[6][5];
            _hess_val[6][6] = RHS._hess_val[6][6];
        }

        ~quadlog()
        {
            delete _val;

            delete _grad_val[0];
            delete _grad_val[1];
            delete _grad_val[2];
            delete _grad_val[3];
            delete _grad_val[4];
            delete _grad_val[5];
            delete _grad_val[6];
            delete [] _grad_val;

            delete _hess_val[0][0];
            delete _hess_val[0][1];
            delete _hess_val[0][2];
            delete _hess_val[0][3];
            delete _hess_val[0][4];
            delete _hess_val[0][5];
            delete _hess_val[0][6];
            delete _hess_val[1][0];
            delete _hess_val[1][1];
            delete _hess_val[1][2];
            delete _hess_val[1][3];
            delete _hess_val[1][4];
            delete _hess_val[1][5];
            delete _hess_val[1][6];
            delete _hess_val[2][0];
            delete _hess_val[2][1];
            delete _hess_val[2][2];
            delete _hess_val[2][3];
            delete _hess_val[2][4];
            delete _hess_val[2][5];
            delete _hess_val[2][6];
            delete _hess_val[3][0];
            delete _hess_val[3][1];
            delete _hess_val[3][2];
            delete _hess_val[3][3];
            delete _hess_val[3][4];
            delete _hess_val[3][5];
            delete _hess_val[3][6];
            delete _hess_val[4][0];
            delete _hess_val[4][1];
            delete _hess_val[4][2];
            delete _hess_val[4][3];
            delete _hess_val[4][4];
            delete _hess_val[4][5];
            delete _hess_val[4][6];
            delete _hess_val[5][0];
            delete _hess_val[5][1];
            delete _hess_val[5][2];
            delete _hess_val[5][3];
            delete _hess_val[5][4];
            delete _hess_val[5][5];
            delete _hess_val[5][6];
            delete _hess_val[6][0];
            delete _hess_val[6][1];
            delete _hess_val[6][2];
            delete _hess_val[6][3];
            delete _hess_val[6][4];
            delete _hess_val[6][5];
            delete _hess_val[6][6];
            delete [] _hess_val[0];
            delete [] _hess_val[1];
            delete [] _hess_val[2];
            delete [] _hess_val[3];
            delete [] _hess_val[4];
            delete [] _hess_val[5];
            delete [] _hess_val[6];
            delete [] _hess_val;
        };

        quadlog<VarContainer>* clone() const
        {
            return new quadlog<VarContainer>(*this);
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
            (*_grad_val[5])(var);
            (*_grad_val[6])(var);
        };

        void eval_hess(const VarContainer &var)
        {
            (*_hess_val[0][0])(var);
            (*_hess_val[0][1])(var);
            (*_hess_val[0][2])(var);
            (*_hess_val[0][3])(var);
            (*_hess_val[0][4])(var);
            (*_hess_val[0][5])(var);
            (*_hess_val[0][6])(var);
            (*_hess_val[1][0])(var);
            (*_hess_val[1][1])(var);
            (*_hess_val[1][2])(var);
            (*_hess_val[1][3])(var);
            (*_hess_val[1][4])(var);
            (*_hess_val[1][5])(var);
            (*_hess_val[1][6])(var);
            (*_hess_val[2][0])(var);
            (*_hess_val[2][1])(var);
            (*_hess_val[2][2])(var);
            (*_hess_val[2][3])(var);
            (*_hess_val[2][4])(var);
            (*_hess_val[2][5])(var);
            (*_hess_val[2][6])(var);
            (*_hess_val[3][0])(var);
            (*_hess_val[3][1])(var);
            (*_hess_val[3][2])(var);
            (*_hess_val[3][3])(var);
            (*_hess_val[3][4])(var);
            (*_hess_val[3][5])(var);
            (*_hess_val[3][6])(var);
            (*_hess_val[4][0])(var);
            (*_hess_val[4][1])(var);
            (*_hess_val[4][2])(var);
            (*_hess_val[4][3])(var);
            (*_hess_val[4][4])(var);
            (*_hess_val[4][5])(var);
            (*_hess_val[4][6])(var);
            (*_hess_val[5][0])(var);
            (*_hess_val[5][1])(var);
            (*_hess_val[5][2])(var);
            (*_hess_val[5][3])(var);
            (*_hess_val[5][4])(var);
            (*_hess_val[5][5])(var);
            (*_hess_val[5][6])(var);
            (*_hess_val[6][0])(var);
            (*_hess_val[6][1])(var);
            (*_hess_val[6][2])(var);
            (*_hess_val[6][3])(var);
            (*_hess_val[6][4])(var);
            (*_hess_val[6][5])(var);
            (*_hess_val[6][6])(var);
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
            _name = "quadlog";
            _var_name.clear();
            _var_name.push_back("lambda");
            _var_name.push_back("mu");
            _var_name.push_back("lambda1");
            _var_name.push_back("lambda2");
            _var_name.push_back("lambda3");
            _var_name.push_back("K");
            _var_name.push_back("alpha");
            _var_description.clear();
            _var_description.push_back("First Lame parameter");
            _var_description.push_back("Second Lame parameter");
            _var_description.push_back("First principle stretch");
            _var_description.push_back("Second principle stretch");
            _var_description.push_back("Third principle stretch");
            _var_description.push_back("Strain hardening coefficient");
            _var_description.push_back("Equivalent plastic strain");
            
            _val = new quadlog_f<VarContainer>();
            
            _grad_val = new PSimpleBase< VarContainer, double>*[7];
            _grad_val[0] = new quadlog_grad_0<VarContainer>();
            _grad_val[1] = new quadlog_grad_1<VarContainer>();
            _grad_val[2] = new quadlog_grad_2<VarContainer>();
            _grad_val[3] = new quadlog_grad_3<VarContainer>();
            _grad_val[4] = new quadlog_grad_4<VarContainer>();
            _grad_val[5] = new quadlog_grad_5<VarContainer>();
            _grad_val[6] = new quadlog_grad_6<VarContainer>();
            
            _hess_val = new PSimpleBase< VarContainer, double>**[7];
            _hess_val[0] = new PSimpleBase< VarContainer, double>*[7];
            _hess_val[1] = new PSimpleBase< VarContainer, double>*[7];
            _hess_val[2] = new PSimpleBase< VarContainer, double>*[7];
            _hess_val[3] = new PSimpleBase< VarContainer, double>*[7];
            _hess_val[4] = new PSimpleBase< VarContainer, double>*[7];
            _hess_val[5] = new PSimpleBase< VarContainer, double>*[7];
            _hess_val[6] = new PSimpleBase< VarContainer, double>*[7];
            _hess_val[0][0] = new quadlog_hess_0_0<VarContainer>();
            _hess_val[0][1] = new quadlog_hess_0_1<VarContainer>();
            _hess_val[0][2] = new quadlog_hess_0_2<VarContainer>();
            _hess_val[0][3] = new quadlog_hess_0_3<VarContainer>();
            _hess_val[0][4] = new quadlog_hess_0_4<VarContainer>();
            _hess_val[0][5] = new quadlog_hess_0_5<VarContainer>();
            _hess_val[0][6] = new quadlog_hess_0_6<VarContainer>();
            _hess_val[1][0] = new quadlog_hess_1_0<VarContainer>();
            _hess_val[1][1] = new quadlog_hess_1_1<VarContainer>();
            _hess_val[1][2] = new quadlog_hess_1_2<VarContainer>();
            _hess_val[1][3] = new quadlog_hess_1_3<VarContainer>();
            _hess_val[1][4] = new quadlog_hess_1_4<VarContainer>();
            _hess_val[1][5] = new quadlog_hess_1_5<VarContainer>();
            _hess_val[1][6] = new quadlog_hess_1_6<VarContainer>();
            _hess_val[2][0] = new quadlog_hess_2_0<VarContainer>();
            _hess_val[2][1] = new quadlog_hess_2_1<VarContainer>();
            _hess_val[2][2] = new quadlog_hess_2_2<VarContainer>();
            _hess_val[2][3] = new quadlog_hess_2_3<VarContainer>();
            _hess_val[2][4] = new quadlog_hess_2_4<VarContainer>();
            _hess_val[2][5] = new quadlog_hess_2_5<VarContainer>();
            _hess_val[2][6] = new quadlog_hess_2_6<VarContainer>();
            _hess_val[3][0] = new quadlog_hess_3_0<VarContainer>();
            _hess_val[3][1] = new quadlog_hess_3_1<VarContainer>();
            _hess_val[3][2] = new quadlog_hess_3_2<VarContainer>();
            _hess_val[3][3] = new quadlog_hess_3_3<VarContainer>();
            _hess_val[3][4] = new quadlog_hess_3_4<VarContainer>();
            _hess_val[3][5] = new quadlog_hess_3_5<VarContainer>();
            _hess_val[3][6] = new quadlog_hess_3_6<VarContainer>();
            _hess_val[4][0] = new quadlog_hess_4_0<VarContainer>();
            _hess_val[4][1] = new quadlog_hess_4_1<VarContainer>();
            _hess_val[4][2] = new quadlog_hess_4_2<VarContainer>();
            _hess_val[4][3] = new quadlog_hess_4_3<VarContainer>();
            _hess_val[4][4] = new quadlog_hess_4_4<VarContainer>();
            _hess_val[4][5] = new quadlog_hess_4_5<VarContainer>();
            _hess_val[4][6] = new quadlog_hess_4_6<VarContainer>();
            _hess_val[5][0] = new quadlog_hess_5_0<VarContainer>();
            _hess_val[5][1] = new quadlog_hess_5_1<VarContainer>();
            _hess_val[5][2] = new quadlog_hess_5_2<VarContainer>();
            _hess_val[5][3] = new quadlog_hess_5_3<VarContainer>();
            _hess_val[5][4] = new quadlog_hess_5_4<VarContainer>();
            _hess_val[5][5] = new quadlog_hess_5_5<VarContainer>();
            _hess_val[5][6] = new quadlog_hess_5_6<VarContainer>();
            _hess_val[6][0] = new quadlog_hess_6_0<VarContainer>();
            _hess_val[6][1] = new quadlog_hess_6_1<VarContainer>();
            _hess_val[6][2] = new quadlog_hess_6_2<VarContainer>();
            _hess_val[6][3] = new quadlog_hess_6_3<VarContainer>();
            _hess_val[6][4] = new quadlog_hess_6_4<VarContainer>();
            _hess_val[6][5] = new quadlog_hess_6_5<VarContainer>();
            _hess_val[6][6] = new quadlog_hess_6_6<VarContainer>();
        };

    };
}
#endif
