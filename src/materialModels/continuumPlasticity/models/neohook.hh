// created: 2016-8-5 14:00:36
// version: master
// url: https://github.com/prisms-center/IntegrationToolsWriter.git
// commit: 8a15adf67355fad30bd75ce9ba6b1f8d24b9a537

#ifndef neohook_HH
#define neohook_HH

#include <cmath>
#include <cstdlib>
#include "../../../../utils/IntegrationTools/PFunction.hh"


namespace PRISMS
{
    template< class VarContainer>
    class neohook_f : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  5.0000000000000000e-01*var[5]*pow(var[6],2.0000000000000000e+00)-log(var[2]*var[3]*var[4])*( 5.0000000000000000e-01*var[0]+var[1])+2.5000000000000000e-01*( pow(var[2],2.0000000000000000e+00)*pow(var[3],2.0000000000000000e+00)*pow(var[4],2.0000000000000000e+00)-1.0000000000000000e+00)*var[0]+5.0000000000000000e-01*( pow(var[2],2.0000000000000000e+00)+pow(var[4],2.0000000000000000e+00)+pow(var[3],2.0000000000000000e+00)-3.0000000000000000e+00)*var[1];
        }

    public:

        neohook_f()
        {
            this->_name = "neohook_f";
        }

        std::string csrc() const
        {
            return " 5.0000000000000000e-01*var[5]*pow(var[6],2.0000000000000000e+00)-log(var[2]*var[3]*var[4])*( 5.0000000000000000e-01*var[0]+var[1])+2.5000000000000000e-01*( pow(var[2],2.0000000000000000e+00)*pow(var[3],2.0000000000000000e+00)*pow(var[4],2.0000000000000000e+00)-1.0000000000000000e+00)*var[0]+5.0000000000000000e-01*( pow(var[2],2.0000000000000000e+00)+pow(var[4],2.0000000000000000e+00)+pow(var[3],2.0000000000000000e+00)-3.0000000000000000e+00)*var[1]";
        }

        std::string sym() const
        {
            return "-log(lambda3*lambda1*lambda2)*((0.5)*lambda+mu)+(0.25)*lambda*(-1.0+lambda3^(2.0)*lambda1^(2.0)*lambda2^(2.0))+(0.5)*K*alpha^(2.0)+(0.5)*(-3.0+lambda1^(2.0)+lambda3^(2.0)+lambda2^(2.0))*mu";
        }

        std::string latex() const
        {
            return "{(0.5)}  \\mu {(-3.0+lambda1^{{(2.0)}}+lambda3^{{(2.0)}}+lambda2^{{(2.0)}})}+{(0.25)}  {(-1.0+ lambda3^{{(2.0)}} lambda1^{{(2.0)}} lambda2^{{(2.0)}})} \\lambda+{(0.5)}  K \\alpha^{{(2.0)}}- \\ln( lambda3 lambda1 lambda2) {(\\mu+{(0.5)} \\lambda)}";
        }

        neohook_f* clone() const
        {
            return new neohook_f(*this);
        }
    };

    template< class VarContainer>
    class neohook_grad_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  2.5000000000000000e-01*pow(var[2],2.0000000000000000e+00)*pow(var[3],2.0000000000000000e+00)*pow(var[4],2.0000000000000000e+00)+-5.0000000000000000e-01*log(var[2]*var[3]*var[4])-2.5000000000000000e-01;
        }

    public:

        neohook_grad_0()
        {
            this->_name = "neohook_grad_0";
        }

        std::string csrc() const
        {
            return " 2.5000000000000000e-01*pow(var[2],2.0000000000000000e+00)*pow(var[3],2.0000000000000000e+00)*pow(var[4],2.0000000000000000e+00)+-5.0000000000000000e-01*log(var[2]*var[3]*var[4])-2.5000000000000000e-01";
        }

        std::string sym() const
        {
            return "-0.25-(0.5)*log(lambda3*lambda1*lambda2)+(0.25)*lambda3^(2.0)*lambda1^(2.0)*lambda2^(2.0)";
        }

        std::string latex() const
        {
            return "-0.25+{(0.25)}  lambda1^{{(2.0)}} lambda2^{{(2.0)}} lambda3^{{(2.0)}}-{(0.5)} \\ln( lambda1 lambda2 lambda3)";
        }

        neohook_grad_0* clone() const
        {
            return new neohook_grad_0(*this);
        }
    };

    template< class VarContainer>
    class neohook_grad_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  5.0000000000000000e-01*pow(var[4],2.0000000000000000e+00)+5.0000000000000000e-01*pow(var[2],2.0000000000000000e+00)-log(var[4]*var[2]*var[3])+5.0000000000000000e-01*pow(var[3],2.0000000000000000e+00)-1.5000000000000000e+00;
        }

    public:

        neohook_grad_1()
        {
            this->_name = "neohook_grad_1";
        }

        std::string csrc() const
        {
            return " 5.0000000000000000e-01*pow(var[4],2.0000000000000000e+00)+5.0000000000000000e-01*pow(var[2],2.0000000000000000e+00)-log(var[4]*var[2]*var[3])+5.0000000000000000e-01*pow(var[3],2.0000000000000000e+00)-1.5000000000000000e+00";
        }

        std::string sym() const
        {
            return "-1.5+(0.5)*lambda3^(2.0)-log(lambda2*lambda3*lambda1)+(0.5)*lambda2^(2.0)+(0.5)*lambda1^(2.0)";
        }

        std::string latex() const
        {
            return "-1.5+{(0.5)} lambda3^{{(2.0)}}-\\ln( lambda3 lambda1 lambda2)+{(0.5)} lambda2^{{(2.0)}}+{(0.5)} lambda1^{{(2.0)}}";
        }

        neohook_grad_1* clone() const
        {
            return new neohook_grad_1(*this);
        }
    };

    template< class VarContainer>
    class neohook_grad_2 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return -( 5.0000000000000000e-01*var[0]+var[1])/var[2]+var[1]*var[2]+5.0000000000000000e-01*pow(var[3],2.0000000000000000e+00)*var[0]*pow(var[4],2.0000000000000000e+00)*var[2];
        }

    public:

        neohook_grad_2()
        {
            this->_name = "neohook_grad_2";
        }

        std::string csrc() const
        {
            return "-( 5.0000000000000000e-01*var[0]+var[1])/var[2]+var[1]*var[2]+5.0000000000000000e-01*pow(var[3],2.0000000000000000e+00)*var[0]*pow(var[4],2.0000000000000000e+00)*var[2]";
        }

        std::string sym() const
        {
            return "mu*lambda1+(0.5)*lambda3^(2.0)*lambda1*lambda2^(2.0)*lambda-(mu+(0.5)*lambda)*lambda1^(-1)";
        }

        std::string latex() const
        {
            return " \\mu lambda1+{(0.5)}  lambda2^{{(2.0)}} \\lambda lambda3^{{(2.0)}} lambda1-\\frac{{(0.5)} \\lambda+\\mu}{lambda1}";
        }

        neohook_grad_2* clone() const
        {
            return new neohook_grad_2(*this);
        }
    };

    template< class VarContainer>
    class neohook_grad_3 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return -1.0/var[3]*( var[1]+5.0000000000000000e-01*var[0])+5.0000000000000000e-01*pow(var[2],2.0000000000000000e+00)*var[3]*var[0]*pow(var[4],2.0000000000000000e+00)+var[1]*var[3];
        }

    public:

        neohook_grad_3()
        {
            this->_name = "neohook_grad_3";
        }

        std::string csrc() const
        {
            return "-1.0/var[3]*( var[1]+5.0000000000000000e-01*var[0])+5.0000000000000000e-01*pow(var[2],2.0000000000000000e+00)*var[3]*var[0]*pow(var[4],2.0000000000000000e+00)+var[1]*var[3]";
        }

        std::string sym() const
        {
            return "lambda2*mu-lambda2^(-1)*((0.5)*lambda+mu)+(0.5)*lambda2*lambda*lambda3^(2.0)*lambda1^(2.0)";
        }

        std::string latex() const
        {
            return "-\\frac{{(0.5)} \\lambda+\\mu}{lambda2}+ lambda2 \\mu+{(0.5)}  lambda1^{{(2.0)}} lambda2 \\lambda lambda3^{{(2.0)}}";
        }

        neohook_grad_3* clone() const
        {
            return new neohook_grad_3(*this);
        }
    };

    template< class VarContainer>
    class neohook_grad_4 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  var[4]*var[1]-1.0/var[4]*( 5.0000000000000000e-01*var[0]+var[1])+5.0000000000000000e-01*pow(var[3],2.0000000000000000e+00)*var[0]*var[4]*pow(var[2],2.0000000000000000e+00);
        }

    public:

        neohook_grad_4()
        {
            this->_name = "neohook_grad_4";
        }

        std::string csrc() const
        {
            return " var[4]*var[1]-1.0/var[4]*( 5.0000000000000000e-01*var[0]+var[1])+5.0000000000000000e-01*pow(var[3],2.0000000000000000e+00)*var[0]*var[4]*pow(var[2],2.0000000000000000e+00)";
        }

        std::string sym() const
        {
            return "lambda3*mu+(0.5)*lambda1^(2.0)*lambda2^(2.0)*lambda*lambda3-((0.5)*lambda+mu)*lambda3^(-1)";
        }

        std::string latex() const
        {
            return "{(0.5)}  lambda2^{{(2.0)}} \\lambda lambda3 lambda1^{{(2.0)}}+ lambda3 \\mu-\\frac{{(0.5)} \\lambda+\\mu}{lambda3}";
        }

        neohook_grad_4* clone() const
        {
            return new neohook_grad_4(*this);
        }
    };

    template< class VarContainer>
    class neohook_grad_5 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 5.0000000000000000e-01*pow(var[6],2.0000000000000000e+00);
        }

    public:

        neohook_grad_5()
        {
            this->_name = "neohook_grad_5";
        }

        std::string csrc() const
        {
            return "5.0000000000000000e-01*pow(var[6],2.0000000000000000e+00)";
        }

        std::string sym() const
        {
            return "(0.5)*alpha^(2.0)";
        }

        std::string latex() const
        {
            return "{(0.5)}  \\alpha^{{(2.0)}}";
        }

        neohook_grad_5* clone() const
        {
            return new neohook_grad_5(*this);
        }
    };

    template< class VarContainer>
    class neohook_grad_6 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return var[5]*var[6];
        }

    public:

        neohook_grad_6()
        {
            this->_name = "neohook_grad_6";
        }

        std::string csrc() const
        {
            return "var[5]*var[6]";
        }

        std::string sym() const
        {
            return "alpha*K";
        }

        std::string latex() const
        {
            return " K \\alpha";
        }

        neohook_grad_6* clone() const
        {
            return new neohook_grad_6(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_0_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        neohook_hess_0_0()
        {
            this->_name = "neohook_hess_0_0";
        }

        std::string csrc() const
        {
            return "0.0";
        }

        std::string sym() const
        {
            return "0";
        }

        std::string latex() const
        {
            return "0";
        }

        neohook_hess_0_0* clone() const
        {
            return new neohook_hess_0_0(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_0_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        neohook_hess_0_1()
        {
            this->_name = "neohook_hess_0_1";
        }

        std::string csrc() const
        {
            return "0.0";
        }

        std::string sym() const
        {
            return "0";
        }

        std::string latex() const
        {
            return "0";
        }

        neohook_hess_0_1* clone() const
        {
            return new neohook_hess_0_1(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_0_2 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  -5.0000000000000000e-01*1.0/(var[2])+5.0000000000000000e-01*pow(var[3],2.0000000000000000e+00)*pow(var[4],2.0000000000000000e+00)*var[2];
        }

    public:

        neohook_hess_0_2()
        {
            this->_name = "neohook_hess_0_2";
        }

        std::string csrc() const
        {
            return " -5.0000000000000000e-01*1.0/(var[2])+5.0000000000000000e-01*pow(var[3],2.0000000000000000e+00)*pow(var[4],2.0000000000000000e+00)*var[2]";
        }

        std::string sym() const
        {
            return "-(0.5)*lambda1^(-1)+(0.5)*lambda1*lambda2^(2.0)*lambda3^(2.0)";
        }

        std::string latex() const
        {
            return "-{(0.5)} \\frac{1}{lambda1}+{(0.5)}  lambda2^{{(2.0)}} lambda3^{{(2.0)}} lambda1";
        }

        neohook_hess_0_2* clone() const
        {
            return new neohook_hess_0_2(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_0_3 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  -5.0000000000000000e-01*1.0/(var[3])+5.0000000000000000e-01*pow(var[2],2.0000000000000000e+00)*var[3]*pow(var[4],2.0000000000000000e+00);
        }

    public:

        neohook_hess_0_3()
        {
            this->_name = "neohook_hess_0_3";
        }

        std::string csrc() const
        {
            return " -5.0000000000000000e-01*1.0/(var[3])+5.0000000000000000e-01*pow(var[2],2.0000000000000000e+00)*var[3]*pow(var[4],2.0000000000000000e+00)";
        }

        std::string sym() const
        {
            return "(0.5)*lambda2*lambda3^(2.0)*lambda1^(2.0)-(0.5)*lambda2^(-1)";
        }

        std::string latex() const
        {
            return "-{(0.5)} \\frac{1}{lambda2}+{(0.5)}  lambda1^{{(2.0)}} lambda2 lambda3^{{(2.0)}}";
        }

        neohook_hess_0_3* clone() const
        {
            return new neohook_hess_0_3(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_0_4 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  -5.0000000000000000e-01*1.0/(var[4])+5.0000000000000000e-01*var[4]*pow(var[2],2.0000000000000000e+00)*pow(var[3],2.0000000000000000e+00);
        }

    public:

        neohook_hess_0_4()
        {
            this->_name = "neohook_hess_0_4";
        }

        std::string csrc() const
        {
            return " -5.0000000000000000e-01*1.0/(var[4])+5.0000000000000000e-01*var[4]*pow(var[2],2.0000000000000000e+00)*pow(var[3],2.0000000000000000e+00)";
        }

        std::string sym() const
        {
            return "-(0.5)*lambda3^(-1)+(0.5)*lambda1^(2.0)*lambda2^(2.0)*lambda3";
        }

        std::string latex() const
        {
            return "-{(0.5)} \\frac{1}{lambda3}+{(0.5)}  lambda3 lambda1^{{(2.0)}} lambda2^{{(2.0)}}";
        }

        neohook_hess_0_4* clone() const
        {
            return new neohook_hess_0_4(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_0_5 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        neohook_hess_0_5()
        {
            this->_name = "neohook_hess_0_5";
        }

        std::string csrc() const
        {
            return "0.0";
        }

        std::string sym() const
        {
            return "0";
        }

        std::string latex() const
        {
            return "0";
        }

        neohook_hess_0_5* clone() const
        {
            return new neohook_hess_0_5(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_0_6 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        neohook_hess_0_6()
        {
            this->_name = "neohook_hess_0_6";
        }

        std::string csrc() const
        {
            return "0.0";
        }

        std::string sym() const
        {
            return "0";
        }

        std::string latex() const
        {
            return "0";
        }

        neohook_hess_0_6* clone() const
        {
            return new neohook_hess_0_6(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_1_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        neohook_hess_1_0()
        {
            this->_name = "neohook_hess_1_0";
        }

        std::string csrc() const
        {
            return "0.0";
        }

        std::string sym() const
        {
            return "0";
        }

        std::string latex() const
        {
            return "0";
        }

        neohook_hess_1_0* clone() const
        {
            return new neohook_hess_1_0(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_1_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        neohook_hess_1_1()
        {
            this->_name = "neohook_hess_1_1";
        }

        std::string csrc() const
        {
            return "0.0";
        }

        std::string sym() const
        {
            return "0";
        }

        std::string latex() const
        {
            return "0";
        }

        neohook_hess_1_1* clone() const
        {
            return new neohook_hess_1_1(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_1_2 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  var[2]-1.0/(var[2]);
        }

    public:

        neohook_hess_1_2()
        {
            this->_name = "neohook_hess_1_2";
        }

        std::string csrc() const
        {
            return " var[2]-1.0/(var[2])";
        }

        std::string sym() const
        {
            return "lambda1-lambda1^(-1)";
        }

        std::string latex() const
        {
            return "lambda1-\\frac{1}{lambda1}";
        }

        neohook_hess_1_2* clone() const
        {
            return new neohook_hess_1_2(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_1_3 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return -1.0/(var[3])+var[3];
        }

    public:

        neohook_hess_1_3()
        {
            this->_name = "neohook_hess_1_3";
        }

        std::string csrc() const
        {
            return "-1.0/(var[3])+var[3]";
        }

        std::string sym() const
        {
            return "-lambda2^(-1)+lambda2";
        }

        std::string latex() const
        {
            return "-\\frac{1}{lambda2}+lambda2";
        }

        neohook_hess_1_3* clone() const
        {
            return new neohook_hess_1_3(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_1_4 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  var[4]-1.0/(var[4]);
        }

    public:

        neohook_hess_1_4()
        {
            this->_name = "neohook_hess_1_4";
        }

        std::string csrc() const
        {
            return " var[4]-1.0/(var[4])";
        }

        std::string sym() const
        {
            return "lambda3-lambda3^(-1)";
        }

        std::string latex() const
        {
            return "lambda3-\\frac{1}{lambda3}";
        }

        neohook_hess_1_4* clone() const
        {
            return new neohook_hess_1_4(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_1_5 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        neohook_hess_1_5()
        {
            this->_name = "neohook_hess_1_5";
        }

        std::string csrc() const
        {
            return "0.0";
        }

        std::string sym() const
        {
            return "0";
        }

        std::string latex() const
        {
            return "0";
        }

        neohook_hess_1_5* clone() const
        {
            return new neohook_hess_1_5(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_1_6 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        neohook_hess_1_6()
        {
            this->_name = "neohook_hess_1_6";
        }

        std::string csrc() const
        {
            return "0.0";
        }

        std::string sym() const
        {
            return "0";
        }

        std::string latex() const
        {
            return "0";
        }

        neohook_hess_1_6* clone() const
        {
            return new neohook_hess_1_6(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_2_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  -5.0000000000000000e-01*1.0/(var[2])+5.0000000000000000e-01*pow(var[3],2.0000000000000000e+00)*pow(var[4],2.0000000000000000e+00)*var[2];
        }

    public:

        neohook_hess_2_0()
        {
            this->_name = "neohook_hess_2_0";
        }

        std::string csrc() const
        {
            return " -5.0000000000000000e-01*1.0/(var[2])+5.0000000000000000e-01*pow(var[3],2.0000000000000000e+00)*pow(var[4],2.0000000000000000e+00)*var[2]";
        }

        std::string sym() const
        {
            return "-(0.5)*lambda1^(-1)+(0.5)*lambda1*lambda2^(2.0)*lambda3^(2.0)";
        }

        std::string latex() const
        {
            return "-{(0.5)} \\frac{1}{lambda1}+{(0.5)}  lambda2^{{(2.0)}} lambda3^{{(2.0)}} lambda1";
        }

        neohook_hess_2_0* clone() const
        {
            return new neohook_hess_2_0(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_2_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  var[2]-1.0/(var[2]);
        }

    public:

        neohook_hess_2_1()
        {
            this->_name = "neohook_hess_2_1";
        }

        std::string csrc() const
        {
            return " var[2]-1.0/(var[2])";
        }

        std::string sym() const
        {
            return "lambda1-lambda1^(-1)";
        }

        std::string latex() const
        {
            return "lambda1-\\frac{1}{lambda1}";
        }

        neohook_hess_2_1* clone() const
        {
            return new neohook_hess_2_1(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_2_2 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  ( 5.0000000000000000e-01*var[0]+var[1])/(var[2]*var[2])+5.0000000000000000e-01*var[0]*pow(var[4],2.0000000000000000e+00)*pow(var[3],2.0000000000000000e+00)+var[1];
        }

    public:

        neohook_hess_2_2()
        {
            this->_name = "neohook_hess_2_2";
        }

        std::string csrc() const
        {
            return " ( 5.0000000000000000e-01*var[0]+var[1])/(var[2]*var[2])+5.0000000000000000e-01*var[0]*pow(var[4],2.0000000000000000e+00)*pow(var[3],2.0000000000000000e+00)+var[1]";
        }

        std::string sym() const
        {
            return "(0.5)*lambda2^(2.0)*lambda*lambda3^(2.0)+mu+lambda1^(-2)*((0.5)*lambda+mu)";
        }

        std::string latex() const
        {
            return "\\mu+{(0.5)}  lambda3^{{(2.0)}} lambda2^{{(2.0)}} \\lambda+\\frac{\\mu+{(0.5)} \\lambda}{lambda1^{2}}";
        }

        neohook_hess_2_2* clone() const
        {
            return new neohook_hess_2_2(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_2_3 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return var[2]*var[3]*var[0]*pow(var[4],2.0000000000000000e+00);
        }

    public:

        neohook_hess_2_3()
        {
            this->_name = "neohook_hess_2_3";
        }

        std::string csrc() const
        {
            return "var[2]*var[3]*var[0]*pow(var[4],2.0000000000000000e+00)";
        }

        std::string sym() const
        {
            return "lambda3^(2.0)*lambda1*lambda2*lambda";
        }

        std::string latex() const
        {
            return " lambda1 lambda2 \\lambda lambda3^{{(2.0)}}";
        }

        neohook_hess_2_3* clone() const
        {
            return new neohook_hess_2_3(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_2_4 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return var[4]*var[2]*pow(var[3],2.0000000000000000e+00)*var[0];
        }

    public:

        neohook_hess_2_4()
        {
            this->_name = "neohook_hess_2_4";
        }

        std::string csrc() const
        {
            return "var[4]*var[2]*pow(var[3],2.0000000000000000e+00)*var[0]";
        }

        std::string sym() const
        {
            return "lambda2^(2.0)*lambda*lambda3*lambda1";
        }

        std::string latex() const
        {
            return " lambda3 lambda1 lambda2^{{(2.0)}} \\lambda";
        }

        neohook_hess_2_4* clone() const
        {
            return new neohook_hess_2_4(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_2_5 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        neohook_hess_2_5()
        {
            this->_name = "neohook_hess_2_5";
        }

        std::string csrc() const
        {
            return "0.0";
        }

        std::string sym() const
        {
            return "0";
        }

        std::string latex() const
        {
            return "0";
        }

        neohook_hess_2_5* clone() const
        {
            return new neohook_hess_2_5(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_2_6 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        neohook_hess_2_6()
        {
            this->_name = "neohook_hess_2_6";
        }

        std::string csrc() const
        {
            return "0.0";
        }

        std::string sym() const
        {
            return "0";
        }

        std::string latex() const
        {
            return "0";
        }

        neohook_hess_2_6* clone() const
        {
            return new neohook_hess_2_6(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_3_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  -5.0000000000000000e-01*1.0/(var[3])+5.0000000000000000e-01*var[3]*pow(var[4],2.0000000000000000e+00)*pow(var[2],2.0000000000000000e+00);
        }

    public:

        neohook_hess_3_0()
        {
            this->_name = "neohook_hess_3_0";
        }

        std::string csrc() const
        {
            return " -5.0000000000000000e-01*1.0/(var[3])+5.0000000000000000e-01*var[3]*pow(var[4],2.0000000000000000e+00)*pow(var[2],2.0000000000000000e+00)";
        }

        std::string sym() const
        {
            return "-(0.5)*lambda2^(-1)+(0.5)*lambda1^(2.0)*lambda2*lambda3^(2.0)";
        }

        std::string latex() const
        {
            return "-{(0.5)} \\frac{1}{lambda2}+{(0.5)}  lambda3^{{(2.0)}} lambda1^{{(2.0)}} lambda2";
        }

        neohook_hess_3_0* clone() const
        {
            return new neohook_hess_3_0(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_3_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return -1.0/(var[3])+var[3];
        }

    public:

        neohook_hess_3_1()
        {
            this->_name = "neohook_hess_3_1";
        }

        std::string csrc() const
        {
            return "-1.0/(var[3])+var[3]";
        }

        std::string sym() const
        {
            return "-lambda2^(-1)+lambda2";
        }

        std::string latex() const
        {
            return "-\\frac{1}{lambda2}+lambda2";
        }

        neohook_hess_3_1* clone() const
        {
            return new neohook_hess_3_1(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_3_2 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return pow(var[4],2.0000000000000000e+00)*var[2]*var[3]*var[0];
        }

    public:

        neohook_hess_3_2()
        {
            this->_name = "neohook_hess_3_2";
        }

        std::string csrc() const
        {
            return "pow(var[4],2.0000000000000000e+00)*var[2]*var[3]*var[0]";
        }

        std::string sym() const
        {
            return "lambda1*lambda2*lambda*lambda3^(2.0)";
        }

        std::string latex() const
        {
            return " lambda3^{{(2.0)}} lambda1 lambda2 \\lambda";
        }

        neohook_hess_3_2* clone() const
        {
            return new neohook_hess_3_2(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_3_3 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  5.0000000000000000e-01*var[0]*pow(var[4],2.0000000000000000e+00)*pow(var[2],2.0000000000000000e+00)+var[1]+( 5.0000000000000000e-01*var[0]+var[1])/(var[3]*var[3]);
        }

    public:

        neohook_hess_3_3()
        {
            this->_name = "neohook_hess_3_3";
        }

        std::string csrc() const
        {
            return " 5.0000000000000000e-01*var[0]*pow(var[4],2.0000000000000000e+00)*pow(var[2],2.0000000000000000e+00)+var[1]+( 5.0000000000000000e-01*var[0]+var[1])/(var[3]*var[3])";
        }

        std::string sym() const
        {
            return "(mu+(0.5)*lambda)*lambda2^(-2)+mu+(0.5)*lambda3^(2.0)*lambda1^(2.0)*lambda";
        }

        std::string latex() const
        {
            return "\\mu+\\frac{{(0.5)} \\lambda+\\mu}{lambda2^{2}}+{(0.5)}  \\lambda lambda3^{{(2.0)}} lambda1^{{(2.0)}}";
        }

        neohook_hess_3_3* clone() const
        {
            return new neohook_hess_3_3(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_3_4 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return var[4]*pow(var[2],2.0000000000000000e+00)*var[3]*var[0];
        }

    public:

        neohook_hess_3_4()
        {
            this->_name = "neohook_hess_3_4";
        }

        std::string csrc() const
        {
            return "var[4]*pow(var[2],2.0000000000000000e+00)*var[3]*var[0]";
        }

        std::string sym() const
        {
            return "lambda2*lambda*lambda3*lambda1^(2.0)";
        }

        std::string latex() const
        {
            return " lambda1^{{(2.0)}} lambda2 \\lambda lambda3";
        }

        neohook_hess_3_4* clone() const
        {
            return new neohook_hess_3_4(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_3_5 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        neohook_hess_3_5()
        {
            this->_name = "neohook_hess_3_5";
        }

        std::string csrc() const
        {
            return "0.0";
        }

        std::string sym() const
        {
            return "0";
        }

        std::string latex() const
        {
            return "0";
        }

        neohook_hess_3_5* clone() const
        {
            return new neohook_hess_3_5(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_3_6 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        neohook_hess_3_6()
        {
            this->_name = "neohook_hess_3_6";
        }

        std::string csrc() const
        {
            return "0.0";
        }

        std::string sym() const
        {
            return "0";
        }

        std::string latex() const
        {
            return "0";
        }

        neohook_hess_3_6* clone() const
        {
            return new neohook_hess_3_6(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_4_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  -5.0000000000000000e-01*1.0/(var[4])+5.0000000000000000e-01*var[4]*pow(var[2],2.0000000000000000e+00)*pow(var[3],2.0000000000000000e+00);
        }

    public:

        neohook_hess_4_0()
        {
            this->_name = "neohook_hess_4_0";
        }

        std::string csrc() const
        {
            return " -5.0000000000000000e-01*1.0/(var[4])+5.0000000000000000e-01*var[4]*pow(var[2],2.0000000000000000e+00)*pow(var[3],2.0000000000000000e+00)";
        }

        std::string sym() const
        {
            return "(0.5)*lambda1^(2.0)*lambda2^(2.0)*lambda3-(0.5)*lambda3^(-1)";
        }

        std::string latex() const
        {
            return "-{(0.5)} \\frac{1}{lambda3}+{(0.5)}  lambda3 lambda1^{{(2.0)}} lambda2^{{(2.0)}}";
        }

        neohook_hess_4_0* clone() const
        {
            return new neohook_hess_4_0(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_4_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  var[4]-1.0/(var[4]);
        }

    public:

        neohook_hess_4_1()
        {
            this->_name = "neohook_hess_4_1";
        }

        std::string csrc() const
        {
            return " var[4]-1.0/(var[4])";
        }

        std::string sym() const
        {
            return "lambda3-lambda3^(-1)";
        }

        std::string latex() const
        {
            return "lambda3-\\frac{1}{lambda3}";
        }

        neohook_hess_4_1* clone() const
        {
            return new neohook_hess_4_1(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_4_2 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return var[4]*var[2]*pow(var[3],2.0000000000000000e+00)*var[0];
        }

    public:

        neohook_hess_4_2()
        {
            this->_name = "neohook_hess_4_2";
        }

        std::string csrc() const
        {
            return "var[4]*var[2]*pow(var[3],2.0000000000000000e+00)*var[0]";
        }

        std::string sym() const
        {
            return "lambda2^(2.0)*lambda*lambda3*lambda1";
        }

        std::string latex() const
        {
            return " lambda3 lambda1 lambda2^{{(2.0)}} \\lambda";
        }

        neohook_hess_4_2* clone() const
        {
            return new neohook_hess_4_2(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_4_3 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return var[3]*var[0]*var[4]*pow(var[2],2.0000000000000000e+00);
        }

    public:

        neohook_hess_4_3()
        {
            this->_name = "neohook_hess_4_3";
        }

        std::string csrc() const
        {
            return "var[3]*var[0]*var[4]*pow(var[2],2.0000000000000000e+00)";
        }

        std::string sym() const
        {
            return "lambda1^(2.0)*lambda2*lambda*lambda3";
        }

        std::string latex() const
        {
            return " lambda2 \\lambda lambda3 lambda1^{{(2.0)}}";
        }

        neohook_hess_4_3* clone() const
        {
            return new neohook_hess_4_3(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_4_4 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  var[1]+5.0000000000000000e-01*pow(var[2],2.0000000000000000e+00)*pow(var[3],2.0000000000000000e+00)*var[0]+1.0/(var[4]*var[4])*( var[1]+5.0000000000000000e-01*var[0]);
        }

    public:

        neohook_hess_4_4()
        {
            this->_name = "neohook_hess_4_4";
        }

        std::string csrc() const
        {
            return " var[1]+5.0000000000000000e-01*pow(var[2],2.0000000000000000e+00)*pow(var[3],2.0000000000000000e+00)*var[0]+1.0/(var[4]*var[4])*( var[1]+5.0000000000000000e-01*var[0])";
        }

        std::string sym() const
        {
            return "lambda3^(-2)*((0.5)*lambda+mu)+(0.5)*lambda2^(2.0)*lambda*lambda1^(2.0)+mu";
        }

        std::string latex() const
        {
            return "{(0.5)}  lambda1^{{(2.0)}} lambda2^{{(2.0)}} \\lambda+\\frac{{(0.5)} \\lambda+\\mu}{lambda3^{2}}+\\mu";
        }

        neohook_hess_4_4* clone() const
        {
            return new neohook_hess_4_4(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_4_5 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        neohook_hess_4_5()
        {
            this->_name = "neohook_hess_4_5";
        }

        std::string csrc() const
        {
            return "0.0";
        }

        std::string sym() const
        {
            return "0";
        }

        std::string latex() const
        {
            return "0";
        }

        neohook_hess_4_5* clone() const
        {
            return new neohook_hess_4_5(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_4_6 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        neohook_hess_4_6()
        {
            this->_name = "neohook_hess_4_6";
        }

        std::string csrc() const
        {
            return "0.0";
        }

        std::string sym() const
        {
            return "0";
        }

        std::string latex() const
        {
            return "0";
        }

        neohook_hess_4_6* clone() const
        {
            return new neohook_hess_4_6(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_5_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        neohook_hess_5_0()
        {
            this->_name = "neohook_hess_5_0";
        }

        std::string csrc() const
        {
            return "0.0";
        }

        std::string sym() const
        {
            return "0";
        }

        std::string latex() const
        {
            return "0";
        }

        neohook_hess_5_0* clone() const
        {
            return new neohook_hess_5_0(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_5_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        neohook_hess_5_1()
        {
            this->_name = "neohook_hess_5_1";
        }

        std::string csrc() const
        {
            return "0.0";
        }

        std::string sym() const
        {
            return "0";
        }

        std::string latex() const
        {
            return "0";
        }

        neohook_hess_5_1* clone() const
        {
            return new neohook_hess_5_1(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_5_2 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        neohook_hess_5_2()
        {
            this->_name = "neohook_hess_5_2";
        }

        std::string csrc() const
        {
            return "0.0";
        }

        std::string sym() const
        {
            return "0";
        }

        std::string latex() const
        {
            return "0";
        }

        neohook_hess_5_2* clone() const
        {
            return new neohook_hess_5_2(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_5_3 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        neohook_hess_5_3()
        {
            this->_name = "neohook_hess_5_3";
        }

        std::string csrc() const
        {
            return "0.0";
        }

        std::string sym() const
        {
            return "0";
        }

        std::string latex() const
        {
            return "0";
        }

        neohook_hess_5_3* clone() const
        {
            return new neohook_hess_5_3(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_5_4 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        neohook_hess_5_4()
        {
            this->_name = "neohook_hess_5_4";
        }

        std::string csrc() const
        {
            return "0.0";
        }

        std::string sym() const
        {
            return "0";
        }

        std::string latex() const
        {
            return "0";
        }

        neohook_hess_5_4* clone() const
        {
            return new neohook_hess_5_4(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_5_5 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        neohook_hess_5_5()
        {
            this->_name = "neohook_hess_5_5";
        }

        std::string csrc() const
        {
            return "0.0";
        }

        std::string sym() const
        {
            return "0";
        }

        std::string latex() const
        {
            return "0";
        }

        neohook_hess_5_5* clone() const
        {
            return new neohook_hess_5_5(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_5_6 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return var[6];
        }

    public:

        neohook_hess_5_6()
        {
            this->_name = "neohook_hess_5_6";
        }

        std::string csrc() const
        {
            return "var[6]";
        }

        std::string sym() const
        {
            return "alpha";
        }

        std::string latex() const
        {
            return "\\alpha";
        }

        neohook_hess_5_6* clone() const
        {
            return new neohook_hess_5_6(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_6_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        neohook_hess_6_0()
        {
            this->_name = "neohook_hess_6_0";
        }

        std::string csrc() const
        {
            return "0.0";
        }

        std::string sym() const
        {
            return "0";
        }

        std::string latex() const
        {
            return "0";
        }

        neohook_hess_6_0* clone() const
        {
            return new neohook_hess_6_0(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_6_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        neohook_hess_6_1()
        {
            this->_name = "neohook_hess_6_1";
        }

        std::string csrc() const
        {
            return "0.0";
        }

        std::string sym() const
        {
            return "0";
        }

        std::string latex() const
        {
            return "0";
        }

        neohook_hess_6_1* clone() const
        {
            return new neohook_hess_6_1(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_6_2 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        neohook_hess_6_2()
        {
            this->_name = "neohook_hess_6_2";
        }

        std::string csrc() const
        {
            return "0.0";
        }

        std::string sym() const
        {
            return "0";
        }

        std::string latex() const
        {
            return "0";
        }

        neohook_hess_6_2* clone() const
        {
            return new neohook_hess_6_2(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_6_3 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        neohook_hess_6_3()
        {
            this->_name = "neohook_hess_6_3";
        }

        std::string csrc() const
        {
            return "0.0";
        }

        std::string sym() const
        {
            return "0";
        }

        std::string latex() const
        {
            return "0";
        }

        neohook_hess_6_3* clone() const
        {
            return new neohook_hess_6_3(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_6_4 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        neohook_hess_6_4()
        {
            this->_name = "neohook_hess_6_4";
        }

        std::string csrc() const
        {
            return "0.0";
        }

        std::string sym() const
        {
            return "0";
        }

        std::string latex() const
        {
            return "0";
        }

        neohook_hess_6_4* clone() const
        {
            return new neohook_hess_6_4(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_6_5 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return var[6];
        }

    public:

        neohook_hess_6_5()
        {
            this->_name = "neohook_hess_6_5";
        }

        std::string csrc() const
        {
            return "var[6]";
        }

        std::string sym() const
        {
            return "alpha";
        }

        std::string latex() const
        {
            return "\\alpha";
        }

        neohook_hess_6_5* clone() const
        {
            return new neohook_hess_6_5(*this);
        }
    };

    template< class VarContainer>
    class neohook_hess_6_6 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return var[5];
        }

    public:

        neohook_hess_6_6()
        {
            this->_name = "neohook_hess_6_6";
        }

        std::string csrc() const
        {
            return "var[5]";
        }

        std::string sym() const
        {
            return "K";
        }

        std::string latex() const
        {
            return "K";
        }

        neohook_hess_6_6* clone() const
        {
            return new neohook_hess_6_6(*this);
        }
    };

    template<class VarContainer>
    class neohook : public PFuncBase< VarContainer, double>
    {
    public:
        
        typedef typename PFuncBase< VarContainer, double>::size_type size_type;

        PSimpleBase< VarContainer, double> *_val;
        PSimpleBase< VarContainer, double> **_grad_val;
        PSimpleBase< VarContainer, double> ***_hess_val;
        
        neohook()
        {
            construct();
        }

        neohook(const neohook &RHS )
        {
            construct(false);
            
            _val = RHS._val->clone();
            _grad_val[0] = RHS._grad_val[0]->clone();
            _grad_val[1] = RHS._grad_val[1]->clone();
            _grad_val[2] = RHS._grad_val[2]->clone();
            _grad_val[3] = RHS._grad_val[3]->clone();
            _grad_val[4] = RHS._grad_val[4]->clone();
            _grad_val[5] = RHS._grad_val[5]->clone();
            _grad_val[6] = RHS._grad_val[6]->clone();
            _hess_val[0][0] = RHS._hess_val[0][0]->clone();
            _hess_val[0][1] = RHS._hess_val[0][1]->clone();
            _hess_val[0][2] = RHS._hess_val[0][2]->clone();
            _hess_val[0][3] = RHS._hess_val[0][3]->clone();
            _hess_val[0][4] = RHS._hess_val[0][4]->clone();
            _hess_val[0][5] = RHS._hess_val[0][5]->clone();
            _hess_val[0][6] = RHS._hess_val[0][6]->clone();
            _hess_val[1][0] = RHS._hess_val[1][0]->clone();
            _hess_val[1][1] = RHS._hess_val[1][1]->clone();
            _hess_val[1][2] = RHS._hess_val[1][2]->clone();
            _hess_val[1][3] = RHS._hess_val[1][3]->clone();
            _hess_val[1][4] = RHS._hess_val[1][4]->clone();
            _hess_val[1][5] = RHS._hess_val[1][5]->clone();
            _hess_val[1][6] = RHS._hess_val[1][6]->clone();
            _hess_val[2][0] = RHS._hess_val[2][0]->clone();
            _hess_val[2][1] = RHS._hess_val[2][1]->clone();
            _hess_val[2][2] = RHS._hess_val[2][2]->clone();
            _hess_val[2][3] = RHS._hess_val[2][3]->clone();
            _hess_val[2][4] = RHS._hess_val[2][4]->clone();
            _hess_val[2][5] = RHS._hess_val[2][5]->clone();
            _hess_val[2][6] = RHS._hess_val[2][6]->clone();
            _hess_val[3][0] = RHS._hess_val[3][0]->clone();
            _hess_val[3][1] = RHS._hess_val[3][1]->clone();
            _hess_val[3][2] = RHS._hess_val[3][2]->clone();
            _hess_val[3][3] = RHS._hess_val[3][3]->clone();
            _hess_val[3][4] = RHS._hess_val[3][4]->clone();
            _hess_val[3][5] = RHS._hess_val[3][5]->clone();
            _hess_val[3][6] = RHS._hess_val[3][6]->clone();
            _hess_val[4][0] = RHS._hess_val[4][0]->clone();
            _hess_val[4][1] = RHS._hess_val[4][1]->clone();
            _hess_val[4][2] = RHS._hess_val[4][2]->clone();
            _hess_val[4][3] = RHS._hess_val[4][3]->clone();
            _hess_val[4][4] = RHS._hess_val[4][4]->clone();
            _hess_val[4][5] = RHS._hess_val[4][5]->clone();
            _hess_val[4][6] = RHS._hess_val[4][6]->clone();
            _hess_val[5][0] = RHS._hess_val[5][0]->clone();
            _hess_val[5][1] = RHS._hess_val[5][1]->clone();
            _hess_val[5][2] = RHS._hess_val[5][2]->clone();
            _hess_val[5][3] = RHS._hess_val[5][3]->clone();
            _hess_val[5][4] = RHS._hess_val[5][4]->clone();
            _hess_val[5][5] = RHS._hess_val[5][5]->clone();
            _hess_val[5][6] = RHS._hess_val[5][6]->clone();
            _hess_val[6][0] = RHS._hess_val[6][0]->clone();
            _hess_val[6][1] = RHS._hess_val[6][1]->clone();
            _hess_val[6][2] = RHS._hess_val[6][2]->clone();
            _hess_val[6][3] = RHS._hess_val[6][3]->clone();
            _hess_val[6][4] = RHS._hess_val[6][4]->clone();
            _hess_val[6][5] = RHS._hess_val[6][5]->clone();
            _hess_val[6][6] = RHS._hess_val[6][6]->clone();
            
        }

        neohook& operator=( neohook RHS )
        {
            using std::swap;
            
            swap(_val, RHS._val);
            swap(_grad_val[0], RHS._grad_val[0]);
            swap(_grad_val[1], RHS._grad_val[1]);
            swap(_grad_val[2], RHS._grad_val[2]);
            swap(_grad_val[3], RHS._grad_val[3]);
            swap(_grad_val[4], RHS._grad_val[4]);
            swap(_grad_val[5], RHS._grad_val[5]);
            swap(_grad_val[6], RHS._grad_val[6]);
            swap(_hess_val[0][0], RHS._hess_val[0][0]);
            swap(_hess_val[0][1], RHS._hess_val[0][1]);
            swap(_hess_val[0][2], RHS._hess_val[0][2]);
            swap(_hess_val[0][3], RHS._hess_val[0][3]);
            swap(_hess_val[0][4], RHS._hess_val[0][4]);
            swap(_hess_val[0][5], RHS._hess_val[0][5]);
            swap(_hess_val[0][6], RHS._hess_val[0][6]);
            swap(_hess_val[1][0], RHS._hess_val[1][0]);
            swap(_hess_val[1][1], RHS._hess_val[1][1]);
            swap(_hess_val[1][2], RHS._hess_val[1][2]);
            swap(_hess_val[1][3], RHS._hess_val[1][3]);
            swap(_hess_val[1][4], RHS._hess_val[1][4]);
            swap(_hess_val[1][5], RHS._hess_val[1][5]);
            swap(_hess_val[1][6], RHS._hess_val[1][6]);
            swap(_hess_val[2][0], RHS._hess_val[2][0]);
            swap(_hess_val[2][1], RHS._hess_val[2][1]);
            swap(_hess_val[2][2], RHS._hess_val[2][2]);
            swap(_hess_val[2][3], RHS._hess_val[2][3]);
            swap(_hess_val[2][4], RHS._hess_val[2][4]);
            swap(_hess_val[2][5], RHS._hess_val[2][5]);
            swap(_hess_val[2][6], RHS._hess_val[2][6]);
            swap(_hess_val[3][0], RHS._hess_val[3][0]);
            swap(_hess_val[3][1], RHS._hess_val[3][1]);
            swap(_hess_val[3][2], RHS._hess_val[3][2]);
            swap(_hess_val[3][3], RHS._hess_val[3][3]);
            swap(_hess_val[3][4], RHS._hess_val[3][4]);
            swap(_hess_val[3][5], RHS._hess_val[3][5]);
            swap(_hess_val[3][6], RHS._hess_val[3][6]);
            swap(_hess_val[4][0], RHS._hess_val[4][0]);
            swap(_hess_val[4][1], RHS._hess_val[4][1]);
            swap(_hess_val[4][2], RHS._hess_val[4][2]);
            swap(_hess_val[4][3], RHS._hess_val[4][3]);
            swap(_hess_val[4][4], RHS._hess_val[4][4]);
            swap(_hess_val[4][5], RHS._hess_val[4][5]);
            swap(_hess_val[4][6], RHS._hess_val[4][6]);
            swap(_hess_val[5][0], RHS._hess_val[5][0]);
            swap(_hess_val[5][1], RHS._hess_val[5][1]);
            swap(_hess_val[5][2], RHS._hess_val[5][2]);
            swap(_hess_val[5][3], RHS._hess_val[5][3]);
            swap(_hess_val[5][4], RHS._hess_val[5][4]);
            swap(_hess_val[5][5], RHS._hess_val[5][5]);
            swap(_hess_val[5][6], RHS._hess_val[5][6]);
            swap(_hess_val[6][0], RHS._hess_val[6][0]);
            swap(_hess_val[6][1], RHS._hess_val[6][1]);
            swap(_hess_val[6][2], RHS._hess_val[6][2]);
            swap(_hess_val[6][3], RHS._hess_val[6][3]);
            swap(_hess_val[6][4], RHS._hess_val[6][4]);
            swap(_hess_val[6][5], RHS._hess_val[6][5]);
            swap(_hess_val[6][6], RHS._hess_val[6][6]);
            
            return *this;
        }

        ~neohook()
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
        }

        neohook<VarContainer>* clone() const
        {
            return new neohook<VarContainer>(*this);
        }

        PSimpleFunction< VarContainer, double> simplefunction() const
        {
            return PSimpleFunction< VarContainer, double>( *_val );
        }

        PSimpleFunction< VarContainer, double> grad_simplefunction(size_type di) const
        {
            return PSimpleFunction< VarContainer, double>( *_grad_val[di] );
        }

        PSimpleFunction< VarContainer, double> hess_simplefunction(size_type di, size_type dj) const
        {
            return PSimpleFunction< VarContainer, double>( *_hess_val[di][dj] );
        }

        double operator()(const VarContainer &var)
        {
            return (*_val)(var);
        }

        double grad(const VarContainer &var, size_type di)
        {
            return (*_grad_val[di])(var);
        }

        double hess(const VarContainer &var, size_type di, size_type dj)
        {
            return (*_hess_val[di][dj])(var);
        }

        void eval(const VarContainer &var)
        {
            (*_val)(var);
        }

        void eval_grad(const VarContainer &var)
        {
            (*_grad_val[0])(var);
            (*_grad_val[1])(var);
            (*_grad_val[2])(var);
            (*_grad_val[3])(var);
            (*_grad_val[4])(var);
            (*_grad_val[5])(var);
            (*_grad_val[6])(var);
        }

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
        }

        double operator()() const
        {
            return (*_val)();
        }

        double grad(size_type di) const
        {
            return (*_grad_val[di])();
        }

        double hess(size_type di, size_type dj) const
        {
            return (*_hess_val[di][dj])();
        }

    private:
        void construct(bool allocate = true)
        {
            this->_name = "neohook";
            this->_var_name.clear();
            this->_var_name.push_back("lambda");
            this->_var_name.push_back("mu");
            this->_var_name.push_back("lambda1");
            this->_var_name.push_back("lambda2");
            this->_var_name.push_back("lambda3");
            this->_var_name.push_back("K");
            this->_var_name.push_back("alpha");
            this->_var_description.clear();
            this->_var_description.push_back("First Lame parameter");
            this->_var_description.push_back("Second Lame parameter");
            this->_var_description.push_back("First principle stretch");
            this->_var_description.push_back("Second principle stretch");
            this->_var_description.push_back("Third principle stretch");
            this->_var_description.push_back("Strain hardening coefficient");
            this->_var_description.push_back("Equivalent plastic strain");
            
            _grad_val = new PSimpleBase< VarContainer, double>*[7];
            
            _hess_val = new PSimpleBase< VarContainer, double>**[7];
            _hess_val[0] = new PSimpleBase< VarContainer, double>*[7];
            _hess_val[1] = new PSimpleBase< VarContainer, double>*[7];
            _hess_val[2] = new PSimpleBase< VarContainer, double>*[7];
            _hess_val[3] = new PSimpleBase< VarContainer, double>*[7];
            _hess_val[4] = new PSimpleBase< VarContainer, double>*[7];
            _hess_val[5] = new PSimpleBase< VarContainer, double>*[7];
            _hess_val[6] = new PSimpleBase< VarContainer, double>*[7];
            
            if(!allocate) return;
            
            _val = new neohook_f<VarContainer>();
            
            _grad_val[0] = new neohook_grad_0<VarContainer>();
            _grad_val[1] = new neohook_grad_1<VarContainer>();
            _grad_val[2] = new neohook_grad_2<VarContainer>();
            _grad_val[3] = new neohook_grad_3<VarContainer>();
            _grad_val[4] = new neohook_grad_4<VarContainer>();
            _grad_val[5] = new neohook_grad_5<VarContainer>();
            _grad_val[6] = new neohook_grad_6<VarContainer>();
            
            _hess_val[0][0] = new neohook_hess_0_0<VarContainer>();
            _hess_val[0][1] = new neohook_hess_0_1<VarContainer>();
            _hess_val[0][2] = new neohook_hess_0_2<VarContainer>();
            _hess_val[0][3] = new neohook_hess_0_3<VarContainer>();
            _hess_val[0][4] = new neohook_hess_0_4<VarContainer>();
            _hess_val[0][5] = new neohook_hess_0_5<VarContainer>();
            _hess_val[0][6] = new neohook_hess_0_6<VarContainer>();
            _hess_val[1][0] = new neohook_hess_1_0<VarContainer>();
            _hess_val[1][1] = new neohook_hess_1_1<VarContainer>();
            _hess_val[1][2] = new neohook_hess_1_2<VarContainer>();
            _hess_val[1][3] = new neohook_hess_1_3<VarContainer>();
            _hess_val[1][4] = new neohook_hess_1_4<VarContainer>();
            _hess_val[1][5] = new neohook_hess_1_5<VarContainer>();
            _hess_val[1][6] = new neohook_hess_1_6<VarContainer>();
            _hess_val[2][0] = new neohook_hess_2_0<VarContainer>();
            _hess_val[2][1] = new neohook_hess_2_1<VarContainer>();
            _hess_val[2][2] = new neohook_hess_2_2<VarContainer>();
            _hess_val[2][3] = new neohook_hess_2_3<VarContainer>();
            _hess_val[2][4] = new neohook_hess_2_4<VarContainer>();
            _hess_val[2][5] = new neohook_hess_2_5<VarContainer>();
            _hess_val[2][6] = new neohook_hess_2_6<VarContainer>();
            _hess_val[3][0] = new neohook_hess_3_0<VarContainer>();
            _hess_val[3][1] = new neohook_hess_3_1<VarContainer>();
            _hess_val[3][2] = new neohook_hess_3_2<VarContainer>();
            _hess_val[3][3] = new neohook_hess_3_3<VarContainer>();
            _hess_val[3][4] = new neohook_hess_3_4<VarContainer>();
            _hess_val[3][5] = new neohook_hess_3_5<VarContainer>();
            _hess_val[3][6] = new neohook_hess_3_6<VarContainer>();
            _hess_val[4][0] = new neohook_hess_4_0<VarContainer>();
            _hess_val[4][1] = new neohook_hess_4_1<VarContainer>();
            _hess_val[4][2] = new neohook_hess_4_2<VarContainer>();
            _hess_val[4][3] = new neohook_hess_4_3<VarContainer>();
            _hess_val[4][4] = new neohook_hess_4_4<VarContainer>();
            _hess_val[4][5] = new neohook_hess_4_5<VarContainer>();
            _hess_val[4][6] = new neohook_hess_4_6<VarContainer>();
            _hess_val[5][0] = new neohook_hess_5_0<VarContainer>();
            _hess_val[5][1] = new neohook_hess_5_1<VarContainer>();
            _hess_val[5][2] = new neohook_hess_5_2<VarContainer>();
            _hess_val[5][3] = new neohook_hess_5_3<VarContainer>();
            _hess_val[5][4] = new neohook_hess_5_4<VarContainer>();
            _hess_val[5][5] = new neohook_hess_5_5<VarContainer>();
            _hess_val[5][6] = new neohook_hess_5_6<VarContainer>();
            _hess_val[6][0] = new neohook_hess_6_0<VarContainer>();
            _hess_val[6][1] = new neohook_hess_6_1<VarContainer>();
            _hess_val[6][2] = new neohook_hess_6_2<VarContainer>();
            _hess_val[6][3] = new neohook_hess_6_3<VarContainer>();
            _hess_val[6][4] = new neohook_hess_6_4<VarContainer>();
            _hess_val[6][5] = new neohook_hess_6_5<VarContainer>();
            _hess_val[6][6] = new neohook_hess_6_6<VarContainer>();
        }

    };


}
#endif
