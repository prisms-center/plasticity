// created: 2016-8-5 14:00:36
// version: master
// url: https://github.com/prisms-center/IntegrationToolsWriter.git
// commit: 8a15adf67355fad30bd75ce9ba6b1f8d24b9a537

#ifndef Cu_hardening_HH
#define Cu_hardening_HH

#include <cmath>
#include <cstdlib>
#include "../../../../utils/IntegrationTools/PFunction.hh"

namespace PRISMS
{
    template< class VarContainer>
    class Cu_hardening_f : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  3.2986102779887700e+08*pow(var[0],3.3333333333333331e-01)+3.2738656842295301e+08*var[0]+-1.1216965746072850e+09*pow(var[0],5.0000000000000000e-01);
        }

    public:

        Cu_hardening_f()
        {
            this->_name = "Cu_hardening_f";
        }

        std::string csrc() const
        {
            return " 3.2986102779887700e+08*pow(var[0],3.3333333333333331e-01)+3.2738656842295301e+08*var[0]+-1.1216965746072850e+09*pow(var[0],5.0000000000000000e-01)";
        }

        std::string sym() const
        {
            return "(3.2738656842295301E8)*alpha+(3.29861027798877E8)*alpha^(0.33333333333333333334)-(1.121696574607285E9)*sqrt(alpha)";
        }

        std::string latex() const
        {
            return "-{(1.121696574607285E9)} \\sqrt{\\alpha}+{(3.29861027798877E8)} \\alpha^{{(0.33333333333333333334)}}+{(3.2738656842295301E8)} \\alpha";
        }

        Cu_hardening_f* clone() const
        {
            return new Cu_hardening_f(*this);
        }
    };

    template< class VarContainer>
    class Cu_hardening_grad_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  -5.6084828730364251e+08*pow(var[0],-5.0000000000000000e-01)+1.0995367593295901e+08*pow(var[0],-6.6666666666666663e-01)+3.2738656842295301e+08;
        }

    public:

        Cu_hardening_grad_0()
        {
            this->_name = "Cu_hardening_grad_0";
        }

        std::string csrc() const
        {
            return " -5.6084828730364251e+08*pow(var[0],-5.0000000000000000e-01)+1.0995367593295901e+08*pow(var[0],-6.6666666666666663e-01)+3.2738656842295301e+08";
        }

        std::string sym() const
        {
            return "3.2738656842295301E8+(1.09953675932959E8)*alpha^(-0.66666666666666666663)-(5.608482873036425E8)*alpha^(-0.5)";
        }

        std::string latex() const
        {
            return "3.2738656842295301E8+{(1.09953675932959E8)} \\frac{1}{\\alpha^{{(0.66666666666666666663)}}}-{(5.608482873036425E8)} \\frac{1}{\\sqrt{\\alpha}}";
        }

        Cu_hardening_grad_0* clone() const
        {
            return new Cu_hardening_grad_0(*this);
        }
    };

    template< class VarContainer>
    class Cu_hardening_grad_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        Cu_hardening_grad_1()
        {
            this->_name = "Cu_hardening_grad_1";
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

        Cu_hardening_grad_1* clone() const
        {
            return new Cu_hardening_grad_1(*this);
        }
    };

    template< class VarContainer>
    class Cu_hardening_hess_0_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  2.8042414365182126e+08*pow(var[0],-1.5000000000000000e+00)+-7.3302450621972665e+07*pow(var[0],-1.6666666666666667e+00);
        }

    public:

        Cu_hardening_hess_0_0()
        {
            this->_name = "Cu_hardening_hess_0_0";
        }

        std::string csrc() const
        {
            return " 2.8042414365182126e+08*pow(var[0],-1.5000000000000000e+00)+-7.3302450621972665e+07*pow(var[0],-1.6666666666666667e+00)";
        }

        std::string sym() const
        {
            return "(2.8042414365182125E8)*alpha^(-1.5)-(7.330245062197266666E7)*alpha^(-1.6666666666666666666)";
        }

        std::string latex() const
        {
            return "-{(7.330245062197266666E7)} \\frac{1}{\\alpha^{{(1.6666666666666666666)}}}+{(2.8042414365182125E8)} \\frac{1}{\\alpha^{{(1.5)}}}";
        }

        Cu_hardening_hess_0_0* clone() const
        {
            return new Cu_hardening_hess_0_0(*this);
        }
    };

    template< class VarContainer>
    class Cu_hardening_hess_0_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        Cu_hardening_hess_0_1()
        {
            this->_name = "Cu_hardening_hess_0_1";
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

        Cu_hardening_hess_0_1* clone() const
        {
            return new Cu_hardening_hess_0_1(*this);
        }
    };

    template< class VarContainer>
    class Cu_hardening_hess_1_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        Cu_hardening_hess_1_0()
        {
            this->_name = "Cu_hardening_hess_1_0";
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

        Cu_hardening_hess_1_0* clone() const
        {
            return new Cu_hardening_hess_1_0(*this);
        }
    };

    template< class VarContainer>
    class Cu_hardening_hess_1_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        Cu_hardening_hess_1_1()
        {
            this->_name = "Cu_hardening_hess_1_1";
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

        Cu_hardening_hess_1_1* clone() const
        {
            return new Cu_hardening_hess_1_1(*this);
        }
    };

    template<class VarContainer>
    class Cu_hardening : public PFuncBase< VarContainer, double>
    {
    public:
        
        typedef typename PFuncBase< VarContainer, double>::size_type size_type;

        PSimpleBase< VarContainer, double> *_val;
        PSimpleBase< VarContainer, double> **_grad_val;
        PSimpleBase< VarContainer, double> ***_hess_val;
        
        Cu_hardening()
        {
            construct();
        }

        Cu_hardening(const Cu_hardening &RHS )
        {
            construct(false);
            
            _val = RHS._val->clone();
            _grad_val[0] = RHS._grad_val[0]->clone();
            _grad_val[1] = RHS._grad_val[1]->clone();
            _hess_val[0][0] = RHS._hess_val[0][0]->clone();
            _hess_val[0][1] = RHS._hess_val[0][1]->clone();
            _hess_val[1][0] = RHS._hess_val[1][0]->clone();
            _hess_val[1][1] = RHS._hess_val[1][1]->clone();
            
        }

        Cu_hardening& operator=( Cu_hardening RHS )
        {
            using std::swap;
            
            swap(_val, RHS._val);
            swap(_grad_val[0], RHS._grad_val[0]);
            swap(_grad_val[1], RHS._grad_val[1]);
            swap(_hess_val[0][0], RHS._hess_val[0][0]);
            swap(_hess_val[0][1], RHS._hess_val[0][1]);
            swap(_hess_val[1][0], RHS._hess_val[1][0]);
            swap(_hess_val[1][1], RHS._hess_val[1][1]);
            
            return *this;
        }

        ~Cu_hardening()
        {
            delete _val;

            delete _grad_val[0];
            delete _grad_val[1];
            delete [] _grad_val;

            delete _hess_val[0][0];
            delete _hess_val[0][1];
            delete _hess_val[1][0];
            delete _hess_val[1][1];
            delete [] _hess_val[0];
            delete [] _hess_val[1];
            delete [] _hess_val;
        }

        Cu_hardening<VarContainer>* clone() const
        {
            return new Cu_hardening<VarContainer>(*this);
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
        }

        void eval_hess(const VarContainer &var)
        {
            (*_hess_val[0][0])(var);
            (*_hess_val[0][1])(var);
            (*_hess_val[1][0])(var);
            (*_hess_val[1][1])(var);
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
            this->_name = "Cu_hardening";
            this->_var_name.clear();
            this->_var_name.push_back("alpha");
            this->_var_name.push_back("K");
            this->_var_description.clear();
            this->_var_description.push_back("Equivalent plastic strain");
            this->_var_description.push_back("Hardening parameter");
            
            _grad_val = new PSimpleBase< VarContainer, double>*[2];
            
            _hess_val = new PSimpleBase< VarContainer, double>**[2];
            _hess_val[0] = new PSimpleBase< VarContainer, double>*[2];
            _hess_val[1] = new PSimpleBase< VarContainer, double>*[2];
            
            if(!allocate) return;
            
            _val = new Cu_hardening_f<VarContainer>();
            
            _grad_val[0] = new Cu_hardening_grad_0<VarContainer>();
            _grad_val[1] = new Cu_hardening_grad_1<VarContainer>();
            
            _hess_val[0][0] = new Cu_hardening_hess_0_0<VarContainer>();
            _hess_val[0][1] = new Cu_hardening_hess_0_1<VarContainer>();
            _hess_val[1][0] = new Cu_hardening_hess_1_0<VarContainer>();
            _hess_val[1][1] = new Cu_hardening_hess_1_1<VarContainer>();
        }

    };


}
#endif
