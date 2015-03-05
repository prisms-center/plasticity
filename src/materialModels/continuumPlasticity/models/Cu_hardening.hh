// created: 2015-3-5 16:54:09
// version: 0.1.0
// url: git@github.com:prisms-center/IntegrationTools.git
// commit: 947b873eb0296ab1631408af59ecf7768231ff63

#ifndef Cu_hardening_HH
#define Cu_hardening_HH

#include <cmath>
#include <cstdlib>
#include "PFunction.hh"

namespace PRISMS
{
    template< class VarContainer>
    class Cu_hardening_f : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  3.2986102779887700e+08*pow(var[0],3.3333333333333331e-01)+3.2738656842295301e+08*var[0]+-1.1216965746072850e+09*pow(var[0],5.0000000000000000e-01);
        };

    public:

        Cu_hardening_f()
        {
            this->_name = "Cu_hardening_f";
        };

        Cu_hardening_f* clone() const
        {
            return new Cu_hardening_f(*this);
        };
    };

    template< class VarContainer>
    class Cu_hardening_grad_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  -5.6084828730364251e+08*pow(var[0],-5.0000000000000000e-01)+1.0995367593295901e+08*pow(var[0],-6.6666666666666663e-01)+3.2738656842295301e+08;
        };

    public:

        Cu_hardening_grad_0()
        {
            this->_name = "Cu_hardening_grad_0";
        };

        Cu_hardening_grad_0* clone() const
        {
            return new Cu_hardening_grad_0(*this);
        };
    };

    template< class VarContainer>
    class Cu_hardening_grad_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        Cu_hardening_grad_1()
        {
            this->_name = "Cu_hardening_grad_1";
        };

        Cu_hardening_grad_1* clone() const
        {
            return new Cu_hardening_grad_1(*this);
        };
    };

    template< class VarContainer>
    class Cu_hardening_hess_0_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  2.8042414365182126e+08*pow(var[0],-1.5000000000000000e+00)+-7.3302450621972665e+07*pow(var[0],-1.6666666666666667e+00);
        };

    public:

        Cu_hardening_hess_0_0()
        {
            this->_name = "Cu_hardening_hess_0_0";
        };

        Cu_hardening_hess_0_0* clone() const
        {
            return new Cu_hardening_hess_0_0(*this);
        };
    };

    template< class VarContainer>
    class Cu_hardening_hess_0_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        Cu_hardening_hess_0_1()
        {
            this->_name = "Cu_hardening_hess_0_1";
        };

        Cu_hardening_hess_0_1* clone() const
        {
            return new Cu_hardening_hess_0_1(*this);
        };
    };

    template< class VarContainer>
    class Cu_hardening_hess_1_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        Cu_hardening_hess_1_0()
        {
            this->_name = "Cu_hardening_hess_1_0";
        };

        Cu_hardening_hess_1_0* clone() const
        {
            return new Cu_hardening_hess_1_0(*this);
        };
    };

    template< class VarContainer>
    class Cu_hardening_hess_1_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        };

    public:

        Cu_hardening_hess_1_1()
        {
            this->_name = "Cu_hardening_hess_1_1";
        };

        Cu_hardening_hess_1_1* clone() const
        {
            return new Cu_hardening_hess_1_1(*this);
        };
    };

    template<class VarContainer>
    class Cu_hardening : public PFuncBase< VarContainer, double>
    {
    public:
        using PFuncBase< VarContainer, double>::_name;
        using PFuncBase< VarContainer, double>::_var_name;
        using PFuncBase< VarContainer, double>::_var_description;
        
        PSimpleBase< VarContainer, double> *_val;
        PSimpleBase< VarContainer, double> **_grad_val;
        PSimpleBase< VarContainer, double> ***_hess_val;
        
        Cu_hardening()
        {
            construct();
        }

        Cu_hardening(const Cu_hardening &RHS )
        {
            construct();
        }

        Cu_hardening& operator=(const Cu_hardening &RHS )
        {
            _val = RHS._val;
            
            _grad_val[0] = RHS._grad_val[0];
            _grad_val[1] = RHS._grad_val[1];
            _hess_val[0][0] = RHS._hess_val[0][0];
            _hess_val[0][1] = RHS._hess_val[0][1];
            _hess_val[1][0] = RHS._hess_val[1][0];
            _hess_val[1][1] = RHS._hess_val[1][1];
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
        };

        Cu_hardening<VarContainer>* clone() const
        {
            return new Cu_hardening<VarContainer>(*this);
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
        };

        void eval_hess(const VarContainer &var)
        {
            (*_hess_val[0][0])(var);
            (*_hess_val[0][1])(var);
            (*_hess_val[1][0])(var);
            (*_hess_val[1][1])(var);
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
            _name = "Cu_hardening";
            _var_name.clear();
            _var_name.push_back("alpha");
            _var_name.push_back("K");
            _var_description.clear();
            _var_description.push_back("Equivalent plastic strain");
            _var_description.push_back("Hardening parameter");
            
            _val = new Cu_hardening_f<VarContainer>();
            
            _grad_val = new PSimpleBase< VarContainer, double>*[2];
            _grad_val[0] = new Cu_hardening_grad_0<VarContainer>();
            _grad_val[1] = new Cu_hardening_grad_1<VarContainer>();
            
            _hess_val = new PSimpleBase< VarContainer, double>**[2];
            _hess_val[0] = new PSimpleBase< VarContainer, double>*[2];
            _hess_val[1] = new PSimpleBase< VarContainer, double>*[2];
            _hess_val[0][0] = new Cu_hardening_hess_0_0<VarContainer>();
            _hess_val[0][1] = new Cu_hardening_hess_0_1<VarContainer>();
            _hess_val[1][0] = new Cu_hardening_hess_1_0<VarContainer>();
            _hess_val[1][1] = new Cu_hardening_hess_1_1<VarContainer>();
        };

    };
}
#endif
