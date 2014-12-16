// created: 2014-7-29 9:30:37
// version: develop
// url: git@github.com:prisms-center/IntegrationTools.git
// commit: b9581986f73c383e5c32629a49d7077746004a54

#ifndef hardening_HH
#define hardening_HH

#include <cmath>
#include <cstdlib>
#include "PFunction.hh"

namespace PRISMS
{
    template< class VarContainer>
    class hardening_f : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  -1.1216965746072850e+09*pow(var[0],5.0000000000000000e-01)+3.2738656842295301e+08*var[0]+3.2986102779887700e+08*pow(var[0],3.3333333333333331e-01);
        }

    public:

        hardening_f()
        {
            this->_name = "hardening_f";
        }

        hardening_f* clone() const
        {
            return new hardening_f(*this);
        }
    };

    template< class VarContainer>
    class hardening_grad_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  -5.6084828730364251e+08*pow(var[0],-5.0000000000000000e-01)+1.0995367593295901e+08*pow(var[0],-6.6666666666666663e-01)+3.2738656842295301e+08;
        }

    public:

        hardening_grad_0()
        {
            this->_name = "hardening_grad_0";
        }

        hardening_grad_0* clone() const
        {
            return new hardening_grad_0(*this);
        }
    };

    template< class VarContainer>
    class hardening_hess_0_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  2.8042414365182126e+08*pow(var[0],-1.5000000000000000e+00)+-7.3302450621972665e+07*pow(var[0],-1.6666666666666667e+00);
        }

    public:

        hardening_hess_0_0()
        {
            this->_name = "hardening_hess_0_0";
        }

        hardening_hess_0_0* clone() const
        {
            return new hardening_hess_0_0(*this);
        }
    };

    template<class VarContainer>
    class hardening : public PFuncBase< VarContainer, double>
    {
    public:
        
        PSimpleBase< VarContainer, double> *_val;
        PSimpleBase< VarContainer, double> **_grad_val;
        PSimpleBase< VarContainer, double> ***_hess_val;
        
        hardening()
        {
            construct();
        }

        hardening(const hardening &RHS )
        {
            construct();
        }

        hardening& operator=(const hardening &RHS )
        {
            _val = RHS._val;
            
            _grad_val[0] = RHS._grad_val[0];
            _hess_val[0][0] = RHS._hess_val[0][0];
        }

        ~hardening()
        {
            delete _val;

            delete _grad_val[0];
            delete [] _grad_val;

            delete _hess_val[0][0];
            delete [] _hess_val[0];
            delete [] _hess_val;
        }

        hardening<VarContainer>* clone() const
        {
            return new hardening<VarContainer>(*this);
        }

        PSimpleFunction< VarContainer, double> simplefunction() const
        {
            return PSimpleFunction< VarContainer, double>( *_val );
        }

        PSimpleFunction< VarContainer, double> grad_simplefunction(int di) const
        {
            return PSimpleFunction< VarContainer, double>( *_grad_val[di] );
        }

        PSimpleFunction< VarContainer, double> hess_simplefunction(int di, int dj) const
        {
            return PSimpleFunction< VarContainer, double>( *_hess_val[di][dj] );
        }

        double operator()(const VarContainer &var)
        {
            return (*_val)(var);
        }

        double grad(const VarContainer &var, int di)
        {
            return (*_grad_val[di])(var);
        }

        double hess(const VarContainer &var, int di, int dj)
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
        }

        void eval_hess(const VarContainer &var)
        {
            (*_hess_val[0][0])(var);
        }

        double operator()() const
        {
            return (*_val)();
        }

        double grad(int di) const
        {
            return (*_grad_val[di])();
        }

        double hess(int di, int dj) const
        {
            return (*_hess_val[di][dj])();
        }

    private:
        void construct()
        {
            this->_name = "hardening";
            this->_var_name.clear();
            this->_var_name.push_back("alpha");
            this->_var_description.clear();
            this->_var_description.push_back("Equivalent plastic strain");
            
            _val = new hardening_f<VarContainer>();
            
            _grad_val = new PSimpleBase< VarContainer, double>*[1];
            _grad_val[0] = new hardening_grad_0<VarContainer>();
            
            _hess_val = new PSimpleBase< VarContainer, double>**[1];
            _hess_val[0] = new PSimpleBase< VarContainer, double>*[1];
            _hess_val[0][0] = new hardening_hess_0_0<VarContainer>();
        }

    };


}
#endif
