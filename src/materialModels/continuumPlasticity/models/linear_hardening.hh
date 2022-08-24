// created: 2016-8-5 14:00:36
// version: master
// url: https://github.com/prisms-center/IntegrationToolsWriter.git
// commit: 8a15adf67355fad30bd75ce9ba6b1f8d24b9a537

#ifndef linear_hardening_HH
#define linear_hardening_HH

#include <cmath>
#include <cstdlib>
#include "../../../../utils/IntegrationTools/PFunction.hh"

namespace PRISMS
{
    template< class VarContainer>
    class linear_hardening_f : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return -var[0]*var[1];
        }

    public:

        linear_hardening_f()
        {
            this->_name = "linear_hardening_f";
        }

        std::string csrc() const
        {
            return "-var[0]*var[1]";
        }

        std::string sym() const
        {
            return "-alpha*K";
        }

        std::string latex() const
        {
            return "- \\alpha K";
        }

        linear_hardening_f* clone() const
        {
            return new linear_hardening_f(*this);
        }
    };

    template< class VarContainer>
    class linear_hardening_grad_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return -var[1];
        }

    public:

        linear_hardening_grad_0()
        {
            this->_name = "linear_hardening_grad_0";
        }

        std::string csrc() const
        {
            return "-var[1]";
        }

        std::string sym() const
        {
            return "-K";
        }

        std::string latex() const
        {
            return "- K";
        }

        linear_hardening_grad_0* clone() const
        {
            return new linear_hardening_grad_0(*this);
        }
    };

    template< class VarContainer>
    class linear_hardening_grad_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return -var[0];
        }

    public:

        linear_hardening_grad_1()
        {
            this->_name = "linear_hardening_grad_1";
        }

        std::string csrc() const
        {
            return "-var[0]";
        }

        std::string sym() const
        {
            return "-alpha";
        }

        std::string latex() const
        {
            return "- \\alpha";
        }

        linear_hardening_grad_1* clone() const
        {
            return new linear_hardening_grad_1(*this);
        }
    };

    template< class VarContainer>
    class linear_hardening_hess_0_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        linear_hardening_hess_0_0()
        {
            this->_name = "linear_hardening_hess_0_0";
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

        linear_hardening_hess_0_0* clone() const
        {
            return new linear_hardening_hess_0_0(*this);
        }
    };

    template< class VarContainer>
    class linear_hardening_hess_0_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return -1.0;
        }

    public:

        linear_hardening_hess_0_1()
        {
            this->_name = "linear_hardening_hess_0_1";
        }

        std::string csrc() const
        {
            return "-1.0";
        }

        std::string sym() const
        {
            return "-1";
        }

        std::string latex() const
        {
            return "-1";
        }

        linear_hardening_hess_0_1* clone() const
        {
            return new linear_hardening_hess_0_1(*this);
        }
    };

    template< class VarContainer>
    class linear_hardening_hess_1_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return -1.0;
        }

    public:

        linear_hardening_hess_1_0()
        {
            this->_name = "linear_hardening_hess_1_0";
        }

        std::string csrc() const
        {
            return "-1.0";
        }

        std::string sym() const
        {
            return "-1";
        }

        std::string latex() const
        {
            return "-1";
        }

        linear_hardening_hess_1_0* clone() const
        {
            return new linear_hardening_hess_1_0(*this);
        }
    };

    template< class VarContainer>
    class linear_hardening_hess_1_1 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 0.0;
        }

    public:

        linear_hardening_hess_1_1()
        {
            this->_name = "linear_hardening_hess_1_1";
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

        linear_hardening_hess_1_1* clone() const
        {
            return new linear_hardening_hess_1_1(*this);
        }
    };

    template<class VarContainer>
    class linear_hardening : public PFuncBase< VarContainer, double>
    {
    public:
        
        typedef typename PFuncBase< VarContainer, double>::size_type size_type;

        PSimpleBase< VarContainer, double> *_val;
        PSimpleBase< VarContainer, double> **_grad_val;
        PSimpleBase< VarContainer, double> ***_hess_val;
        
        linear_hardening()
        {
            construct();
        }

        linear_hardening(const linear_hardening &RHS )
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

        linear_hardening& operator=( linear_hardening RHS )
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

        ~linear_hardening()
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

        linear_hardening<VarContainer>* clone() const
        {
            return new linear_hardening<VarContainer>(*this);
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
            this->_name = "linear_hardening";
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
            
            _val = new linear_hardening_f<VarContainer>();
            
            _grad_val[0] = new linear_hardening_grad_0<VarContainer>();
            _grad_val[1] = new linear_hardening_grad_1<VarContainer>();
            
            _hess_val[0][0] = new linear_hardening_hess_0_0<VarContainer>();
            _hess_val[0][1] = new linear_hardening_hess_0_1<VarContainer>();
            _hess_val[1][0] = new linear_hardening_hess_1_0<VarContainer>();
            _hess_val[1][1] = new linear_hardening_hess_1_1<VarContainer>();
        }

    };


}
#endif
