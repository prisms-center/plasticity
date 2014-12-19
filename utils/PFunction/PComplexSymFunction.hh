
#ifndef PComplexSymFunction_HH
#define PComplexSymFunction_HH

#include "PFunction.hh"
#include <ginac/ginac.h>    // compile with: -lcln -lginac
#include <complex>

namespace PRISMS
{

    /// Real valued symbolic functions

    template< class VarContainer>
    class PComplexSymFunction : public PFuncBase< VarContainer, std::complex<double> >
    {
        GiNaC::ex _e;
        std::vector< GiNaC::symbol> _sym;

        std::complex<double> _val;
        std::vector< std::complex<double> > _grad;
        std::vector< std::vector< std::complex<double> > > _hess;

    public:

        // ----------------------------------------------------------
        //   Non-Inherited:

        void set(const std::string &name, const std::vector<GiNaC::symbol> &sym, const GiNaC::ex &e)
        {
            _sym.clear();
            _sym.resize(_sym.size());

            _var_name.clear();
            _var_name.resize(_sym.size());

            for(int i = 0; i < sym.size(); i++)
            {
                _sym.push_back(sym[i]);
                _sym[i] = sym[i];
                _var_name.push_back(_sym[i].get_name());
            }

            _name = name;

            _e = e;
        }


        // ----------------------------------------------------------
        //   Inherited:

        using PFuncBase< VarContainer, std::complex<double> >::_name;
        using PFuncBase< VarContainer, std::complex<double> >::_var_name;

        virtual PComplexSymFunction<VarContainer> *clone() const
        {
            return new PComplexSymFunction<VarContainer>(*this);
        };

        // ----------------------------------------------------------
        // Use these functions if you want to evaluate a single value
        std::complex<double> operator()(const VarContainer &var);
        std::complex<double> grad(const VarContainer &var, int di);
        std::complex<double> hess(const VarContainer &var, int di, int dj);

    };

    template< class VarContainer>
    std::complex<double> PComplexSymFunction<VarContainer>::operator()(const VarContainer &var)
    {
        GiNaC::exmap m;
        for(int i = 0; i < var.size(); i++)
            m[_sym[i]] = var[i].real() + var[i].imag() * GiNaC::I;
        
        return std::complex<double>(GiNaC::ex_to<GiNaC::numeric>(GiNaC::evalf(_e.real_part().subs(m))).to_double(),
                               GiNaC::ex_to<GiNaC::numeric>(GiNaC::evalf(_e.imag_part().subs(m))).to_double());
    }


    template< class VarContainer>
    std::complex<double> PComplexSymFunction<VarContainer>::grad(const VarContainer &var, int di)
    {
        GiNaC::ex de = _e.diff(_sym[di]);

        GiNaC::exmap m;
        for(int i = 0; i < var.size(); i++)
            m[_sym[i]] = var[i].real() + var[i].imag() * GiNaC::I;
        
        return std::complex<double>(GiNaC::ex_to<GiNaC::numeric>(GiNaC::evalf(de.real_part().subs(m))).to_double(),
                               GiNaC::ex_to<GiNaC::numeric>(GiNaC::evalf(de.imag_part().subs(m))).to_double());
    };

    template< class VarContainer>
    std::complex<double> PComplexSymFunction<VarContainer>::hess(const VarContainer &var, int di, int dj)
    {
        GiNaC::ex de = _e.diff(_sym[di]).diff(_sym[dj]);

        GiNaC::exmap m;
        for(int i = 0; i < var.size(); i++)
            m[_sym[i]] = var[i].real() + var[i].imag() * GiNaC::I;
        
        return std::complex<double>(GiNaC::ex_to<GiNaC::numeric>(GiNaC::evalf(de.real_part().subs(m))).to_double(),
                               GiNaC::ex_to<GiNaC::numeric>(GiNaC::evalf(de.imag_part().subs(m))).to_double());
    };

}


#endif