
#ifndef PBasisSetWriter_HH
#define PBasisSetWriter_HH

#include<cstring>
#include<iostream>
#include<vector>

namespace PRISMS
{

    /// Class used to write PBasisSet classes

    class PBasisSetWriter
    {
    public:

        std::string _name;
        std::string _description;
        
        std::string _basic_indent;
        
        std::string _intype;
        std::string _outtype;
        
        bool _write_phi;
        std::vector< std::string> _phi;
        
        bool _write_grad;
        std::vector< std::string> _grad;
        
        bool _write_hess;
        std::vector< std::string> _hess;
        
        
        // Constructor initializes strings to call 'undefined' message
        PBasisSetWriter(const std::string &name, const std::string &description);
        
        // After construction, need to set things
        void set_basic_indent(std::string basic_indent);
        void set_intype(std::string intype);
        void set_outtype(std::string outtype);
        void phi_on();
        void phi_off();
        void grad_on();
        void grad_off();
        void hess_on();
        void hess_off();
        void set_types(std::string intype, std::string outtype);
        
        // Write the PBasisSet
        
        // 1) Non-recursive generating expression (example: phi_i = x^i;)
        void sym2code( std::string f, 
                              std::string var, 
                              std::string index,
                              int N,
                              std::ostream &sout);
        
        // 2) Recursive generating expression (example: phi_2 = 2*x*phi_1 - phi_0; phi_0 = 1; phi_1 = x;)
        void sym2code( std::string f, 
                              std::string var, 
                              const std::vector<std::string> &phi_init,
                              const std::vector<std::string> &phi_sym,
                              int N,
                              std::ostream &sout);
        
        void code(     const std::vector<std::string> &phi, 
                              const std::vector<std::string> &grad,
                              const std::vector<std::string> &hess,
                              std::ostream &sout);
        
        
        // TODO:
        //void code(     const std::string &json_str, std::ostream &sout);
        
        void head( std::ostream &sout) const;
        
        void foot( std::ostream &sout) const;
        

    private:
        
        std::string indent(int step) const;
        
        void write_basis_function(int I, const std::string &name, const std::string &f, std::ostream &sout) const;
        
        void code( std::ostream &sout) const;
        
        
        
    };

}


#endif