
#ifndef PFUNCTION_HH
#define PFUNCTION_HH

#include<cstring>
#include<iostream>
#include<vector>
#include<cstdlib>

namespace PRISMS
{

    /// Base classes for functions that can be hard-coded,
    ///   then shared and used elsewhere
    
    /// A simple expression evaluator
    ///
    template< class VarContainer, class OutType>
    class PSimpleBase
    {
        public:
        std::string _name;
        OutType _val;
        
        std::string name() const { return _name;};
        OutType operator()( const VarContainer &var){ return _val = eval(var);};
        OutType operator()() const { return _val;};
        
        void is_derived_from_PSimpleBase() const
        {
            return;
        };
        
        virtual PSimpleBase<VarContainer, OutType>* clone() const
        {
            return new PSimpleBase<VarContainer, OutType>(*this);
        }
        
        private:
        virtual OutType eval( const VarContainer &var) const { undefined("OutType eval( const VarContainer &var)"); return OutType();};
        
        void undefined(std::string fname) const
        {
            std::cout << "Error in PSimpleBase '" << _name << "'." << std::endl;
            std::cout << "   The member function '" << fname << "' has not been defined." << std::endl;
            exit(1);
        }
    };
    
    
    template< class VarContainer, class OutType>
    class PSimpleFunction 
    {
    private:
        PSimpleBase<VarContainer,OutType> *ptr;

    public:

        std::string name() const
        {
            return (*ptr).name();
        };
        

        // ----------------------------------------------------------
        // Use this function if you want to evaluate,
        //   return and store result
        OutType operator()(const VarContainer &var)
        {
            return (*ptr)(var);
        };

        // ----------------------------------------------------------
        // Then use 'get' methods to access results later
        void eval(const VarContainer &var)
        {
            (*ptr)(var);
        };

        OutType operator()() const
        {
            return (*ptr)();
        };

        // PFunction unique members ------------------------------------------

        PSimpleFunction& operator=(const PSimpleFunction &RHS)
        {
            if(ptr != NULL)
                delete ptr;
            ptr = RHS.ptr->clone();
            return *this;
        };

        template<class T> 
        PSimpleFunction& operator=(const T &RHS)
        {
            RHS.is_derived_from_PSimpleBase();

            if(ptr != NULL)
                delete ptr;
            ptr = RHS.clone();
            return *this;
        };

        // If you use this, PSimpleFunction becomes the 'owner' of the function RHS points to
        //    and it will delete it
        PSimpleFunction& set( PSimpleBase<VarContainer,OutType> *RHS)
        {
            if(RHS == NULL)
            {
                std::cout << "Error in PSimpleFunction::set. RHS == NULL." << std::endl;
                exit(1);
            }
            if(ptr != NULL)
                delete ptr;
            ptr = RHS;
            return *this;
        }
        
        PSimpleFunction()
        {
            ptr = NULL;
        }

        PSimpleFunction(const PSimpleFunction &RHS)
        {
            if( RHS.ptr != NULL)
                ptr = RHS.ptr->clone();
            else 
                ptr = NULL;
        }
        
        template<class T> PSimpleFunction(const T &RHS)
        {
            RHS.is_derived_from_PSimpleBase();

            ptr = RHS.clone();
            
        };

        ~PSimpleFunction()
        {
            if(ptr != NULL)
                delete ptr;
        };

    
    };

    
    /// A Base class for a function, including grad & hess
    /// 
    template< class VarContainer, class OutType>
    class PFuncBase
    {
    public:

        std::string _name;
        std::vector<std::string> _var_name;
        std::vector<std::string> _var_description;
        
        virtual ~PFuncBase(){};

        std::string name()
        {
            return _name;
        };
        int size() const
        {
            return _var_name.size();
        };
        std::string var_name(int i)
        {
            return _var_name[i];
        };
        std::string var_description(int i)
        {
            return _var_description[i];
        };

        void is_derived_from_PFuncBase() const
        {
            return;
        };

        virtual PFuncBase<VarContainer, OutType> *clone() const
        {
            return new PFuncBase<VarContainer, OutType>(*this);
        };
        
        virtual PSimpleFunction<VarContainer, OutType> simplefunction() const
        {
            undefined("PSimpleFunction<VarContainer, OutType> simplefunction() const");
            return PSimpleFunction<VarContainer, OutType>();
        };
        
        virtual PSimpleFunction<VarContainer, OutType> grad_simplefunction(int di) const
        {
            undefined("PSimpleFunction<VarContainer, OutType> grad_simplefunction() const");
            return PSimpleFunction<VarContainer, OutType>();
        };
        
        virtual PSimpleFunction<VarContainer, OutType> hess_simplefunction(int di, int dj) const
        {
            undefined("PSimpleFunction<VarContainer, OutType> hess_simplefunction(int di, int dj) const");
            return PSimpleFunction<VarContainer, OutType>();
        };

        // ----------------------------------------------------------
        // Use these functions if you want to evaluate a single value
        virtual OutType operator()(const VarContainer &var)
        {
            undefined("OutType operator()(const VarContainer &var)");
            return OutType();
        };
        virtual OutType grad(const VarContainer &var, int di)
        {
            undefined("OutType grad(const VarContainer &var, int di)");
            return OutType();
        };
        virtual OutType hess(const VarContainer &var, int di, int dj)
        {
            undefined("OutType hess(const VarContainer &var, int di, int dj)");
            return OutType();
        };

        // ----------------------------------------------------------
        // Use these functions to evaluate several values, then use 'get' methods to access results
        virtual void eval(const VarContainer &var)
        {
            undefined("void eval_grad( const VarContainer &var)");
        };
        virtual void eval_grad(const VarContainer &var)
        {
            undefined("void eval_grad( const VarContainer &var)");
        };
        virtual void eval_hess(const VarContainer &var)
        {
            undefined("void eval_hess( const VarContainer &var)");
        };

        virtual OutType operator()() const
        {
            undefined("OutType operator()");
            return OutType();
        };
        virtual OutType grad(int di) const
        {
            undefined("OutType grad(int di)");
            return OutType();
        };
        virtual OutType hess(int di, int dj) const
        {
            undefined("OutType hess(int di, int dj)");
            return OutType();
        };

    private:
        void undefined(std::string fname) const
        {
            std::cout << "Error in PFuncBase '" << _name << "'." << std::endl;
            std::cout << "   The member function '" << fname << "' has not been defined." << std::endl;
            exit(1);
        }

    };
    
    
    template<class VarContainer, class OutType>
    class PFlexFunction : public PFuncBase< VarContainer, OutType>
    {
    public:
        using PFuncBase< VarContainer, OutType>::_name;
        using PFuncBase< VarContainer, OutType>::_var_name;
        using PFuncBase< VarContainer, OutType>::_var_description;
        
        PSimpleFunction< VarContainer, OutType> _val;
        std::vector< PSimpleFunction< VarContainer, OutType> > _grad_val;
        std::vector< std::vector< PSimpleFunction< VarContainer, OutType> > > _hess_val;
        
        PFlexFunction()
        {
            
        }

        PFlexFunction(const PFlexFunction &RHS )
        {
            _name = RHS._name;
            _var_name = RHS._var_name;
            _var_description = RHS._var_description;
            
            _val = RHS._val;
            _grad_val = RHS._grad_val;
            _hess_val = RHS._hess_val;
        }

        PFlexFunction& operator=(const PFlexFunction &RHS )
        {
            _name = RHS._name;
            _var_name = RHS._var_name;
            _var_description = RHS._var_description;
            
            _val = RHS._val;
            _grad_val = RHS._grad_val;
            _hess_val = RHS._hess_val;
            
        }

        ~PFlexFunction()
        {
            
        };

        PFlexFunction<VarContainer, OutType>* clone() const
        {
            return new PFlexFunction<VarContainer, OutType>(*this);
        };

        PSimpleFunction< VarContainer, OutType> simplefunction() const
        {
            return  _val;
        };

        PSimpleFunction< VarContainer, OutType> grad_simplefunction(int di) const
        {
            return _grad_val[di];
        };

        PSimpleFunction< VarContainer, OutType> hess_simplefunction(int di, int dj) const
        {
            return _hess_val[di][dj];
        };

        OutType operator()(const VarContainer &var)
        {
            return _val(var);
        };

        OutType grad(const VarContainer &var, int di)
        {
            return _grad_val[di](var);
        };

        OutType hess(const VarContainer &var, int di, int dj)
        {
            return _hess_val[di][dj](var);
        };
        
        void eval(const VarContainer &var)
        {
            (*this)(var);
        };

        void eval_grad(const VarContainer &var)
        {
            for( int i=0; i<_grad_val.size(); i++)
                _grad_val[i](var);
        };

        void eval_hess(const VarContainer &var)
        {
            for( int i=0; i<_hess_val.size(); i++)
                for( int j=0; j<_hess_val[i].size(); j++)
                    _hess_val[i][j](var);
        };

        OutType operator()() const
        {
            return _val();
        };

        OutType grad(int di) const
        {
            return _grad_val[di]();
        };

        OutType hess(int di, int dj) const
        {
            return _hess_val[di][dj]();
        };

        

    };
    
    
    /// A class that contains a ptr to a PFuncBase object
    ///   - like a smart ptr class
    ///   - same interface as PFuncBase
    ///   - allows for using PFuncBase objects polymorphically 
    ///     without dereferencing and without worrying about new/delete
    ///    
    ///   example: MyFuncA, MyFuncB, MyFuncC, etc. are defined:
    ///     template< class VarContainer>
    ///     MyFuncX : public PFuncBase<VarContainer, double>
    ///
    ///   // Then you can do things like this:
    ///   
    ///   MyFuncA<std::vector<double> > my_func_a;
    ///   MyFuncB<std::vector<double> > my_func_b;
    ///
    ///   PFuncBase<std::vector<double>, double>* my_func_c_ptr;
    ///   my_func_c_ptr = new MyFuncC<std::vector<double>, double >();
    ///
    ///   PFunction<std::vector<double>, double > f, g, h;
    ///
    ///   f = my_func_a;
    ///   f = my_func_b;
    ///   g.set(my_func_c_ptr->clone());
    ///   h.set(my_func_c_ptr);
    ///   double result = f(3.0) + g(4.0) + h(5.0);
    ///   
    ///   - No deletions are used in this example.  
    ///     PFunction::set makes PFunction the 'owner' of the MyFuncC object and it will delete it.
    ///
    template< class VarContainer, class OutType>
    class PFunction 
    {
    public:

        std::string name() const
        {
            return (*ptr).name();
        };
        int size() const
        {
            return (*ptr).size();
        };
        std::string var_name(int i)
        {
            return (*ptr).var_name(i);
        };
        std::string var_description(int i)
        {
            return (*ptr).var_description(i);
        };
        
        PSimpleFunction<VarContainer, OutType> simplefunction() const
        {
            return (*ptr).simplefunction();
        };
        
        PSimpleFunction<VarContainer, OutType> grad_simplefunction(int di) const
        {
            return (*ptr).grad_simplefunction(di);
        };
        
        PSimpleFunction<VarContainer, OutType> hess_simplefunction(int di, int dj) const
        {
            return (*ptr).hess_simplefunction(di, dj);
        };

        // ----------------------------------------------------------
        // Use these functions if you want to evaluate a single value
        OutType operator()(const VarContainer &var)
        {
            return (*ptr)(var);
        };
        OutType grad(const VarContainer &var, int di)
        {
            return (*ptr).grad(var, di);
        };
        OutType hess(const VarContainer &var, int di, int dj)
        {
            return (*ptr).hess(var, di, dj);
        };

        // ----------------------------------------------------------
        // Use these functions to evaluate several values, then use 'get' methods to access results
        void eval(const VarContainer &var)
        {
            return (*ptr).eval(var);
        };
        void eval_grad(const VarContainer &var)
        {
            return (*ptr).eval_grad(var);
        };
        void eval_hess(const VarContainer &var)
        {
            return (*ptr).eval_hess(var);
        };

        OutType operator()() const
        {
            return (*ptr)();
        };
        OutType grad(int di) const
        {
            return (*ptr).grad(di);
        };
        OutType hess(int di, int dj) const
        {
            return (*ptr).hess(di, dj);
        };


        // PFunction unique members ------------------------------------------

        PFunction &operator=(const PFunction &RHS)
        {
            if(ptr != NULL)
                delete ptr;
            ptr = RHS.ptr->clone();
            return *this;
        };

        template<class T> PFunction &operator=(const T &RHS)
        {
            RHS.is_derived_from_PFuncBase();

            if(ptr != NULL)
                delete ptr;
            ptr = RHS.clone();
            return *this;
        };

        // If you use this, PFunction becomes the 'owner' of the function RHS points to
        //    and it will delete it
        PFunction &set( PFuncBase<VarContainer,OutType> *RHS)
        {
            if(RHS == NULL)
            {
                std::cout << "Error in PFunction::set. RHS == NULL." << std::endl;
                exit(1);
            }
            if(ptr != NULL)
                delete ptr;
            ptr = RHS;
            return *this;
        };

        PFunction()
        {
            ptr = NULL;
        };

        PFunction(const PFunction &RHS)
        {
            if( RHS.ptr != NULL)
                ptr = RHS.ptr->clone();
            else 
                ptr = NULL;
        };
        
        template<class T> PFunction(const T &RHS)
        {
            RHS.is_derived_from_PFuncBase();

            ptr = RHS.clone();
            
        };

        ~PFunction()
        {
            if(ptr != NULL)
                delete ptr;
        };

    private:
        PFuncBase<VarContainer,OutType> *ptr;

    };

}


#endif