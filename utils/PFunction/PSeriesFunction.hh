
#ifndef PSERIESFUNCTION_HH
#define PSERIESFUNCTION_HH

#include<cstring>
#include<iostream>
#include<vector>
#include<limits>

#include "PFunction.hh"

namespace PRISMS
{
    /// Evaluate basis functions and their derivatives. Store and access the results.
    ///
    template< class InType, class OutType>
    class PBasisSetBase
    {
        public:
        
        std::string _name;
        std::string _description;
        
        std::vector<OutType> _val;
        std::vector<OutType> _grad_val;
        std::vector<OutType> _hess_val;
        
        std::string name() const
        {
            return _name;
        };
        std::string description() const
        {
            return _description;
        };
        int size() const
        {
            return _val.size();
        };
        
        PBasisSetBase( int N)
        {
            resize(N);
        };
        
        void is_derived_from_PBasisSetBase() const
        {
            return;
        };
        
        virtual void resize( int N)
        {
            _val.resize(N);
            _grad_val.resize(N);
            _hess_val.resize(N);
        };
        
        virtual ~PBasisSetBase()
        {
            
        };
        
        virtual int max_size()
        {
            // default to (essentially) no limit
            return std::numeric_limits<int>::max();
        };
        
        virtual PBasisSetBase<InType,OutType>* clone() const
        {
            return new PBasisSetBase<InType, OutType>(*this);
        };
        
        virtual PFunction<InType, OutType> basis_function(int term) const
        {
            undefined("const PFunction<InType, OutType>& basis_function(int term) const");
            return PFunction<InType, OutType>();
        };
        
        // ----------------------------------------------------------
        // Use these functions if you want to evaluate a single value
        OutType operator()(int term, const InType &var)
        {
            return _val[term] = eval(term, var);
        };
        OutType grad(int term, const InType &var)
        {
            return _grad_val[term] = eval_grad(term, var);
        };
        OutType hess(int term, const InType &var)
        {
            return _hess_val[term] = eval_hess(term, var);
        };
        

        // ----------------------------------------------------------
        // Use these functions to evaluate several values, then use 'get' methods to access results
        
        // By default, evaluate each individual term one-by-one. 
        // These are virtual so derived classes may implement more efficient methods.
        //    Returns vector containing results
        virtual const std::vector<OutType>& eval(const InType &var)
        {
            for( int i=0; i<_val.size(); i++)
                (*this)(i,var);
            return _val;
        };
        virtual const std::vector<OutType>& eval_grad(const InType &var)
        {
            for( int i=0; i<_val.size(); i++)
                (*this).grad(i,var);
            return _grad_val;
        };
        virtual const std::vector<OutType>& eval_hess(const InType &var)
        {
            for( int i=0; i<_val.size(); i++)
                (*this).hess(i,var);
            return _hess_val;
        };

        // Getters for individual terms
        OutType operator()(int term) const
        {
            return _val[term];
        };
        OutType grad(int term) const
        {
            return _grad_val[term];
        };
        OutType hess(int term) const
        {
            return _hess_val[term];
        };
        
        // Getters returning vector containing all terms
        const std::vector<OutType>& operator()() const
        {
            return _val;
        };
        const std::vector<OutType>& grad() const
        {
            return _grad_val;
        };
        const std::vector<OutType>& hess() const
        {
            return _hess_val;
        };

    private:
        
        /// ----------------------------------------------------------
        /// !!! Derived classes must define these functions !!!
        ///   Usually this is done using PBasisSetWriter
        virtual OutType eval(int term, const InType &var)
        {
            undefined("OutType PBasisSetBase::eval(int term, const InType &var)");
            return OutType();
        };
        virtual OutType eval_grad(int term, const InType &var)
        {
            undefined("OutType PBasisSetBase::eval_grad(int term, const InType &var)");
            return OutType();
        };
        virtual OutType eval_hess(int term, const InType &var)
        {
            undefined("OutType PBasisSetBase::eval_hess(int term, const InType &var)");
            return OutType();
        };
    
        void undefined(std::string fname) const
        {
            std::cout << "Error. In PBasisSetBase '" << _name << "'." << std::endl;
            std::cout << "   The member function '" << fname << "' has not been defined." << std::endl;
            exit(1);
        }
    };
    
    
    
    
    /// Evaluate basis functions and their derivatives. Store and access the results.
    ///
    template< class InType, class OutType>
    class PBasisSet
    {
        public:
        
        
        std::string name() const
        {
            return (*ptr).name();
        };
        std::string description() const
        {
            return (*ptr)._description();
        };
        int size() const
        {
            return (*ptr).size();
        };
        
        
        virtual void resize( int N)
        {
            (*ptr).resize(N);
        };
        
        virtual int max_size()
        {
            // default to (essentially) no limit
            return (*ptr).max_size();
        };
        
        
        // Output PFunction for individual basis functions:
        
        
        PFunction<InType, OutType> basis_function(int term)
        {
            return (*ptr).basis_function(term);
        };
        
        
        // ----------------------------------------------------------
        // Use these functions if you want to evaluate a single value
        OutType operator()(int term, const InType &var)
        {
            return (*ptr)(term,var);
        };
        OutType grad(int term, const InType &var)
        {
            return (*ptr).grad(term, var);
        };
        OutType hess(int term, const InType &var)
        {
            return (*ptr).hess(term, var);
        };
        

        // ----------------------------------------------------------
        // Use these functions to evaluate several values, then use 'get' methods to access results
        
        // By default, evaluate each individual term one-by-one. 
        // These are virtual so derived classes may implement more efficient methods.
        //    Returns vector containing results
        virtual const std::vector<OutType>& eval(const InType &var)
        {
            return (*ptr).eval(var);
        };
        virtual const std::vector<OutType>& eval_grad(const InType &var)
        {
            return (*ptr).eval_grad(var);
        };
        virtual const std::vector<OutType>& eval_hess(const InType &var)
        {
            return (*ptr).eval_hess(var);
        };

        // Getters for individual terms
        OutType operator()(int term) const
        {
            return (*ptr)(term);
        };
        OutType grad(int term) const
        {
            return (*ptr).grad(term);
        };
        OutType hess(int term) const
        {
            return (*ptr).hess(term);
        };
        
        // Getters returning vector containing all terms
        const std::vector<OutType>& operator()() const
        {
            return (*ptr)();
        };
        const std::vector<OutType>& grad() const
        {
            return (*ptr).grad();
        };
        const std::vector<OutType>& hess() const
        {
            return (*ptr).hess();
        };
        
        // PFunction unique members ------------------------------------------

        PBasisSet &operator=(const PBasisSet &RHS)
        {
            if(ptr != NULL)
                delete ptr;
            ptr = RHS.ptr->clone();
            return *this;
        }

        template<class T> PBasisSet &operator=(const T &RHS)
        {
            RHS.is_derived_from_PBasisSetBase();

            if(ptr != NULL)
                delete ptr;
            ptr = RHS.clone();
            return *this;
        }

        // If you use this, PBasisSet becomes the 'owner' of the function RHS points to
        //    and it will delete it
        PBasisSet &set( PBasisSet *RHS)
        {
            if(RHS == NULL)
            {
                std::cout << "Error in PBasisSet::set. RHS == NULL." << std::endl;
                exit(1);
            }
            if(ptr != NULL)
                delete ptr;
            ptr = RHS;
            return *this;
        }

        PBasisSet()
        {
            ptr = NULL;
        }

        PBasisSet(const PBasisSet &RHS)
        {
            if( RHS.ptr != NULL)
                ptr = RHS.ptr->clone();
            else 
                ptr = NULL;
        }
        
        template<class T> PBasisSet(const T &RHS)
        {
            RHS.is_derived_from_PBasisSetBase();

            ptr = RHS.clone();
        };

        ~PBasisSet()
        {
            if(ptr != NULL)
                delete ptr;
        };


    private:
    
        PBasisSetBase<InType, OutType> *ptr;
        
    };
    
    
    
    
    /// This class is for a series function that is the tensor product of indepedent bases
    ///
    ///   Example for order-3 tensor:
    ///
    ///    f(x,y,z) = sum_i sum_j sum_k T_ijk * phix_i(x) * phiy_j(y) * phiz_k(z)
    ///
    ///    summations run from i to _dim[i]
    ///    x, y, z are independent variables (dx/dy = 0, etc.)
    ///
    /// Template parameters:
    ///   - InType: input type for BasisSets (must be the same for all BasisSets)
    ///   - OutType: output type for SeriesFunction and BasisSets
    ///   - VarContainer: container for input
    ///   - IndexContainer: container for tensor indices
    ///
    /// Example template parameters are: <double, double, std::vector<double>, std::vector<int> >
    ///
    template< class InType, class OutType, class VarContainer, class IndexContainer>
    class PSeriesFunction : public PFuncBase<VarContainer, OutType>
    {
    public:
        using PFuncBase<VarContainer, OutType>::_name;
        using PFuncBase<VarContainer, OutType>::_var_name;
        using PFuncBase<VarContainer, OutType>::_var_description;
        
        // tensor info
        std::vector< int > _dim;                        // number of basis functions in each basis set
        std::vector< int > _unroll;                     // used to translate from tensor indices to linear index
        std::vector< OutType > _coeff_tensor;           // unrolled tensor of coefficients (first index is outer loop)
        int _order;                                     // _dim.size()
        int _volume;                                    // _coeff_tensor.size() = product_i(_dim[i]) 
        
        OutType _identity_val;
        OutType _zero_val;
        
        // basis sets
        std::vector< PBasisSet<InType, OutType> > _basis_set;
        
        // evaluated values
        OutType _val;
        std::vector<OutType> _grad_val;
        std::vector< std::vector< OutType> > _hess_val;
        
        
        PSeriesFunction(OutType zero, OutType identity)
        {
            _zero_val = zero;
            _identity_val = identity;
        };
        
        PSeriesFunction( OutType zero, OutType identity, const std::vector<PBasisSet<InType,OutType> > &basis_set) 
        {
            _zero_val = zero;
            _identity_val = identity;

            set(basis_set);
        };
        
        virtual ~PSeriesFunction()
        {
            clear();
        };
        
        void clear()
        {
            _basis_set.clear();
            _coeff_tensor.clear();
            _dim.clear();
            _unroll.clear();
            _grad_val.clear();
            _hess_val.clear();
            
        }
        
        /// generate the PSeriesFunction, by cloning 'basis_set' and resizing everything to match
        ///
        void set( const std::vector<PBasisSet<InType,OutType> > &basis_set)
        {
            //std::cout << "begin PSeriesFunction::set()" << std::endl;
            clear();
            
            //std::cout << "  basis_set.size(): " << basis_set.size() << std::endl;
            _order = basis_set.size();
            //std::cout << "  resize: " << basis_set.size() << std::endl;
            _basis_set.resize(basis_set.size());
            
            //std::cout << "  _order: " << _order << std::endl;
            _dim.resize( _order );
            
            _grad_val.resize( _order );
            _hess_val.resize( _order );
            
            for( int i=0; i<_order; i++)
            {
                //std::cout << "  i: " << i << std::endl;
                _dim[i] = basis_set[i].size();
                _hess_val[i].resize( _order );
                //std::cout << "  copy basis set" << std::endl;
                _basis_set[i] = basis_set[i];
            }
            
            //std::cout << "  generate_unroll" << std::endl;
                
            generate_unroll();
            
            _volume = _dim[0]*_unroll[0];
            _coeff_tensor.resize(_volume);
            //std::cout << "finish PSeriesFunction::set()" << std::endl;
            
        };
        
        int order() const
        {
            return _order;
        }
        
        int volume() const
        {
            return _volume;
        }
        
        int dim(int i) const
        {
            return _basis_set[i].size();
        }
        
        OutType& coeff( int i)
        {
            return _coeff_tensor[i];
        }
        
        OutType& coeff( const IndexContainer &term)
        {
            return _coeff_tensor[linear_index(term)];
        }
        
        int linear_index( const IndexContainer &term) const
        {
            int lindex = 0;
            for( int i=0; i<_unroll.size(); i++)
                lindex += term[i]*_unroll[i];
            return lindex;
        };
        
        template <class Container> void tensor_indices( int lindex, Container &term) const
        {   // assumes term.size() == basis_dim.size()  (the tensor dimensions)
            //   not sure if this is how we want to do it, but it avoids assuming push_back()
            
            for( int i=0; i<_unroll.size(); i++)
            {
                term[i] = lindex/_unroll[i];
                lindex -= term[i]*_unroll[i];
            }
        }

        virtual PSeriesFunction<InType, OutType, VarContainer, IndexContainer> *clone() const
        {
            return new PSeriesFunction<InType, OutType, VarContainer, IndexContainer>(*this);
        };

        // ----------------------------------------------------------
        // Use these functions if you want to evaluate a single value
        virtual OutType operator()(const VarContainer &var)
        {
            // evaluate basis functions needed
            eval_basis(var);
            
            // evaluate basis function products, multiply by _coeff_tensor & sum
            return calc_val();
        };
        
        virtual OutType grad(const VarContainer &var, int di)
        {
            // evaluate basis functions needed
            for( int i=0; i<_order; i++)
            {
                if( i == di)
                {   
                    _basis_set[i].eval_grad(var[i]);
                }
                else
                {
                    _basis_set[i](var[i]);
                }
            }
            
            // evaluate basis function products, multiply by _coeff_tensor & sum
            return calc_grad_val(di);
        };
        
        virtual OutType hess(const VarContainer &var, int di, int dj)
        {
            // evaluate basis functions needed
            for( int i=0; i<_order; i++)
            {
                if( i == di || i == dj)
                {   
                    if( di == dj)
                    {
                        _basis_set[i].eval_hess(var[i]);
                    }
                    else
                    {
                        _basis_set[i].eval_grad(var[i]);
                    }
                }
                else
                {
                    _basis_set[i](var[i]);
                }
            }
            
            // evaluate basis function products, multiply by _coeff_tensor & sum
            return calc_hess_val(di,dj);
        };
        
        // ----------------------------------------------------------
        // Use these functions to evaluate several values, then use 'get' methods to access results
        virtual void eval(const VarContainer &var)
        {
            (*this)(var);
        };
        virtual void eval_grad(const VarContainer &var)
        {
            eval_basis(var);
            
            eval_basis_grad(var);
            
            for( int i=0; i<_order; i++)
                (*this).calc_grad_val(i);
        };
        virtual void eval_hess(const VarContainer &var)
        {
            eval_basis(var);
            
            eval_basis_grad(var);
            
            eval_basis_hess(var);
            
            
            for( int i=0; i<_order; i++)
            for( int j=0; j<_order; j++)
                (*this).calc_hess_val(i,j);
        };

        virtual OutType operator()() const
        {
            return _val;
        };
        virtual OutType grad(int di) const
        {
            return _grad_val[di];
        };
        virtual OutType hess(int di, int dj) const
        {
            return _hess_val[di][dj];
        };


        // Functions for evaluating basis functions & their derivatives:

        // ----------------------------------------------------------
        // Use these functions if you want to evaluate a single value

        //   use basis index and term index for individual basis function
        virtual OutType basis(int bindex, int term, const VarContainer &var)
        {
            return _basis_set[bindex](term,var[bindex]);
        };
        virtual OutType basis_grad(int bindex, int term, const VarContainer &var)
        {
            return _basis_set[bindex].grad(term, var[bindex]);
        };
        virtual OutType basis_hess(int bindex, int term, const VarContainer &var)
        {
            return _basis_set[bindex].hess(term, var[bindex]);
        };

        //   or use tensor indices to evaluate basis function product
        virtual OutType basis(const IndexContainer &term, const VarContainer &var)
        {
            // evaluate basis functions needed
            OutType tmp = _identity_val;
            for( int i=0; i<_order; i++)
                tmp *= _basis_set[i](term[i],var[i]);
                
            return tmp;
        };
        virtual OutType basis_grad(const IndexContainer &term, const VarContainer &var, int di)
        {
            // evaluate basis functions needed
            OutType tmp = _identity_val;
            for( int i=0; i<_order; i++)
            {
                if( i == di)
                {   
                    tmp *= _basis_set[i].grad(term[i], var[i]);
                }
                else
                {
                    tmp *= _basis_set[i](term[i], var[i]);
                }
            }
            
            return tmp;
        };
        virtual OutType basis_hess(const IndexContainer &term, const VarContainer &var, int di, int dj)
        {
            // evaluate basis functions needed
            OutType tmp = _identity_val;
            for( int i=0; i<_order; i++)
            {
                if( i == di || i == dj)
                {   
                    if( di == dj)
                    {
                        tmp *= _basis_set[i].hess(term[i], var[i]);
                    }
                    else
                    {
                        tmp *= _basis_set[i].grad(term[i], var[i]);
                    }
                }
                else
                {
                    tmp *= _basis_set[i](term[i], var[i]);
                }
            }
            
            return tmp;
        };

        // ----------------------------------------------------------
        // Use these functions to evaluate all basis functions,
        //   then use following methods to access results.

        virtual void eval_basis(const VarContainer &var)
        {
            // evaluate basis functions
            for( int i=0; i<_order; i++)
                _basis_set[i].eval(var[i]);
            
        };
        virtual void eval_basis(const VarContainer &var, int i)
        {
            // evaluate basis functions
            _basis_set[i].eval(var[i]);
            
        };
        virtual void eval_basis_grad(const VarContainer &var)
        {
            // evaluate basis grad functions
            for( int i=0; i<_order; i++)
                _basis_set[i].eval_grad(var[i]);
        };
        virtual void eval_basis_grad(const VarContainer &var, int i)
        {
            // evaluate basis grad functions
            _basis_set[i].eval_grad(var[i]);
        };
        virtual void eval_basis_hess(const VarContainer &var)
        {
            // evaluate basis hess functions
            for( int i=0; i<_order; i++)
                _basis_set[i].eval_hess(var[i]);
        };
        virtual void eval_basis_hess(const VarContainer &var, int i)
        {
            // evaluate basis hess functions
            _basis_set[i].eval_hess(var[i]);
        };

        //   use basis index and term index for individual basis function
        virtual OutType basis(int bindex, int term) const
        {
            return _basis_set[bindex](term);
        };
        virtual OutType basis_grad(int bindex, int term) const
        {
            return _basis_set[bindex].grad(term);
        };
        virtual OutType basis_hess(int bindex, int term) const
        {
            return _basis_set[bindex].hess(term);
        };

        //   or use tensor indices to evaluate basis function product
        virtual OutType basis(const IndexContainer &term) const
        {
            // evaluate basis function product
            OutType tmp = _identity_val;
            for( int i=0; i<_order; i++)
                tmp *= _basis_set[i](term[i]);
                
            return tmp;
        };
        virtual OutType basis_grad(const IndexContainer &term, int di) const
        {
            // evaluate basis function product
            OutType tmp = _identity_val;
            for( int i=0; i<_order; i++)
            {
                if( i == di)
                {   
                    tmp *= _basis_set[i].grad(term[i]);
                }
                else
                {
                    tmp *= _basis_set[i](term[i]);
                }
            }
            
            return tmp;
        };
        virtual OutType basis_hess(const IndexContainer &term, int di, int dj) const
        {
            // evaluate basis function product
            OutType tmp = _identity_val;
            for( int i=0; i<_order; i++)
            {
                if( i == di || i == dj)
                {   
                    if( di == dj)
                    {
                        tmp *= _basis_set[i].hess(term[i]);
                    }
                    else
                    {
                        tmp *= _basis_set[i].grad(term[i]);
                    }
                }
                else
                {
                    tmp *= _basis_set[i](term[i]);
                }
            }
            
            return tmp;
        };
        
        // ----------------------------------------------------------
        // Read and write routines: 
        //   Write routines only uses 'getters', they assume you have 
        //   already evaluated the necessary basis functions
        //
        // These are not included in PExtern
        
        void print_basis(std::ostream &sout)
        {
            std::vector<int> term = std::vector<int>(_order,0);
            for( int i=0; i<_volume; i++)
            {
                tensor_indices(i,term);
                for( int j=0; j<term.size(); j++)
                    sout << term[j] << " ";
                sout << basis(term) << "\n";
            }
                
        };
        void print_basis_grad(std::ostream &sout, int di)
        {
            std::vector<int> term = std::vector<int>(_order,0);
            for( int i=0; i<_volume; i++)
            {
                tensor_indices(i,term);
                for( int j=0; j<term.size(); j++)
                    sout << term[j] << " ";
                sout << basis_grad(term,di) << "\n";
            }
                
        };
        void print_basis_hess(std::ostream &sout, int di, int dj)
        {
            std::vector<int> term = std::vector<int>(_order,0);
            for( int i=0; i<_volume; i++)
            {
                tensor_indices(i,term);
                for( int j=0; j<term.size(); j++)
                    sout << term[j] << " ";
                sout << basis_hess(term,di,dj) << "\n";
            }
                
        };
        void print_coeff(std::ostream &sout)
        {
            std::vector<int> term = std::vector<int>(_order,0);
            for( int i=0; i<_volume; i++)
            {
                tensor_indices(i,term);
                for( int j=0; j<term.size(); j++)
                    sout << term[j] << " ";
                sout << coeff(i) << "\n";
            }
        };
        void read_coeff(std::istream &sin)
        {
            // no error checking, 
            //   assumes format identical to print_coeff output
            std::vector<int> term = std::vector<int>(_order,0);
            for( int i=0; i<_volume; i++)
            {
                // ignores tensor indices
                for( int j=0; j<term.size(); j++)
                {    
                    sin >> term[j];
                }
                sin >> coeff(i);
            }
        };
        
        
        ///
        private:
        
        void generate_unroll()
        {
            _unroll.resize(_dim.size());
            _unroll[_dim.size()-1] = 1;
            for( int i=_dim.size()-2; i>=0; i--)
                _unroll[i] = _unroll[i+1]*_dim[i+1];
        };
        
        // ----------------------------------------------------------
        // These assume you've evaluated the necessary _basis, _basis_grad, and _basis_hess functions
        //   and just need to take products & sum
        
        // evaluate basis function products, multiply by _coeff_tensor & sum
        OutType calc_val()
        {
            std::vector<int> tindex(_order);
            OutType tmp;
            _val = _zero_val;
            for( int i=0; i<_volume; i++)
            {
                tmp = _coeff_tensor[i];
                tensor_indices(i, tindex);
                
                for( int j=0; j<_order; j++)
                {
                    tmp *= _basis_set[j](tindex[j]);
                }
                _val += tmp;
            }
            
            return _val;
        };
        
        // evaluate basis function products, multiply by _coeff_tensor & sum
        OutType calc_grad_val(int di)
        {
            std::vector<int> tindex(_order);
            OutType tmp;
            _grad_val[di] = _zero_val;
            for( int i=0; i<_volume; i++)
            {
                tmp = _coeff_tensor[i];
                tensor_indices(i, tindex);
                for( int j=0; j<_order; j++)
                {
                    if( j == di)
                    {
                        tmp *= _basis_set[j].grad(tindex[j]);
                    }
                    else
                    {
                        tmp *= _basis_set[j](tindex[j]);
                    }
                }
                
                _grad_val[di] += tmp;
                
            }
            return _grad_val[di];
        };
        
        // evaluate basis function products, multiply by _coeff_tensor & sum
        OutType calc_hess_val(int di, int dj)
        {
            std::vector<int> tindex(_order);
            OutType tmp;
            _hess_val[di][dj] = _zero_val;
            for( int i=0; i<_volume; i++)
            {
                tmp = _coeff_tensor[i];
                tensor_indices(i, tindex);
                for( int j=0; j<_order; j++)
                {
                    if( j == di || j == dj)
                    {
                        if( di == dj)
                        {
                            tmp *= _basis_set[j].hess(tindex[j]);
                        }
                        else
                        {
                            tmp *= _basis_set[j].grad(tindex[j]);
                        }
                    }
                    else
                    {
                        tmp *= _basis_set[j](tindex[j]);
                    }
                }
                _hess_val[di][dj] += tmp;
            }
            return _hess_val[di][dj];
        }
    };
    

}


#endif