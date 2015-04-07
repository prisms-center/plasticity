#ifndef PExtern_CC
#define PExtern_CC

#include<cstring>
#include<iostream>
#include<vector>
#include<cstdlib>

#include "../../../../utils/IntegrationTools/PExtern.hh"
#include "../../../../utils/IntegrationTools/PFunction.hh"
#include "../../../../utils/IntegrationTools/PField.hh"

// In future, might have more complicated OutType, 
//   so make all have 'void' return and pass everything by reference


extern "C"
{
    // Functions for using a PSimpleBase externally (say Python or Fortran)
    //   written for VarContainer=double*, OutType=double, hence 'dsd' in function names
    
    void PSimpleFunction_dsd_new(char* name, PRISMS::PSimpleBase<double*,double>* &f)
    { 
        PRISMS::PLibrary::checkout(std::string(name), f);
    }
    
    void PSimpleFunction_dsd_delete(PRISMS::PSimpleBase<double*,double>* &f)
    { 
        delete f;
        f = NULL;
    }
    
    void PSimpleFunction_dsd_name(PRISMS::PSimpleBase<double*,double> *f, char* name)
    {
        std::strcpy(name, f->name().c_str());
    }
    
    void PSimpleFunction_dsd_calc( PRISMS::PSimpleBase<double*,double> *f, double* var, double &val)
    {
        val = (*f)(var);
    }
    
    void PSimpleFunction_dsd_get( PRISMS::PSimpleBase<double*,double> *f, double  &val)
    {
        val = (*f)();
    }
    
    
    
    
    
    // Functions for using a PSimpleBase externally (say Python or Fortran)
    //   written for VarContainer=double, OutType=double, hence 'dd' in function names
    
    void PSimpleFunction_dd_new(char* name, PRISMS::PSimpleBase<double,double>* &f)
    { 
        PRISMS::PLibrary::checkout(std::string(name), f);
    }
    
    void PSimpleFunction_dd_delete(PRISMS::PSimpleBase<double,double>* &f)
    { 
        delete f;
        f = NULL;
    }
    
    void PSimpleFunction_dd_name(PRISMS::PSimpleBase<double,double> *f, char* name)
    {
        std::strcpy(name, f->name().c_str());
    }
    
    void PSimpleFunction_dd_calc( PRISMS::PSimpleBase<double,double> *f, double var, double &val)
    {
        val = (*f)(var);
    }
    
    void PSimpleFunction_dd_get( PRISMS::PSimpleBase<double,double> *f, double  &val)
    {
        val = (*f)();
    }
    
    
    
    
    
    // Functions for using a PFuncBase externally (say Python or Fortran)
    //   written for VarContainer=double*, OutType=double, hence 'dsd' in function names
    
    void PFunction_dsd_new(char* name, PRISMS::PFuncBase<double*,double>* &f)
    { 
        PRISMS::PLibrary::checkout(std::string(name), f);
    }
    
    void PFunction_dsd_delete(PRISMS::PFuncBase<double*,double>* &f)
    { 
        delete f;
        f = NULL;
    }
    
    void PFunction_dsd_name(PRISMS::PFuncBase<double*,double>* f, char* name)
    {
        std::strcpy(name, f->name().c_str());
    }
    
    void PFunction_dsd_size(PRISMS::PFuncBase<double*,double>* f, int &size)
    {
        size = f->size();
    }
    
    void PFunction_dsd_var_name(PRISMS::PFuncBase<double*,double>* f, int i, char* var_name)
    {
        std::strcpy(var_name, f->var_name(i).c_str());
    }
    
    void PFunction_dsd_var_description(PRISMS::PFuncBase<double*,double>* f, int i, char* var_description)
    {
        std::strcpy(var_description, f->var_description(i).c_str());
    }
    
    //void PFunction_dsd_simplefunc(PRISMS::PFuncBase<double*,double> *f, PSimpleBase<double*, double> *simplefunc);
    //void PFunction_dsd_grad_simplefunc(PRISMS::PFuncBase<double*,double> *f, int *di, PSimpleBase<double*, double> *simplefunc);
    //void PFunction_dsd_hess_simplefunc(PRISMS::PFuncBase<double*,double> *f, int *di, int *dj, PSimpleBase<double*, double> *simplefunc);
    
    void PFunction_dsd_calc(PRISMS::PFuncBase<double*,double>* f, double* var, double &val)
    {
        val = (*f)(var);
    }
    
    void PFunction_dsd_calc_grad(PRISMS::PFuncBase<double*,double>* f, double* var, int di, double &val)
    {
        val = (*f).grad(var, di);
    }
    
    void PFunction_dsd_calc_hess(PRISMS::PFuncBase<double*,double>* f, double* var, int di, int dj, double &val)
    {
        val = (*f).hess(var, di, dj);
    }
    
    void PFunction_dsd_eval(PRISMS::PFuncBase<double*,double>* f, double* var)
    {
        (*f)(var);
    }
    
    void PFunction_dsd_eval_grad(PRISMS::PFuncBase<double*,double>* f, double* var, int di)
    {
        (*f).grad(var, di);
    }
    
    void PFunction_dsd_eval_hess(PRISMS::PFuncBase<double*,double>* f, double* var, int di, int dj)
    {
        (*f).hess(var, di, dj);
    }
    
    void PFunction_dsd_get(PRISMS::PFuncBase<double*,double>* f, double &val)
    {
        val = (*f)();
    }
    
    void PFunction_dsd_get_grad(PRISMS::PFuncBase<double*,double>* f, int di, double &val)
    {
        val = (*f).grad(di);
    }
    
    void PFunction_dsd_get_hess(PRISMS::PFuncBase<double*,double>* f, int di, int dj, double &val)
    {
        val = (*f).hess(di, dj);
    }
    
    
    
    
    
    // Functions for using a PBasisSetBase externally (say Python or Fortran)
    //   written for InType=double, OutType=double, hence 'dd' in function names
    
    void PBasisSet_dd_new(char* name, PRISMS::PBasisSetBase<double,double>* &b, int N)
    { 
        PRISMS::PLibrary::checkout(std::string(name), b, N);
    }
    
    void PBasisSet_dd_delete(PRISMS::PBasisSetBase<double,double>* &b)
    { 
        delete b;
        b = NULL;
    }
    
    void PBasisSet_dd_name(PRISMS::PBasisSetBase<double,double>* b, char* name)
    {
        std::strcpy(name, b->name().c_str());
    }
    
    void PBasisSet_dd_description(PRISMS::PBasisSetBase<double,double>* b, char* description)
    {
        std::strcpy(description, b->description().c_str());
    }
    
    void PBasisSet_dd_size(PRISMS::PBasisSetBase<double,double>* b, int &size)
    {
        size = b->size();
    }
    
    void PBasisSet_dd_resize(PRISMS::PBasisSetBase<double,double>* b, int N)
    {
        b->resize(N);
    }
    
    void PBasisSet_dd_max_size(PRISMS::PBasisSetBase<double,double>* b, int &max_size)
    {
        max_size = b->max_size();
    }
    
    //void PBasisSet_dd_basis_function(PRISMS::PFuncBase<double*,double> *b, int* term, PFuncBase<double,double> *f);
    
    void PBasisSet_dd_calc(PRISMS::PBasisSetBase<double,double>* b, int term, double var, double &val)
    {
        val = (*b)(term, var);
    }
    
    void PBasisSet_dd_calc_grad(PRISMS::PBasisSetBase<double,double>* b, int term, double var, double &val)
    {
        val = (*b).grad(term, var);
    }
    
    void PBasisSet_dd_calc_hess(PRISMS::PBasisSetBase<double,double>* b, int term, double var, double &val)
    {
        val = (*b).hess(term, var);
    }
        
    void PBasisSet_dd_eval(PRISMS::PBasisSetBase<double,double>* b, double var)
    {
        (*b).eval(var);
    }
    
    void PBasisSet_dd_eval_grad(PRISMS::PBasisSetBase<double,double>* b, double var)
    {
        (*b).eval_grad(var);
    }
    
    void PBasisSet_dd_eval_hess(PRISMS::PBasisSetBase<double,double>* b, double var)
    {
        (*b).eval_hess(var);
    }
    
    void PBasisSet_dd_get(PRISMS::PBasisSetBase<double,double>* b, int term, double &val)
    {
        val = (*b)(term);
    }
    
    void PBasisSet_dd_get_grad(PRISMS::PBasisSetBase<double,double>* b, int term, double &val)
    {
        val = (*b).grad(term);
    }
    
    void PBasisSet_dd_get_hess(PRISMS::PBasisSetBase<double,double>* b, int term, double &val)
    {
        val = (*b).hess(term);
    }
    
    //void PBasisSet_dd_getall(PRISMS::PBasisSetBase<double,double>* b, const double* &val);
    //void PBasisSet_dd_getall_grad(PRISMS::PBasisSetBase<double,double>* b, const double* &val);
    //void PBasisSet_dd_getall_hess(PRISMS::PBasisSetBase<double,double>* b, const double* &val);
    
    
    
    
    
    // Functions for using a PSeriesFunction externally (say Python or Fortran)
    //   written for InType=double, OutType=double, VarContainer=double*, IndexContainer=int*, hence 'dsis'
    
    void PSeriesFunction_dsis_new(PRISMS::PSeriesFunction<double,double,double*,int*>* &f)
    {
        f = new PRISMS::PSeriesFunction<double,double,double*,int*>(0.0, 1.0);
    }
    
    void PSeriesFunction_dsis_setnew(PRISMS::PSeriesFunction<double,double,double*,int*>* &f, PRISMS::PBasisSetBase<double,double>** basis_set, int order)
    {
        std::vector<PRISMS::PBasisSet<double, double> > basis_set_vec;
        for( int i=0; i<order; i++)
            basis_set_vec.push_back( *basis_set[i]);
        
        f = new PRISMS::PSeriesFunction<double,double,double*,int*>(0.0, 1.0, basis_set_vec);
    }
        
    void PSeriesFunction_dsis_delete(PRISMS::PSeriesFunction<double,double,double*,int*>* &f)
    {
        delete f;
        f = NULL;
    }
    
    void PSeriesFunction_dsis_clear(PRISMS::PSeriesFunction<double,double,double*,int*>* f)
    {
        (*f).clear();
    }
    
    void PSeriesFunction_dsis_set(PRISMS::PSeriesFunction<double,double,double*,int*>* f, PRISMS::PBasisSetBase<double,double>** basis_set, int order)
    {
        std::vector<PRISMS::PBasisSet<double,double> > basis_set_vec;
        PRISMS::PBasisSet<double,double> tmp;
        for( int i=0; i<order; i++)
        {
            PRISMS::PLibrary::checkout( (*(basis_set[i])).name(), tmp, (*(basis_set[i])).size());
            basis_set_vec.push_back( tmp);
        }
        
        (*f).set(basis_set_vec);
    }
    
    void PSeriesFunction_dsis_order(PRISMS::PSeriesFunction<double,double,double*,int*>* f, int &order)
    {
        order = (*f).coeff().order();
    }
    
    void PSeriesFunction_dsis_volume(PRISMS::PSeriesFunction<double,double,double*,int*>* f, int &volume)
    {
        volume = (*f).coeff().volume();
    }
    
    void PSeriesFunction_dsis_dim(PRISMS::PSeriesFunction<double,double,double*,int*>* f, int i, int &dim)
    {
        dim = (*f).coeff().dim(i);
    }
    
    void PSeriesFunction_dsis_get_linear_coeff(PRISMS::PSeriesFunction<double,double,double*,int*>* f, int i, double &coeff)
    {
        coeff = (*f).coeff()(i);
    }
    
    void PSeriesFunction_dsis_get_tensor_coeff(PRISMS::PSeriesFunction<double,double,double*,int*>* f, int* term, double &coeff)
    {
        coeff = (*f).coeff()(term);
    }
    
    void PSeriesFunction_dsis_set_linear_coeff(PRISMS::PSeriesFunction<double,double,double*,int*>* f, int i, double coeff)
    {
        (*f).coeff()(i) = coeff;
    }
    
    void PSeriesFunction_dsis_set_tensor_coeff(PRISMS::PSeriesFunction<double,double,double*,int*>* f, int* term, double coeff)
    {
        (*f).coeff()(term) = coeff;
    }
    
    void PSeriesFunction_dsis_linear_index(PRISMS::PSeriesFunction<double,double,double*,int*>* f, int* term, int &linear_index)
    {
        linear_index = (*f).coeff().linear_index(term);
    }
    
    void PSeriesFunction_dsis_tensor_indices(PRISMS::PSeriesFunction<double,double,double*,int*>* f, int linear_index, int* term)
    {
        // assumes term.size() == order()  (the tensor dimensions)
        //   not sure if this is how we want to do it, but it avoids assuming push_back()
        (*f).coeff().tensor_indices(linear_index, term);
    }
    
    // ----------------------------------------------------------
    // Use these functions if you want to evaluate a single value
    
    void PSeriesFunction_dsis_calc(PRISMS::PSeriesFunction<double,double,double*,int*>* f, double* var, double &val)
    {
        val = (*f)(var);
    }
    
    void PSeriesFunction_dsis_calc_grad(PRISMS::PSeriesFunction<double,double,double*,int*>* f, double* var, int di, double &val)
    {
        val = (*f).grad(var,di);
    }
    
    void PSeriesFunction_dsis_calc_hess(PRISMS::PSeriesFunction<double,double,double*,int*>* f, double* var, int di, int dj, double &val)
    {
        val = (*f).hess(var,di,dj);
    }
    
    // ----------------------------------------------------------
    // Use these functions to evaluate several values, then use 'get' methods to access results
    
    void PSeriesFunction_dsis_eval(PRISMS::PSeriesFunction<double,double,double*,int*>* f, double* var)
    {
        (*f).eval(var);
    }
    
    void PSeriesFunction_dsis_eval_grad(PRISMS::PSeriesFunction<double,double,double*,int*>* f, double* var)
    {
        (*f).eval_grad(var);
    }
    
    void PSeriesFunction_dsis_eval_hess(PRISMS::PSeriesFunction<double,double,double*,int*>* f, double* var)
    {
        (*f).eval_hess(var);
    }
    
    void PSeriesFunction_dsis_get(PRISMS::PSeriesFunction<double,double,double*,int*>* f, double &val)
    {
        val = (*f)();
    }
    
    void PSeriesFunction_dsis_get_grad(PRISMS::PSeriesFunction<double,double,double*,int*>* f, int di, double &val)
    {
        val = (*f).grad(di);
    }
    
    void PSeriesFunction_dsis_get_hess(PRISMS::PSeriesFunction<double,double,double*,int*>* f, int di, int dj, double &val)
    {
        val = (*f).hess(di,dj);
    }
    
    
    // ----------------------------------------------------------
    // Functions for evaluating basis functions & their derivatives:

    // Use these functions if you want to evaluate a single value

    //   use basis index and term index for individual basis function
    
    void PSeriesFunction_dsis_calc_basis(PRISMS::PSeriesFunction<double,double,double*,int*>* f, int bindex, int term, double* var, double &val)
    {
        val = (*f).basis(bindex,term,var);
    }
    
    void PSeriesFunction_dsis_calc_basis_grad(PRISMS::PSeriesFunction<double,double,double*,int*>* f, int bindex, int term, double* var, double &val)
    {
        val = (*f).basis_grad(bindex,term,var);
    }
    
    void PSeriesFunction_dsis_calc_basis_hess(PRISMS::PSeriesFunction<double,double,double*,int*>* f, int bindex, int term, double* var, double &val)
    {
        val = (*f).basis_hess(bindex,term,var);
    }
    
    //   or use tensor indices to evaluate basis function product
    void PSeriesFunction_dsis_calc_tensor_basis(PRISMS::PSeriesFunction<double,double,double*,int*>* f, int* term, double* var, double &val)
    {
        val = (*f).basis(term,var);
    }
    
    void PSeriesFunction_dsis_calc_tensor_basis_grad(PRISMS::PSeriesFunction<double,double,double*,int*>* f, int* term, double* var, int di, double &val)
    {
        val = (*f).basis_grad(term,var,di);
    }
    
    void PSeriesFunction_dsis_calc_tensor_basis_hess(PRISMS::PSeriesFunction<double,double,double*,int*>* f, int* term, double* var, int di, int dj, double &val)
    {
        val = (*f).basis_hess(term,var,di,dj);
    }
    
    // ----------------------------------------------------------
    // Use these functions to evaluate all basis functions,
    //   then use following methods to access results.
    
    void PSeriesFunction_dsis_eval_basis_all(PRISMS::PSeriesFunction<double,double,double*,int*>* f, double* var)
    {
        (*f).eval_basis(var);
    }
    
    void PSeriesFunction_dsis_eval_basis(PRISMS::PSeriesFunction<double,double,double*,int*>* f, double* var, int i)
    {
        (*f).eval_basis(var,i);
    }
    
    void PSeriesFunction_dsis_eval_basis_grad_all(PRISMS::PSeriesFunction<double,double,double*,int*>* f, double* var)
    {
        (*f).eval_basis_grad(var);
    }
    
    void PSeriesFunction_dsis_eval_basis_grad(PRISMS::PSeriesFunction<double,double,double*,int*>* f, double* var, int i)
    {
        (*f).eval_basis_grad(var,i);
    }
    
    void PSeriesFunction_dsis_eval_basis_hess_all(PRISMS::PSeriesFunction<double,double,double*,int*>* f, double* var)
    {
        (*f).eval_basis_hess(var);
    }
    
    void PSeriesFunction_dsis_eval_basis_hess(PRISMS::PSeriesFunction<double,double,double*,int*>* f, double* var, int i)
    {
        (*f).eval_basis_hess(var,i);
    }
    

    //   use basis index and term index for individual basis function
    void PSeriesFunction_dsis_get_basis(PRISMS::PSeriesFunction<double,double,double*,int*>* f, int bindex, int term, double &val)
    {
        val = (*f).basis(bindex,term);
    }
    
    void PSeriesFunction_dsis_get_basis_grad(PRISMS::PSeriesFunction<double,double,double*,int*>* f, int bindex, int term, double &val)
    {
        val = (*f).basis_grad(bindex,term);
    }
    
    void PSeriesFunction_dsis_get_basis_hess(PRISMS::PSeriesFunction<double,double,double*,int*>* f, int bindex, int term, double &val)
    {
        val = (*f).basis_hess(bindex,term);
    }
    
    //   or use tensor indices to evaluate basis function product
    void PSeriesFunction_dsis_get_tensor_basis(PRISMS::PSeriesFunction<double,double,double*,int*>* f, int* term, double &val)
    {
        val = (*f).basis(term);
    }
    
    void PSeriesFunction_dsis_get_tensor_basis_grad(PRISMS::PSeriesFunction<double,double,double*,int*>* f, int* term, int di, double &val)
    {
        val = (*f).basis_grad(term,di);
    }
    
    void PSeriesFunction_dsis_get_tensor_basis_hess(PRISMS::PSeriesFunction<double,double,double*,int*>* f, int* term, int di, int dj, double &val)
    {
        val = (*f).basis_hess(term,di,dj);
    }
    
    // Functions for using constructing a 2D PRISMS::Body externally (say Python or Fortran),
    //   allowing access to PFields
    //   written for Coordinate=double*, OutType=double, DIM=2
    
    void Body2D_new(char* vtkfile, PRISMS::Body<double*,2>* &b)
    {
        b = new PRISMS::Body<double*,2>();
        (*b).read_vtk(std::string(vtkfile));
    };
    
    void Body2D_delete(PRISMS::Body<double*,2>* &b)
    {
        delete b;
        b = NULL;
    };
    
    
    // Functions for using a 2D scalar PField externally (say Python or Fortran), as a PFunction.
    //   From a Body pointer, returns a pointer to a PFuncBase
    //   written for Coordinate=double*, OutType=double, DIM=2
    //   don't delete this! it will be deleted by deleting the Body
    
    void ScalarField2D(char* name, PRISMS::Body<double*,2>* b, PRISMS::PFuncBase<double*,double>* &f)
    {
        f = &((*b).find_scalar_field(std::string(name)));
    };
    
    
    // Functions for using constructing a 3D PRISMS::Body externally (say Python or Fortran),
    //   allowing access to PFields
    //   written for Coordinate=double*, OutType=double, DIM=3
    
    void Body3D_new(char* vtkfile, PRISMS::Body<double*,3>* &b)
    {
        b = new PRISMS::Body<double*,3>();
        (*b).read_vtk(std::string(vtkfile));
    };
    
    void Body3D_delete(PRISMS::Body<double*,3>* &b)
    {
        delete b;
        b = NULL;
    };
    
    
    // Functions for using a 3D scalar PField externally (say Python or Fortran), as a PFunction.
    //   From a Body pointer, returns a pointer to a PFuncBase
    //   written for Coordinate=double*, OutType=double, DIM=2
    //   don't delete this! it will be deleted by deleting the Body
    
    void ScalarField3D(char* name, PRISMS::Body<double*,3>* b, PRISMS::PFuncBase<double*,double>* &f)
    {
        f = &((*b).find_scalar_field(std::string(name)));
    };
    
    
}


#endif
