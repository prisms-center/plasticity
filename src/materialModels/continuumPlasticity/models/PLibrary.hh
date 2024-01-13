// created: 2016-8-5 14:00:38
// version: HEAD
// url: git@github.com:prisms-center/IntegrationTools.git
// commit: 010dc90537a86ed3458c4b12292b6af9545b3c1f

#ifndef PLIBRARY_HH
#define PLIBRARY_HH

#include<cstring>
#include "../../../../utils/IntegrationTools/PFunction.hh"
#include "../../../../utils/IntegrationTools/PPieceWise.hh"
//#include "../../../../src/ellipticBVP/ellipticBVP.cc"

namespace PRISMS
{

    /// Library where you can find functions and basis sets
    ///
    namespace PLibrary
    {
        // Use these functions to checkout objects which manage their own memory

        void checkout( std::string name, PSimpleFunction< std::vector<double>, double > &simplefunc);
        void checkout( std::string name, PSimpleFunction< double*, double > &simplefunc);

        void checkout( std::string name, PFunction< std::vector<double>, double > &func);
        void checkout( std::string name, PFunction< double*, double > &func);




        // Use these functions to checkout new 'Base' objects which the user must delete

        void checkout( std::string name, PSimpleBase< std::vector<double>, double > *&simplefunc);
        void checkout( std::string name, PSimpleBase< double*, double > *&simplefunc);

        void checkout( std::string name, PFuncBase< std::vector<double>, double > *&func);
        void checkout( std::string name, PFuncBase< double*, double > *&func);



    }

}


#endif
