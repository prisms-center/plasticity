// created: 2014-12-20 0:32:00
// version: develop
// url: git@github.com:prisms-center/IntegrationTools.git
// commit: 3e86b7184d2ec2de450d3668f7dd88dcd0396838

#ifndef PLIBRARY_CC
#define PLIBRARY_CC

#include<cstring>
#include<stdexcept>
#include<vector>
<<<<<<< HEAD
#include "neohook.hh"
#include "stvenkir.hh"
#include "linear_hardening.hh"
#include "Cu_hardening.hh"
#include "von_mises.hh"
=======
#include "stvenkir.hh"
>>>>>>> 1853fc05524906a2a2d829468a5c2e263978073c
#include "quadlog.hh"
#include "hardening.hh"
#include "von_mises.hh"
#include "neohook.hh"
#include "PLibrary.hh"

namespace PRISMS
{

        void PLibrary::checkout( std::string name, PSimpleFunction< std::vector<double>, double > &simplefunc)
        {
<<<<<<< HEAD
            if( name == "neohook_f") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_f< std::vector<double> >() );
            if( name == "neohook_grad_0") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_grad_0< std::vector<double> >() );
            if( name == "neohook_grad_1") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_grad_1< std::vector<double> >() );
            if( name == "neohook_grad_2") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_grad_2< std::vector<double> >() );
            if( name == "neohook_grad_3") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_grad_3< std::vector<double> >() );
            if( name == "neohook_grad_4") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_grad_4< std::vector<double> >() );
            if( name == "neohook_grad_5") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_grad_5< std::vector<double> >() );
            if( name == "neohook_grad_6") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_grad_6< std::vector<double> >() );
            if( name == "neohook_hess_0_0") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_0_0< std::vector<double> >() );
            if( name == "neohook_hess_0_1") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_0_1< std::vector<double> >() );
            if( name == "neohook_hess_0_2") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_0_2< std::vector<double> >() );
            if( name == "neohook_hess_0_3") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_0_3< std::vector<double> >() );
            if( name == "neohook_hess_0_4") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_0_4< std::vector<double> >() );
            if( name == "neohook_hess_0_5") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_0_5< std::vector<double> >() );
            if( name == "neohook_hess_0_6") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_0_6< std::vector<double> >() );
            if( name == "neohook_hess_1_0") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_1_0< std::vector<double> >() );
            if( name == "neohook_hess_1_1") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_1_1< std::vector<double> >() );
            if( name == "neohook_hess_1_2") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_1_2< std::vector<double> >() );
            if( name == "neohook_hess_1_3") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_1_3< std::vector<double> >() );
            if( name == "neohook_hess_1_4") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_1_4< std::vector<double> >() );
            if( name == "neohook_hess_1_5") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_1_5< std::vector<double> >() );
            if( name == "neohook_hess_1_6") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_1_6< std::vector<double> >() );
            if( name == "neohook_hess_2_0") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_2_0< std::vector<double> >() );
            if( name == "neohook_hess_2_1") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_2_1< std::vector<double> >() );
            if( name == "neohook_hess_2_2") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_2_2< std::vector<double> >() );
            if( name == "neohook_hess_2_3") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_2_3< std::vector<double> >() );
            if( name == "neohook_hess_2_4") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_2_4< std::vector<double> >() );
            if( name == "neohook_hess_2_5") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_2_5< std::vector<double> >() );
            if( name == "neohook_hess_2_6") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_2_6< std::vector<double> >() );
            if( name == "neohook_hess_3_0") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_3_0< std::vector<double> >() );
            if( name == "neohook_hess_3_1") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_3_1< std::vector<double> >() );
            if( name == "neohook_hess_3_2") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_3_2< std::vector<double> >() );
            if( name == "neohook_hess_3_3") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_3_3< std::vector<double> >() );
            if( name == "neohook_hess_3_4") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_3_4< std::vector<double> >() );
            if( name == "neohook_hess_3_5") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_3_5< std::vector<double> >() );
            if( name == "neohook_hess_3_6") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_3_6< std::vector<double> >() );
            if( name == "neohook_hess_4_0") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_4_0< std::vector<double> >() );
            if( name == "neohook_hess_4_1") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_4_1< std::vector<double> >() );
            if( name == "neohook_hess_4_2") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_4_2< std::vector<double> >() );
            if( name == "neohook_hess_4_3") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_4_3< std::vector<double> >() );
            if( name == "neohook_hess_4_4") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_4_4< std::vector<double> >() );
            if( name == "neohook_hess_4_5") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_4_5< std::vector<double> >() );
            if( name == "neohook_hess_4_6") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_4_6< std::vector<double> >() );
            if( name == "neohook_hess_5_0") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_5_0< std::vector<double> >() );
            if( name == "neohook_hess_5_1") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_5_1< std::vector<double> >() );
            if( name == "neohook_hess_5_2") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_5_2< std::vector<double> >() );
            if( name == "neohook_hess_5_3") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_5_3< std::vector<double> >() );
            if( name == "neohook_hess_5_4") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_5_4< std::vector<double> >() );
            if( name == "neohook_hess_5_5") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_5_5< std::vector<double> >() );
            if( name == "neohook_hess_5_6") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_5_6< std::vector<double> >() );
            if( name == "neohook_hess_6_0") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_6_0< std::vector<double> >() );
            if( name == "neohook_hess_6_1") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_6_1< std::vector<double> >() );
            if( name == "neohook_hess_6_2") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_6_2< std::vector<double> >() );
            if( name == "neohook_hess_6_3") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_6_3< std::vector<double> >() );
            if( name == "neohook_hess_6_4") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_6_4< std::vector<double> >() );
            if( name == "neohook_hess_6_5") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_6_5< std::vector<double> >() );
            if( name == "neohook_hess_6_6") simplefunc = PSimpleFunction< std::vector<double>, double >( neohook_hess_6_6< std::vector<double> >() );
            if( name == "stvenkir_f") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_f< std::vector<double> >() );
            if( name == "stvenkir_grad_0") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_grad_0< std::vector<double> >() );
            if( name == "stvenkir_grad_1") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_grad_1< std::vector<double> >() );
            if( name == "stvenkir_grad_2") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_grad_2< std::vector<double> >() );
            if( name == "stvenkir_grad_3") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_grad_3< std::vector<double> >() );
            if( name == "stvenkir_grad_4") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_grad_4< std::vector<double> >() );
            if( name == "stvenkir_grad_5") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_grad_5< std::vector<double> >() );
            if( name == "stvenkir_grad_6") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_grad_6< std::vector<double> >() );
            if( name == "stvenkir_hess_0_0") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_0_0< std::vector<double> >() );
            if( name == "stvenkir_hess_0_1") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_0_1< std::vector<double> >() );
            if( name == "stvenkir_hess_0_2") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_0_2< std::vector<double> >() );
            if( name == "stvenkir_hess_0_3") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_0_3< std::vector<double> >() );
            if( name == "stvenkir_hess_0_4") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_0_4< std::vector<double> >() );
            if( name == "stvenkir_hess_0_5") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_0_5< std::vector<double> >() );
            if( name == "stvenkir_hess_0_6") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_0_6< std::vector<double> >() );
            if( name == "stvenkir_hess_1_0") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_1_0< std::vector<double> >() );
            if( name == "stvenkir_hess_1_1") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_1_1< std::vector<double> >() );
            if( name == "stvenkir_hess_1_2") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_1_2< std::vector<double> >() );
            if( name == "stvenkir_hess_1_3") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_1_3< std::vector<double> >() );
            if( name == "stvenkir_hess_1_4") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_1_4< std::vector<double> >() );
            if( name == "stvenkir_hess_1_5") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_1_5< std::vector<double> >() );
            if( name == "stvenkir_hess_1_6") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_1_6< std::vector<double> >() );
            if( name == "stvenkir_hess_2_0") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_2_0< std::vector<double> >() );
            if( name == "stvenkir_hess_2_1") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_2_1< std::vector<double> >() );
            if( name == "stvenkir_hess_2_2") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_2_2< std::vector<double> >() );
            if( name == "stvenkir_hess_2_3") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_2_3< std::vector<double> >() );
            if( name == "stvenkir_hess_2_4") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_2_4< std::vector<double> >() );
            if( name == "stvenkir_hess_2_5") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_2_5< std::vector<double> >() );
            if( name == "stvenkir_hess_2_6") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_2_6< std::vector<double> >() );
            if( name == "stvenkir_hess_3_0") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_3_0< std::vector<double> >() );
            if( name == "stvenkir_hess_3_1") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_3_1< std::vector<double> >() );
            if( name == "stvenkir_hess_3_2") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_3_2< std::vector<double> >() );
            if( name == "stvenkir_hess_3_3") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_3_3< std::vector<double> >() );
            if( name == "stvenkir_hess_3_4") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_3_4< std::vector<double> >() );
            if( name == "stvenkir_hess_3_5") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_3_5< std::vector<double> >() );
            if( name == "stvenkir_hess_3_6") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_3_6< std::vector<double> >() );
            if( name == "stvenkir_hess_4_0") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_4_0< std::vector<double> >() );
            if( name == "stvenkir_hess_4_1") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_4_1< std::vector<double> >() );
            if( name == "stvenkir_hess_4_2") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_4_2< std::vector<double> >() );
            if( name == "stvenkir_hess_4_3") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_4_3< std::vector<double> >() );
            if( name == "stvenkir_hess_4_4") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_4_4< std::vector<double> >() );
            if( name == "stvenkir_hess_4_5") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_4_5< std::vector<double> >() );
            if( name == "stvenkir_hess_4_6") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_4_6< std::vector<double> >() );
            if( name == "stvenkir_hess_5_0") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_5_0< std::vector<double> >() );
            if( name == "stvenkir_hess_5_1") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_5_1< std::vector<double> >() );
            if( name == "stvenkir_hess_5_2") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_5_2< std::vector<double> >() );
            if( name == "stvenkir_hess_5_3") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_5_3< std::vector<double> >() );
            if( name == "stvenkir_hess_5_4") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_5_4< std::vector<double> >() );
            if( name == "stvenkir_hess_5_5") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_5_5< std::vector<double> >() );
            if( name == "stvenkir_hess_5_6") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_5_6< std::vector<double> >() );
            if( name == "stvenkir_hess_6_0") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_6_0< std::vector<double> >() );
            if( name == "stvenkir_hess_6_1") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_6_1< std::vector<double> >() );
            if( name == "stvenkir_hess_6_2") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_6_2< std::vector<double> >() );
            if( name == "stvenkir_hess_6_3") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_6_3< std::vector<double> >() );
            if( name == "stvenkir_hess_6_4") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_6_4< std::vector<double> >() );
            if( name == "stvenkir_hess_6_5") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_6_5< std::vector<double> >() );
            if( name == "stvenkir_hess_6_6") simplefunc = PSimpleFunction< std::vector<double>, double >( stvenkir_hess_6_6< std::vector<double> >() );
            if( name == "linear_hardening_f") simplefunc = PSimpleFunction< std::vector<double>, double >( linear_hardening_f< std::vector<double> >() );
            if( name == "linear_hardening_grad_0") simplefunc = PSimpleFunction< std::vector<double>, double >( linear_hardening_grad_0< std::vector<double> >() );
            if( name == "linear_hardening_grad_1") simplefunc = PSimpleFunction< std::vector<double>, double >( linear_hardening_grad_1< std::vector<double> >() );
            if( name == "linear_hardening_hess_0_0") simplefunc = PSimpleFunction< std::vector<double>, double >( linear_hardening_hess_0_0< std::vector<double> >() );
            if( name == "linear_hardening_hess_0_1") simplefunc = PSimpleFunction< std::vector<double>, double >( linear_hardening_hess_0_1< std::vector<double> >() );
            if( name == "linear_hardening_hess_1_0") simplefunc = PSimpleFunction< std::vector<double>, double >( linear_hardening_hess_1_0< std::vector<double> >() );
            if( name == "linear_hardening_hess_1_1") simplefunc = PSimpleFunction< std::vector<double>, double >( linear_hardening_hess_1_1< std::vector<double> >() );
            if( name == "Cu_hardening_f") simplefunc = PSimpleFunction< std::vector<double>, double >( Cu_hardening_f< std::vector<double> >() );
            if( name == "Cu_hardening_grad_0") simplefunc = PSimpleFunction< std::vector<double>, double >( Cu_hardening_grad_0< std::vector<double> >() );
            if( name == "Cu_hardening_grad_1") simplefunc = PSimpleFunction< std::vector<double>, double >( Cu_hardening_grad_1< std::vector<double> >() );
            if( name == "Cu_hardening_hess_0_0") simplefunc = PSimpleFunction< std::vector<double>, double >( Cu_hardening_hess_0_0< std::vector<double> >() );
            if( name == "Cu_hardening_hess_0_1") simplefunc = PSimpleFunction< std::vector<double>, double >( Cu_hardening_hess_0_1< std::vector<double> >() );
            if( name == "Cu_hardening_hess_1_0") simplefunc = PSimpleFunction< std::vector<double>, double >( Cu_hardening_hess_1_0< std::vector<double> >() );
            if( name == "Cu_hardening_hess_1_1") simplefunc = PSimpleFunction< std::vector<double>, double >( Cu_hardening_hess_1_1< std::vector<double> >() );
            if( name == "von_mises_f") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_f< std::vector<double> >() );
            if( name == "von_mises_grad_0") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_grad_0< std::vector<double> >() );
            if( name == "von_mises_grad_1") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_grad_1< std::vector<double> >() );
            if( name == "von_mises_grad_2") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_grad_2< std::vector<double> >() );
            if( name == "von_mises_grad_3") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_grad_3< std::vector<double> >() );
            if( name == "von_mises_grad_4") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_grad_4< std::vector<double> >() );
            if( name == "von_mises_hess_0_0") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_hess_0_0< std::vector<double> >() );
            if( name == "von_mises_hess_0_1") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_hess_0_1< std::vector<double> >() );
            if( name == "von_mises_hess_0_2") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_hess_0_2< std::vector<double> >() );
            if( name == "von_mises_hess_0_3") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_hess_0_3< std::vector<double> >() );
            if( name == "von_mises_hess_0_4") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_hess_0_4< std::vector<double> >() );
            if( name == "von_mises_hess_1_0") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_hess_1_0< std::vector<double> >() );
            if( name == "von_mises_hess_1_1") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_hess_1_1< std::vector<double> >() );
            if( name == "von_mises_hess_1_2") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_hess_1_2< std::vector<double> >() );
            if( name == "von_mises_hess_1_3") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_hess_1_3< std::vector<double> >() );
            if( name == "von_mises_hess_1_4") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_hess_1_4< std::vector<double> >() );
            if( name == "von_mises_hess_2_0") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_hess_2_0< std::vector<double> >() );
            if( name == "von_mises_hess_2_1") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_hess_2_1< std::vector<double> >() );
            if( name == "von_mises_hess_2_2") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_hess_2_2< std::vector<double> >() );
            if( name == "von_mises_hess_2_3") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_hess_2_3< std::vector<double> >() );
            if( name == "von_mises_hess_2_4") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_hess_2_4< std::vector<double> >() );
            if( name == "von_mises_hess_3_0") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_hess_3_0< std::vector<double> >() );
            if( name == "von_mises_hess_3_1") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_hess_3_1< std::vector<double> >() );
            if( name == "von_mises_hess_3_2") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_hess_3_2< std::vector<double> >() );
            if( name == "von_mises_hess_3_3") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_hess_3_3< std::vector<double> >() );
            if( name == "von_mises_hess_3_4") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_hess_3_4< std::vector<double> >() );
            if( name == "von_mises_hess_4_0") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_hess_4_0< std::vector<double> >() );
            if( name == "von_mises_hess_4_1") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_hess_4_1< std::vector<double> >() );
            if( name == "von_mises_hess_4_2") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_hess_4_2< std::vector<double> >() );
            if( name == "von_mises_hess_4_3") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_hess_4_3< std::vector<double> >() );
            if( name == "von_mises_hess_4_4") simplefunc = PSimpleFunction< std::vector<double>, double >( von_mises_hess_4_4< std::vector<double> >() );
            if( name == "quadlog_f") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_f< std::vector<double> >() );
            if( name == "quadlog_grad_0") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_grad_0< std::vector<double> >() );
            if( name == "quadlog_grad_1") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_grad_1< std::vector<double> >() );
            if( name == "quadlog_grad_2") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_grad_2< std::vector<double> >() );
            if( name == "quadlog_grad_3") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_grad_3< std::vector<double> >() );
            if( name == "quadlog_grad_4") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_grad_4< std::vector<double> >() );
            if( name == "quadlog_grad_5") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_grad_5< std::vector<double> >() );
            if( name == "quadlog_grad_6") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_grad_6< std::vector<double> >() );
            if( name == "quadlog_hess_0_0") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_0_0< std::vector<double> >() );
            if( name == "quadlog_hess_0_1") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_0_1< std::vector<double> >() );
            if( name == "quadlog_hess_0_2") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_0_2< std::vector<double> >() );
            if( name == "quadlog_hess_0_3") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_0_3< std::vector<double> >() );
            if( name == "quadlog_hess_0_4") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_0_4< std::vector<double> >() );
            if( name == "quadlog_hess_0_5") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_0_5< std::vector<double> >() );
            if( name == "quadlog_hess_0_6") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_0_6< std::vector<double> >() );
            if( name == "quadlog_hess_1_0") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_1_0< std::vector<double> >() );
            if( name == "quadlog_hess_1_1") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_1_1< std::vector<double> >() );
            if( name == "quadlog_hess_1_2") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_1_2< std::vector<double> >() );
            if( name == "quadlog_hess_1_3") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_1_3< std::vector<double> >() );
            if( name == "quadlog_hess_1_4") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_1_4< std::vector<double> >() );
            if( name == "quadlog_hess_1_5") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_1_5< std::vector<double> >() );
            if( name == "quadlog_hess_1_6") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_1_6< std::vector<double> >() );
            if( name == "quadlog_hess_2_0") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_2_0< std::vector<double> >() );
            if( name == "quadlog_hess_2_1") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_2_1< std::vector<double> >() );
            if( name == "quadlog_hess_2_2") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_2_2< std::vector<double> >() );
            if( name == "quadlog_hess_2_3") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_2_3< std::vector<double> >() );
            if( name == "quadlog_hess_2_4") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_2_4< std::vector<double> >() );
            if( name == "quadlog_hess_2_5") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_2_5< std::vector<double> >() );
            if( name == "quadlog_hess_2_6") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_2_6< std::vector<double> >() );
            if( name == "quadlog_hess_3_0") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_3_0< std::vector<double> >() );
            if( name == "quadlog_hess_3_1") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_3_1< std::vector<double> >() );
            if( name == "quadlog_hess_3_2") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_3_2< std::vector<double> >() );
            if( name == "quadlog_hess_3_3") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_3_3< std::vector<double> >() );
            if( name == "quadlog_hess_3_4") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_3_4< std::vector<double> >() );
            if( name == "quadlog_hess_3_5") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_3_5< std::vector<double> >() );
            if( name == "quadlog_hess_3_6") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_3_6< std::vector<double> >() );
            if( name == "quadlog_hess_4_0") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_4_0< std::vector<double> >() );
            if( name == "quadlog_hess_4_1") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_4_1< std::vector<double> >() );
            if( name == "quadlog_hess_4_2") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_4_2< std::vector<double> >() );
            if( name == "quadlog_hess_4_3") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_4_3< std::vector<double> >() );
            if( name == "quadlog_hess_4_4") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_4_4< std::vector<double> >() );
            if( name == "quadlog_hess_4_5") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_4_5< std::vector<double> >() );
            if( name == "quadlog_hess_4_6") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_4_6< std::vector<double> >() );
            if( name == "quadlog_hess_5_0") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_5_0< std::vector<double> >() );
            if( name == "quadlog_hess_5_1") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_5_1< std::vector<double> >() );
            if( name == "quadlog_hess_5_2") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_5_2< std::vector<double> >() );
            if( name == "quadlog_hess_5_3") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_5_3< std::vector<double> >() );
            if( name == "quadlog_hess_5_4") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_5_4< std::vector<double> >() );
            if( name == "quadlog_hess_5_5") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_5_5< std::vector<double> >() );
            if( name == "quadlog_hess_5_6") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_5_6< std::vector<double> >() );
            if( name == "quadlog_hess_6_0") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_6_0< std::vector<double> >() );
            if( name == "quadlog_hess_6_1") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_6_1< std::vector<double> >() );
            if( name == "quadlog_hess_6_2") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_6_2< std::vector<double> >() );
            if( name == "quadlog_hess_6_3") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_6_3< std::vector<double> >() );
            if( name == "quadlog_hess_6_4") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_6_4< std::vector<double> >() );
            if( name == "quadlog_hess_6_5") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_6_5< std::vector<double> >() );
            if( name == "quadlog_hess_6_6") simplefunc = PSimpleFunction< std::vector<double>, double >( quadlog_hess_6_6< std::vector<double> >() );
        }
        void PLibrary::checkout( std::string name, PSimpleFunction< double*, double > &simplefunc)
        {
            if( name == "neohook_f") simplefunc = PSimpleFunction< double*, double >( neohook_f< double* >() );
            if( name == "neohook_grad_0") simplefunc = PSimpleFunction< double*, double >( neohook_grad_0< double* >() );
            if( name == "neohook_grad_1") simplefunc = PSimpleFunction< double*, double >( neohook_grad_1< double* >() );
            if( name == "neohook_grad_2") simplefunc = PSimpleFunction< double*, double >( neohook_grad_2< double* >() );
            if( name == "neohook_grad_3") simplefunc = PSimpleFunction< double*, double >( neohook_grad_3< double* >() );
            if( name == "neohook_grad_4") simplefunc = PSimpleFunction< double*, double >( neohook_grad_4< double* >() );
            if( name == "neohook_grad_5") simplefunc = PSimpleFunction< double*, double >( neohook_grad_5< double* >() );
            if( name == "neohook_grad_6") simplefunc = PSimpleFunction< double*, double >( neohook_grad_6< double* >() );
            if( name == "neohook_hess_0_0") simplefunc = PSimpleFunction< double*, double >( neohook_hess_0_0< double* >() );
            if( name == "neohook_hess_0_1") simplefunc = PSimpleFunction< double*, double >( neohook_hess_0_1< double* >() );
            if( name == "neohook_hess_0_2") simplefunc = PSimpleFunction< double*, double >( neohook_hess_0_2< double* >() );
            if( name == "neohook_hess_0_3") simplefunc = PSimpleFunction< double*, double >( neohook_hess_0_3< double* >() );
            if( name == "neohook_hess_0_4") simplefunc = PSimpleFunction< double*, double >( neohook_hess_0_4< double* >() );
            if( name == "neohook_hess_0_5") simplefunc = PSimpleFunction< double*, double >( neohook_hess_0_5< double* >() );
            if( name == "neohook_hess_0_6") simplefunc = PSimpleFunction< double*, double >( neohook_hess_0_6< double* >() );
            if( name == "neohook_hess_1_0") simplefunc = PSimpleFunction< double*, double >( neohook_hess_1_0< double* >() );
            if( name == "neohook_hess_1_1") simplefunc = PSimpleFunction< double*, double >( neohook_hess_1_1< double* >() );
            if( name == "neohook_hess_1_2") simplefunc = PSimpleFunction< double*, double >( neohook_hess_1_2< double* >() );
            if( name == "neohook_hess_1_3") simplefunc = PSimpleFunction< double*, double >( neohook_hess_1_3< double* >() );
            if( name == "neohook_hess_1_4") simplefunc = PSimpleFunction< double*, double >( neohook_hess_1_4< double* >() );
            if( name == "neohook_hess_1_5") simplefunc = PSimpleFunction< double*, double >( neohook_hess_1_5< double* >() );
            if( name == "neohook_hess_1_6") simplefunc = PSimpleFunction< double*, double >( neohook_hess_1_6< double* >() );
            if( name == "neohook_hess_2_0") simplefunc = PSimpleFunction< double*, double >( neohook_hess_2_0< double* >() );
            if( name == "neohook_hess_2_1") simplefunc = PSimpleFunction< double*, double >( neohook_hess_2_1< double* >() );
            if( name == "neohook_hess_2_2") simplefunc = PSimpleFunction< double*, double >( neohook_hess_2_2< double* >() );
            if( name == "neohook_hess_2_3") simplefunc = PSimpleFunction< double*, double >( neohook_hess_2_3< double* >() );
            if( name == "neohook_hess_2_4") simplefunc = PSimpleFunction< double*, double >( neohook_hess_2_4< double* >() );
            if( name == "neohook_hess_2_5") simplefunc = PSimpleFunction< double*, double >( neohook_hess_2_5< double* >() );
            if( name == "neohook_hess_2_6") simplefunc = PSimpleFunction< double*, double >( neohook_hess_2_6< double* >() );
            if( name == "neohook_hess_3_0") simplefunc = PSimpleFunction< double*, double >( neohook_hess_3_0< double* >() );
            if( name == "neohook_hess_3_1") simplefunc = PSimpleFunction< double*, double >( neohook_hess_3_1< double* >() );
            if( name == "neohook_hess_3_2") simplefunc = PSimpleFunction< double*, double >( neohook_hess_3_2< double* >() );
            if( name == "neohook_hess_3_3") simplefunc = PSimpleFunction< double*, double >( neohook_hess_3_3< double* >() );
            if( name == "neohook_hess_3_4") simplefunc = PSimpleFunction< double*, double >( neohook_hess_3_4< double* >() );
            if( name == "neohook_hess_3_5") simplefunc = PSimpleFunction< double*, double >( neohook_hess_3_5< double* >() );
            if( name == "neohook_hess_3_6") simplefunc = PSimpleFunction< double*, double >( neohook_hess_3_6< double* >() );
            if( name == "neohook_hess_4_0") simplefunc = PSimpleFunction< double*, double >( neohook_hess_4_0< double* >() );
            if( name == "neohook_hess_4_1") simplefunc = PSimpleFunction< double*, double >( neohook_hess_4_1< double* >() );
            if( name == "neohook_hess_4_2") simplefunc = PSimpleFunction< double*, double >( neohook_hess_4_2< double* >() );
            if( name == "neohook_hess_4_3") simplefunc = PSimpleFunction< double*, double >( neohook_hess_4_3< double* >() );
            if( name == "neohook_hess_4_4") simplefunc = PSimpleFunction< double*, double >( neohook_hess_4_4< double* >() );
            if( name == "neohook_hess_4_5") simplefunc = PSimpleFunction< double*, double >( neohook_hess_4_5< double* >() );
            if( name == "neohook_hess_4_6") simplefunc = PSimpleFunction< double*, double >( neohook_hess_4_6< double* >() );
            if( name == "neohook_hess_5_0") simplefunc = PSimpleFunction< double*, double >( neohook_hess_5_0< double* >() );
            if( name == "neohook_hess_5_1") simplefunc = PSimpleFunction< double*, double >( neohook_hess_5_1< double* >() );
            if( name == "neohook_hess_5_2") simplefunc = PSimpleFunction< double*, double >( neohook_hess_5_2< double* >() );
            if( name == "neohook_hess_5_3") simplefunc = PSimpleFunction< double*, double >( neohook_hess_5_3< double* >() );
            if( name == "neohook_hess_5_4") simplefunc = PSimpleFunction< double*, double >( neohook_hess_5_4< double* >() );
            if( name == "neohook_hess_5_5") simplefunc = PSimpleFunction< double*, double >( neohook_hess_5_5< double* >() );
            if( name == "neohook_hess_5_6") simplefunc = PSimpleFunction< double*, double >( neohook_hess_5_6< double* >() );
            if( name == "neohook_hess_6_0") simplefunc = PSimpleFunction< double*, double >( neohook_hess_6_0< double* >() );
            if( name == "neohook_hess_6_1") simplefunc = PSimpleFunction< double*, double >( neohook_hess_6_1< double* >() );
            if( name == "neohook_hess_6_2") simplefunc = PSimpleFunction< double*, double >( neohook_hess_6_2< double* >() );
            if( name == "neohook_hess_6_3") simplefunc = PSimpleFunction< double*, double >( neohook_hess_6_3< double* >() );
            if( name == "neohook_hess_6_4") simplefunc = PSimpleFunction< double*, double >( neohook_hess_6_4< double* >() );
            if( name == "neohook_hess_6_5") simplefunc = PSimpleFunction< double*, double >( neohook_hess_6_5< double* >() );
            if( name == "neohook_hess_6_6") simplefunc = PSimpleFunction< double*, double >( neohook_hess_6_6< double* >() );
            if( name == "stvenkir_f") simplefunc = PSimpleFunction< double*, double >( stvenkir_f< double* >() );
            if( name == "stvenkir_grad_0") simplefunc = PSimpleFunction< double*, double >( stvenkir_grad_0< double* >() );
            if( name == "stvenkir_grad_1") simplefunc = PSimpleFunction< double*, double >( stvenkir_grad_1< double* >() );
            if( name == "stvenkir_grad_2") simplefunc = PSimpleFunction< double*, double >( stvenkir_grad_2< double* >() );
            if( name == "stvenkir_grad_3") simplefunc = PSimpleFunction< double*, double >( stvenkir_grad_3< double* >() );
            if( name == "stvenkir_grad_4") simplefunc = PSimpleFunction< double*, double >( stvenkir_grad_4< double* >() );
            if( name == "stvenkir_grad_5") simplefunc = PSimpleFunction< double*, double >( stvenkir_grad_5< double* >() );
            if( name == "stvenkir_grad_6") simplefunc = PSimpleFunction< double*, double >( stvenkir_grad_6< double* >() );
            if( name == "stvenkir_hess_0_0") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_0_0< double* >() );
            if( name == "stvenkir_hess_0_1") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_0_1< double* >() );
            if( name == "stvenkir_hess_0_2") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_0_2< double* >() );
            if( name == "stvenkir_hess_0_3") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_0_3< double* >() );
            if( name == "stvenkir_hess_0_4") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_0_4< double* >() );
            if( name == "stvenkir_hess_0_5") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_0_5< double* >() );
            if( name == "stvenkir_hess_0_6") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_0_6< double* >() );
            if( name == "stvenkir_hess_1_0") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_1_0< double* >() );
            if( name == "stvenkir_hess_1_1") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_1_1< double* >() );
            if( name == "stvenkir_hess_1_2") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_1_2< double* >() );
            if( name == "stvenkir_hess_1_3") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_1_3< double* >() );
            if( name == "stvenkir_hess_1_4") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_1_4< double* >() );
            if( name == "stvenkir_hess_1_5") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_1_5< double* >() );
            if( name == "stvenkir_hess_1_6") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_1_6< double* >() );
            if( name == "stvenkir_hess_2_0") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_2_0< double* >() );
            if( name == "stvenkir_hess_2_1") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_2_1< double* >() );
            if( name == "stvenkir_hess_2_2") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_2_2< double* >() );
            if( name == "stvenkir_hess_2_3") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_2_3< double* >() );
            if( name == "stvenkir_hess_2_4") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_2_4< double* >() );
            if( name == "stvenkir_hess_2_5") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_2_5< double* >() );
            if( name == "stvenkir_hess_2_6") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_2_6< double* >() );
            if( name == "stvenkir_hess_3_0") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_3_0< double* >() );
            if( name == "stvenkir_hess_3_1") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_3_1< double* >() );
            if( name == "stvenkir_hess_3_2") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_3_2< double* >() );
            if( name == "stvenkir_hess_3_3") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_3_3< double* >() );
            if( name == "stvenkir_hess_3_4") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_3_4< double* >() );
            if( name == "stvenkir_hess_3_5") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_3_5< double* >() );
            if( name == "stvenkir_hess_3_6") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_3_6< double* >() );
            if( name == "stvenkir_hess_4_0") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_4_0< double* >() );
            if( name == "stvenkir_hess_4_1") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_4_1< double* >() );
            if( name == "stvenkir_hess_4_2") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_4_2< double* >() );
            if( name == "stvenkir_hess_4_3") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_4_3< double* >() );
            if( name == "stvenkir_hess_4_4") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_4_4< double* >() );
            if( name == "stvenkir_hess_4_5") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_4_5< double* >() );
            if( name == "stvenkir_hess_4_6") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_4_6< double* >() );
            if( name == "stvenkir_hess_5_0") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_5_0< double* >() );
            if( name == "stvenkir_hess_5_1") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_5_1< double* >() );
            if( name == "stvenkir_hess_5_2") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_5_2< double* >() );
            if( name == "stvenkir_hess_5_3") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_5_3< double* >() );
            if( name == "stvenkir_hess_5_4") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_5_4< double* >() );
            if( name == "stvenkir_hess_5_5") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_5_5< double* >() );
            if( name == "stvenkir_hess_5_6") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_5_6< double* >() );
            if( name == "stvenkir_hess_6_0") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_6_0< double* >() );
            if( name == "stvenkir_hess_6_1") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_6_1< double* >() );
            if( name == "stvenkir_hess_6_2") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_6_2< double* >() );
            if( name == "stvenkir_hess_6_3") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_6_3< double* >() );
            if( name == "stvenkir_hess_6_4") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_6_4< double* >() );
            if( name == "stvenkir_hess_6_5") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_6_5< double* >() );
            if( name == "stvenkir_hess_6_6") simplefunc = PSimpleFunction< double*, double >( stvenkir_hess_6_6< double* >() );
            if( name == "linear_hardening_f") simplefunc = PSimpleFunction< double*, double >( linear_hardening_f< double* >() );
            if( name == "linear_hardening_grad_0") simplefunc = PSimpleFunction< double*, double >( linear_hardening_grad_0< double* >() );
            if( name == "linear_hardening_grad_1") simplefunc = PSimpleFunction< double*, double >( linear_hardening_grad_1< double* >() );
            if( name == "linear_hardening_hess_0_0") simplefunc = PSimpleFunction< double*, double >( linear_hardening_hess_0_0< double* >() );
            if( name == "linear_hardening_hess_0_1") simplefunc = PSimpleFunction< double*, double >( linear_hardening_hess_0_1< double* >() );
            if( name == "linear_hardening_hess_1_0") simplefunc = PSimpleFunction< double*, double >( linear_hardening_hess_1_0< double* >() );
            if( name == "linear_hardening_hess_1_1") simplefunc = PSimpleFunction< double*, double >( linear_hardening_hess_1_1< double* >() );
            if( name == "Cu_hardening_f") simplefunc = PSimpleFunction< double*, double >( Cu_hardening_f< double* >() );
            if( name == "Cu_hardening_grad_0") simplefunc = PSimpleFunction< double*, double >( Cu_hardening_grad_0< double* >() );
            if( name == "Cu_hardening_grad_1") simplefunc = PSimpleFunction< double*, double >( Cu_hardening_grad_1< double* >() );
            if( name == "Cu_hardening_hess_0_0") simplefunc = PSimpleFunction< double*, double >( Cu_hardening_hess_0_0< double* >() );
            if( name == "Cu_hardening_hess_0_1") simplefunc = PSimpleFunction< double*, double >( Cu_hardening_hess_0_1< double* >() );
            if( name == "Cu_hardening_hess_1_0") simplefunc = PSimpleFunction< double*, double >( Cu_hardening_hess_1_0< double* >() );
            if( name == "Cu_hardening_hess_1_1") simplefunc = PSimpleFunction< double*, double >( Cu_hardening_hess_1_1< double* >() );
            if( name == "von_mises_f") simplefunc = PSimpleFunction< double*, double >( von_mises_f< double* >() );
            if( name == "von_mises_grad_0") simplefunc = PSimpleFunction< double*, double >( von_mises_grad_0< double* >() );
            if( name == "von_mises_grad_1") simplefunc = PSimpleFunction< double*, double >( von_mises_grad_1< double* >() );
            if( name == "von_mises_grad_2") simplefunc = PSimpleFunction< double*, double >( von_mises_grad_2< double* >() );
            if( name == "von_mises_grad_3") simplefunc = PSimpleFunction< double*, double >( von_mises_grad_3< double* >() );
            if( name == "von_mises_grad_4") simplefunc = PSimpleFunction< double*, double >( von_mises_grad_4< double* >() );
            if( name == "von_mises_hess_0_0") simplefunc = PSimpleFunction< double*, double >( von_mises_hess_0_0< double* >() );
            if( name == "von_mises_hess_0_1") simplefunc = PSimpleFunction< double*, double >( von_mises_hess_0_1< double* >() );
            if( name == "von_mises_hess_0_2") simplefunc = PSimpleFunction< double*, double >( von_mises_hess_0_2< double* >() );
            if( name == "von_mises_hess_0_3") simplefunc = PSimpleFunction< double*, double >( von_mises_hess_0_3< double* >() );
            if( name == "von_mises_hess_0_4") simplefunc = PSimpleFunction< double*, double >( von_mises_hess_0_4< double* >() );
            if( name == "von_mises_hess_1_0") simplefunc = PSimpleFunction< double*, double >( von_mises_hess_1_0< double* >() );
            if( name == "von_mises_hess_1_1") simplefunc = PSimpleFunction< double*, double >( von_mises_hess_1_1< double* >() );
            if( name == "von_mises_hess_1_2") simplefunc = PSimpleFunction< double*, double >( von_mises_hess_1_2< double* >() );
            if( name == "von_mises_hess_1_3") simplefunc = PSimpleFunction< double*, double >( von_mises_hess_1_3< double* >() );
            if( name == "von_mises_hess_1_4") simplefunc = PSimpleFunction< double*, double >( von_mises_hess_1_4< double* >() );
            if( name == "von_mises_hess_2_0") simplefunc = PSimpleFunction< double*, double >( von_mises_hess_2_0< double* >() );
            if( name == "von_mises_hess_2_1") simplefunc = PSimpleFunction< double*, double >( von_mises_hess_2_1< double* >() );
            if( name == "von_mises_hess_2_2") simplefunc = PSimpleFunction< double*, double >( von_mises_hess_2_2< double* >() );
            if( name == "von_mises_hess_2_3") simplefunc = PSimpleFunction< double*, double >( von_mises_hess_2_3< double* >() );
            if( name == "von_mises_hess_2_4") simplefunc = PSimpleFunction< double*, double >( von_mises_hess_2_4< double* >() );
            if( name == "von_mises_hess_3_0") simplefunc = PSimpleFunction< double*, double >( von_mises_hess_3_0< double* >() );
            if( name == "von_mises_hess_3_1") simplefunc = PSimpleFunction< double*, double >( von_mises_hess_3_1< double* >() );
            if( name == "von_mises_hess_3_2") simplefunc = PSimpleFunction< double*, double >( von_mises_hess_3_2< double* >() );
            if( name == "von_mises_hess_3_3") simplefunc = PSimpleFunction< double*, double >( von_mises_hess_3_3< double* >() );
            if( name == "von_mises_hess_3_4") simplefunc = PSimpleFunction< double*, double >( von_mises_hess_3_4< double* >() );
            if( name == "von_mises_hess_4_0") simplefunc = PSimpleFunction< double*, double >( von_mises_hess_4_0< double* >() );
            if( name == "von_mises_hess_4_1") simplefunc = PSimpleFunction< double*, double >( von_mises_hess_4_1< double* >() );
            if( name == "von_mises_hess_4_2") simplefunc = PSimpleFunction< double*, double >( von_mises_hess_4_2< double* >() );
            if( name == "von_mises_hess_4_3") simplefunc = PSimpleFunction< double*, double >( von_mises_hess_4_3< double* >() );
            if( name == "von_mises_hess_4_4") simplefunc = PSimpleFunction< double*, double >( von_mises_hess_4_4< double* >() );
            if( name == "quadlog_f") simplefunc = PSimpleFunction< double*, double >( quadlog_f< double* >() );
            if( name == "quadlog_grad_0") simplefunc = PSimpleFunction< double*, double >( quadlog_grad_0< double* >() );
            if( name == "quadlog_grad_1") simplefunc = PSimpleFunction< double*, double >( quadlog_grad_1< double* >() );
            if( name == "quadlog_grad_2") simplefunc = PSimpleFunction< double*, double >( quadlog_grad_2< double* >() );
            if( name == "quadlog_grad_3") simplefunc = PSimpleFunction< double*, double >( quadlog_grad_3< double* >() );
            if( name == "quadlog_grad_4") simplefunc = PSimpleFunction< double*, double >( quadlog_grad_4< double* >() );
            if( name == "quadlog_grad_5") simplefunc = PSimpleFunction< double*, double >( quadlog_grad_5< double* >() );
            if( name == "quadlog_grad_6") simplefunc = PSimpleFunction< double*, double >( quadlog_grad_6< double* >() );
            if( name == "quadlog_hess_0_0") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_0_0< double* >() );
            if( name == "quadlog_hess_0_1") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_0_1< double* >() );
            if( name == "quadlog_hess_0_2") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_0_2< double* >() );
            if( name == "quadlog_hess_0_3") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_0_3< double* >() );
            if( name == "quadlog_hess_0_4") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_0_4< double* >() );
            if( name == "quadlog_hess_0_5") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_0_5< double* >() );
            if( name == "quadlog_hess_0_6") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_0_6< double* >() );
            if( name == "quadlog_hess_1_0") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_1_0< double* >() );
            if( name == "quadlog_hess_1_1") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_1_1< double* >() );
            if( name == "quadlog_hess_1_2") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_1_2< double* >() );
            if( name == "quadlog_hess_1_3") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_1_3< double* >() );
            if( name == "quadlog_hess_1_4") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_1_4< double* >() );
            if( name == "quadlog_hess_1_5") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_1_5< double* >() );
            if( name == "quadlog_hess_1_6") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_1_6< double* >() );
            if( name == "quadlog_hess_2_0") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_2_0< double* >() );
            if( name == "quadlog_hess_2_1") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_2_1< double* >() );
            if( name == "quadlog_hess_2_2") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_2_2< double* >() );
            if( name == "quadlog_hess_2_3") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_2_3< double* >() );
            if( name == "quadlog_hess_2_4") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_2_4< double* >() );
            if( name == "quadlog_hess_2_5") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_2_5< double* >() );
            if( name == "quadlog_hess_2_6") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_2_6< double* >() );
            if( name == "quadlog_hess_3_0") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_3_0< double* >() );
            if( name == "quadlog_hess_3_1") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_3_1< double* >() );
            if( name == "quadlog_hess_3_2") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_3_2< double* >() );
            if( name == "quadlog_hess_3_3") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_3_3< double* >() );
            if( name == "quadlog_hess_3_4") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_3_4< double* >() );
            if( name == "quadlog_hess_3_5") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_3_5< double* >() );
            if( name == "quadlog_hess_3_6") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_3_6< double* >() );
            if( name == "quadlog_hess_4_0") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_4_0< double* >() );
            if( name == "quadlog_hess_4_1") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_4_1< double* >() );
            if( name == "quadlog_hess_4_2") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_4_2< double* >() );
            if( name == "quadlog_hess_4_3") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_4_3< double* >() );
            if( name == "quadlog_hess_4_4") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_4_4< double* >() );
            if( name == "quadlog_hess_4_5") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_4_5< double* >() );
            if( name == "quadlog_hess_4_6") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_4_6< double* >() );
            if( name == "quadlog_hess_5_0") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_5_0< double* >() );
            if( name == "quadlog_hess_5_1") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_5_1< double* >() );
            if( name == "quadlog_hess_5_2") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_5_2< double* >() );
            if( name == "quadlog_hess_5_3") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_5_3< double* >() );
            if( name == "quadlog_hess_5_4") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_5_4< double* >() );
            if( name == "quadlog_hess_5_5") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_5_5< double* >() );
            if( name == "quadlog_hess_5_6") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_5_6< double* >() );
            if( name == "quadlog_hess_6_0") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_6_0< double* >() );
            if( name == "quadlog_hess_6_1") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_6_1< double* >() );
            if( name == "quadlog_hess_6_2") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_6_2< double* >() );
            if( name == "quadlog_hess_6_3") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_6_3< double* >() );
            if( name == "quadlog_hess_6_4") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_6_4< double* >() );
            if( name == "quadlog_hess_6_5") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_6_5< double* >() );
            if( name == "quadlog_hess_6_6") simplefunc = PSimpleFunction< double*, double >( quadlog_hess_6_6< double* >() );
=======
            typedef PSimpleFunction< std::vector<double>, double > psf;
            if( name == "stvenkir_f") { simplefunc = psf( stvenkir_f< std::vector<double> >() ); return;}
            if( name == "stvenkir_grad_0") { simplefunc = psf( stvenkir_grad_0< std::vector<double> >() ); return;}
            if( name == "stvenkir_grad_1") { simplefunc = psf( stvenkir_grad_1< std::vector<double> >() ); return;}
            if( name == "stvenkir_grad_2") { simplefunc = psf( stvenkir_grad_2< std::vector<double> >() ); return;}
            if( name == "stvenkir_grad_3") { simplefunc = psf( stvenkir_grad_3< std::vector<double> >() ); return;}
            if( name == "stvenkir_grad_4") { simplefunc = psf( stvenkir_grad_4< std::vector<double> >() ); return;}
            if( name == "stvenkir_grad_5") { simplefunc = psf( stvenkir_grad_5< std::vector<double> >() ); return;}
            if( name == "stvenkir_grad_6") { simplefunc = psf( stvenkir_grad_6< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_0_0") { simplefunc = psf( stvenkir_hess_0_0< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_0_1") { simplefunc = psf( stvenkir_hess_0_1< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_0_2") { simplefunc = psf( stvenkir_hess_0_2< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_0_3") { simplefunc = psf( stvenkir_hess_0_3< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_0_4") { simplefunc = psf( stvenkir_hess_0_4< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_0_5") { simplefunc = psf( stvenkir_hess_0_5< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_0_6") { simplefunc = psf( stvenkir_hess_0_6< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_1_0") { simplefunc = psf( stvenkir_hess_1_0< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_1_1") { simplefunc = psf( stvenkir_hess_1_1< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_1_2") { simplefunc = psf( stvenkir_hess_1_2< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_1_3") { simplefunc = psf( stvenkir_hess_1_3< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_1_4") { simplefunc = psf( stvenkir_hess_1_4< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_1_5") { simplefunc = psf( stvenkir_hess_1_5< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_1_6") { simplefunc = psf( stvenkir_hess_1_6< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_2_0") { simplefunc = psf( stvenkir_hess_2_0< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_2_1") { simplefunc = psf( stvenkir_hess_2_1< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_2_2") { simplefunc = psf( stvenkir_hess_2_2< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_2_3") { simplefunc = psf( stvenkir_hess_2_3< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_2_4") { simplefunc = psf( stvenkir_hess_2_4< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_2_5") { simplefunc = psf( stvenkir_hess_2_5< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_2_6") { simplefunc = psf( stvenkir_hess_2_6< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_3_0") { simplefunc = psf( stvenkir_hess_3_0< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_3_1") { simplefunc = psf( stvenkir_hess_3_1< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_3_2") { simplefunc = psf( stvenkir_hess_3_2< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_3_3") { simplefunc = psf( stvenkir_hess_3_3< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_3_4") { simplefunc = psf( stvenkir_hess_3_4< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_3_5") { simplefunc = psf( stvenkir_hess_3_5< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_3_6") { simplefunc = psf( stvenkir_hess_3_6< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_4_0") { simplefunc = psf( stvenkir_hess_4_0< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_4_1") { simplefunc = psf( stvenkir_hess_4_1< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_4_2") { simplefunc = psf( stvenkir_hess_4_2< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_4_3") { simplefunc = psf( stvenkir_hess_4_3< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_4_4") { simplefunc = psf( stvenkir_hess_4_4< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_4_5") { simplefunc = psf( stvenkir_hess_4_5< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_4_6") { simplefunc = psf( stvenkir_hess_4_6< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_5_0") { simplefunc = psf( stvenkir_hess_5_0< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_5_1") { simplefunc = psf( stvenkir_hess_5_1< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_5_2") { simplefunc = psf( stvenkir_hess_5_2< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_5_3") { simplefunc = psf( stvenkir_hess_5_3< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_5_4") { simplefunc = psf( stvenkir_hess_5_4< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_5_5") { simplefunc = psf( stvenkir_hess_5_5< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_5_6") { simplefunc = psf( stvenkir_hess_5_6< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_6_0") { simplefunc = psf( stvenkir_hess_6_0< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_6_1") { simplefunc = psf( stvenkir_hess_6_1< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_6_2") { simplefunc = psf( stvenkir_hess_6_2< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_6_3") { simplefunc = psf( stvenkir_hess_6_3< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_6_4") { simplefunc = psf( stvenkir_hess_6_4< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_6_5") { simplefunc = psf( stvenkir_hess_6_5< std::vector<double> >() ); return;}
            if( name == "stvenkir_hess_6_6") { simplefunc = psf( stvenkir_hess_6_6< std::vector<double> >() ); return;}
            if( name == "quadlog_f") { simplefunc = psf( quadlog_f< std::vector<double> >() ); return;}
            if( name == "quadlog_grad_0") { simplefunc = psf( quadlog_grad_0< std::vector<double> >() ); return;}
            if( name == "quadlog_grad_1") { simplefunc = psf( quadlog_grad_1< std::vector<double> >() ); return;}
            if( name == "quadlog_grad_2") { simplefunc = psf( quadlog_grad_2< std::vector<double> >() ); return;}
            if( name == "quadlog_grad_3") { simplefunc = psf( quadlog_grad_3< std::vector<double> >() ); return;}
            if( name == "quadlog_grad_4") { simplefunc = psf( quadlog_grad_4< std::vector<double> >() ); return;}
            if( name == "quadlog_grad_5") { simplefunc = psf( quadlog_grad_5< std::vector<double> >() ); return;}
            if( name == "quadlog_grad_6") { simplefunc = psf( quadlog_grad_6< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_0_0") { simplefunc = psf( quadlog_hess_0_0< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_0_1") { simplefunc = psf( quadlog_hess_0_1< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_0_2") { simplefunc = psf( quadlog_hess_0_2< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_0_3") { simplefunc = psf( quadlog_hess_0_3< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_0_4") { simplefunc = psf( quadlog_hess_0_4< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_0_5") { simplefunc = psf( quadlog_hess_0_5< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_0_6") { simplefunc = psf( quadlog_hess_0_6< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_1_0") { simplefunc = psf( quadlog_hess_1_0< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_1_1") { simplefunc = psf( quadlog_hess_1_1< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_1_2") { simplefunc = psf( quadlog_hess_1_2< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_1_3") { simplefunc = psf( quadlog_hess_1_3< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_1_4") { simplefunc = psf( quadlog_hess_1_4< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_1_5") { simplefunc = psf( quadlog_hess_1_5< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_1_6") { simplefunc = psf( quadlog_hess_1_6< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_2_0") { simplefunc = psf( quadlog_hess_2_0< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_2_1") { simplefunc = psf( quadlog_hess_2_1< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_2_2") { simplefunc = psf( quadlog_hess_2_2< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_2_3") { simplefunc = psf( quadlog_hess_2_3< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_2_4") { simplefunc = psf( quadlog_hess_2_4< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_2_5") { simplefunc = psf( quadlog_hess_2_5< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_2_6") { simplefunc = psf( quadlog_hess_2_6< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_3_0") { simplefunc = psf( quadlog_hess_3_0< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_3_1") { simplefunc = psf( quadlog_hess_3_1< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_3_2") { simplefunc = psf( quadlog_hess_3_2< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_3_3") { simplefunc = psf( quadlog_hess_3_3< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_3_4") { simplefunc = psf( quadlog_hess_3_4< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_3_5") { simplefunc = psf( quadlog_hess_3_5< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_3_6") { simplefunc = psf( quadlog_hess_3_6< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_4_0") { simplefunc = psf( quadlog_hess_4_0< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_4_1") { simplefunc = psf( quadlog_hess_4_1< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_4_2") { simplefunc = psf( quadlog_hess_4_2< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_4_3") { simplefunc = psf( quadlog_hess_4_3< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_4_4") { simplefunc = psf( quadlog_hess_4_4< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_4_5") { simplefunc = psf( quadlog_hess_4_5< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_4_6") { simplefunc = psf( quadlog_hess_4_6< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_5_0") { simplefunc = psf( quadlog_hess_5_0< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_5_1") { simplefunc = psf( quadlog_hess_5_1< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_5_2") { simplefunc = psf( quadlog_hess_5_2< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_5_3") { simplefunc = psf( quadlog_hess_5_3< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_5_4") { simplefunc = psf( quadlog_hess_5_4< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_5_5") { simplefunc = psf( quadlog_hess_5_5< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_5_6") { simplefunc = psf( quadlog_hess_5_6< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_6_0") { simplefunc = psf( quadlog_hess_6_0< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_6_1") { simplefunc = psf( quadlog_hess_6_1< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_6_2") { simplefunc = psf( quadlog_hess_6_2< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_6_3") { simplefunc = psf( quadlog_hess_6_3< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_6_4") { simplefunc = psf( quadlog_hess_6_4< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_6_5") { simplefunc = psf( quadlog_hess_6_5< std::vector<double> >() ); return;}
            if( name == "quadlog_hess_6_6") { simplefunc = psf( quadlog_hess_6_6< std::vector<double> >() ); return;}
            if( name == "hardening_f") { simplefunc = psf( hardening_f< std::vector<double> >() ); return;}
            if( name == "hardening_grad_0") { simplefunc = psf( hardening_grad_0< std::vector<double> >() ); return;}
            if( name == "hardening_hess_0_0") { simplefunc = psf( hardening_hess_0_0< std::vector<double> >() ); return;}
            if( name == "von_mises_f") { simplefunc = psf( von_mises_f< std::vector<double> >() ); return;}
            if( name == "von_mises_grad_0") { simplefunc = psf( von_mises_grad_0< std::vector<double> >() ); return;}
            if( name == "von_mises_grad_1") { simplefunc = psf( von_mises_grad_1< std::vector<double> >() ); return;}
            if( name == "von_mises_grad_2") { simplefunc = psf( von_mises_grad_2< std::vector<double> >() ); return;}
            if( name == "von_mises_grad_3") { simplefunc = psf( von_mises_grad_3< std::vector<double> >() ); return;}
            if( name == "von_mises_grad_4") { simplefunc = psf( von_mises_grad_4< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_0_0") { simplefunc = psf( von_mises_hess_0_0< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_0_1") { simplefunc = psf( von_mises_hess_0_1< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_0_2") { simplefunc = psf( von_mises_hess_0_2< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_0_3") { simplefunc = psf( von_mises_hess_0_3< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_0_4") { simplefunc = psf( von_mises_hess_0_4< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_1_0") { simplefunc = psf( von_mises_hess_1_0< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_1_1") { simplefunc = psf( von_mises_hess_1_1< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_1_2") { simplefunc = psf( von_mises_hess_1_2< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_1_3") { simplefunc = psf( von_mises_hess_1_3< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_1_4") { simplefunc = psf( von_mises_hess_1_4< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_2_0") { simplefunc = psf( von_mises_hess_2_0< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_2_1") { simplefunc = psf( von_mises_hess_2_1< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_2_2") { simplefunc = psf( von_mises_hess_2_2< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_2_3") { simplefunc = psf( von_mises_hess_2_3< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_2_4") { simplefunc = psf( von_mises_hess_2_4< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_3_0") { simplefunc = psf( von_mises_hess_3_0< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_3_1") { simplefunc = psf( von_mises_hess_3_1< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_3_2") { simplefunc = psf( von_mises_hess_3_2< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_3_3") { simplefunc = psf( von_mises_hess_3_3< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_3_4") { simplefunc = psf( von_mises_hess_3_4< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_4_0") { simplefunc = psf( von_mises_hess_4_0< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_4_1") { simplefunc = psf( von_mises_hess_4_1< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_4_2") { simplefunc = psf( von_mises_hess_4_2< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_4_3") { simplefunc = psf( von_mises_hess_4_3< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_4_4") { simplefunc = psf( von_mises_hess_4_4< std::vector<double> >() ); return;}
            if( name == "neohook_f") { simplefunc = psf( neohook_f< std::vector<double> >() ); return;}
            if( name == "neohook_grad_0") { simplefunc = psf( neohook_grad_0< std::vector<double> >() ); return;}
            if( name == "neohook_grad_1") { simplefunc = psf( neohook_grad_1< std::vector<double> >() ); return;}
            if( name == "neohook_grad_2") { simplefunc = psf( neohook_grad_2< std::vector<double> >() ); return;}
            if( name == "neohook_grad_3") { simplefunc = psf( neohook_grad_3< std::vector<double> >() ); return;}
            if( name == "neohook_grad_4") { simplefunc = psf( neohook_grad_4< std::vector<double> >() ); return;}
            if( name == "neohook_grad_5") { simplefunc = psf( neohook_grad_5< std::vector<double> >() ); return;}
            if( name == "neohook_grad_6") { simplefunc = psf( neohook_grad_6< std::vector<double> >() ); return;}
            if( name == "neohook_hess_0_0") { simplefunc = psf( neohook_hess_0_0< std::vector<double> >() ); return;}
            if( name == "neohook_hess_0_1") { simplefunc = psf( neohook_hess_0_1< std::vector<double> >() ); return;}
            if( name == "neohook_hess_0_2") { simplefunc = psf( neohook_hess_0_2< std::vector<double> >() ); return;}
            if( name == "neohook_hess_0_3") { simplefunc = psf( neohook_hess_0_3< std::vector<double> >() ); return;}
            if( name == "neohook_hess_0_4") { simplefunc = psf( neohook_hess_0_4< std::vector<double> >() ); return;}
            if( name == "neohook_hess_0_5") { simplefunc = psf( neohook_hess_0_5< std::vector<double> >() ); return;}
            if( name == "neohook_hess_0_6") { simplefunc = psf( neohook_hess_0_6< std::vector<double> >() ); return;}
            if( name == "neohook_hess_1_0") { simplefunc = psf( neohook_hess_1_0< std::vector<double> >() ); return;}
            if( name == "neohook_hess_1_1") { simplefunc = psf( neohook_hess_1_1< std::vector<double> >() ); return;}
            if( name == "neohook_hess_1_2") { simplefunc = psf( neohook_hess_1_2< std::vector<double> >() ); return;}
            if( name == "neohook_hess_1_3") { simplefunc = psf( neohook_hess_1_3< std::vector<double> >() ); return;}
            if( name == "neohook_hess_1_4") { simplefunc = psf( neohook_hess_1_4< std::vector<double> >() ); return;}
            if( name == "neohook_hess_1_5") { simplefunc = psf( neohook_hess_1_5< std::vector<double> >() ); return;}
            if( name == "neohook_hess_1_6") { simplefunc = psf( neohook_hess_1_6< std::vector<double> >() ); return;}
            if( name == "neohook_hess_2_0") { simplefunc = psf( neohook_hess_2_0< std::vector<double> >() ); return;}
            if( name == "neohook_hess_2_1") { simplefunc = psf( neohook_hess_2_1< std::vector<double> >() ); return;}
            if( name == "neohook_hess_2_2") { simplefunc = psf( neohook_hess_2_2< std::vector<double> >() ); return;}
            if( name == "neohook_hess_2_3") { simplefunc = psf( neohook_hess_2_3< std::vector<double> >() ); return;}
            if( name == "neohook_hess_2_4") { simplefunc = psf( neohook_hess_2_4< std::vector<double> >() ); return;}
            if( name == "neohook_hess_2_5") { simplefunc = psf( neohook_hess_2_5< std::vector<double> >() ); return;}
            if( name == "neohook_hess_2_6") { simplefunc = psf( neohook_hess_2_6< std::vector<double> >() ); return;}
            if( name == "neohook_hess_3_0") { simplefunc = psf( neohook_hess_3_0< std::vector<double> >() ); return;}
            if( name == "neohook_hess_3_1") { simplefunc = psf( neohook_hess_3_1< std::vector<double> >() ); return;}
            if( name == "neohook_hess_3_2") { simplefunc = psf( neohook_hess_3_2< std::vector<double> >() ); return;}
            if( name == "neohook_hess_3_3") { simplefunc = psf( neohook_hess_3_3< std::vector<double> >() ); return;}
            if( name == "neohook_hess_3_4") { simplefunc = psf( neohook_hess_3_4< std::vector<double> >() ); return;}
            if( name == "neohook_hess_3_5") { simplefunc = psf( neohook_hess_3_5< std::vector<double> >() ); return;}
            if( name == "neohook_hess_3_6") { simplefunc = psf( neohook_hess_3_6< std::vector<double> >() ); return;}
            if( name == "neohook_hess_4_0") { simplefunc = psf( neohook_hess_4_0< std::vector<double> >() ); return;}
            if( name == "neohook_hess_4_1") { simplefunc = psf( neohook_hess_4_1< std::vector<double> >() ); return;}
            if( name == "neohook_hess_4_2") { simplefunc = psf( neohook_hess_4_2< std::vector<double> >() ); return;}
            if( name == "neohook_hess_4_3") { simplefunc = psf( neohook_hess_4_3< std::vector<double> >() ); return;}
            if( name == "neohook_hess_4_4") { simplefunc = psf( neohook_hess_4_4< std::vector<double> >() ); return;}
            if( name == "neohook_hess_4_5") { simplefunc = psf( neohook_hess_4_5< std::vector<double> >() ); return;}
            if( name == "neohook_hess_4_6") { simplefunc = psf( neohook_hess_4_6< std::vector<double> >() ); return;}
            if( name == "neohook_hess_5_0") { simplefunc = psf( neohook_hess_5_0< std::vector<double> >() ); return;}
            if( name == "neohook_hess_5_1") { simplefunc = psf( neohook_hess_5_1< std::vector<double> >() ); return;}
            if( name == "neohook_hess_5_2") { simplefunc = psf( neohook_hess_5_2< std::vector<double> >() ); return;}
            if( name == "neohook_hess_5_3") { simplefunc = psf( neohook_hess_5_3< std::vector<double> >() ); return;}
            if( name == "neohook_hess_5_4") { simplefunc = psf( neohook_hess_5_4< std::vector<double> >() ); return;}
            if( name == "neohook_hess_5_5") { simplefunc = psf( neohook_hess_5_5< std::vector<double> >() ); return;}
            if( name == "neohook_hess_5_6") { simplefunc = psf( neohook_hess_5_6< std::vector<double> >() ); return;}
            if( name == "neohook_hess_6_0") { simplefunc = psf( neohook_hess_6_0< std::vector<double> >() ); return;}
            if( name == "neohook_hess_6_1") { simplefunc = psf( neohook_hess_6_1< std::vector<double> >() ); return;}
            if( name == "neohook_hess_6_2") { simplefunc = psf( neohook_hess_6_2< std::vector<double> >() ); return;}
            if( name == "neohook_hess_6_3") { simplefunc = psf( neohook_hess_6_3< std::vector<double> >() ); return;}
            if( name == "neohook_hess_6_4") { simplefunc = psf( neohook_hess_6_4< std::vector<double> >() ); return;}
            if( name == "neohook_hess_6_5") { simplefunc = psf( neohook_hess_6_5< std::vector<double> >() ); return;}
            if( name == "neohook_hess_6_6") { simplefunc = psf( neohook_hess_6_6< std::vector<double> >() ); return;}
            else throw std::runtime_error( "PSimpleFunction< std::vector<double>, double > " + name + " was not found in the PLibrary");
        }
        void PLibrary::checkout( std::string name, PSimpleFunction< double*, double > &simplefunc)
        {
            typedef PSimpleFunction< double*, double > psf;
            if( name == "stvenkir_f") { simplefunc = psf( stvenkir_f< double* >() ); return;}
            if( name == "stvenkir_grad_0") { simplefunc = psf( stvenkir_grad_0< double* >() ); return;}
            if( name == "stvenkir_grad_1") { simplefunc = psf( stvenkir_grad_1< double* >() ); return;}
            if( name == "stvenkir_grad_2") { simplefunc = psf( stvenkir_grad_2< double* >() ); return;}
            if( name == "stvenkir_grad_3") { simplefunc = psf( stvenkir_grad_3< double* >() ); return;}
            if( name == "stvenkir_grad_4") { simplefunc = psf( stvenkir_grad_4< double* >() ); return;}
            if( name == "stvenkir_grad_5") { simplefunc = psf( stvenkir_grad_5< double* >() ); return;}
            if( name == "stvenkir_grad_6") { simplefunc = psf( stvenkir_grad_6< double* >() ); return;}
            if( name == "stvenkir_hess_0_0") { simplefunc = psf( stvenkir_hess_0_0< double* >() ); return;}
            if( name == "stvenkir_hess_0_1") { simplefunc = psf( stvenkir_hess_0_1< double* >() ); return;}
            if( name == "stvenkir_hess_0_2") { simplefunc = psf( stvenkir_hess_0_2< double* >() ); return;}
            if( name == "stvenkir_hess_0_3") { simplefunc = psf( stvenkir_hess_0_3< double* >() ); return;}
            if( name == "stvenkir_hess_0_4") { simplefunc = psf( stvenkir_hess_0_4< double* >() ); return;}
            if( name == "stvenkir_hess_0_5") { simplefunc = psf( stvenkir_hess_0_5< double* >() ); return;}
            if( name == "stvenkir_hess_0_6") { simplefunc = psf( stvenkir_hess_0_6< double* >() ); return;}
            if( name == "stvenkir_hess_1_0") { simplefunc = psf( stvenkir_hess_1_0< double* >() ); return;}
            if( name == "stvenkir_hess_1_1") { simplefunc = psf( stvenkir_hess_1_1< double* >() ); return;}
            if( name == "stvenkir_hess_1_2") { simplefunc = psf( stvenkir_hess_1_2< double* >() ); return;}
            if( name == "stvenkir_hess_1_3") { simplefunc = psf( stvenkir_hess_1_3< double* >() ); return;}
            if( name == "stvenkir_hess_1_4") { simplefunc = psf( stvenkir_hess_1_4< double* >() ); return;}
            if( name == "stvenkir_hess_1_5") { simplefunc = psf( stvenkir_hess_1_5< double* >() ); return;}
            if( name == "stvenkir_hess_1_6") { simplefunc = psf( stvenkir_hess_1_6< double* >() ); return;}
            if( name == "stvenkir_hess_2_0") { simplefunc = psf( stvenkir_hess_2_0< double* >() ); return;}
            if( name == "stvenkir_hess_2_1") { simplefunc = psf( stvenkir_hess_2_1< double* >() ); return;}
            if( name == "stvenkir_hess_2_2") { simplefunc = psf( stvenkir_hess_2_2< double* >() ); return;}
            if( name == "stvenkir_hess_2_3") { simplefunc = psf( stvenkir_hess_2_3< double* >() ); return;}
            if( name == "stvenkir_hess_2_4") { simplefunc = psf( stvenkir_hess_2_4< double* >() ); return;}
            if( name == "stvenkir_hess_2_5") { simplefunc = psf( stvenkir_hess_2_5< double* >() ); return;}
            if( name == "stvenkir_hess_2_6") { simplefunc = psf( stvenkir_hess_2_6< double* >() ); return;}
            if( name == "stvenkir_hess_3_0") { simplefunc = psf( stvenkir_hess_3_0< double* >() ); return;}
            if( name == "stvenkir_hess_3_1") { simplefunc = psf( stvenkir_hess_3_1< double* >() ); return;}
            if( name == "stvenkir_hess_3_2") { simplefunc = psf( stvenkir_hess_3_2< double* >() ); return;}
            if( name == "stvenkir_hess_3_3") { simplefunc = psf( stvenkir_hess_3_3< double* >() ); return;}
            if( name == "stvenkir_hess_3_4") { simplefunc = psf( stvenkir_hess_3_4< double* >() ); return;}
            if( name == "stvenkir_hess_3_5") { simplefunc = psf( stvenkir_hess_3_5< double* >() ); return;}
            if( name == "stvenkir_hess_3_6") { simplefunc = psf( stvenkir_hess_3_6< double* >() ); return;}
            if( name == "stvenkir_hess_4_0") { simplefunc = psf( stvenkir_hess_4_0< double* >() ); return;}
            if( name == "stvenkir_hess_4_1") { simplefunc = psf( stvenkir_hess_4_1< double* >() ); return;}
            if( name == "stvenkir_hess_4_2") { simplefunc = psf( stvenkir_hess_4_2< double* >() ); return;}
            if( name == "stvenkir_hess_4_3") { simplefunc = psf( stvenkir_hess_4_3< double* >() ); return;}
            if( name == "stvenkir_hess_4_4") { simplefunc = psf( stvenkir_hess_4_4< double* >() ); return;}
            if( name == "stvenkir_hess_4_5") { simplefunc = psf( stvenkir_hess_4_5< double* >() ); return;}
            if( name == "stvenkir_hess_4_6") { simplefunc = psf( stvenkir_hess_4_6< double* >() ); return;}
            if( name == "stvenkir_hess_5_0") { simplefunc = psf( stvenkir_hess_5_0< double* >() ); return;}
            if( name == "stvenkir_hess_5_1") { simplefunc = psf( stvenkir_hess_5_1< double* >() ); return;}
            if( name == "stvenkir_hess_5_2") { simplefunc = psf( stvenkir_hess_5_2< double* >() ); return;}
            if( name == "stvenkir_hess_5_3") { simplefunc = psf( stvenkir_hess_5_3< double* >() ); return;}
            if( name == "stvenkir_hess_5_4") { simplefunc = psf( stvenkir_hess_5_4< double* >() ); return;}
            if( name == "stvenkir_hess_5_5") { simplefunc = psf( stvenkir_hess_5_5< double* >() ); return;}
            if( name == "stvenkir_hess_5_6") { simplefunc = psf( stvenkir_hess_5_6< double* >() ); return;}
            if( name == "stvenkir_hess_6_0") { simplefunc = psf( stvenkir_hess_6_0< double* >() ); return;}
            if( name == "stvenkir_hess_6_1") { simplefunc = psf( stvenkir_hess_6_1< double* >() ); return;}
            if( name == "stvenkir_hess_6_2") { simplefunc = psf( stvenkir_hess_6_2< double* >() ); return;}
            if( name == "stvenkir_hess_6_3") { simplefunc = psf( stvenkir_hess_6_3< double* >() ); return;}
            if( name == "stvenkir_hess_6_4") { simplefunc = psf( stvenkir_hess_6_4< double* >() ); return;}
            if( name == "stvenkir_hess_6_5") { simplefunc = psf( stvenkir_hess_6_5< double* >() ); return;}
            if( name == "stvenkir_hess_6_6") { simplefunc = psf( stvenkir_hess_6_6< double* >() ); return;}
            if( name == "quadlog_f") { simplefunc = psf( quadlog_f< double* >() ); return;}
            if( name == "quadlog_grad_0") { simplefunc = psf( quadlog_grad_0< double* >() ); return;}
            if( name == "quadlog_grad_1") { simplefunc = psf( quadlog_grad_1< double* >() ); return;}
            if( name == "quadlog_grad_2") { simplefunc = psf( quadlog_grad_2< double* >() ); return;}
            if( name == "quadlog_grad_3") { simplefunc = psf( quadlog_grad_3< double* >() ); return;}
            if( name == "quadlog_grad_4") { simplefunc = psf( quadlog_grad_4< double* >() ); return;}
            if( name == "quadlog_grad_5") { simplefunc = psf( quadlog_grad_5< double* >() ); return;}
            if( name == "quadlog_grad_6") { simplefunc = psf( quadlog_grad_6< double* >() ); return;}
            if( name == "quadlog_hess_0_0") { simplefunc = psf( quadlog_hess_0_0< double* >() ); return;}
            if( name == "quadlog_hess_0_1") { simplefunc = psf( quadlog_hess_0_1< double* >() ); return;}
            if( name == "quadlog_hess_0_2") { simplefunc = psf( quadlog_hess_0_2< double* >() ); return;}
            if( name == "quadlog_hess_0_3") { simplefunc = psf( quadlog_hess_0_3< double* >() ); return;}
            if( name == "quadlog_hess_0_4") { simplefunc = psf( quadlog_hess_0_4< double* >() ); return;}
            if( name == "quadlog_hess_0_5") { simplefunc = psf( quadlog_hess_0_5< double* >() ); return;}
            if( name == "quadlog_hess_0_6") { simplefunc = psf( quadlog_hess_0_6< double* >() ); return;}
            if( name == "quadlog_hess_1_0") { simplefunc = psf( quadlog_hess_1_0< double* >() ); return;}
            if( name == "quadlog_hess_1_1") { simplefunc = psf( quadlog_hess_1_1< double* >() ); return;}
            if( name == "quadlog_hess_1_2") { simplefunc = psf( quadlog_hess_1_2< double* >() ); return;}
            if( name == "quadlog_hess_1_3") { simplefunc = psf( quadlog_hess_1_3< double* >() ); return;}
            if( name == "quadlog_hess_1_4") { simplefunc = psf( quadlog_hess_1_4< double* >() ); return;}
            if( name == "quadlog_hess_1_5") { simplefunc = psf( quadlog_hess_1_5< double* >() ); return;}
            if( name == "quadlog_hess_1_6") { simplefunc = psf( quadlog_hess_1_6< double* >() ); return;}
            if( name == "quadlog_hess_2_0") { simplefunc = psf( quadlog_hess_2_0< double* >() ); return;}
            if( name == "quadlog_hess_2_1") { simplefunc = psf( quadlog_hess_2_1< double* >() ); return;}
            if( name == "quadlog_hess_2_2") { simplefunc = psf( quadlog_hess_2_2< double* >() ); return;}
            if( name == "quadlog_hess_2_3") { simplefunc = psf( quadlog_hess_2_3< double* >() ); return;}
            if( name == "quadlog_hess_2_4") { simplefunc = psf( quadlog_hess_2_4< double* >() ); return;}
            if( name == "quadlog_hess_2_5") { simplefunc = psf( quadlog_hess_2_5< double* >() ); return;}
            if( name == "quadlog_hess_2_6") { simplefunc = psf( quadlog_hess_2_6< double* >() ); return;}
            if( name == "quadlog_hess_3_0") { simplefunc = psf( quadlog_hess_3_0< double* >() ); return;}
            if( name == "quadlog_hess_3_1") { simplefunc = psf( quadlog_hess_3_1< double* >() ); return;}
            if( name == "quadlog_hess_3_2") { simplefunc = psf( quadlog_hess_3_2< double* >() ); return;}
            if( name == "quadlog_hess_3_3") { simplefunc = psf( quadlog_hess_3_3< double* >() ); return;}
            if( name == "quadlog_hess_3_4") { simplefunc = psf( quadlog_hess_3_4< double* >() ); return;}
            if( name == "quadlog_hess_3_5") { simplefunc = psf( quadlog_hess_3_5< double* >() ); return;}
            if( name == "quadlog_hess_3_6") { simplefunc = psf( quadlog_hess_3_6< double* >() ); return;}
            if( name == "quadlog_hess_4_0") { simplefunc = psf( quadlog_hess_4_0< double* >() ); return;}
            if( name == "quadlog_hess_4_1") { simplefunc = psf( quadlog_hess_4_1< double* >() ); return;}
            if( name == "quadlog_hess_4_2") { simplefunc = psf( quadlog_hess_4_2< double* >() ); return;}
            if( name == "quadlog_hess_4_3") { simplefunc = psf( quadlog_hess_4_3< double* >() ); return;}
            if( name == "quadlog_hess_4_4") { simplefunc = psf( quadlog_hess_4_4< double* >() ); return;}
            if( name == "quadlog_hess_4_5") { simplefunc = psf( quadlog_hess_4_5< double* >() ); return;}
            if( name == "quadlog_hess_4_6") { simplefunc = psf( quadlog_hess_4_6< double* >() ); return;}
            if( name == "quadlog_hess_5_0") { simplefunc = psf( quadlog_hess_5_0< double* >() ); return;}
            if( name == "quadlog_hess_5_1") { simplefunc = psf( quadlog_hess_5_1< double* >() ); return;}
            if( name == "quadlog_hess_5_2") { simplefunc = psf( quadlog_hess_5_2< double* >() ); return;}
            if( name == "quadlog_hess_5_3") { simplefunc = psf( quadlog_hess_5_3< double* >() ); return;}
            if( name == "quadlog_hess_5_4") { simplefunc = psf( quadlog_hess_5_4< double* >() ); return;}
            if( name == "quadlog_hess_5_5") { simplefunc = psf( quadlog_hess_5_5< double* >() ); return;}
            if( name == "quadlog_hess_5_6") { simplefunc = psf( quadlog_hess_5_6< double* >() ); return;}
            if( name == "quadlog_hess_6_0") { simplefunc = psf( quadlog_hess_6_0< double* >() ); return;}
            if( name == "quadlog_hess_6_1") { simplefunc = psf( quadlog_hess_6_1< double* >() ); return;}
            if( name == "quadlog_hess_6_2") { simplefunc = psf( quadlog_hess_6_2< double* >() ); return;}
            if( name == "quadlog_hess_6_3") { simplefunc = psf( quadlog_hess_6_3< double* >() ); return;}
            if( name == "quadlog_hess_6_4") { simplefunc = psf( quadlog_hess_6_4< double* >() ); return;}
            if( name == "quadlog_hess_6_5") { simplefunc = psf( quadlog_hess_6_5< double* >() ); return;}
            if( name == "quadlog_hess_6_6") { simplefunc = psf( quadlog_hess_6_6< double* >() ); return;}
            if( name == "hardening_f") { simplefunc = psf( hardening_f< double* >() ); return;}
            if( name == "hardening_grad_0") { simplefunc = psf( hardening_grad_0< double* >() ); return;}
            if( name == "hardening_hess_0_0") { simplefunc = psf( hardening_hess_0_0< double* >() ); return;}
            if( name == "von_mises_f") { simplefunc = psf( von_mises_f< double* >() ); return;}
            if( name == "von_mises_grad_0") { simplefunc = psf( von_mises_grad_0< double* >() ); return;}
            if( name == "von_mises_grad_1") { simplefunc = psf( von_mises_grad_1< double* >() ); return;}
            if( name == "von_mises_grad_2") { simplefunc = psf( von_mises_grad_2< double* >() ); return;}
            if( name == "von_mises_grad_3") { simplefunc = psf( von_mises_grad_3< double* >() ); return;}
            if( name == "von_mises_grad_4") { simplefunc = psf( von_mises_grad_4< double* >() ); return;}
            if( name == "von_mises_hess_0_0") { simplefunc = psf( von_mises_hess_0_0< double* >() ); return;}
            if( name == "von_mises_hess_0_1") { simplefunc = psf( von_mises_hess_0_1< double* >() ); return;}
            if( name == "von_mises_hess_0_2") { simplefunc = psf( von_mises_hess_0_2< double* >() ); return;}
            if( name == "von_mises_hess_0_3") { simplefunc = psf( von_mises_hess_0_3< double* >() ); return;}
            if( name == "von_mises_hess_0_4") { simplefunc = psf( von_mises_hess_0_4< double* >() ); return;}
            if( name == "von_mises_hess_1_0") { simplefunc = psf( von_mises_hess_1_0< double* >() ); return;}
            if( name == "von_mises_hess_1_1") { simplefunc = psf( von_mises_hess_1_1< double* >() ); return;}
            if( name == "von_mises_hess_1_2") { simplefunc = psf( von_mises_hess_1_2< double* >() ); return;}
            if( name == "von_mises_hess_1_3") { simplefunc = psf( von_mises_hess_1_3< double* >() ); return;}
            if( name == "von_mises_hess_1_4") { simplefunc = psf( von_mises_hess_1_4< double* >() ); return;}
            if( name == "von_mises_hess_2_0") { simplefunc = psf( von_mises_hess_2_0< double* >() ); return;}
            if( name == "von_mises_hess_2_1") { simplefunc = psf( von_mises_hess_2_1< double* >() ); return;}
            if( name == "von_mises_hess_2_2") { simplefunc = psf( von_mises_hess_2_2< double* >() ); return;}
            if( name == "von_mises_hess_2_3") { simplefunc = psf( von_mises_hess_2_3< double* >() ); return;}
            if( name == "von_mises_hess_2_4") { simplefunc = psf( von_mises_hess_2_4< double* >() ); return;}
            if( name == "von_mises_hess_3_0") { simplefunc = psf( von_mises_hess_3_0< double* >() ); return;}
            if( name == "von_mises_hess_3_1") { simplefunc = psf( von_mises_hess_3_1< double* >() ); return;}
            if( name == "von_mises_hess_3_2") { simplefunc = psf( von_mises_hess_3_2< double* >() ); return;}
            if( name == "von_mises_hess_3_3") { simplefunc = psf( von_mises_hess_3_3< double* >() ); return;}
            if( name == "von_mises_hess_3_4") { simplefunc = psf( von_mises_hess_3_4< double* >() ); return;}
            if( name == "von_mises_hess_4_0") { simplefunc = psf( von_mises_hess_4_0< double* >() ); return;}
            if( name == "von_mises_hess_4_1") { simplefunc = psf( von_mises_hess_4_1< double* >() ); return;}
            if( name == "von_mises_hess_4_2") { simplefunc = psf( von_mises_hess_4_2< double* >() ); return;}
            if( name == "von_mises_hess_4_3") { simplefunc = psf( von_mises_hess_4_3< double* >() ); return;}
            if( name == "von_mises_hess_4_4") { simplefunc = psf( von_mises_hess_4_4< double* >() ); return;}
            if( name == "neohook_f") { simplefunc = psf( neohook_f< double* >() ); return;}
            if( name == "neohook_grad_0") { simplefunc = psf( neohook_grad_0< double* >() ); return;}
            if( name == "neohook_grad_1") { simplefunc = psf( neohook_grad_1< double* >() ); return;}
            if( name == "neohook_grad_2") { simplefunc = psf( neohook_grad_2< double* >() ); return;}
            if( name == "neohook_grad_3") { simplefunc = psf( neohook_grad_3< double* >() ); return;}
            if( name == "neohook_grad_4") { simplefunc = psf( neohook_grad_4< double* >() ); return;}
            if( name == "neohook_grad_5") { simplefunc = psf( neohook_grad_5< double* >() ); return;}
            if( name == "neohook_grad_6") { simplefunc = psf( neohook_grad_6< double* >() ); return;}
            if( name == "neohook_hess_0_0") { simplefunc = psf( neohook_hess_0_0< double* >() ); return;}
            if( name == "neohook_hess_0_1") { simplefunc = psf( neohook_hess_0_1< double* >() ); return;}
            if( name == "neohook_hess_0_2") { simplefunc = psf( neohook_hess_0_2< double* >() ); return;}
            if( name == "neohook_hess_0_3") { simplefunc = psf( neohook_hess_0_3< double* >() ); return;}
            if( name == "neohook_hess_0_4") { simplefunc = psf( neohook_hess_0_4< double* >() ); return;}
            if( name == "neohook_hess_0_5") { simplefunc = psf( neohook_hess_0_5< double* >() ); return;}
            if( name == "neohook_hess_0_6") { simplefunc = psf( neohook_hess_0_6< double* >() ); return;}
            if( name == "neohook_hess_1_0") { simplefunc = psf( neohook_hess_1_0< double* >() ); return;}
            if( name == "neohook_hess_1_1") { simplefunc = psf( neohook_hess_1_1< double* >() ); return;}
            if( name == "neohook_hess_1_2") { simplefunc = psf( neohook_hess_1_2< double* >() ); return;}
            if( name == "neohook_hess_1_3") { simplefunc = psf( neohook_hess_1_3< double* >() ); return;}
            if( name == "neohook_hess_1_4") { simplefunc = psf( neohook_hess_1_4< double* >() ); return;}
            if( name == "neohook_hess_1_5") { simplefunc = psf( neohook_hess_1_5< double* >() ); return;}
            if( name == "neohook_hess_1_6") { simplefunc = psf( neohook_hess_1_6< double* >() ); return;}
            if( name == "neohook_hess_2_0") { simplefunc = psf( neohook_hess_2_0< double* >() ); return;}
            if( name == "neohook_hess_2_1") { simplefunc = psf( neohook_hess_2_1< double* >() ); return;}
            if( name == "neohook_hess_2_2") { simplefunc = psf( neohook_hess_2_2< double* >() ); return;}
            if( name == "neohook_hess_2_3") { simplefunc = psf( neohook_hess_2_3< double* >() ); return;}
            if( name == "neohook_hess_2_4") { simplefunc = psf( neohook_hess_2_4< double* >() ); return;}
            if( name == "neohook_hess_2_5") { simplefunc = psf( neohook_hess_2_5< double* >() ); return;}
            if( name == "neohook_hess_2_6") { simplefunc = psf( neohook_hess_2_6< double* >() ); return;}
            if( name == "neohook_hess_3_0") { simplefunc = psf( neohook_hess_3_0< double* >() ); return;}
            if( name == "neohook_hess_3_1") { simplefunc = psf( neohook_hess_3_1< double* >() ); return;}
            if( name == "neohook_hess_3_2") { simplefunc = psf( neohook_hess_3_2< double* >() ); return;}
            if( name == "neohook_hess_3_3") { simplefunc = psf( neohook_hess_3_3< double* >() ); return;}
            if( name == "neohook_hess_3_4") { simplefunc = psf( neohook_hess_3_4< double* >() ); return;}
            if( name == "neohook_hess_3_5") { simplefunc = psf( neohook_hess_3_5< double* >() ); return;}
            if( name == "neohook_hess_3_6") { simplefunc = psf( neohook_hess_3_6< double* >() ); return;}
            if( name == "neohook_hess_4_0") { simplefunc = psf( neohook_hess_4_0< double* >() ); return;}
            if( name == "neohook_hess_4_1") { simplefunc = psf( neohook_hess_4_1< double* >() ); return;}
            if( name == "neohook_hess_4_2") { simplefunc = psf( neohook_hess_4_2< double* >() ); return;}
            if( name == "neohook_hess_4_3") { simplefunc = psf( neohook_hess_4_3< double* >() ); return;}
            if( name == "neohook_hess_4_4") { simplefunc = psf( neohook_hess_4_4< double* >() ); return;}
            if( name == "neohook_hess_4_5") { simplefunc = psf( neohook_hess_4_5< double* >() ); return;}
            if( name == "neohook_hess_4_6") { simplefunc = psf( neohook_hess_4_6< double* >() ); return;}
            if( name == "neohook_hess_5_0") { simplefunc = psf( neohook_hess_5_0< double* >() ); return;}
            if( name == "neohook_hess_5_1") { simplefunc = psf( neohook_hess_5_1< double* >() ); return;}
            if( name == "neohook_hess_5_2") { simplefunc = psf( neohook_hess_5_2< double* >() ); return;}
            if( name == "neohook_hess_5_3") { simplefunc = psf( neohook_hess_5_3< double* >() ); return;}
            if( name == "neohook_hess_5_4") { simplefunc = psf( neohook_hess_5_4< double* >() ); return;}
            if( name == "neohook_hess_5_5") { simplefunc = psf( neohook_hess_5_5< double* >() ); return;}
            if( name == "neohook_hess_5_6") { simplefunc = psf( neohook_hess_5_6< double* >() ); return;}
            if( name == "neohook_hess_6_0") { simplefunc = psf( neohook_hess_6_0< double* >() ); return;}
            if( name == "neohook_hess_6_1") { simplefunc = psf( neohook_hess_6_1< double* >() ); return;}
            if( name == "neohook_hess_6_2") { simplefunc = psf( neohook_hess_6_2< double* >() ); return;}
            if( name == "neohook_hess_6_3") { simplefunc = psf( neohook_hess_6_3< double* >() ); return;}
            if( name == "neohook_hess_6_4") { simplefunc = psf( neohook_hess_6_4< double* >() ); return;}
            if( name == "neohook_hess_6_5") { simplefunc = psf( neohook_hess_6_5< double* >() ); return;}
            if( name == "neohook_hess_6_6") { simplefunc = psf( neohook_hess_6_6< double* >() ); return;}
            else throw std::runtime_error( "PSimpleFunction< double*, double > " + name + " was not found in the PLibrary");
>>>>>>> 1853fc05524906a2a2d829468a5c2e263978073c
        }

        void PLibrary::checkout( std::string name, PFunction< std::vector<double>, double > &func)
        {
<<<<<<< HEAD
            if( name == "neohook") func = PFunction< std::vector<double>, double >( neohook< std::vector<double> >() );
            if( name == "stvenkir") func = PFunction< std::vector<double>, double >( stvenkir< std::vector<double> >() );
            if( name == "linear_hardening") func = PFunction< std::vector<double>, double >( linear_hardening< std::vector<double> >() );
            if( name == "Cu_hardening") func = PFunction< std::vector<double>, double >( Cu_hardening< std::vector<double> >() );
            if( name == "von_mises") func = PFunction< std::vector<double>, double >( von_mises< std::vector<double> >() );
            if( name == "quadlog") func = PFunction< std::vector<double>, double >( quadlog< std::vector<double> >() );
        }
        void PLibrary::checkout( std::string name, PFunction< double*, double > &func)
        {
            if( name == "neohook") func = PFunction< double*, double >( neohook< double* >() );
            if( name == "stvenkir") func = PFunction< double*, double >( stvenkir< double* >() );
            if( name == "linear_hardening") func = PFunction< double*, double >( linear_hardening< double* >() );
            if( name == "Cu_hardening") func = PFunction< double*, double >( Cu_hardening< double* >() );
            if( name == "von_mises") func = PFunction< double*, double >( von_mises< double* >() );
            if( name == "quadlog") func = PFunction< double*, double >( quadlog< double* >() );
=======
            typedef PFunction< std::vector<double>, double > pf;
            if( name == "stvenkir") { func = pf( stvenkir< std::vector<double> >() ); return;}
            if( name == "quadlog") { func = pf( quadlog< std::vector<double> >() ); return;}
            if( name == "hardening") { func = pf( hardening< std::vector<double> >() ); return;}
            if( name == "von_mises") { func = pf( von_mises< std::vector<double> >() ); return;}
            if( name == "neohook") { func = pf( neohook< std::vector<double> >() ); return;}
            throw std::runtime_error( "PFunction< std::vector<double>, double > " + name + " was not found in the PLibrary");
        }
        void PLibrary::checkout( std::string name, PFunction< double*, double > &func)
        {
            typedef PFunction< double*, double > pf;
            if( name == "stvenkir") { func = pf( stvenkir< double* >() ); return;}
            if( name == "quadlog") { func = pf( quadlog< double* >() ); return;}
            if( name == "hardening") { func = pf( hardening< double* >() ); return;}
            if( name == "von_mises") { func = pf( von_mises< double* >() ); return;}
            if( name == "neohook") { func = pf( neohook< double* >() ); return;}
            throw std::runtime_error( "PFunction< double*, double > " + name + " was not found in the PLibrary");
>>>>>>> 1853fc05524906a2a2d829468a5c2e263978073c
        }




        void PLibrary::checkout( std::string name, PSimpleBase< std::vector<double>, double > *&simplefunc)
        {
<<<<<<< HEAD
            if( name == "neohook_f") simplefunc = new neohook_f< std::vector<double> >();
            if( name == "neohook_grad_0") simplefunc = new neohook_grad_0< std::vector<double> >();
            if( name == "neohook_grad_1") simplefunc = new neohook_grad_1< std::vector<double> >();
            if( name == "neohook_grad_2") simplefunc = new neohook_grad_2< std::vector<double> >();
            if( name == "neohook_grad_3") simplefunc = new neohook_grad_3< std::vector<double> >();
            if( name == "neohook_grad_4") simplefunc = new neohook_grad_4< std::vector<double> >();
            if( name == "neohook_grad_5") simplefunc = new neohook_grad_5< std::vector<double> >();
            if( name == "neohook_grad_6") simplefunc = new neohook_grad_6< std::vector<double> >();
            if( name == "neohook_hess_0_0") simplefunc = new neohook_hess_0_0< std::vector<double> >();
            if( name == "neohook_hess_0_1") simplefunc = new neohook_hess_0_1< std::vector<double> >();
            if( name == "neohook_hess_0_2") simplefunc = new neohook_hess_0_2< std::vector<double> >();
            if( name == "neohook_hess_0_3") simplefunc = new neohook_hess_0_3< std::vector<double> >();
            if( name == "neohook_hess_0_4") simplefunc = new neohook_hess_0_4< std::vector<double> >();
            if( name == "neohook_hess_0_5") simplefunc = new neohook_hess_0_5< std::vector<double> >();
            if( name == "neohook_hess_0_6") simplefunc = new neohook_hess_0_6< std::vector<double> >();
            if( name == "neohook_hess_1_0") simplefunc = new neohook_hess_1_0< std::vector<double> >();
            if( name == "neohook_hess_1_1") simplefunc = new neohook_hess_1_1< std::vector<double> >();
            if( name == "neohook_hess_1_2") simplefunc = new neohook_hess_1_2< std::vector<double> >();
            if( name == "neohook_hess_1_3") simplefunc = new neohook_hess_1_3< std::vector<double> >();
            if( name == "neohook_hess_1_4") simplefunc = new neohook_hess_1_4< std::vector<double> >();
            if( name == "neohook_hess_1_5") simplefunc = new neohook_hess_1_5< std::vector<double> >();
            if( name == "neohook_hess_1_6") simplefunc = new neohook_hess_1_6< std::vector<double> >();
            if( name == "neohook_hess_2_0") simplefunc = new neohook_hess_2_0< std::vector<double> >();
            if( name == "neohook_hess_2_1") simplefunc = new neohook_hess_2_1< std::vector<double> >();
            if( name == "neohook_hess_2_2") simplefunc = new neohook_hess_2_2< std::vector<double> >();
            if( name == "neohook_hess_2_3") simplefunc = new neohook_hess_2_3< std::vector<double> >();
            if( name == "neohook_hess_2_4") simplefunc = new neohook_hess_2_4< std::vector<double> >();
            if( name == "neohook_hess_2_5") simplefunc = new neohook_hess_2_5< std::vector<double> >();
            if( name == "neohook_hess_2_6") simplefunc = new neohook_hess_2_6< std::vector<double> >();
            if( name == "neohook_hess_3_0") simplefunc = new neohook_hess_3_0< std::vector<double> >();
            if( name == "neohook_hess_3_1") simplefunc = new neohook_hess_3_1< std::vector<double> >();
            if( name == "neohook_hess_3_2") simplefunc = new neohook_hess_3_2< std::vector<double> >();
            if( name == "neohook_hess_3_3") simplefunc = new neohook_hess_3_3< std::vector<double> >();
            if( name == "neohook_hess_3_4") simplefunc = new neohook_hess_3_4< std::vector<double> >();
            if( name == "neohook_hess_3_5") simplefunc = new neohook_hess_3_5< std::vector<double> >();
            if( name == "neohook_hess_3_6") simplefunc = new neohook_hess_3_6< std::vector<double> >();
            if( name == "neohook_hess_4_0") simplefunc = new neohook_hess_4_0< std::vector<double> >();
            if( name == "neohook_hess_4_1") simplefunc = new neohook_hess_4_1< std::vector<double> >();
            if( name == "neohook_hess_4_2") simplefunc = new neohook_hess_4_2< std::vector<double> >();
            if( name == "neohook_hess_4_3") simplefunc = new neohook_hess_4_3< std::vector<double> >();
            if( name == "neohook_hess_4_4") simplefunc = new neohook_hess_4_4< std::vector<double> >();
            if( name == "neohook_hess_4_5") simplefunc = new neohook_hess_4_5< std::vector<double> >();
            if( name == "neohook_hess_4_6") simplefunc = new neohook_hess_4_6< std::vector<double> >();
            if( name == "neohook_hess_5_0") simplefunc = new neohook_hess_5_0< std::vector<double> >();
            if( name == "neohook_hess_5_1") simplefunc = new neohook_hess_5_1< std::vector<double> >();
            if( name == "neohook_hess_5_2") simplefunc = new neohook_hess_5_2< std::vector<double> >();
            if( name == "neohook_hess_5_3") simplefunc = new neohook_hess_5_3< std::vector<double> >();
            if( name == "neohook_hess_5_4") simplefunc = new neohook_hess_5_4< std::vector<double> >();
            if( name == "neohook_hess_5_5") simplefunc = new neohook_hess_5_5< std::vector<double> >();
            if( name == "neohook_hess_5_6") simplefunc = new neohook_hess_5_6< std::vector<double> >();
            if( name == "neohook_hess_6_0") simplefunc = new neohook_hess_6_0< std::vector<double> >();
            if( name == "neohook_hess_6_1") simplefunc = new neohook_hess_6_1< std::vector<double> >();
            if( name == "neohook_hess_6_2") simplefunc = new neohook_hess_6_2< std::vector<double> >();
            if( name == "neohook_hess_6_3") simplefunc = new neohook_hess_6_3< std::vector<double> >();
            if( name == "neohook_hess_6_4") simplefunc = new neohook_hess_6_4< std::vector<double> >();
            if( name == "neohook_hess_6_5") simplefunc = new neohook_hess_6_5< std::vector<double> >();
            if( name == "neohook_hess_6_6") simplefunc = new neohook_hess_6_6< std::vector<double> >();
            if( name == "stvenkir_f") simplefunc = new stvenkir_f< std::vector<double> >();
            if( name == "stvenkir_grad_0") simplefunc = new stvenkir_grad_0< std::vector<double> >();
            if( name == "stvenkir_grad_1") simplefunc = new stvenkir_grad_1< std::vector<double> >();
            if( name == "stvenkir_grad_2") simplefunc = new stvenkir_grad_2< std::vector<double> >();
            if( name == "stvenkir_grad_3") simplefunc = new stvenkir_grad_3< std::vector<double> >();
            if( name == "stvenkir_grad_4") simplefunc = new stvenkir_grad_4< std::vector<double> >();
            if( name == "stvenkir_grad_5") simplefunc = new stvenkir_grad_5< std::vector<double> >();
            if( name == "stvenkir_grad_6") simplefunc = new stvenkir_grad_6< std::vector<double> >();
            if( name == "stvenkir_hess_0_0") simplefunc = new stvenkir_hess_0_0< std::vector<double> >();
            if( name == "stvenkir_hess_0_1") simplefunc = new stvenkir_hess_0_1< std::vector<double> >();
            if( name == "stvenkir_hess_0_2") simplefunc = new stvenkir_hess_0_2< std::vector<double> >();
            if( name == "stvenkir_hess_0_3") simplefunc = new stvenkir_hess_0_3< std::vector<double> >();
            if( name == "stvenkir_hess_0_4") simplefunc = new stvenkir_hess_0_4< std::vector<double> >();
            if( name == "stvenkir_hess_0_5") simplefunc = new stvenkir_hess_0_5< std::vector<double> >();
            if( name == "stvenkir_hess_0_6") simplefunc = new stvenkir_hess_0_6< std::vector<double> >();
            if( name == "stvenkir_hess_1_0") simplefunc = new stvenkir_hess_1_0< std::vector<double> >();
            if( name == "stvenkir_hess_1_1") simplefunc = new stvenkir_hess_1_1< std::vector<double> >();
            if( name == "stvenkir_hess_1_2") simplefunc = new stvenkir_hess_1_2< std::vector<double> >();
            if( name == "stvenkir_hess_1_3") simplefunc = new stvenkir_hess_1_3< std::vector<double> >();
            if( name == "stvenkir_hess_1_4") simplefunc = new stvenkir_hess_1_4< std::vector<double> >();
            if( name == "stvenkir_hess_1_5") simplefunc = new stvenkir_hess_1_5< std::vector<double> >();
            if( name == "stvenkir_hess_1_6") simplefunc = new stvenkir_hess_1_6< std::vector<double> >();
            if( name == "stvenkir_hess_2_0") simplefunc = new stvenkir_hess_2_0< std::vector<double> >();
            if( name == "stvenkir_hess_2_1") simplefunc = new stvenkir_hess_2_1< std::vector<double> >();
            if( name == "stvenkir_hess_2_2") simplefunc = new stvenkir_hess_2_2< std::vector<double> >();
            if( name == "stvenkir_hess_2_3") simplefunc = new stvenkir_hess_2_3< std::vector<double> >();
            if( name == "stvenkir_hess_2_4") simplefunc = new stvenkir_hess_2_4< std::vector<double> >();
            if( name == "stvenkir_hess_2_5") simplefunc = new stvenkir_hess_2_5< std::vector<double> >();
            if( name == "stvenkir_hess_2_6") simplefunc = new stvenkir_hess_2_6< std::vector<double> >();
            if( name == "stvenkir_hess_3_0") simplefunc = new stvenkir_hess_3_0< std::vector<double> >();
            if( name == "stvenkir_hess_3_1") simplefunc = new stvenkir_hess_3_1< std::vector<double> >();
            if( name == "stvenkir_hess_3_2") simplefunc = new stvenkir_hess_3_2< std::vector<double> >();
            if( name == "stvenkir_hess_3_3") simplefunc = new stvenkir_hess_3_3< std::vector<double> >();
            if( name == "stvenkir_hess_3_4") simplefunc = new stvenkir_hess_3_4< std::vector<double> >();
            if( name == "stvenkir_hess_3_5") simplefunc = new stvenkir_hess_3_5< std::vector<double> >();
            if( name == "stvenkir_hess_3_6") simplefunc = new stvenkir_hess_3_6< std::vector<double> >();
            if( name == "stvenkir_hess_4_0") simplefunc = new stvenkir_hess_4_0< std::vector<double> >();
            if( name == "stvenkir_hess_4_1") simplefunc = new stvenkir_hess_4_1< std::vector<double> >();
            if( name == "stvenkir_hess_4_2") simplefunc = new stvenkir_hess_4_2< std::vector<double> >();
            if( name == "stvenkir_hess_4_3") simplefunc = new stvenkir_hess_4_3< std::vector<double> >();
            if( name == "stvenkir_hess_4_4") simplefunc = new stvenkir_hess_4_4< std::vector<double> >();
            if( name == "stvenkir_hess_4_5") simplefunc = new stvenkir_hess_4_5< std::vector<double> >();
            if( name == "stvenkir_hess_4_6") simplefunc = new stvenkir_hess_4_6< std::vector<double> >();
            if( name == "stvenkir_hess_5_0") simplefunc = new stvenkir_hess_5_0< std::vector<double> >();
            if( name == "stvenkir_hess_5_1") simplefunc = new stvenkir_hess_5_1< std::vector<double> >();
            if( name == "stvenkir_hess_5_2") simplefunc = new stvenkir_hess_5_2< std::vector<double> >();
            if( name == "stvenkir_hess_5_3") simplefunc = new stvenkir_hess_5_3< std::vector<double> >();
            if( name == "stvenkir_hess_5_4") simplefunc = new stvenkir_hess_5_4< std::vector<double> >();
            if( name == "stvenkir_hess_5_5") simplefunc = new stvenkir_hess_5_5< std::vector<double> >();
            if( name == "stvenkir_hess_5_6") simplefunc = new stvenkir_hess_5_6< std::vector<double> >();
            if( name == "stvenkir_hess_6_0") simplefunc = new stvenkir_hess_6_0< std::vector<double> >();
            if( name == "stvenkir_hess_6_1") simplefunc = new stvenkir_hess_6_1< std::vector<double> >();
            if( name == "stvenkir_hess_6_2") simplefunc = new stvenkir_hess_6_2< std::vector<double> >();
            if( name == "stvenkir_hess_6_3") simplefunc = new stvenkir_hess_6_3< std::vector<double> >();
            if( name == "stvenkir_hess_6_4") simplefunc = new stvenkir_hess_6_4< std::vector<double> >();
            if( name == "stvenkir_hess_6_5") simplefunc = new stvenkir_hess_6_5< std::vector<double> >();
            if( name == "stvenkir_hess_6_6") simplefunc = new stvenkir_hess_6_6< std::vector<double> >();
            if( name == "linear_hardening_f") simplefunc = new linear_hardening_f< std::vector<double> >();
            if( name == "linear_hardening_grad_0") simplefunc = new linear_hardening_grad_0< std::vector<double> >();
            if( name == "linear_hardening_grad_1") simplefunc = new linear_hardening_grad_1< std::vector<double> >();
            if( name == "linear_hardening_hess_0_0") simplefunc = new linear_hardening_hess_0_0< std::vector<double> >();
            if( name == "linear_hardening_hess_0_1") simplefunc = new linear_hardening_hess_0_1< std::vector<double> >();
            if( name == "linear_hardening_hess_1_0") simplefunc = new linear_hardening_hess_1_0< std::vector<double> >();
            if( name == "linear_hardening_hess_1_1") simplefunc = new linear_hardening_hess_1_1< std::vector<double> >();
            if( name == "Cu_hardening_f") simplefunc = new Cu_hardening_f< std::vector<double> >();
            if( name == "Cu_hardening_grad_0") simplefunc = new Cu_hardening_grad_0< std::vector<double> >();
            if( name == "Cu_hardening_grad_1") simplefunc = new Cu_hardening_grad_1< std::vector<double> >();
            if( name == "Cu_hardening_hess_0_0") simplefunc = new Cu_hardening_hess_0_0< std::vector<double> >();
            if( name == "Cu_hardening_hess_0_1") simplefunc = new Cu_hardening_hess_0_1< std::vector<double> >();
            if( name == "Cu_hardening_hess_1_0") simplefunc = new Cu_hardening_hess_1_0< std::vector<double> >();
            if( name == "Cu_hardening_hess_1_1") simplefunc = new Cu_hardening_hess_1_1< std::vector<double> >();
            if( name == "von_mises_f") simplefunc = new von_mises_f< std::vector<double> >();
            if( name == "von_mises_grad_0") simplefunc = new von_mises_grad_0< std::vector<double> >();
            if( name == "von_mises_grad_1") simplefunc = new von_mises_grad_1< std::vector<double> >();
            if( name == "von_mises_grad_2") simplefunc = new von_mises_grad_2< std::vector<double> >();
            if( name == "von_mises_grad_3") simplefunc = new von_mises_grad_3< std::vector<double> >();
            if( name == "von_mises_grad_4") simplefunc = new von_mises_grad_4< std::vector<double> >();
            if( name == "von_mises_hess_0_0") simplefunc = new von_mises_hess_0_0< std::vector<double> >();
            if( name == "von_mises_hess_0_1") simplefunc = new von_mises_hess_0_1< std::vector<double> >();
            if( name == "von_mises_hess_0_2") simplefunc = new von_mises_hess_0_2< std::vector<double> >();
            if( name == "von_mises_hess_0_3") simplefunc = new von_mises_hess_0_3< std::vector<double> >();
            if( name == "von_mises_hess_0_4") simplefunc = new von_mises_hess_0_4< std::vector<double> >();
            if( name == "von_mises_hess_1_0") simplefunc = new von_mises_hess_1_0< std::vector<double> >();
            if( name == "von_mises_hess_1_1") simplefunc = new von_mises_hess_1_1< std::vector<double> >();
            if( name == "von_mises_hess_1_2") simplefunc = new von_mises_hess_1_2< std::vector<double> >();
            if( name == "von_mises_hess_1_3") simplefunc = new von_mises_hess_1_3< std::vector<double> >();
            if( name == "von_mises_hess_1_4") simplefunc = new von_mises_hess_1_4< std::vector<double> >();
            if( name == "von_mises_hess_2_0") simplefunc = new von_mises_hess_2_0< std::vector<double> >();
            if( name == "von_mises_hess_2_1") simplefunc = new von_mises_hess_2_1< std::vector<double> >();
            if( name == "von_mises_hess_2_2") simplefunc = new von_mises_hess_2_2< std::vector<double> >();
            if( name == "von_mises_hess_2_3") simplefunc = new von_mises_hess_2_3< std::vector<double> >();
            if( name == "von_mises_hess_2_4") simplefunc = new von_mises_hess_2_4< std::vector<double> >();
            if( name == "von_mises_hess_3_0") simplefunc = new von_mises_hess_3_0< std::vector<double> >();
            if( name == "von_mises_hess_3_1") simplefunc = new von_mises_hess_3_1< std::vector<double> >();
            if( name == "von_mises_hess_3_2") simplefunc = new von_mises_hess_3_2< std::vector<double> >();
            if( name == "von_mises_hess_3_3") simplefunc = new von_mises_hess_3_3< std::vector<double> >();
            if( name == "von_mises_hess_3_4") simplefunc = new von_mises_hess_3_4< std::vector<double> >();
            if( name == "von_mises_hess_4_0") simplefunc = new von_mises_hess_4_0< std::vector<double> >();
            if( name == "von_mises_hess_4_1") simplefunc = new von_mises_hess_4_1< std::vector<double> >();
            if( name == "von_mises_hess_4_2") simplefunc = new von_mises_hess_4_2< std::vector<double> >();
            if( name == "von_mises_hess_4_3") simplefunc = new von_mises_hess_4_3< std::vector<double> >();
            if( name == "von_mises_hess_4_4") simplefunc = new von_mises_hess_4_4< std::vector<double> >();
            if( name == "quadlog_f") simplefunc = new quadlog_f< std::vector<double> >();
            if( name == "quadlog_grad_0") simplefunc = new quadlog_grad_0< std::vector<double> >();
            if( name == "quadlog_grad_1") simplefunc = new quadlog_grad_1< std::vector<double> >();
            if( name == "quadlog_grad_2") simplefunc = new quadlog_grad_2< std::vector<double> >();
            if( name == "quadlog_grad_3") simplefunc = new quadlog_grad_3< std::vector<double> >();
            if( name == "quadlog_grad_4") simplefunc = new quadlog_grad_4< std::vector<double> >();
            if( name == "quadlog_grad_5") simplefunc = new quadlog_grad_5< std::vector<double> >();
            if( name == "quadlog_grad_6") simplefunc = new quadlog_grad_6< std::vector<double> >();
            if( name == "quadlog_hess_0_0") simplefunc = new quadlog_hess_0_0< std::vector<double> >();
            if( name == "quadlog_hess_0_1") simplefunc = new quadlog_hess_0_1< std::vector<double> >();
            if( name == "quadlog_hess_0_2") simplefunc = new quadlog_hess_0_2< std::vector<double> >();
            if( name == "quadlog_hess_0_3") simplefunc = new quadlog_hess_0_3< std::vector<double> >();
            if( name == "quadlog_hess_0_4") simplefunc = new quadlog_hess_0_4< std::vector<double> >();
            if( name == "quadlog_hess_0_5") simplefunc = new quadlog_hess_0_5< std::vector<double> >();
            if( name == "quadlog_hess_0_6") simplefunc = new quadlog_hess_0_6< std::vector<double> >();
            if( name == "quadlog_hess_1_0") simplefunc = new quadlog_hess_1_0< std::vector<double> >();
            if( name == "quadlog_hess_1_1") simplefunc = new quadlog_hess_1_1< std::vector<double> >();
            if( name == "quadlog_hess_1_2") simplefunc = new quadlog_hess_1_2< std::vector<double> >();
            if( name == "quadlog_hess_1_3") simplefunc = new quadlog_hess_1_3< std::vector<double> >();
            if( name == "quadlog_hess_1_4") simplefunc = new quadlog_hess_1_4< std::vector<double> >();
            if( name == "quadlog_hess_1_5") simplefunc = new quadlog_hess_1_5< std::vector<double> >();
            if( name == "quadlog_hess_1_6") simplefunc = new quadlog_hess_1_6< std::vector<double> >();
            if( name == "quadlog_hess_2_0") simplefunc = new quadlog_hess_2_0< std::vector<double> >();
            if( name == "quadlog_hess_2_1") simplefunc = new quadlog_hess_2_1< std::vector<double> >();
            if( name == "quadlog_hess_2_2") simplefunc = new quadlog_hess_2_2< std::vector<double> >();
            if( name == "quadlog_hess_2_3") simplefunc = new quadlog_hess_2_3< std::vector<double> >();
            if( name == "quadlog_hess_2_4") simplefunc = new quadlog_hess_2_4< std::vector<double> >();
            if( name == "quadlog_hess_2_5") simplefunc = new quadlog_hess_2_5< std::vector<double> >();
            if( name == "quadlog_hess_2_6") simplefunc = new quadlog_hess_2_6< std::vector<double> >();
            if( name == "quadlog_hess_3_0") simplefunc = new quadlog_hess_3_0< std::vector<double> >();
            if( name == "quadlog_hess_3_1") simplefunc = new quadlog_hess_3_1< std::vector<double> >();
            if( name == "quadlog_hess_3_2") simplefunc = new quadlog_hess_3_2< std::vector<double> >();
            if( name == "quadlog_hess_3_3") simplefunc = new quadlog_hess_3_3< std::vector<double> >();
            if( name == "quadlog_hess_3_4") simplefunc = new quadlog_hess_3_4< std::vector<double> >();
            if( name == "quadlog_hess_3_5") simplefunc = new quadlog_hess_3_5< std::vector<double> >();
            if( name == "quadlog_hess_3_6") simplefunc = new quadlog_hess_3_6< std::vector<double> >();
            if( name == "quadlog_hess_4_0") simplefunc = new quadlog_hess_4_0< std::vector<double> >();
            if( name == "quadlog_hess_4_1") simplefunc = new quadlog_hess_4_1< std::vector<double> >();
            if( name == "quadlog_hess_4_2") simplefunc = new quadlog_hess_4_2< std::vector<double> >();
            if( name == "quadlog_hess_4_3") simplefunc = new quadlog_hess_4_3< std::vector<double> >();
            if( name == "quadlog_hess_4_4") simplefunc = new quadlog_hess_4_4< std::vector<double> >();
            if( name == "quadlog_hess_4_5") simplefunc = new quadlog_hess_4_5< std::vector<double> >();
            if( name == "quadlog_hess_4_6") simplefunc = new quadlog_hess_4_6< std::vector<double> >();
            if( name == "quadlog_hess_5_0") simplefunc = new quadlog_hess_5_0< std::vector<double> >();
            if( name == "quadlog_hess_5_1") simplefunc = new quadlog_hess_5_1< std::vector<double> >();
            if( name == "quadlog_hess_5_2") simplefunc = new quadlog_hess_5_2< std::vector<double> >();
            if( name == "quadlog_hess_5_3") simplefunc = new quadlog_hess_5_3< std::vector<double> >();
            if( name == "quadlog_hess_5_4") simplefunc = new quadlog_hess_5_4< std::vector<double> >();
            if( name == "quadlog_hess_5_5") simplefunc = new quadlog_hess_5_5< std::vector<double> >();
            if( name == "quadlog_hess_5_6") simplefunc = new quadlog_hess_5_6< std::vector<double> >();
            if( name == "quadlog_hess_6_0") simplefunc = new quadlog_hess_6_0< std::vector<double> >();
            if( name == "quadlog_hess_6_1") simplefunc = new quadlog_hess_6_1< std::vector<double> >();
            if( name == "quadlog_hess_6_2") simplefunc = new quadlog_hess_6_2< std::vector<double> >();
            if( name == "quadlog_hess_6_3") simplefunc = new quadlog_hess_6_3< std::vector<double> >();
            if( name == "quadlog_hess_6_4") simplefunc = new quadlog_hess_6_4< std::vector<double> >();
            if( name == "quadlog_hess_6_5") simplefunc = new quadlog_hess_6_5< std::vector<double> >();
            if( name == "quadlog_hess_6_6") simplefunc = new quadlog_hess_6_6< std::vector<double> >();
        }
        void PLibrary::checkout( std::string name, PSimpleBase< double*, double > *&simplefunc)
        {
            if( name == "neohook_f") simplefunc = new neohook_f< double* >();
            if( name == "neohook_grad_0") simplefunc = new neohook_grad_0< double* >();
            if( name == "neohook_grad_1") simplefunc = new neohook_grad_1< double* >();
            if( name == "neohook_grad_2") simplefunc = new neohook_grad_2< double* >();
            if( name == "neohook_grad_3") simplefunc = new neohook_grad_3< double* >();
            if( name == "neohook_grad_4") simplefunc = new neohook_grad_4< double* >();
            if( name == "neohook_grad_5") simplefunc = new neohook_grad_5< double* >();
            if( name == "neohook_grad_6") simplefunc = new neohook_grad_6< double* >();
            if( name == "neohook_hess_0_0") simplefunc = new neohook_hess_0_0< double* >();
            if( name == "neohook_hess_0_1") simplefunc = new neohook_hess_0_1< double* >();
            if( name == "neohook_hess_0_2") simplefunc = new neohook_hess_0_2< double* >();
            if( name == "neohook_hess_0_3") simplefunc = new neohook_hess_0_3< double* >();
            if( name == "neohook_hess_0_4") simplefunc = new neohook_hess_0_4< double* >();
            if( name == "neohook_hess_0_5") simplefunc = new neohook_hess_0_5< double* >();
            if( name == "neohook_hess_0_6") simplefunc = new neohook_hess_0_6< double* >();
            if( name == "neohook_hess_1_0") simplefunc = new neohook_hess_1_0< double* >();
            if( name == "neohook_hess_1_1") simplefunc = new neohook_hess_1_1< double* >();
            if( name == "neohook_hess_1_2") simplefunc = new neohook_hess_1_2< double* >();
            if( name == "neohook_hess_1_3") simplefunc = new neohook_hess_1_3< double* >();
            if( name == "neohook_hess_1_4") simplefunc = new neohook_hess_1_4< double* >();
            if( name == "neohook_hess_1_5") simplefunc = new neohook_hess_1_5< double* >();
            if( name == "neohook_hess_1_6") simplefunc = new neohook_hess_1_6< double* >();
            if( name == "neohook_hess_2_0") simplefunc = new neohook_hess_2_0< double* >();
            if( name == "neohook_hess_2_1") simplefunc = new neohook_hess_2_1< double* >();
            if( name == "neohook_hess_2_2") simplefunc = new neohook_hess_2_2< double* >();
            if( name == "neohook_hess_2_3") simplefunc = new neohook_hess_2_3< double* >();
            if( name == "neohook_hess_2_4") simplefunc = new neohook_hess_2_4< double* >();
            if( name == "neohook_hess_2_5") simplefunc = new neohook_hess_2_5< double* >();
            if( name == "neohook_hess_2_6") simplefunc = new neohook_hess_2_6< double* >();
            if( name == "neohook_hess_3_0") simplefunc = new neohook_hess_3_0< double* >();
            if( name == "neohook_hess_3_1") simplefunc = new neohook_hess_3_1< double* >();
            if( name == "neohook_hess_3_2") simplefunc = new neohook_hess_3_2< double* >();
            if( name == "neohook_hess_3_3") simplefunc = new neohook_hess_3_3< double* >();
            if( name == "neohook_hess_3_4") simplefunc = new neohook_hess_3_4< double* >();
            if( name == "neohook_hess_3_5") simplefunc = new neohook_hess_3_5< double* >();
            if( name == "neohook_hess_3_6") simplefunc = new neohook_hess_3_6< double* >();
            if( name == "neohook_hess_4_0") simplefunc = new neohook_hess_4_0< double* >();
            if( name == "neohook_hess_4_1") simplefunc = new neohook_hess_4_1< double* >();
            if( name == "neohook_hess_4_2") simplefunc = new neohook_hess_4_2< double* >();
            if( name == "neohook_hess_4_3") simplefunc = new neohook_hess_4_3< double* >();
            if( name == "neohook_hess_4_4") simplefunc = new neohook_hess_4_4< double* >();
            if( name == "neohook_hess_4_5") simplefunc = new neohook_hess_4_5< double* >();
            if( name == "neohook_hess_4_6") simplefunc = new neohook_hess_4_6< double* >();
            if( name == "neohook_hess_5_0") simplefunc = new neohook_hess_5_0< double* >();
            if( name == "neohook_hess_5_1") simplefunc = new neohook_hess_5_1< double* >();
            if( name == "neohook_hess_5_2") simplefunc = new neohook_hess_5_2< double* >();
            if( name == "neohook_hess_5_3") simplefunc = new neohook_hess_5_3< double* >();
            if( name == "neohook_hess_5_4") simplefunc = new neohook_hess_5_4< double* >();
            if( name == "neohook_hess_5_5") simplefunc = new neohook_hess_5_5< double* >();
            if( name == "neohook_hess_5_6") simplefunc = new neohook_hess_5_6< double* >();
            if( name == "neohook_hess_6_0") simplefunc = new neohook_hess_6_0< double* >();
            if( name == "neohook_hess_6_1") simplefunc = new neohook_hess_6_1< double* >();
            if( name == "neohook_hess_6_2") simplefunc = new neohook_hess_6_2< double* >();
            if( name == "neohook_hess_6_3") simplefunc = new neohook_hess_6_3< double* >();
            if( name == "neohook_hess_6_4") simplefunc = new neohook_hess_6_4< double* >();
            if( name == "neohook_hess_6_5") simplefunc = new neohook_hess_6_5< double* >();
            if( name == "neohook_hess_6_6") simplefunc = new neohook_hess_6_6< double* >();
            if( name == "stvenkir_f") simplefunc = new stvenkir_f< double* >();
            if( name == "stvenkir_grad_0") simplefunc = new stvenkir_grad_0< double* >();
            if( name == "stvenkir_grad_1") simplefunc = new stvenkir_grad_1< double* >();
            if( name == "stvenkir_grad_2") simplefunc = new stvenkir_grad_2< double* >();
            if( name == "stvenkir_grad_3") simplefunc = new stvenkir_grad_3< double* >();
            if( name == "stvenkir_grad_4") simplefunc = new stvenkir_grad_4< double* >();
            if( name == "stvenkir_grad_5") simplefunc = new stvenkir_grad_5< double* >();
            if( name == "stvenkir_grad_6") simplefunc = new stvenkir_grad_6< double* >();
            if( name == "stvenkir_hess_0_0") simplefunc = new stvenkir_hess_0_0< double* >();
            if( name == "stvenkir_hess_0_1") simplefunc = new stvenkir_hess_0_1< double* >();
            if( name == "stvenkir_hess_0_2") simplefunc = new stvenkir_hess_0_2< double* >();
            if( name == "stvenkir_hess_0_3") simplefunc = new stvenkir_hess_0_3< double* >();
            if( name == "stvenkir_hess_0_4") simplefunc = new stvenkir_hess_0_4< double* >();
            if( name == "stvenkir_hess_0_5") simplefunc = new stvenkir_hess_0_5< double* >();
            if( name == "stvenkir_hess_0_6") simplefunc = new stvenkir_hess_0_6< double* >();
            if( name == "stvenkir_hess_1_0") simplefunc = new stvenkir_hess_1_0< double* >();
            if( name == "stvenkir_hess_1_1") simplefunc = new stvenkir_hess_1_1< double* >();
            if( name == "stvenkir_hess_1_2") simplefunc = new stvenkir_hess_1_2< double* >();
            if( name == "stvenkir_hess_1_3") simplefunc = new stvenkir_hess_1_3< double* >();
            if( name == "stvenkir_hess_1_4") simplefunc = new stvenkir_hess_1_4< double* >();
            if( name == "stvenkir_hess_1_5") simplefunc = new stvenkir_hess_1_5< double* >();
            if( name == "stvenkir_hess_1_6") simplefunc = new stvenkir_hess_1_6< double* >();
            if( name == "stvenkir_hess_2_0") simplefunc = new stvenkir_hess_2_0< double* >();
            if( name == "stvenkir_hess_2_1") simplefunc = new stvenkir_hess_2_1< double* >();
            if( name == "stvenkir_hess_2_2") simplefunc = new stvenkir_hess_2_2< double* >();
            if( name == "stvenkir_hess_2_3") simplefunc = new stvenkir_hess_2_3< double* >();
            if( name == "stvenkir_hess_2_4") simplefunc = new stvenkir_hess_2_4< double* >();
            if( name == "stvenkir_hess_2_5") simplefunc = new stvenkir_hess_2_5< double* >();
            if( name == "stvenkir_hess_2_6") simplefunc = new stvenkir_hess_2_6< double* >();
            if( name == "stvenkir_hess_3_0") simplefunc = new stvenkir_hess_3_0< double* >();
            if( name == "stvenkir_hess_3_1") simplefunc = new stvenkir_hess_3_1< double* >();
            if( name == "stvenkir_hess_3_2") simplefunc = new stvenkir_hess_3_2< double* >();
            if( name == "stvenkir_hess_3_3") simplefunc = new stvenkir_hess_3_3< double* >();
            if( name == "stvenkir_hess_3_4") simplefunc = new stvenkir_hess_3_4< double* >();
            if( name == "stvenkir_hess_3_5") simplefunc = new stvenkir_hess_3_5< double* >();
            if( name == "stvenkir_hess_3_6") simplefunc = new stvenkir_hess_3_6< double* >();
            if( name == "stvenkir_hess_4_0") simplefunc = new stvenkir_hess_4_0< double* >();
            if( name == "stvenkir_hess_4_1") simplefunc = new stvenkir_hess_4_1< double* >();
            if( name == "stvenkir_hess_4_2") simplefunc = new stvenkir_hess_4_2< double* >();
            if( name == "stvenkir_hess_4_3") simplefunc = new stvenkir_hess_4_3< double* >();
            if( name == "stvenkir_hess_4_4") simplefunc = new stvenkir_hess_4_4< double* >();
            if( name == "stvenkir_hess_4_5") simplefunc = new stvenkir_hess_4_5< double* >();
            if( name == "stvenkir_hess_4_6") simplefunc = new stvenkir_hess_4_6< double* >();
            if( name == "stvenkir_hess_5_0") simplefunc = new stvenkir_hess_5_0< double* >();
            if( name == "stvenkir_hess_5_1") simplefunc = new stvenkir_hess_5_1< double* >();
            if( name == "stvenkir_hess_5_2") simplefunc = new stvenkir_hess_5_2< double* >();
            if( name == "stvenkir_hess_5_3") simplefunc = new stvenkir_hess_5_3< double* >();
            if( name == "stvenkir_hess_5_4") simplefunc = new stvenkir_hess_5_4< double* >();
            if( name == "stvenkir_hess_5_5") simplefunc = new stvenkir_hess_5_5< double* >();
            if( name == "stvenkir_hess_5_6") simplefunc = new stvenkir_hess_5_6< double* >();
            if( name == "stvenkir_hess_6_0") simplefunc = new stvenkir_hess_6_0< double* >();
            if( name == "stvenkir_hess_6_1") simplefunc = new stvenkir_hess_6_1< double* >();
            if( name == "stvenkir_hess_6_2") simplefunc = new stvenkir_hess_6_2< double* >();
            if( name == "stvenkir_hess_6_3") simplefunc = new stvenkir_hess_6_3< double* >();
            if( name == "stvenkir_hess_6_4") simplefunc = new stvenkir_hess_6_4< double* >();
            if( name == "stvenkir_hess_6_5") simplefunc = new stvenkir_hess_6_5< double* >();
            if( name == "stvenkir_hess_6_6") simplefunc = new stvenkir_hess_6_6< double* >();
            if( name == "linear_hardening_f") simplefunc = new linear_hardening_f< double* >();
            if( name == "linear_hardening_grad_0") simplefunc = new linear_hardening_grad_0< double* >();
            if( name == "linear_hardening_grad_1") simplefunc = new linear_hardening_grad_1< double* >();
            if( name == "linear_hardening_hess_0_0") simplefunc = new linear_hardening_hess_0_0< double* >();
            if( name == "linear_hardening_hess_0_1") simplefunc = new linear_hardening_hess_0_1< double* >();
            if( name == "linear_hardening_hess_1_0") simplefunc = new linear_hardening_hess_1_0< double* >();
            if( name == "linear_hardening_hess_1_1") simplefunc = new linear_hardening_hess_1_1< double* >();
            if( name == "Cu_hardening_f") simplefunc = new Cu_hardening_f< double* >();
            if( name == "Cu_hardening_grad_0") simplefunc = new Cu_hardening_grad_0< double* >();
            if( name == "Cu_hardening_grad_1") simplefunc = new Cu_hardening_grad_1< double* >();
            if( name == "Cu_hardening_hess_0_0") simplefunc = new Cu_hardening_hess_0_0< double* >();
            if( name == "Cu_hardening_hess_0_1") simplefunc = new Cu_hardening_hess_0_1< double* >();
            if( name == "Cu_hardening_hess_1_0") simplefunc = new Cu_hardening_hess_1_0< double* >();
            if( name == "Cu_hardening_hess_1_1") simplefunc = new Cu_hardening_hess_1_1< double* >();
            if( name == "von_mises_f") simplefunc = new von_mises_f< double* >();
            if( name == "von_mises_grad_0") simplefunc = new von_mises_grad_0< double* >();
            if( name == "von_mises_grad_1") simplefunc = new von_mises_grad_1< double* >();
            if( name == "von_mises_grad_2") simplefunc = new von_mises_grad_2< double* >();
            if( name == "von_mises_grad_3") simplefunc = new von_mises_grad_3< double* >();
            if( name == "von_mises_grad_4") simplefunc = new von_mises_grad_4< double* >();
            if( name == "von_mises_hess_0_0") simplefunc = new von_mises_hess_0_0< double* >();
            if( name == "von_mises_hess_0_1") simplefunc = new von_mises_hess_0_1< double* >();
            if( name == "von_mises_hess_0_2") simplefunc = new von_mises_hess_0_2< double* >();
            if( name == "von_mises_hess_0_3") simplefunc = new von_mises_hess_0_3< double* >();
            if( name == "von_mises_hess_0_4") simplefunc = new von_mises_hess_0_4< double* >();
            if( name == "von_mises_hess_1_0") simplefunc = new von_mises_hess_1_0< double* >();
            if( name == "von_mises_hess_1_1") simplefunc = new von_mises_hess_1_1< double* >();
            if( name == "von_mises_hess_1_2") simplefunc = new von_mises_hess_1_2< double* >();
            if( name == "von_mises_hess_1_3") simplefunc = new von_mises_hess_1_3< double* >();
            if( name == "von_mises_hess_1_4") simplefunc = new von_mises_hess_1_4< double* >();
            if( name == "von_mises_hess_2_0") simplefunc = new von_mises_hess_2_0< double* >();
            if( name == "von_mises_hess_2_1") simplefunc = new von_mises_hess_2_1< double* >();
            if( name == "von_mises_hess_2_2") simplefunc = new von_mises_hess_2_2< double* >();
            if( name == "von_mises_hess_2_3") simplefunc = new von_mises_hess_2_3< double* >();
            if( name == "von_mises_hess_2_4") simplefunc = new von_mises_hess_2_4< double* >();
            if( name == "von_mises_hess_3_0") simplefunc = new von_mises_hess_3_0< double* >();
            if( name == "von_mises_hess_3_1") simplefunc = new von_mises_hess_3_1< double* >();
            if( name == "von_mises_hess_3_2") simplefunc = new von_mises_hess_3_2< double* >();
            if( name == "von_mises_hess_3_3") simplefunc = new von_mises_hess_3_3< double* >();
            if( name == "von_mises_hess_3_4") simplefunc = new von_mises_hess_3_4< double* >();
            if( name == "von_mises_hess_4_0") simplefunc = new von_mises_hess_4_0< double* >();
            if( name == "von_mises_hess_4_1") simplefunc = new von_mises_hess_4_1< double* >();
            if( name == "von_mises_hess_4_2") simplefunc = new von_mises_hess_4_2< double* >();
            if( name == "von_mises_hess_4_3") simplefunc = new von_mises_hess_4_3< double* >();
            if( name == "von_mises_hess_4_4") simplefunc = new von_mises_hess_4_4< double* >();
            if( name == "quadlog_f") simplefunc = new quadlog_f< double* >();
            if( name == "quadlog_grad_0") simplefunc = new quadlog_grad_0< double* >();
            if( name == "quadlog_grad_1") simplefunc = new quadlog_grad_1< double* >();
            if( name == "quadlog_grad_2") simplefunc = new quadlog_grad_2< double* >();
            if( name == "quadlog_grad_3") simplefunc = new quadlog_grad_3< double* >();
            if( name == "quadlog_grad_4") simplefunc = new quadlog_grad_4< double* >();
            if( name == "quadlog_grad_5") simplefunc = new quadlog_grad_5< double* >();
            if( name == "quadlog_grad_6") simplefunc = new quadlog_grad_6< double* >();
            if( name == "quadlog_hess_0_0") simplefunc = new quadlog_hess_0_0< double* >();
            if( name == "quadlog_hess_0_1") simplefunc = new quadlog_hess_0_1< double* >();
            if( name == "quadlog_hess_0_2") simplefunc = new quadlog_hess_0_2< double* >();
            if( name == "quadlog_hess_0_3") simplefunc = new quadlog_hess_0_3< double* >();
            if( name == "quadlog_hess_0_4") simplefunc = new quadlog_hess_0_4< double* >();
            if( name == "quadlog_hess_0_5") simplefunc = new quadlog_hess_0_5< double* >();
            if( name == "quadlog_hess_0_6") simplefunc = new quadlog_hess_0_6< double* >();
            if( name == "quadlog_hess_1_0") simplefunc = new quadlog_hess_1_0< double* >();
            if( name == "quadlog_hess_1_1") simplefunc = new quadlog_hess_1_1< double* >();
            if( name == "quadlog_hess_1_2") simplefunc = new quadlog_hess_1_2< double* >();
            if( name == "quadlog_hess_1_3") simplefunc = new quadlog_hess_1_3< double* >();
            if( name == "quadlog_hess_1_4") simplefunc = new quadlog_hess_1_4< double* >();
            if( name == "quadlog_hess_1_5") simplefunc = new quadlog_hess_1_5< double* >();
            if( name == "quadlog_hess_1_6") simplefunc = new quadlog_hess_1_6< double* >();
            if( name == "quadlog_hess_2_0") simplefunc = new quadlog_hess_2_0< double* >();
            if( name == "quadlog_hess_2_1") simplefunc = new quadlog_hess_2_1< double* >();
            if( name == "quadlog_hess_2_2") simplefunc = new quadlog_hess_2_2< double* >();
            if( name == "quadlog_hess_2_3") simplefunc = new quadlog_hess_2_3< double* >();
            if( name == "quadlog_hess_2_4") simplefunc = new quadlog_hess_2_4< double* >();
            if( name == "quadlog_hess_2_5") simplefunc = new quadlog_hess_2_5< double* >();
            if( name == "quadlog_hess_2_6") simplefunc = new quadlog_hess_2_6< double* >();
            if( name == "quadlog_hess_3_0") simplefunc = new quadlog_hess_3_0< double* >();
            if( name == "quadlog_hess_3_1") simplefunc = new quadlog_hess_3_1< double* >();
            if( name == "quadlog_hess_3_2") simplefunc = new quadlog_hess_3_2< double* >();
            if( name == "quadlog_hess_3_3") simplefunc = new quadlog_hess_3_3< double* >();
            if( name == "quadlog_hess_3_4") simplefunc = new quadlog_hess_3_4< double* >();
            if( name == "quadlog_hess_3_5") simplefunc = new quadlog_hess_3_5< double* >();
            if( name == "quadlog_hess_3_6") simplefunc = new quadlog_hess_3_6< double* >();
            if( name == "quadlog_hess_4_0") simplefunc = new quadlog_hess_4_0< double* >();
            if( name == "quadlog_hess_4_1") simplefunc = new quadlog_hess_4_1< double* >();
            if( name == "quadlog_hess_4_2") simplefunc = new quadlog_hess_4_2< double* >();
            if( name == "quadlog_hess_4_3") simplefunc = new quadlog_hess_4_3< double* >();
            if( name == "quadlog_hess_4_4") simplefunc = new quadlog_hess_4_4< double* >();
            if( name == "quadlog_hess_4_5") simplefunc = new quadlog_hess_4_5< double* >();
            if( name == "quadlog_hess_4_6") simplefunc = new quadlog_hess_4_6< double* >();
            if( name == "quadlog_hess_5_0") simplefunc = new quadlog_hess_5_0< double* >();
            if( name == "quadlog_hess_5_1") simplefunc = new quadlog_hess_5_1< double* >();
            if( name == "quadlog_hess_5_2") simplefunc = new quadlog_hess_5_2< double* >();
            if( name == "quadlog_hess_5_3") simplefunc = new quadlog_hess_5_3< double* >();
            if( name == "quadlog_hess_5_4") simplefunc = new quadlog_hess_5_4< double* >();
            if( name == "quadlog_hess_5_5") simplefunc = new quadlog_hess_5_5< double* >();
            if( name == "quadlog_hess_5_6") simplefunc = new quadlog_hess_5_6< double* >();
            if( name == "quadlog_hess_6_0") simplefunc = new quadlog_hess_6_0< double* >();
            if( name == "quadlog_hess_6_1") simplefunc = new quadlog_hess_6_1< double* >();
            if( name == "quadlog_hess_6_2") simplefunc = new quadlog_hess_6_2< double* >();
            if( name == "quadlog_hess_6_3") simplefunc = new quadlog_hess_6_3< double* >();
            if( name == "quadlog_hess_6_4") simplefunc = new quadlog_hess_6_4< double* >();
            if( name == "quadlog_hess_6_5") simplefunc = new quadlog_hess_6_5< double* >();
            if( name == "quadlog_hess_6_6") simplefunc = new quadlog_hess_6_6< double* >();
=======
            if( name == "stvenkir_f") { simplefunc = new stvenkir_f< std::vector<double> >(); return;}
            if( name == "stvenkir_grad_0") { simplefunc = new stvenkir_grad_0< std::vector<double> >(); return;}
            if( name == "stvenkir_grad_1") { simplefunc = new stvenkir_grad_1< std::vector<double> >(); return;}
            if( name == "stvenkir_grad_2") { simplefunc = new stvenkir_grad_2< std::vector<double> >(); return;}
            if( name == "stvenkir_grad_3") { simplefunc = new stvenkir_grad_3< std::vector<double> >(); return;}
            if( name == "stvenkir_grad_4") { simplefunc = new stvenkir_grad_4< std::vector<double> >(); return;}
            if( name == "stvenkir_grad_5") { simplefunc = new stvenkir_grad_5< std::vector<double> >(); return;}
            if( name == "stvenkir_grad_6") { simplefunc = new stvenkir_grad_6< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_0_0") { simplefunc = new stvenkir_hess_0_0< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_0_1") { simplefunc = new stvenkir_hess_0_1< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_0_2") { simplefunc = new stvenkir_hess_0_2< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_0_3") { simplefunc = new stvenkir_hess_0_3< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_0_4") { simplefunc = new stvenkir_hess_0_4< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_0_5") { simplefunc = new stvenkir_hess_0_5< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_0_6") { simplefunc = new stvenkir_hess_0_6< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_1_0") { simplefunc = new stvenkir_hess_1_0< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_1_1") { simplefunc = new stvenkir_hess_1_1< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_1_2") { simplefunc = new stvenkir_hess_1_2< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_1_3") { simplefunc = new stvenkir_hess_1_3< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_1_4") { simplefunc = new stvenkir_hess_1_4< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_1_5") { simplefunc = new stvenkir_hess_1_5< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_1_6") { simplefunc = new stvenkir_hess_1_6< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_2_0") { simplefunc = new stvenkir_hess_2_0< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_2_1") { simplefunc = new stvenkir_hess_2_1< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_2_2") { simplefunc = new stvenkir_hess_2_2< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_2_3") { simplefunc = new stvenkir_hess_2_3< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_2_4") { simplefunc = new stvenkir_hess_2_4< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_2_5") { simplefunc = new stvenkir_hess_2_5< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_2_6") { simplefunc = new stvenkir_hess_2_6< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_3_0") { simplefunc = new stvenkir_hess_3_0< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_3_1") { simplefunc = new stvenkir_hess_3_1< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_3_2") { simplefunc = new stvenkir_hess_3_2< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_3_3") { simplefunc = new stvenkir_hess_3_3< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_3_4") { simplefunc = new stvenkir_hess_3_4< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_3_5") { simplefunc = new stvenkir_hess_3_5< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_3_6") { simplefunc = new stvenkir_hess_3_6< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_4_0") { simplefunc = new stvenkir_hess_4_0< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_4_1") { simplefunc = new stvenkir_hess_4_1< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_4_2") { simplefunc = new stvenkir_hess_4_2< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_4_3") { simplefunc = new stvenkir_hess_4_3< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_4_4") { simplefunc = new stvenkir_hess_4_4< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_4_5") { simplefunc = new stvenkir_hess_4_5< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_4_6") { simplefunc = new stvenkir_hess_4_6< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_5_0") { simplefunc = new stvenkir_hess_5_0< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_5_1") { simplefunc = new stvenkir_hess_5_1< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_5_2") { simplefunc = new stvenkir_hess_5_2< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_5_3") { simplefunc = new stvenkir_hess_5_3< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_5_4") { simplefunc = new stvenkir_hess_5_4< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_5_5") { simplefunc = new stvenkir_hess_5_5< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_5_6") { simplefunc = new stvenkir_hess_5_6< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_6_0") { simplefunc = new stvenkir_hess_6_0< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_6_1") { simplefunc = new stvenkir_hess_6_1< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_6_2") { simplefunc = new stvenkir_hess_6_2< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_6_3") { simplefunc = new stvenkir_hess_6_3< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_6_4") { simplefunc = new stvenkir_hess_6_4< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_6_5") { simplefunc = new stvenkir_hess_6_5< std::vector<double> >(); return;}
            if( name == "stvenkir_hess_6_6") { simplefunc = new stvenkir_hess_6_6< std::vector<double> >(); return;}
            if( name == "quadlog_f") { simplefunc = new quadlog_f< std::vector<double> >(); return;}
            if( name == "quadlog_grad_0") { simplefunc = new quadlog_grad_0< std::vector<double> >(); return;}
            if( name == "quadlog_grad_1") { simplefunc = new quadlog_grad_1< std::vector<double> >(); return;}
            if( name == "quadlog_grad_2") { simplefunc = new quadlog_grad_2< std::vector<double> >(); return;}
            if( name == "quadlog_grad_3") { simplefunc = new quadlog_grad_3< std::vector<double> >(); return;}
            if( name == "quadlog_grad_4") { simplefunc = new quadlog_grad_4< std::vector<double> >(); return;}
            if( name == "quadlog_grad_5") { simplefunc = new quadlog_grad_5< std::vector<double> >(); return;}
            if( name == "quadlog_grad_6") { simplefunc = new quadlog_grad_6< std::vector<double> >(); return;}
            if( name == "quadlog_hess_0_0") { simplefunc = new quadlog_hess_0_0< std::vector<double> >(); return;}
            if( name == "quadlog_hess_0_1") { simplefunc = new quadlog_hess_0_1< std::vector<double> >(); return;}
            if( name == "quadlog_hess_0_2") { simplefunc = new quadlog_hess_0_2< std::vector<double> >(); return;}
            if( name == "quadlog_hess_0_3") { simplefunc = new quadlog_hess_0_3< std::vector<double> >(); return;}
            if( name == "quadlog_hess_0_4") { simplefunc = new quadlog_hess_0_4< std::vector<double> >(); return;}
            if( name == "quadlog_hess_0_5") { simplefunc = new quadlog_hess_0_5< std::vector<double> >(); return;}
            if( name == "quadlog_hess_0_6") { simplefunc = new quadlog_hess_0_6< std::vector<double> >(); return;}
            if( name == "quadlog_hess_1_0") { simplefunc = new quadlog_hess_1_0< std::vector<double> >(); return;}
            if( name == "quadlog_hess_1_1") { simplefunc = new quadlog_hess_1_1< std::vector<double> >(); return;}
            if( name == "quadlog_hess_1_2") { simplefunc = new quadlog_hess_1_2< std::vector<double> >(); return;}
            if( name == "quadlog_hess_1_3") { simplefunc = new quadlog_hess_1_3< std::vector<double> >(); return;}
            if( name == "quadlog_hess_1_4") { simplefunc = new quadlog_hess_1_4< std::vector<double> >(); return;}
            if( name == "quadlog_hess_1_5") { simplefunc = new quadlog_hess_1_5< std::vector<double> >(); return;}
            if( name == "quadlog_hess_1_6") { simplefunc = new quadlog_hess_1_6< std::vector<double> >(); return;}
            if( name == "quadlog_hess_2_0") { simplefunc = new quadlog_hess_2_0< std::vector<double> >(); return;}
            if( name == "quadlog_hess_2_1") { simplefunc = new quadlog_hess_2_1< std::vector<double> >(); return;}
            if( name == "quadlog_hess_2_2") { simplefunc = new quadlog_hess_2_2< std::vector<double> >(); return;}
            if( name == "quadlog_hess_2_3") { simplefunc = new quadlog_hess_2_3< std::vector<double> >(); return;}
            if( name == "quadlog_hess_2_4") { simplefunc = new quadlog_hess_2_4< std::vector<double> >(); return;}
            if( name == "quadlog_hess_2_5") { simplefunc = new quadlog_hess_2_5< std::vector<double> >(); return;}
            if( name == "quadlog_hess_2_6") { simplefunc = new quadlog_hess_2_6< std::vector<double> >(); return;}
            if( name == "quadlog_hess_3_0") { simplefunc = new quadlog_hess_3_0< std::vector<double> >(); return;}
            if( name == "quadlog_hess_3_1") { simplefunc = new quadlog_hess_3_1< std::vector<double> >(); return;}
            if( name == "quadlog_hess_3_2") { simplefunc = new quadlog_hess_3_2< std::vector<double> >(); return;}
            if( name == "quadlog_hess_3_3") { simplefunc = new quadlog_hess_3_3< std::vector<double> >(); return;}
            if( name == "quadlog_hess_3_4") { simplefunc = new quadlog_hess_3_4< std::vector<double> >(); return;}
            if( name == "quadlog_hess_3_5") { simplefunc = new quadlog_hess_3_5< std::vector<double> >(); return;}
            if( name == "quadlog_hess_3_6") { simplefunc = new quadlog_hess_3_6< std::vector<double> >(); return;}
            if( name == "quadlog_hess_4_0") { simplefunc = new quadlog_hess_4_0< std::vector<double> >(); return;}
            if( name == "quadlog_hess_4_1") { simplefunc = new quadlog_hess_4_1< std::vector<double> >(); return;}
            if( name == "quadlog_hess_4_2") { simplefunc = new quadlog_hess_4_2< std::vector<double> >(); return;}
            if( name == "quadlog_hess_4_3") { simplefunc = new quadlog_hess_4_3< std::vector<double> >(); return;}
            if( name == "quadlog_hess_4_4") { simplefunc = new quadlog_hess_4_4< std::vector<double> >(); return;}
            if( name == "quadlog_hess_4_5") { simplefunc = new quadlog_hess_4_5< std::vector<double> >(); return;}
            if( name == "quadlog_hess_4_6") { simplefunc = new quadlog_hess_4_6< std::vector<double> >(); return;}
            if( name == "quadlog_hess_5_0") { simplefunc = new quadlog_hess_5_0< std::vector<double> >(); return;}
            if( name == "quadlog_hess_5_1") { simplefunc = new quadlog_hess_5_1< std::vector<double> >(); return;}
            if( name == "quadlog_hess_5_2") { simplefunc = new quadlog_hess_5_2< std::vector<double> >(); return;}
            if( name == "quadlog_hess_5_3") { simplefunc = new quadlog_hess_5_3< std::vector<double> >(); return;}
            if( name == "quadlog_hess_5_4") { simplefunc = new quadlog_hess_5_4< std::vector<double> >(); return;}
            if( name == "quadlog_hess_5_5") { simplefunc = new quadlog_hess_5_5< std::vector<double> >(); return;}
            if( name == "quadlog_hess_5_6") { simplefunc = new quadlog_hess_5_6< std::vector<double> >(); return;}
            if( name == "quadlog_hess_6_0") { simplefunc = new quadlog_hess_6_0< std::vector<double> >(); return;}
            if( name == "quadlog_hess_6_1") { simplefunc = new quadlog_hess_6_1< std::vector<double> >(); return;}
            if( name == "quadlog_hess_6_2") { simplefunc = new quadlog_hess_6_2< std::vector<double> >(); return;}
            if( name == "quadlog_hess_6_3") { simplefunc = new quadlog_hess_6_3< std::vector<double> >(); return;}
            if( name == "quadlog_hess_6_4") { simplefunc = new quadlog_hess_6_4< std::vector<double> >(); return;}
            if( name == "quadlog_hess_6_5") { simplefunc = new quadlog_hess_6_5< std::vector<double> >(); return;}
            if( name == "quadlog_hess_6_6") { simplefunc = new quadlog_hess_6_6< std::vector<double> >(); return;}
            if( name == "hardening_f") { simplefunc = new hardening_f< std::vector<double> >(); return;}
            if( name == "hardening_grad_0") { simplefunc = new hardening_grad_0< std::vector<double> >(); return;}
            if( name == "hardening_hess_0_0") { simplefunc = new hardening_hess_0_0< std::vector<double> >(); return;}
            if( name == "von_mises_f") { simplefunc = new von_mises_f< std::vector<double> >(); return;}
            if( name == "von_mises_grad_0") { simplefunc = new von_mises_grad_0< std::vector<double> >(); return;}
            if( name == "von_mises_grad_1") { simplefunc = new von_mises_grad_1< std::vector<double> >(); return;}
            if( name == "von_mises_grad_2") { simplefunc = new von_mises_grad_2< std::vector<double> >(); return;}
            if( name == "von_mises_grad_3") { simplefunc = new von_mises_grad_3< std::vector<double> >(); return;}
            if( name == "von_mises_grad_4") { simplefunc = new von_mises_grad_4< std::vector<double> >(); return;}
            if( name == "von_mises_hess_0_0") { simplefunc = new von_mises_hess_0_0< std::vector<double> >(); return;}
            if( name == "von_mises_hess_0_1") { simplefunc = new von_mises_hess_0_1< std::vector<double> >(); return;}
            if( name == "von_mises_hess_0_2") { simplefunc = new von_mises_hess_0_2< std::vector<double> >(); return;}
            if( name == "von_mises_hess_0_3") { simplefunc = new von_mises_hess_0_3< std::vector<double> >(); return;}
            if( name == "von_mises_hess_0_4") { simplefunc = new von_mises_hess_0_4< std::vector<double> >(); return;}
            if( name == "von_mises_hess_1_0") { simplefunc = new von_mises_hess_1_0< std::vector<double> >(); return;}
            if( name == "von_mises_hess_1_1") { simplefunc = new von_mises_hess_1_1< std::vector<double> >(); return;}
            if( name == "von_mises_hess_1_2") { simplefunc = new von_mises_hess_1_2< std::vector<double> >(); return;}
            if( name == "von_mises_hess_1_3") { simplefunc = new von_mises_hess_1_3< std::vector<double> >(); return;}
            if( name == "von_mises_hess_1_4") { simplefunc = new von_mises_hess_1_4< std::vector<double> >(); return;}
            if( name == "von_mises_hess_2_0") { simplefunc = new von_mises_hess_2_0< std::vector<double> >(); return;}
            if( name == "von_mises_hess_2_1") { simplefunc = new von_mises_hess_2_1< std::vector<double> >(); return;}
            if( name == "von_mises_hess_2_2") { simplefunc = new von_mises_hess_2_2< std::vector<double> >(); return;}
            if( name == "von_mises_hess_2_3") { simplefunc = new von_mises_hess_2_3< std::vector<double> >(); return;}
            if( name == "von_mises_hess_2_4") { simplefunc = new von_mises_hess_2_4< std::vector<double> >(); return;}
            if( name == "von_mises_hess_3_0") { simplefunc = new von_mises_hess_3_0< std::vector<double> >(); return;}
            if( name == "von_mises_hess_3_1") { simplefunc = new von_mises_hess_3_1< std::vector<double> >(); return;}
            if( name == "von_mises_hess_3_2") { simplefunc = new von_mises_hess_3_2< std::vector<double> >(); return;}
            if( name == "von_mises_hess_3_3") { simplefunc = new von_mises_hess_3_3< std::vector<double> >(); return;}
            if( name == "von_mises_hess_3_4") { simplefunc = new von_mises_hess_3_4< std::vector<double> >(); return;}
            if( name == "von_mises_hess_4_0") { simplefunc = new von_mises_hess_4_0< std::vector<double> >(); return;}
            if( name == "von_mises_hess_4_1") { simplefunc = new von_mises_hess_4_1< std::vector<double> >(); return;}
            if( name == "von_mises_hess_4_2") { simplefunc = new von_mises_hess_4_2< std::vector<double> >(); return;}
            if( name == "von_mises_hess_4_3") { simplefunc = new von_mises_hess_4_3< std::vector<double> >(); return;}
            if( name == "von_mises_hess_4_4") { simplefunc = new von_mises_hess_4_4< std::vector<double> >(); return;}
            if( name == "neohook_f") { simplefunc = new neohook_f< std::vector<double> >(); return;}
            if( name == "neohook_grad_0") { simplefunc = new neohook_grad_0< std::vector<double> >(); return;}
            if( name == "neohook_grad_1") { simplefunc = new neohook_grad_1< std::vector<double> >(); return;}
            if( name == "neohook_grad_2") { simplefunc = new neohook_grad_2< std::vector<double> >(); return;}
            if( name == "neohook_grad_3") { simplefunc = new neohook_grad_3< std::vector<double> >(); return;}
            if( name == "neohook_grad_4") { simplefunc = new neohook_grad_4< std::vector<double> >(); return;}
            if( name == "neohook_grad_5") { simplefunc = new neohook_grad_5< std::vector<double> >(); return;}
            if( name == "neohook_grad_6") { simplefunc = new neohook_grad_6< std::vector<double> >(); return;}
            if( name == "neohook_hess_0_0") { simplefunc = new neohook_hess_0_0< std::vector<double> >(); return;}
            if( name == "neohook_hess_0_1") { simplefunc = new neohook_hess_0_1< std::vector<double> >(); return;}
            if( name == "neohook_hess_0_2") { simplefunc = new neohook_hess_0_2< std::vector<double> >(); return;}
            if( name == "neohook_hess_0_3") { simplefunc = new neohook_hess_0_3< std::vector<double> >(); return;}
            if( name == "neohook_hess_0_4") { simplefunc = new neohook_hess_0_4< std::vector<double> >(); return;}
            if( name == "neohook_hess_0_5") { simplefunc = new neohook_hess_0_5< std::vector<double> >(); return;}
            if( name == "neohook_hess_0_6") { simplefunc = new neohook_hess_0_6< std::vector<double> >(); return;}
            if( name == "neohook_hess_1_0") { simplefunc = new neohook_hess_1_0< std::vector<double> >(); return;}
            if( name == "neohook_hess_1_1") { simplefunc = new neohook_hess_1_1< std::vector<double> >(); return;}
            if( name == "neohook_hess_1_2") { simplefunc = new neohook_hess_1_2< std::vector<double> >(); return;}
            if( name == "neohook_hess_1_3") { simplefunc = new neohook_hess_1_3< std::vector<double> >(); return;}
            if( name == "neohook_hess_1_4") { simplefunc = new neohook_hess_1_4< std::vector<double> >(); return;}
            if( name == "neohook_hess_1_5") { simplefunc = new neohook_hess_1_5< std::vector<double> >(); return;}
            if( name == "neohook_hess_1_6") { simplefunc = new neohook_hess_1_6< std::vector<double> >(); return;}
            if( name == "neohook_hess_2_0") { simplefunc = new neohook_hess_2_0< std::vector<double> >(); return;}
            if( name == "neohook_hess_2_1") { simplefunc = new neohook_hess_2_1< std::vector<double> >(); return;}
            if( name == "neohook_hess_2_2") { simplefunc = new neohook_hess_2_2< std::vector<double> >(); return;}
            if( name == "neohook_hess_2_3") { simplefunc = new neohook_hess_2_3< std::vector<double> >(); return;}
            if( name == "neohook_hess_2_4") { simplefunc = new neohook_hess_2_4< std::vector<double> >(); return;}
            if( name == "neohook_hess_2_5") { simplefunc = new neohook_hess_2_5< std::vector<double> >(); return;}
            if( name == "neohook_hess_2_6") { simplefunc = new neohook_hess_2_6< std::vector<double> >(); return;}
            if( name == "neohook_hess_3_0") { simplefunc = new neohook_hess_3_0< std::vector<double> >(); return;}
            if( name == "neohook_hess_3_1") { simplefunc = new neohook_hess_3_1< std::vector<double> >(); return;}
            if( name == "neohook_hess_3_2") { simplefunc = new neohook_hess_3_2< std::vector<double> >(); return;}
            if( name == "neohook_hess_3_3") { simplefunc = new neohook_hess_3_3< std::vector<double> >(); return;}
            if( name == "neohook_hess_3_4") { simplefunc = new neohook_hess_3_4< std::vector<double> >(); return;}
            if( name == "neohook_hess_3_5") { simplefunc = new neohook_hess_3_5< std::vector<double> >(); return;}
            if( name == "neohook_hess_3_6") { simplefunc = new neohook_hess_3_6< std::vector<double> >(); return;}
            if( name == "neohook_hess_4_0") { simplefunc = new neohook_hess_4_0< std::vector<double> >(); return;}
            if( name == "neohook_hess_4_1") { simplefunc = new neohook_hess_4_1< std::vector<double> >(); return;}
            if( name == "neohook_hess_4_2") { simplefunc = new neohook_hess_4_2< std::vector<double> >(); return;}
            if( name == "neohook_hess_4_3") { simplefunc = new neohook_hess_4_3< std::vector<double> >(); return;}
            if( name == "neohook_hess_4_4") { simplefunc = new neohook_hess_4_4< std::vector<double> >(); return;}
            if( name == "neohook_hess_4_5") { simplefunc = new neohook_hess_4_5< std::vector<double> >(); return;}
            if( name == "neohook_hess_4_6") { simplefunc = new neohook_hess_4_6< std::vector<double> >(); return;}
            if( name == "neohook_hess_5_0") { simplefunc = new neohook_hess_5_0< std::vector<double> >(); return;}
            if( name == "neohook_hess_5_1") { simplefunc = new neohook_hess_5_1< std::vector<double> >(); return;}
            if( name == "neohook_hess_5_2") { simplefunc = new neohook_hess_5_2< std::vector<double> >(); return;}
            if( name == "neohook_hess_5_3") { simplefunc = new neohook_hess_5_3< std::vector<double> >(); return;}
            if( name == "neohook_hess_5_4") { simplefunc = new neohook_hess_5_4< std::vector<double> >(); return;}
            if( name == "neohook_hess_5_5") { simplefunc = new neohook_hess_5_5< std::vector<double> >(); return;}
            if( name == "neohook_hess_5_6") { simplefunc = new neohook_hess_5_6< std::vector<double> >(); return;}
            if( name == "neohook_hess_6_0") { simplefunc = new neohook_hess_6_0< std::vector<double> >(); return;}
            if( name == "neohook_hess_6_1") { simplefunc = new neohook_hess_6_1< std::vector<double> >(); return;}
            if( name == "neohook_hess_6_2") { simplefunc = new neohook_hess_6_2< std::vector<double> >(); return;}
            if( name == "neohook_hess_6_3") { simplefunc = new neohook_hess_6_3< std::vector<double> >(); return;}
            if( name == "neohook_hess_6_4") { simplefunc = new neohook_hess_6_4< std::vector<double> >(); return;}
            if( name == "neohook_hess_6_5") { simplefunc = new neohook_hess_6_5< std::vector<double> >(); return;}
            if( name == "neohook_hess_6_6") { simplefunc = new neohook_hess_6_6< std::vector<double> >(); return;}
            throw std::runtime_error( "PSimpleBase< std::vector<double>, double > " + name + " was not found in the PLibrary");
        }
        void PLibrary::checkout( std::string name, PSimpleBase< double*, double > *&simplefunc)
        {
            if( name == "stvenkir_f") { simplefunc = new stvenkir_f< double* >(); return;}
            if( name == "stvenkir_grad_0") { simplefunc = new stvenkir_grad_0< double* >(); return;}
            if( name == "stvenkir_grad_1") { simplefunc = new stvenkir_grad_1< double* >(); return;}
            if( name == "stvenkir_grad_2") { simplefunc = new stvenkir_grad_2< double* >(); return;}
            if( name == "stvenkir_grad_3") { simplefunc = new stvenkir_grad_3< double* >(); return;}
            if( name == "stvenkir_grad_4") { simplefunc = new stvenkir_grad_4< double* >(); return;}
            if( name == "stvenkir_grad_5") { simplefunc = new stvenkir_grad_5< double* >(); return;}
            if( name == "stvenkir_grad_6") { simplefunc = new stvenkir_grad_6< double* >(); return;}
            if( name == "stvenkir_hess_0_0") { simplefunc = new stvenkir_hess_0_0< double* >(); return;}
            if( name == "stvenkir_hess_0_1") { simplefunc = new stvenkir_hess_0_1< double* >(); return;}
            if( name == "stvenkir_hess_0_2") { simplefunc = new stvenkir_hess_0_2< double* >(); return;}
            if( name == "stvenkir_hess_0_3") { simplefunc = new stvenkir_hess_0_3< double* >(); return;}
            if( name == "stvenkir_hess_0_4") { simplefunc = new stvenkir_hess_0_4< double* >(); return;}
            if( name == "stvenkir_hess_0_5") { simplefunc = new stvenkir_hess_0_5< double* >(); return;}
            if( name == "stvenkir_hess_0_6") { simplefunc = new stvenkir_hess_0_6< double* >(); return;}
            if( name == "stvenkir_hess_1_0") { simplefunc = new stvenkir_hess_1_0< double* >(); return;}
            if( name == "stvenkir_hess_1_1") { simplefunc = new stvenkir_hess_1_1< double* >(); return;}
            if( name == "stvenkir_hess_1_2") { simplefunc = new stvenkir_hess_1_2< double* >(); return;}
            if( name == "stvenkir_hess_1_3") { simplefunc = new stvenkir_hess_1_3< double* >(); return;}
            if( name == "stvenkir_hess_1_4") { simplefunc = new stvenkir_hess_1_4< double* >(); return;}
            if( name == "stvenkir_hess_1_5") { simplefunc = new stvenkir_hess_1_5< double* >(); return;}
            if( name == "stvenkir_hess_1_6") { simplefunc = new stvenkir_hess_1_6< double* >(); return;}
            if( name == "stvenkir_hess_2_0") { simplefunc = new stvenkir_hess_2_0< double* >(); return;}
            if( name == "stvenkir_hess_2_1") { simplefunc = new stvenkir_hess_2_1< double* >(); return;}
            if( name == "stvenkir_hess_2_2") { simplefunc = new stvenkir_hess_2_2< double* >(); return;}
            if( name == "stvenkir_hess_2_3") { simplefunc = new stvenkir_hess_2_3< double* >(); return;}
            if( name == "stvenkir_hess_2_4") { simplefunc = new stvenkir_hess_2_4< double* >(); return;}
            if( name == "stvenkir_hess_2_5") { simplefunc = new stvenkir_hess_2_5< double* >(); return;}
            if( name == "stvenkir_hess_2_6") { simplefunc = new stvenkir_hess_2_6< double* >(); return;}
            if( name == "stvenkir_hess_3_0") { simplefunc = new stvenkir_hess_3_0< double* >(); return;}
            if( name == "stvenkir_hess_3_1") { simplefunc = new stvenkir_hess_3_1< double* >(); return;}
            if( name == "stvenkir_hess_3_2") { simplefunc = new stvenkir_hess_3_2< double* >(); return;}
            if( name == "stvenkir_hess_3_3") { simplefunc = new stvenkir_hess_3_3< double* >(); return;}
            if( name == "stvenkir_hess_3_4") { simplefunc = new stvenkir_hess_3_4< double* >(); return;}
            if( name == "stvenkir_hess_3_5") { simplefunc = new stvenkir_hess_3_5< double* >(); return;}
            if( name == "stvenkir_hess_3_6") { simplefunc = new stvenkir_hess_3_6< double* >(); return;}
            if( name == "stvenkir_hess_4_0") { simplefunc = new stvenkir_hess_4_0< double* >(); return;}
            if( name == "stvenkir_hess_4_1") { simplefunc = new stvenkir_hess_4_1< double* >(); return;}
            if( name == "stvenkir_hess_4_2") { simplefunc = new stvenkir_hess_4_2< double* >(); return;}
            if( name == "stvenkir_hess_4_3") { simplefunc = new stvenkir_hess_4_3< double* >(); return;}
            if( name == "stvenkir_hess_4_4") { simplefunc = new stvenkir_hess_4_4< double* >(); return;}
            if( name == "stvenkir_hess_4_5") { simplefunc = new stvenkir_hess_4_5< double* >(); return;}
            if( name == "stvenkir_hess_4_6") { simplefunc = new stvenkir_hess_4_6< double* >(); return;}
            if( name == "stvenkir_hess_5_0") { simplefunc = new stvenkir_hess_5_0< double* >(); return;}
            if( name == "stvenkir_hess_5_1") { simplefunc = new stvenkir_hess_5_1< double* >(); return;}
            if( name == "stvenkir_hess_5_2") { simplefunc = new stvenkir_hess_5_2< double* >(); return;}
            if( name == "stvenkir_hess_5_3") { simplefunc = new stvenkir_hess_5_3< double* >(); return;}
            if( name == "stvenkir_hess_5_4") { simplefunc = new stvenkir_hess_5_4< double* >(); return;}
            if( name == "stvenkir_hess_5_5") { simplefunc = new stvenkir_hess_5_5< double* >(); return;}
            if( name == "stvenkir_hess_5_6") { simplefunc = new stvenkir_hess_5_6< double* >(); return;}
            if( name == "stvenkir_hess_6_0") { simplefunc = new stvenkir_hess_6_0< double* >(); return;}
            if( name == "stvenkir_hess_6_1") { simplefunc = new stvenkir_hess_6_1< double* >(); return;}
            if( name == "stvenkir_hess_6_2") { simplefunc = new stvenkir_hess_6_2< double* >(); return;}
            if( name == "stvenkir_hess_6_3") { simplefunc = new stvenkir_hess_6_3< double* >(); return;}
            if( name == "stvenkir_hess_6_4") { simplefunc = new stvenkir_hess_6_4< double* >(); return;}
            if( name == "stvenkir_hess_6_5") { simplefunc = new stvenkir_hess_6_5< double* >(); return;}
            if( name == "stvenkir_hess_6_6") { simplefunc = new stvenkir_hess_6_6< double* >(); return;}
            if( name == "quadlog_f") { simplefunc = new quadlog_f< double* >(); return;}
            if( name == "quadlog_grad_0") { simplefunc = new quadlog_grad_0< double* >(); return;}
            if( name == "quadlog_grad_1") { simplefunc = new quadlog_grad_1< double* >(); return;}
            if( name == "quadlog_grad_2") { simplefunc = new quadlog_grad_2< double* >(); return;}
            if( name == "quadlog_grad_3") { simplefunc = new quadlog_grad_3< double* >(); return;}
            if( name == "quadlog_grad_4") { simplefunc = new quadlog_grad_4< double* >(); return;}
            if( name == "quadlog_grad_5") { simplefunc = new quadlog_grad_5< double* >(); return;}
            if( name == "quadlog_grad_6") { simplefunc = new quadlog_grad_6< double* >(); return;}
            if( name == "quadlog_hess_0_0") { simplefunc = new quadlog_hess_0_0< double* >(); return;}
            if( name == "quadlog_hess_0_1") { simplefunc = new quadlog_hess_0_1< double* >(); return;}
            if( name == "quadlog_hess_0_2") { simplefunc = new quadlog_hess_0_2< double* >(); return;}
            if( name == "quadlog_hess_0_3") { simplefunc = new quadlog_hess_0_3< double* >(); return;}
            if( name == "quadlog_hess_0_4") { simplefunc = new quadlog_hess_0_4< double* >(); return;}
            if( name == "quadlog_hess_0_5") { simplefunc = new quadlog_hess_0_5< double* >(); return;}
            if( name == "quadlog_hess_0_6") { simplefunc = new quadlog_hess_0_6< double* >(); return;}
            if( name == "quadlog_hess_1_0") { simplefunc = new quadlog_hess_1_0< double* >(); return;}
            if( name == "quadlog_hess_1_1") { simplefunc = new quadlog_hess_1_1< double* >(); return;}
            if( name == "quadlog_hess_1_2") { simplefunc = new quadlog_hess_1_2< double* >(); return;}
            if( name == "quadlog_hess_1_3") { simplefunc = new quadlog_hess_1_3< double* >(); return;}
            if( name == "quadlog_hess_1_4") { simplefunc = new quadlog_hess_1_4< double* >(); return;}
            if( name == "quadlog_hess_1_5") { simplefunc = new quadlog_hess_1_5< double* >(); return;}
            if( name == "quadlog_hess_1_6") { simplefunc = new quadlog_hess_1_6< double* >(); return;}
            if( name == "quadlog_hess_2_0") { simplefunc = new quadlog_hess_2_0< double* >(); return;}
            if( name == "quadlog_hess_2_1") { simplefunc = new quadlog_hess_2_1< double* >(); return;}
            if( name == "quadlog_hess_2_2") { simplefunc = new quadlog_hess_2_2< double* >(); return;}
            if( name == "quadlog_hess_2_3") { simplefunc = new quadlog_hess_2_3< double* >(); return;}
            if( name == "quadlog_hess_2_4") { simplefunc = new quadlog_hess_2_4< double* >(); return;}
            if( name == "quadlog_hess_2_5") { simplefunc = new quadlog_hess_2_5< double* >(); return;}
            if( name == "quadlog_hess_2_6") { simplefunc = new quadlog_hess_2_6< double* >(); return;}
            if( name == "quadlog_hess_3_0") { simplefunc = new quadlog_hess_3_0< double* >(); return;}
            if( name == "quadlog_hess_3_1") { simplefunc = new quadlog_hess_3_1< double* >(); return;}
            if( name == "quadlog_hess_3_2") { simplefunc = new quadlog_hess_3_2< double* >(); return;}
            if( name == "quadlog_hess_3_3") { simplefunc = new quadlog_hess_3_3< double* >(); return;}
            if( name == "quadlog_hess_3_4") { simplefunc = new quadlog_hess_3_4< double* >(); return;}
            if( name == "quadlog_hess_3_5") { simplefunc = new quadlog_hess_3_5< double* >(); return;}
            if( name == "quadlog_hess_3_6") { simplefunc = new quadlog_hess_3_6< double* >(); return;}
            if( name == "quadlog_hess_4_0") { simplefunc = new quadlog_hess_4_0< double* >(); return;}
            if( name == "quadlog_hess_4_1") { simplefunc = new quadlog_hess_4_1< double* >(); return;}
            if( name == "quadlog_hess_4_2") { simplefunc = new quadlog_hess_4_2< double* >(); return;}
            if( name == "quadlog_hess_4_3") { simplefunc = new quadlog_hess_4_3< double* >(); return;}
            if( name == "quadlog_hess_4_4") { simplefunc = new quadlog_hess_4_4< double* >(); return;}
            if( name == "quadlog_hess_4_5") { simplefunc = new quadlog_hess_4_5< double* >(); return;}
            if( name == "quadlog_hess_4_6") { simplefunc = new quadlog_hess_4_6< double* >(); return;}
            if( name == "quadlog_hess_5_0") { simplefunc = new quadlog_hess_5_0< double* >(); return;}
            if( name == "quadlog_hess_5_1") { simplefunc = new quadlog_hess_5_1< double* >(); return;}
            if( name == "quadlog_hess_5_2") { simplefunc = new quadlog_hess_5_2< double* >(); return;}
            if( name == "quadlog_hess_5_3") { simplefunc = new quadlog_hess_5_3< double* >(); return;}
            if( name == "quadlog_hess_5_4") { simplefunc = new quadlog_hess_5_4< double* >(); return;}
            if( name == "quadlog_hess_5_5") { simplefunc = new quadlog_hess_5_5< double* >(); return;}
            if( name == "quadlog_hess_5_6") { simplefunc = new quadlog_hess_5_6< double* >(); return;}
            if( name == "quadlog_hess_6_0") { simplefunc = new quadlog_hess_6_0< double* >(); return;}
            if( name == "quadlog_hess_6_1") { simplefunc = new quadlog_hess_6_1< double* >(); return;}
            if( name == "quadlog_hess_6_2") { simplefunc = new quadlog_hess_6_2< double* >(); return;}
            if( name == "quadlog_hess_6_3") { simplefunc = new quadlog_hess_6_3< double* >(); return;}
            if( name == "quadlog_hess_6_4") { simplefunc = new quadlog_hess_6_4< double* >(); return;}
            if( name == "quadlog_hess_6_5") { simplefunc = new quadlog_hess_6_5< double* >(); return;}
            if( name == "quadlog_hess_6_6") { simplefunc = new quadlog_hess_6_6< double* >(); return;}
            if( name == "hardening_f") { simplefunc = new hardening_f< double* >(); return;}
            if( name == "hardening_grad_0") { simplefunc = new hardening_grad_0< double* >(); return;}
            if( name == "hardening_hess_0_0") { simplefunc = new hardening_hess_0_0< double* >(); return;}
            if( name == "von_mises_f") { simplefunc = new von_mises_f< double* >(); return;}
            if( name == "von_mises_grad_0") { simplefunc = new von_mises_grad_0< double* >(); return;}
            if( name == "von_mises_grad_1") { simplefunc = new von_mises_grad_1< double* >(); return;}
            if( name == "von_mises_grad_2") { simplefunc = new von_mises_grad_2< double* >(); return;}
            if( name == "von_mises_grad_3") { simplefunc = new von_mises_grad_3< double* >(); return;}
            if( name == "von_mises_grad_4") { simplefunc = new von_mises_grad_4< double* >(); return;}
            if( name == "von_mises_hess_0_0") { simplefunc = new von_mises_hess_0_0< double* >(); return;}
            if( name == "von_mises_hess_0_1") { simplefunc = new von_mises_hess_0_1< double* >(); return;}
            if( name == "von_mises_hess_0_2") { simplefunc = new von_mises_hess_0_2< double* >(); return;}
            if( name == "von_mises_hess_0_3") { simplefunc = new von_mises_hess_0_3< double* >(); return;}
            if( name == "von_mises_hess_0_4") { simplefunc = new von_mises_hess_0_4< double* >(); return;}
            if( name == "von_mises_hess_1_0") { simplefunc = new von_mises_hess_1_0< double* >(); return;}
            if( name == "von_mises_hess_1_1") { simplefunc = new von_mises_hess_1_1< double* >(); return;}
            if( name == "von_mises_hess_1_2") { simplefunc = new von_mises_hess_1_2< double* >(); return;}
            if( name == "von_mises_hess_1_3") { simplefunc = new von_mises_hess_1_3< double* >(); return;}
            if( name == "von_mises_hess_1_4") { simplefunc = new von_mises_hess_1_4< double* >(); return;}
            if( name == "von_mises_hess_2_0") { simplefunc = new von_mises_hess_2_0< double* >(); return;}
            if( name == "von_mises_hess_2_1") { simplefunc = new von_mises_hess_2_1< double* >(); return;}
            if( name == "von_mises_hess_2_2") { simplefunc = new von_mises_hess_2_2< double* >(); return;}
            if( name == "von_mises_hess_2_3") { simplefunc = new von_mises_hess_2_3< double* >(); return;}
            if( name == "von_mises_hess_2_4") { simplefunc = new von_mises_hess_2_4< double* >(); return;}
            if( name == "von_mises_hess_3_0") { simplefunc = new von_mises_hess_3_0< double* >(); return;}
            if( name == "von_mises_hess_3_1") { simplefunc = new von_mises_hess_3_1< double* >(); return;}
            if( name == "von_mises_hess_3_2") { simplefunc = new von_mises_hess_3_2< double* >(); return;}
            if( name == "von_mises_hess_3_3") { simplefunc = new von_mises_hess_3_3< double* >(); return;}
            if( name == "von_mises_hess_3_4") { simplefunc = new von_mises_hess_3_4< double* >(); return;}
            if( name == "von_mises_hess_4_0") { simplefunc = new von_mises_hess_4_0< double* >(); return;}
            if( name == "von_mises_hess_4_1") { simplefunc = new von_mises_hess_4_1< double* >(); return;}
            if( name == "von_mises_hess_4_2") { simplefunc = new von_mises_hess_4_2< double* >(); return;}
            if( name == "von_mises_hess_4_3") { simplefunc = new von_mises_hess_4_3< double* >(); return;}
            if( name == "von_mises_hess_4_4") { simplefunc = new von_mises_hess_4_4< double* >(); return;}
            if( name == "neohook_f") { simplefunc = new neohook_f< double* >(); return;}
            if( name == "neohook_grad_0") { simplefunc = new neohook_grad_0< double* >(); return;}
            if( name == "neohook_grad_1") { simplefunc = new neohook_grad_1< double* >(); return;}
            if( name == "neohook_grad_2") { simplefunc = new neohook_grad_2< double* >(); return;}
            if( name == "neohook_grad_3") { simplefunc = new neohook_grad_3< double* >(); return;}
            if( name == "neohook_grad_4") { simplefunc = new neohook_grad_4< double* >(); return;}
            if( name == "neohook_grad_5") { simplefunc = new neohook_grad_5< double* >(); return;}
            if( name == "neohook_grad_6") { simplefunc = new neohook_grad_6< double* >(); return;}
            if( name == "neohook_hess_0_0") { simplefunc = new neohook_hess_0_0< double* >(); return;}
            if( name == "neohook_hess_0_1") { simplefunc = new neohook_hess_0_1< double* >(); return;}
            if( name == "neohook_hess_0_2") { simplefunc = new neohook_hess_0_2< double* >(); return;}
            if( name == "neohook_hess_0_3") { simplefunc = new neohook_hess_0_3< double* >(); return;}
            if( name == "neohook_hess_0_4") { simplefunc = new neohook_hess_0_4< double* >(); return;}
            if( name == "neohook_hess_0_5") { simplefunc = new neohook_hess_0_5< double* >(); return;}
            if( name == "neohook_hess_0_6") { simplefunc = new neohook_hess_0_6< double* >(); return;}
            if( name == "neohook_hess_1_0") { simplefunc = new neohook_hess_1_0< double* >(); return;}
            if( name == "neohook_hess_1_1") { simplefunc = new neohook_hess_1_1< double* >(); return;}
            if( name == "neohook_hess_1_2") { simplefunc = new neohook_hess_1_2< double* >(); return;}
            if( name == "neohook_hess_1_3") { simplefunc = new neohook_hess_1_3< double* >(); return;}
            if( name == "neohook_hess_1_4") { simplefunc = new neohook_hess_1_4< double* >(); return;}
            if( name == "neohook_hess_1_5") { simplefunc = new neohook_hess_1_5< double* >(); return;}
            if( name == "neohook_hess_1_6") { simplefunc = new neohook_hess_1_6< double* >(); return;}
            if( name == "neohook_hess_2_0") { simplefunc = new neohook_hess_2_0< double* >(); return;}
            if( name == "neohook_hess_2_1") { simplefunc = new neohook_hess_2_1< double* >(); return;}
            if( name == "neohook_hess_2_2") { simplefunc = new neohook_hess_2_2< double* >(); return;}
            if( name == "neohook_hess_2_3") { simplefunc = new neohook_hess_2_3< double* >(); return;}
            if( name == "neohook_hess_2_4") { simplefunc = new neohook_hess_2_4< double* >(); return;}
            if( name == "neohook_hess_2_5") { simplefunc = new neohook_hess_2_5< double* >(); return;}
            if( name == "neohook_hess_2_6") { simplefunc = new neohook_hess_2_6< double* >(); return;}
            if( name == "neohook_hess_3_0") { simplefunc = new neohook_hess_3_0< double* >(); return;}
            if( name == "neohook_hess_3_1") { simplefunc = new neohook_hess_3_1< double* >(); return;}
            if( name == "neohook_hess_3_2") { simplefunc = new neohook_hess_3_2< double* >(); return;}
            if( name == "neohook_hess_3_3") { simplefunc = new neohook_hess_3_3< double* >(); return;}
            if( name == "neohook_hess_3_4") { simplefunc = new neohook_hess_3_4< double* >(); return;}
            if( name == "neohook_hess_3_5") { simplefunc = new neohook_hess_3_5< double* >(); return;}
            if( name == "neohook_hess_3_6") { simplefunc = new neohook_hess_3_6< double* >(); return;}
            if( name == "neohook_hess_4_0") { simplefunc = new neohook_hess_4_0< double* >(); return;}
            if( name == "neohook_hess_4_1") { simplefunc = new neohook_hess_4_1< double* >(); return;}
            if( name == "neohook_hess_4_2") { simplefunc = new neohook_hess_4_2< double* >(); return;}
            if( name == "neohook_hess_4_3") { simplefunc = new neohook_hess_4_3< double* >(); return;}
            if( name == "neohook_hess_4_4") { simplefunc = new neohook_hess_4_4< double* >(); return;}
            if( name == "neohook_hess_4_5") { simplefunc = new neohook_hess_4_5< double* >(); return;}
            if( name == "neohook_hess_4_6") { simplefunc = new neohook_hess_4_6< double* >(); return;}
            if( name == "neohook_hess_5_0") { simplefunc = new neohook_hess_5_0< double* >(); return;}
            if( name == "neohook_hess_5_1") { simplefunc = new neohook_hess_5_1< double* >(); return;}
            if( name == "neohook_hess_5_2") { simplefunc = new neohook_hess_5_2< double* >(); return;}
            if( name == "neohook_hess_5_3") { simplefunc = new neohook_hess_5_3< double* >(); return;}
            if( name == "neohook_hess_5_4") { simplefunc = new neohook_hess_5_4< double* >(); return;}
            if( name == "neohook_hess_5_5") { simplefunc = new neohook_hess_5_5< double* >(); return;}
            if( name == "neohook_hess_5_6") { simplefunc = new neohook_hess_5_6< double* >(); return;}
            if( name == "neohook_hess_6_0") { simplefunc = new neohook_hess_6_0< double* >(); return;}
            if( name == "neohook_hess_6_1") { simplefunc = new neohook_hess_6_1< double* >(); return;}
            if( name == "neohook_hess_6_2") { simplefunc = new neohook_hess_6_2< double* >(); return;}
            if( name == "neohook_hess_6_3") { simplefunc = new neohook_hess_6_3< double* >(); return;}
            if( name == "neohook_hess_6_4") { simplefunc = new neohook_hess_6_4< double* >(); return;}
            if( name == "neohook_hess_6_5") { simplefunc = new neohook_hess_6_5< double* >(); return;}
            if( name == "neohook_hess_6_6") { simplefunc = new neohook_hess_6_6< double* >(); return;}
            throw std::runtime_error( "PSimpleBase< double*, double > " + name + " was not found in the PLibrary");
>>>>>>> 1853fc05524906a2a2d829468a5c2e263978073c
        }

        void PLibrary::checkout( std::string name, PFuncBase< std::vector<double>, double > *&func)
        {
<<<<<<< HEAD
            if( name == "neohook") func = new neohook< std::vector<double> >();
            if( name == "stvenkir") func = new stvenkir< std::vector<double> >();
            if( name == "linear_hardening") func = new linear_hardening< std::vector<double> >();
            if( name == "Cu_hardening") func = new Cu_hardening< std::vector<double> >();
            if( name == "von_mises") func = new von_mises< std::vector<double> >();
            if( name == "quadlog") func = new quadlog< std::vector<double> >();
        }
        void PLibrary::checkout( std::string name, PFuncBase< double*, double > *&func)
        {
            if( name == "neohook") func = new neohook< double* >();
            if( name == "stvenkir") func = new stvenkir< double* >();
            if( name == "linear_hardening") func = new linear_hardening< double* >();
            if( name == "Cu_hardening") func = new Cu_hardening< double* >();
            if( name == "von_mises") func = new von_mises< double* >();
            if( name == "quadlog") func = new quadlog< double* >();
=======
            if( name == "stvenkir") { func = new stvenkir< std::vector<double> >(); return;}
            if( name == "quadlog") { func = new quadlog< std::vector<double> >(); return;}
            if( name == "hardening") { func = new hardening< std::vector<double> >(); return;}
            if( name == "von_mises") { func = new von_mises< std::vector<double> >(); return;}
            if( name == "neohook") { func = new neohook< std::vector<double> >(); return;}
            throw std::runtime_error( "PFuncBase< std::vector<double>, double > " + name + " was not found in the PLibrary");
        }
        void PLibrary::checkout( std::string name, PFuncBase< double*, double > *&func)
        {
            if( name == "stvenkir") { func = new stvenkir< double* >(); return;}
            if( name == "quadlog") { func = new quadlog< double* >(); return;}
            if( name == "hardening") { func = new hardening< double* >(); return;}
            if( name == "von_mises") { func = new von_mises< double* >(); return;}
            if( name == "neohook") { func = new neohook< double* >(); return;}
            throw std::runtime_error( "PFuncBase< double*, double > " + name + " was not found in the PLibrary");
>>>>>>> 1853fc05524906a2a2d829468a5c2e263978073c
        }




}


#endif
