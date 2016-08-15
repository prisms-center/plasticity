// created: 2016-8-5 14:00:38
// version: HEAD
// url: git@github.com:prisms-center/IntegrationTools.git
// commit: 010dc90537a86ed3458c4b12292b6af9545b3c1f

#ifndef PLIBRARY_CC
#define PLIBRARY_CC

#include<cstring>
#include<stdexcept>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wextra"

#include "neohook.hh"
#include "stvenkir.hh"
#include "linear_hardening.hh"
#include "Cu_hardening.hh"
#include "von_mises.hh"
#include "quadlog.hh"
#include "PLibrary.hh"

#pragma GCC diagnostic pop

namespace PRISMS
{

        void PLibrary::checkout( std::string name, PSimpleFunction< std::vector<double>, double > &simplefunc)
        {
            typedef PSimpleFunction< std::vector<double>, double > psf;
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
            if( name == "linear_hardening_f") { simplefunc = psf( linear_hardening_f< std::vector<double> >() ); return;}
            if( name == "linear_hardening_grad_0") { simplefunc = psf( linear_hardening_grad_0< std::vector<double> >() ); return;}
            if( name == "linear_hardening_grad_1") { simplefunc = psf( linear_hardening_grad_1< std::vector<double> >() ); return;}
            if( name == "linear_hardening_hess_0_0") { simplefunc = psf( linear_hardening_hess_0_0< std::vector<double> >() ); return;}
            if( name == "linear_hardening_hess_0_1") { simplefunc = psf( linear_hardening_hess_0_1< std::vector<double> >() ); return;}
            if( name == "linear_hardening_hess_1_0") { simplefunc = psf( linear_hardening_hess_1_0< std::vector<double> >() ); return;}
            if( name == "linear_hardening_hess_1_1") { simplefunc = psf( linear_hardening_hess_1_1< std::vector<double> >() ); return;}
            if( name == "Cu_hardening_f") { simplefunc = psf( Cu_hardening_f< std::vector<double> >() ); return;}
            if( name == "Cu_hardening_grad_0") { simplefunc = psf( Cu_hardening_grad_0< std::vector<double> >() ); return;}
            if( name == "Cu_hardening_grad_1") { simplefunc = psf( Cu_hardening_grad_1< std::vector<double> >() ); return;}
            if( name == "Cu_hardening_hess_0_0") { simplefunc = psf( Cu_hardening_hess_0_0< std::vector<double> >() ); return;}
            if( name == "Cu_hardening_hess_0_1") { simplefunc = psf( Cu_hardening_hess_0_1< std::vector<double> >() ); return;}
            if( name == "Cu_hardening_hess_1_0") { simplefunc = psf( Cu_hardening_hess_1_0< std::vector<double> >() ); return;}
            if( name == "Cu_hardening_hess_1_1") { simplefunc = psf( Cu_hardening_hess_1_1< std::vector<double> >() ); return;}
            if( name == "von_mises_f") { simplefunc = psf( von_mises_f< std::vector<double> >() ); return;}
            if( name == "von_mises_grad_0") { simplefunc = psf( von_mises_grad_0< std::vector<double> >() ); return;}
            if( name == "von_mises_grad_1") { simplefunc = psf( von_mises_grad_1< std::vector<double> >() ); return;}
            if( name == "von_mises_grad_2") { simplefunc = psf( von_mises_grad_2< std::vector<double> >() ); return;}
            if( name == "von_mises_grad_3") { simplefunc = psf( von_mises_grad_3< std::vector<double> >() ); return;}
            if( name == "von_mises_grad_4") { simplefunc = psf( von_mises_grad_4< std::vector<double> >() ); return;}
            if( name == "von_mises_grad_5") { simplefunc = psf( von_mises_grad_5< std::vector<double> >() ); return;}
            if( name == "von_mises_grad_6") { simplefunc = psf( von_mises_grad_6< std::vector<double> >() ); return;}
            if( name == "von_mises_grad_7") { simplefunc = psf( von_mises_grad_7< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_0_0") { simplefunc = psf( von_mises_hess_0_0< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_0_1") { simplefunc = psf( von_mises_hess_0_1< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_0_2") { simplefunc = psf( von_mises_hess_0_2< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_0_3") { simplefunc = psf( von_mises_hess_0_3< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_0_4") { simplefunc = psf( von_mises_hess_0_4< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_0_5") { simplefunc = psf( von_mises_hess_0_5< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_0_6") { simplefunc = psf( von_mises_hess_0_6< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_0_7") { simplefunc = psf( von_mises_hess_0_7< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_1_0") { simplefunc = psf( von_mises_hess_1_0< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_1_1") { simplefunc = psf( von_mises_hess_1_1< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_1_2") { simplefunc = psf( von_mises_hess_1_2< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_1_3") { simplefunc = psf( von_mises_hess_1_3< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_1_4") { simplefunc = psf( von_mises_hess_1_4< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_1_5") { simplefunc = psf( von_mises_hess_1_5< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_1_6") { simplefunc = psf( von_mises_hess_1_6< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_1_7") { simplefunc = psf( von_mises_hess_1_7< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_2_0") { simplefunc = psf( von_mises_hess_2_0< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_2_1") { simplefunc = psf( von_mises_hess_2_1< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_2_2") { simplefunc = psf( von_mises_hess_2_2< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_2_3") { simplefunc = psf( von_mises_hess_2_3< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_2_4") { simplefunc = psf( von_mises_hess_2_4< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_2_5") { simplefunc = psf( von_mises_hess_2_5< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_2_6") { simplefunc = psf( von_mises_hess_2_6< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_2_7") { simplefunc = psf( von_mises_hess_2_7< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_3_0") { simplefunc = psf( von_mises_hess_3_0< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_3_1") { simplefunc = psf( von_mises_hess_3_1< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_3_2") { simplefunc = psf( von_mises_hess_3_2< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_3_3") { simplefunc = psf( von_mises_hess_3_3< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_3_4") { simplefunc = psf( von_mises_hess_3_4< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_3_5") { simplefunc = psf( von_mises_hess_3_5< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_3_6") { simplefunc = psf( von_mises_hess_3_6< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_3_7") { simplefunc = psf( von_mises_hess_3_7< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_4_0") { simplefunc = psf( von_mises_hess_4_0< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_4_1") { simplefunc = psf( von_mises_hess_4_1< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_4_2") { simplefunc = psf( von_mises_hess_4_2< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_4_3") { simplefunc = psf( von_mises_hess_4_3< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_4_4") { simplefunc = psf( von_mises_hess_4_4< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_4_5") { simplefunc = psf( von_mises_hess_4_5< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_4_6") { simplefunc = psf( von_mises_hess_4_6< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_4_7") { simplefunc = psf( von_mises_hess_4_7< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_5_0") { simplefunc = psf( von_mises_hess_5_0< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_5_1") { simplefunc = psf( von_mises_hess_5_1< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_5_2") { simplefunc = psf( von_mises_hess_5_2< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_5_3") { simplefunc = psf( von_mises_hess_5_3< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_5_4") { simplefunc = psf( von_mises_hess_5_4< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_5_5") { simplefunc = psf( von_mises_hess_5_5< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_5_6") { simplefunc = psf( von_mises_hess_5_6< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_5_7") { simplefunc = psf( von_mises_hess_5_7< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_6_0") { simplefunc = psf( von_mises_hess_6_0< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_6_1") { simplefunc = psf( von_mises_hess_6_1< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_6_2") { simplefunc = psf( von_mises_hess_6_2< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_6_3") { simplefunc = psf( von_mises_hess_6_3< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_6_4") { simplefunc = psf( von_mises_hess_6_4< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_6_5") { simplefunc = psf( von_mises_hess_6_5< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_6_6") { simplefunc = psf( von_mises_hess_6_6< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_6_7") { simplefunc = psf( von_mises_hess_6_7< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_7_0") { simplefunc = psf( von_mises_hess_7_0< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_7_1") { simplefunc = psf( von_mises_hess_7_1< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_7_2") { simplefunc = psf( von_mises_hess_7_2< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_7_3") { simplefunc = psf( von_mises_hess_7_3< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_7_4") { simplefunc = psf( von_mises_hess_7_4< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_7_5") { simplefunc = psf( von_mises_hess_7_5< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_7_6") { simplefunc = psf( von_mises_hess_7_6< std::vector<double> >() ); return;}
            if( name == "von_mises_hess_7_7") { simplefunc = psf( von_mises_hess_7_7< std::vector<double> >() ); return;}
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
            else throw std::runtime_error( "PSimpleFunction< std::vector<double>, double > " + name + " was not found in the PLibrary");
        }
        void PLibrary::checkout( std::string name, PSimpleFunction< double*, double > &simplefunc)
        {
            typedef PSimpleFunction< double*, double > psf;
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
            if( name == "linear_hardening_f") { simplefunc = psf( linear_hardening_f< double* >() ); return;}
            if( name == "linear_hardening_grad_0") { simplefunc = psf( linear_hardening_grad_0< double* >() ); return;}
            if( name == "linear_hardening_grad_1") { simplefunc = psf( linear_hardening_grad_1< double* >() ); return;}
            if( name == "linear_hardening_hess_0_0") { simplefunc = psf( linear_hardening_hess_0_0< double* >() ); return;}
            if( name == "linear_hardening_hess_0_1") { simplefunc = psf( linear_hardening_hess_0_1< double* >() ); return;}
            if( name == "linear_hardening_hess_1_0") { simplefunc = psf( linear_hardening_hess_1_0< double* >() ); return;}
            if( name == "linear_hardening_hess_1_1") { simplefunc = psf( linear_hardening_hess_1_1< double* >() ); return;}
            if( name == "Cu_hardening_f") { simplefunc = psf( Cu_hardening_f< double* >() ); return;}
            if( name == "Cu_hardening_grad_0") { simplefunc = psf( Cu_hardening_grad_0< double* >() ); return;}
            if( name == "Cu_hardening_grad_1") { simplefunc = psf( Cu_hardening_grad_1< double* >() ); return;}
            if( name == "Cu_hardening_hess_0_0") { simplefunc = psf( Cu_hardening_hess_0_0< double* >() ); return;}
            if( name == "Cu_hardening_hess_0_1") { simplefunc = psf( Cu_hardening_hess_0_1< double* >() ); return;}
            if( name == "Cu_hardening_hess_1_0") { simplefunc = psf( Cu_hardening_hess_1_0< double* >() ); return;}
            if( name == "Cu_hardening_hess_1_1") { simplefunc = psf( Cu_hardening_hess_1_1< double* >() ); return;}
            if( name == "von_mises_f") { simplefunc = psf( von_mises_f< double* >() ); return;}
            if( name == "von_mises_grad_0") { simplefunc = psf( von_mises_grad_0< double* >() ); return;}
            if( name == "von_mises_grad_1") { simplefunc = psf( von_mises_grad_1< double* >() ); return;}
            if( name == "von_mises_grad_2") { simplefunc = psf( von_mises_grad_2< double* >() ); return;}
            if( name == "von_mises_grad_3") { simplefunc = psf( von_mises_grad_3< double* >() ); return;}
            if( name == "von_mises_grad_4") { simplefunc = psf( von_mises_grad_4< double* >() ); return;}
            if( name == "von_mises_grad_5") { simplefunc = psf( von_mises_grad_5< double* >() ); return;}
            if( name == "von_mises_grad_6") { simplefunc = psf( von_mises_grad_6< double* >() ); return;}
            if( name == "von_mises_grad_7") { simplefunc = psf( von_mises_grad_7< double* >() ); return;}
            if( name == "von_mises_hess_0_0") { simplefunc = psf( von_mises_hess_0_0< double* >() ); return;}
            if( name == "von_mises_hess_0_1") { simplefunc = psf( von_mises_hess_0_1< double* >() ); return;}
            if( name == "von_mises_hess_0_2") { simplefunc = psf( von_mises_hess_0_2< double* >() ); return;}
            if( name == "von_mises_hess_0_3") { simplefunc = psf( von_mises_hess_0_3< double* >() ); return;}
            if( name == "von_mises_hess_0_4") { simplefunc = psf( von_mises_hess_0_4< double* >() ); return;}
            if( name == "von_mises_hess_0_5") { simplefunc = psf( von_mises_hess_0_5< double* >() ); return;}
            if( name == "von_mises_hess_0_6") { simplefunc = psf( von_mises_hess_0_6< double* >() ); return;}
            if( name == "von_mises_hess_0_7") { simplefunc = psf( von_mises_hess_0_7< double* >() ); return;}
            if( name == "von_mises_hess_1_0") { simplefunc = psf( von_mises_hess_1_0< double* >() ); return;}
            if( name == "von_mises_hess_1_1") { simplefunc = psf( von_mises_hess_1_1< double* >() ); return;}
            if( name == "von_mises_hess_1_2") { simplefunc = psf( von_mises_hess_1_2< double* >() ); return;}
            if( name == "von_mises_hess_1_3") { simplefunc = psf( von_mises_hess_1_3< double* >() ); return;}
            if( name == "von_mises_hess_1_4") { simplefunc = psf( von_mises_hess_1_4< double* >() ); return;}
            if( name == "von_mises_hess_1_5") { simplefunc = psf( von_mises_hess_1_5< double* >() ); return;}
            if( name == "von_mises_hess_1_6") { simplefunc = psf( von_mises_hess_1_6< double* >() ); return;}
            if( name == "von_mises_hess_1_7") { simplefunc = psf( von_mises_hess_1_7< double* >() ); return;}
            if( name == "von_mises_hess_2_0") { simplefunc = psf( von_mises_hess_2_0< double* >() ); return;}
            if( name == "von_mises_hess_2_1") { simplefunc = psf( von_mises_hess_2_1< double* >() ); return;}
            if( name == "von_mises_hess_2_2") { simplefunc = psf( von_mises_hess_2_2< double* >() ); return;}
            if( name == "von_mises_hess_2_3") { simplefunc = psf( von_mises_hess_2_3< double* >() ); return;}
            if( name == "von_mises_hess_2_4") { simplefunc = psf( von_mises_hess_2_4< double* >() ); return;}
            if( name == "von_mises_hess_2_5") { simplefunc = psf( von_mises_hess_2_5< double* >() ); return;}
            if( name == "von_mises_hess_2_6") { simplefunc = psf( von_mises_hess_2_6< double* >() ); return;}
            if( name == "von_mises_hess_2_7") { simplefunc = psf( von_mises_hess_2_7< double* >() ); return;}
            if( name == "von_mises_hess_3_0") { simplefunc = psf( von_mises_hess_3_0< double* >() ); return;}
            if( name == "von_mises_hess_3_1") { simplefunc = psf( von_mises_hess_3_1< double* >() ); return;}
            if( name == "von_mises_hess_3_2") { simplefunc = psf( von_mises_hess_3_2< double* >() ); return;}
            if( name == "von_mises_hess_3_3") { simplefunc = psf( von_mises_hess_3_3< double* >() ); return;}
            if( name == "von_mises_hess_3_4") { simplefunc = psf( von_mises_hess_3_4< double* >() ); return;}
            if( name == "von_mises_hess_3_5") { simplefunc = psf( von_mises_hess_3_5< double* >() ); return;}
            if( name == "von_mises_hess_3_6") { simplefunc = psf( von_mises_hess_3_6< double* >() ); return;}
            if( name == "von_mises_hess_3_7") { simplefunc = psf( von_mises_hess_3_7< double* >() ); return;}
            if( name == "von_mises_hess_4_0") { simplefunc = psf( von_mises_hess_4_0< double* >() ); return;}
            if( name == "von_mises_hess_4_1") { simplefunc = psf( von_mises_hess_4_1< double* >() ); return;}
            if( name == "von_mises_hess_4_2") { simplefunc = psf( von_mises_hess_4_2< double* >() ); return;}
            if( name == "von_mises_hess_4_3") { simplefunc = psf( von_mises_hess_4_3< double* >() ); return;}
            if( name == "von_mises_hess_4_4") { simplefunc = psf( von_mises_hess_4_4< double* >() ); return;}
            if( name == "von_mises_hess_4_5") { simplefunc = psf( von_mises_hess_4_5< double* >() ); return;}
            if( name == "von_mises_hess_4_6") { simplefunc = psf( von_mises_hess_4_6< double* >() ); return;}
            if( name == "von_mises_hess_4_7") { simplefunc = psf( von_mises_hess_4_7< double* >() ); return;}
            if( name == "von_mises_hess_5_0") { simplefunc = psf( von_mises_hess_5_0< double* >() ); return;}
            if( name == "von_mises_hess_5_1") { simplefunc = psf( von_mises_hess_5_1< double* >() ); return;}
            if( name == "von_mises_hess_5_2") { simplefunc = psf( von_mises_hess_5_2< double* >() ); return;}
            if( name == "von_mises_hess_5_3") { simplefunc = psf( von_mises_hess_5_3< double* >() ); return;}
            if( name == "von_mises_hess_5_4") { simplefunc = psf( von_mises_hess_5_4< double* >() ); return;}
            if( name == "von_mises_hess_5_5") { simplefunc = psf( von_mises_hess_5_5< double* >() ); return;}
            if( name == "von_mises_hess_5_6") { simplefunc = psf( von_mises_hess_5_6< double* >() ); return;}
            if( name == "von_mises_hess_5_7") { simplefunc = psf( von_mises_hess_5_7< double* >() ); return;}
            if( name == "von_mises_hess_6_0") { simplefunc = psf( von_mises_hess_6_0< double* >() ); return;}
            if( name == "von_mises_hess_6_1") { simplefunc = psf( von_mises_hess_6_1< double* >() ); return;}
            if( name == "von_mises_hess_6_2") { simplefunc = psf( von_mises_hess_6_2< double* >() ); return;}
            if( name == "von_mises_hess_6_3") { simplefunc = psf( von_mises_hess_6_3< double* >() ); return;}
            if( name == "von_mises_hess_6_4") { simplefunc = psf( von_mises_hess_6_4< double* >() ); return;}
            if( name == "von_mises_hess_6_5") { simplefunc = psf( von_mises_hess_6_5< double* >() ); return;}
            if( name == "von_mises_hess_6_6") { simplefunc = psf( von_mises_hess_6_6< double* >() ); return;}
            if( name == "von_mises_hess_6_7") { simplefunc = psf( von_mises_hess_6_7< double* >() ); return;}
            if( name == "von_mises_hess_7_0") { simplefunc = psf( von_mises_hess_7_0< double* >() ); return;}
            if( name == "von_mises_hess_7_1") { simplefunc = psf( von_mises_hess_7_1< double* >() ); return;}
            if( name == "von_mises_hess_7_2") { simplefunc = psf( von_mises_hess_7_2< double* >() ); return;}
            if( name == "von_mises_hess_7_3") { simplefunc = psf( von_mises_hess_7_3< double* >() ); return;}
            if( name == "von_mises_hess_7_4") { simplefunc = psf( von_mises_hess_7_4< double* >() ); return;}
            if( name == "von_mises_hess_7_5") { simplefunc = psf( von_mises_hess_7_5< double* >() ); return;}
            if( name == "von_mises_hess_7_6") { simplefunc = psf( von_mises_hess_7_6< double* >() ); return;}
            if( name == "von_mises_hess_7_7") { simplefunc = psf( von_mises_hess_7_7< double* >() ); return;}
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
            else throw std::runtime_error( "PSimpleFunction< double*, double > " + name + " was not found in the PLibrary");
        }

        void PLibrary::checkout( std::string name, PFunction< std::vector<double>, double > &func)
        {
            typedef PFunction< std::vector<double>, double > pf;
            if( name == "neohook") { func = pf( neohook< std::vector<double> >() ); return;}
            if( name == "stvenkir") { func = pf( stvenkir< std::vector<double> >() ); return;}
            if( name == "linear_hardening") { func = pf( linear_hardening< std::vector<double> >() ); return;}
            if( name == "Cu_hardening") { func = pf( Cu_hardening< std::vector<double> >() ); return;}
            if( name == "von_mises") { func = pf( von_mises< std::vector<double> >() ); return;}
            if( name == "quadlog") { func = pf( quadlog< std::vector<double> >() ); return;}
            throw std::runtime_error( "PFunction< std::vector<double>, double > " + name + " was not found in the PLibrary");
        }
        void PLibrary::checkout( std::string name, PFunction< double*, double > &func)
        {
            typedef PFunction< double*, double > pf;
            if( name == "neohook") { func = pf( neohook< double* >() ); return;}
            if( name == "stvenkir") { func = pf( stvenkir< double* >() ); return;}
            if( name == "linear_hardening") { func = pf( linear_hardening< double* >() ); return;}
            if( name == "Cu_hardening") { func = pf( Cu_hardening< double* >() ); return;}
            if( name == "von_mises") { func = pf( von_mises< double* >() ); return;}
            if( name == "quadlog") { func = pf( quadlog< double* >() ); return;}
            throw std::runtime_error( "PFunction< double*, double > " + name + " was not found in the PLibrary");
        }




        void PLibrary::checkout( std::string name, PSimpleBase< std::vector<double>, double > *&simplefunc)
        {
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
            if( name == "linear_hardening_f") { simplefunc = new linear_hardening_f< std::vector<double> >(); return;}
            if( name == "linear_hardening_grad_0") { simplefunc = new linear_hardening_grad_0< std::vector<double> >(); return;}
            if( name == "linear_hardening_grad_1") { simplefunc = new linear_hardening_grad_1< std::vector<double> >(); return;}
            if( name == "linear_hardening_hess_0_0") { simplefunc = new linear_hardening_hess_0_0< std::vector<double> >(); return;}
            if( name == "linear_hardening_hess_0_1") { simplefunc = new linear_hardening_hess_0_1< std::vector<double> >(); return;}
            if( name == "linear_hardening_hess_1_0") { simplefunc = new linear_hardening_hess_1_0< std::vector<double> >(); return;}
            if( name == "linear_hardening_hess_1_1") { simplefunc = new linear_hardening_hess_1_1< std::vector<double> >(); return;}
            if( name == "Cu_hardening_f") { simplefunc = new Cu_hardening_f< std::vector<double> >(); return;}
            if( name == "Cu_hardening_grad_0") { simplefunc = new Cu_hardening_grad_0< std::vector<double> >(); return;}
            if( name == "Cu_hardening_grad_1") { simplefunc = new Cu_hardening_grad_1< std::vector<double> >(); return;}
            if( name == "Cu_hardening_hess_0_0") { simplefunc = new Cu_hardening_hess_0_0< std::vector<double> >(); return;}
            if( name == "Cu_hardening_hess_0_1") { simplefunc = new Cu_hardening_hess_0_1< std::vector<double> >(); return;}
            if( name == "Cu_hardening_hess_1_0") { simplefunc = new Cu_hardening_hess_1_0< std::vector<double> >(); return;}
            if( name == "Cu_hardening_hess_1_1") { simplefunc = new Cu_hardening_hess_1_1< std::vector<double> >(); return;}
            if( name == "von_mises_f") { simplefunc = new von_mises_f< std::vector<double> >(); return;}
            if( name == "von_mises_grad_0") { simplefunc = new von_mises_grad_0< std::vector<double> >(); return;}
            if( name == "von_mises_grad_1") { simplefunc = new von_mises_grad_1< std::vector<double> >(); return;}
            if( name == "von_mises_grad_2") { simplefunc = new von_mises_grad_2< std::vector<double> >(); return;}
            if( name == "von_mises_grad_3") { simplefunc = new von_mises_grad_3< std::vector<double> >(); return;}
            if( name == "von_mises_grad_4") { simplefunc = new von_mises_grad_4< std::vector<double> >(); return;}
            if( name == "von_mises_grad_5") { simplefunc = new von_mises_grad_5< std::vector<double> >(); return;}
            if( name == "von_mises_grad_6") { simplefunc = new von_mises_grad_6< std::vector<double> >(); return;}
            if( name == "von_mises_grad_7") { simplefunc = new von_mises_grad_7< std::vector<double> >(); return;}
            if( name == "von_mises_hess_0_0") { simplefunc = new von_mises_hess_0_0< std::vector<double> >(); return;}
            if( name == "von_mises_hess_0_1") { simplefunc = new von_mises_hess_0_1< std::vector<double> >(); return;}
            if( name == "von_mises_hess_0_2") { simplefunc = new von_mises_hess_0_2< std::vector<double> >(); return;}
            if( name == "von_mises_hess_0_3") { simplefunc = new von_mises_hess_0_3< std::vector<double> >(); return;}
            if( name == "von_mises_hess_0_4") { simplefunc = new von_mises_hess_0_4< std::vector<double> >(); return;}
            if( name == "von_mises_hess_0_5") { simplefunc = new von_mises_hess_0_5< std::vector<double> >(); return;}
            if( name == "von_mises_hess_0_6") { simplefunc = new von_mises_hess_0_6< std::vector<double> >(); return;}
            if( name == "von_mises_hess_0_7") { simplefunc = new von_mises_hess_0_7< std::vector<double> >(); return;}
            if( name == "von_mises_hess_1_0") { simplefunc = new von_mises_hess_1_0< std::vector<double> >(); return;}
            if( name == "von_mises_hess_1_1") { simplefunc = new von_mises_hess_1_1< std::vector<double> >(); return;}
            if( name == "von_mises_hess_1_2") { simplefunc = new von_mises_hess_1_2< std::vector<double> >(); return;}
            if( name == "von_mises_hess_1_3") { simplefunc = new von_mises_hess_1_3< std::vector<double> >(); return;}
            if( name == "von_mises_hess_1_4") { simplefunc = new von_mises_hess_1_4< std::vector<double> >(); return;}
            if( name == "von_mises_hess_1_5") { simplefunc = new von_mises_hess_1_5< std::vector<double> >(); return;}
            if( name == "von_mises_hess_1_6") { simplefunc = new von_mises_hess_1_6< std::vector<double> >(); return;}
            if( name == "von_mises_hess_1_7") { simplefunc = new von_mises_hess_1_7< std::vector<double> >(); return;}
            if( name == "von_mises_hess_2_0") { simplefunc = new von_mises_hess_2_0< std::vector<double> >(); return;}
            if( name == "von_mises_hess_2_1") { simplefunc = new von_mises_hess_2_1< std::vector<double> >(); return;}
            if( name == "von_mises_hess_2_2") { simplefunc = new von_mises_hess_2_2< std::vector<double> >(); return;}
            if( name == "von_mises_hess_2_3") { simplefunc = new von_mises_hess_2_3< std::vector<double> >(); return;}
            if( name == "von_mises_hess_2_4") { simplefunc = new von_mises_hess_2_4< std::vector<double> >(); return;}
            if( name == "von_mises_hess_2_5") { simplefunc = new von_mises_hess_2_5< std::vector<double> >(); return;}
            if( name == "von_mises_hess_2_6") { simplefunc = new von_mises_hess_2_6< std::vector<double> >(); return;}
            if( name == "von_mises_hess_2_7") { simplefunc = new von_mises_hess_2_7< std::vector<double> >(); return;}
            if( name == "von_mises_hess_3_0") { simplefunc = new von_mises_hess_3_0< std::vector<double> >(); return;}
            if( name == "von_mises_hess_3_1") { simplefunc = new von_mises_hess_3_1< std::vector<double> >(); return;}
            if( name == "von_mises_hess_3_2") { simplefunc = new von_mises_hess_3_2< std::vector<double> >(); return;}
            if( name == "von_mises_hess_3_3") { simplefunc = new von_mises_hess_3_3< std::vector<double> >(); return;}
            if( name == "von_mises_hess_3_4") { simplefunc = new von_mises_hess_3_4< std::vector<double> >(); return;}
            if( name == "von_mises_hess_3_5") { simplefunc = new von_mises_hess_3_5< std::vector<double> >(); return;}
            if( name == "von_mises_hess_3_6") { simplefunc = new von_mises_hess_3_6< std::vector<double> >(); return;}
            if( name == "von_mises_hess_3_7") { simplefunc = new von_mises_hess_3_7< std::vector<double> >(); return;}
            if( name == "von_mises_hess_4_0") { simplefunc = new von_mises_hess_4_0< std::vector<double> >(); return;}
            if( name == "von_mises_hess_4_1") { simplefunc = new von_mises_hess_4_1< std::vector<double> >(); return;}
            if( name == "von_mises_hess_4_2") { simplefunc = new von_mises_hess_4_2< std::vector<double> >(); return;}
            if( name == "von_mises_hess_4_3") { simplefunc = new von_mises_hess_4_3< std::vector<double> >(); return;}
            if( name == "von_mises_hess_4_4") { simplefunc = new von_mises_hess_4_4< std::vector<double> >(); return;}
            if( name == "von_mises_hess_4_5") { simplefunc = new von_mises_hess_4_5< std::vector<double> >(); return;}
            if( name == "von_mises_hess_4_6") { simplefunc = new von_mises_hess_4_6< std::vector<double> >(); return;}
            if( name == "von_mises_hess_4_7") { simplefunc = new von_mises_hess_4_7< std::vector<double> >(); return;}
            if( name == "von_mises_hess_5_0") { simplefunc = new von_mises_hess_5_0< std::vector<double> >(); return;}
            if( name == "von_mises_hess_5_1") { simplefunc = new von_mises_hess_5_1< std::vector<double> >(); return;}
            if( name == "von_mises_hess_5_2") { simplefunc = new von_mises_hess_5_2< std::vector<double> >(); return;}
            if( name == "von_mises_hess_5_3") { simplefunc = new von_mises_hess_5_3< std::vector<double> >(); return;}
            if( name == "von_mises_hess_5_4") { simplefunc = new von_mises_hess_5_4< std::vector<double> >(); return;}
            if( name == "von_mises_hess_5_5") { simplefunc = new von_mises_hess_5_5< std::vector<double> >(); return;}
            if( name == "von_mises_hess_5_6") { simplefunc = new von_mises_hess_5_6< std::vector<double> >(); return;}
            if( name == "von_mises_hess_5_7") { simplefunc = new von_mises_hess_5_7< std::vector<double> >(); return;}
            if( name == "von_mises_hess_6_0") { simplefunc = new von_mises_hess_6_0< std::vector<double> >(); return;}
            if( name == "von_mises_hess_6_1") { simplefunc = new von_mises_hess_6_1< std::vector<double> >(); return;}
            if( name == "von_mises_hess_6_2") { simplefunc = new von_mises_hess_6_2< std::vector<double> >(); return;}
            if( name == "von_mises_hess_6_3") { simplefunc = new von_mises_hess_6_3< std::vector<double> >(); return;}
            if( name == "von_mises_hess_6_4") { simplefunc = new von_mises_hess_6_4< std::vector<double> >(); return;}
            if( name == "von_mises_hess_6_5") { simplefunc = new von_mises_hess_6_5< std::vector<double> >(); return;}
            if( name == "von_mises_hess_6_6") { simplefunc = new von_mises_hess_6_6< std::vector<double> >(); return;}
            if( name == "von_mises_hess_6_7") { simplefunc = new von_mises_hess_6_7< std::vector<double> >(); return;}
            if( name == "von_mises_hess_7_0") { simplefunc = new von_mises_hess_7_0< std::vector<double> >(); return;}
            if( name == "von_mises_hess_7_1") { simplefunc = new von_mises_hess_7_1< std::vector<double> >(); return;}
            if( name == "von_mises_hess_7_2") { simplefunc = new von_mises_hess_7_2< std::vector<double> >(); return;}
            if( name == "von_mises_hess_7_3") { simplefunc = new von_mises_hess_7_3< std::vector<double> >(); return;}
            if( name == "von_mises_hess_7_4") { simplefunc = new von_mises_hess_7_4< std::vector<double> >(); return;}
            if( name == "von_mises_hess_7_5") { simplefunc = new von_mises_hess_7_5< std::vector<double> >(); return;}
            if( name == "von_mises_hess_7_6") { simplefunc = new von_mises_hess_7_6< std::vector<double> >(); return;}
            if( name == "von_mises_hess_7_7") { simplefunc = new von_mises_hess_7_7< std::vector<double> >(); return;}
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
            throw std::runtime_error( "PSimpleBase< std::vector<double>, double > " + name + " was not found in the PLibrary");
        }
        void PLibrary::checkout( std::string name, PSimpleBase< double*, double > *&simplefunc)
        {
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
            if( name == "linear_hardening_f") { simplefunc = new linear_hardening_f< double* >(); return;}
            if( name == "linear_hardening_grad_0") { simplefunc = new linear_hardening_grad_0< double* >(); return;}
            if( name == "linear_hardening_grad_1") { simplefunc = new linear_hardening_grad_1< double* >(); return;}
            if( name == "linear_hardening_hess_0_0") { simplefunc = new linear_hardening_hess_0_0< double* >(); return;}
            if( name == "linear_hardening_hess_0_1") { simplefunc = new linear_hardening_hess_0_1< double* >(); return;}
            if( name == "linear_hardening_hess_1_0") { simplefunc = new linear_hardening_hess_1_0< double* >(); return;}
            if( name == "linear_hardening_hess_1_1") { simplefunc = new linear_hardening_hess_1_1< double* >(); return;}
            if( name == "Cu_hardening_f") { simplefunc = new Cu_hardening_f< double* >(); return;}
            if( name == "Cu_hardening_grad_0") { simplefunc = new Cu_hardening_grad_0< double* >(); return;}
            if( name == "Cu_hardening_grad_1") { simplefunc = new Cu_hardening_grad_1< double* >(); return;}
            if( name == "Cu_hardening_hess_0_0") { simplefunc = new Cu_hardening_hess_0_0< double* >(); return;}
            if( name == "Cu_hardening_hess_0_1") { simplefunc = new Cu_hardening_hess_0_1< double* >(); return;}
            if( name == "Cu_hardening_hess_1_0") { simplefunc = new Cu_hardening_hess_1_0< double* >(); return;}
            if( name == "Cu_hardening_hess_1_1") { simplefunc = new Cu_hardening_hess_1_1< double* >(); return;}
            if( name == "von_mises_f") { simplefunc = new von_mises_f< double* >(); return;}
            if( name == "von_mises_grad_0") { simplefunc = new von_mises_grad_0< double* >(); return;}
            if( name == "von_mises_grad_1") { simplefunc = new von_mises_grad_1< double* >(); return;}
            if( name == "von_mises_grad_2") { simplefunc = new von_mises_grad_2< double* >(); return;}
            if( name == "von_mises_grad_3") { simplefunc = new von_mises_grad_3< double* >(); return;}
            if( name == "von_mises_grad_4") { simplefunc = new von_mises_grad_4< double* >(); return;}
            if( name == "von_mises_grad_5") { simplefunc = new von_mises_grad_5< double* >(); return;}
            if( name == "von_mises_grad_6") { simplefunc = new von_mises_grad_6< double* >(); return;}
            if( name == "von_mises_grad_7") { simplefunc = new von_mises_grad_7< double* >(); return;}
            if( name == "von_mises_hess_0_0") { simplefunc = new von_mises_hess_0_0< double* >(); return;}
            if( name == "von_mises_hess_0_1") { simplefunc = new von_mises_hess_0_1< double* >(); return;}
            if( name == "von_mises_hess_0_2") { simplefunc = new von_mises_hess_0_2< double* >(); return;}
            if( name == "von_mises_hess_0_3") { simplefunc = new von_mises_hess_0_3< double* >(); return;}
            if( name == "von_mises_hess_0_4") { simplefunc = new von_mises_hess_0_4< double* >(); return;}
            if( name == "von_mises_hess_0_5") { simplefunc = new von_mises_hess_0_5< double* >(); return;}
            if( name == "von_mises_hess_0_6") { simplefunc = new von_mises_hess_0_6< double* >(); return;}
            if( name == "von_mises_hess_0_7") { simplefunc = new von_mises_hess_0_7< double* >(); return;}
            if( name == "von_mises_hess_1_0") { simplefunc = new von_mises_hess_1_0< double* >(); return;}
            if( name == "von_mises_hess_1_1") { simplefunc = new von_mises_hess_1_1< double* >(); return;}
            if( name == "von_mises_hess_1_2") { simplefunc = new von_mises_hess_1_2< double* >(); return;}
            if( name == "von_mises_hess_1_3") { simplefunc = new von_mises_hess_1_3< double* >(); return;}
            if( name == "von_mises_hess_1_4") { simplefunc = new von_mises_hess_1_4< double* >(); return;}
            if( name == "von_mises_hess_1_5") { simplefunc = new von_mises_hess_1_5< double* >(); return;}
            if( name == "von_mises_hess_1_6") { simplefunc = new von_mises_hess_1_6< double* >(); return;}
            if( name == "von_mises_hess_1_7") { simplefunc = new von_mises_hess_1_7< double* >(); return;}
            if( name == "von_mises_hess_2_0") { simplefunc = new von_mises_hess_2_0< double* >(); return;}
            if( name == "von_mises_hess_2_1") { simplefunc = new von_mises_hess_2_1< double* >(); return;}
            if( name == "von_mises_hess_2_2") { simplefunc = new von_mises_hess_2_2< double* >(); return;}
            if( name == "von_mises_hess_2_3") { simplefunc = new von_mises_hess_2_3< double* >(); return;}
            if( name == "von_mises_hess_2_4") { simplefunc = new von_mises_hess_2_4< double* >(); return;}
            if( name == "von_mises_hess_2_5") { simplefunc = new von_mises_hess_2_5< double* >(); return;}
            if( name == "von_mises_hess_2_6") { simplefunc = new von_mises_hess_2_6< double* >(); return;}
            if( name == "von_mises_hess_2_7") { simplefunc = new von_mises_hess_2_7< double* >(); return;}
            if( name == "von_mises_hess_3_0") { simplefunc = new von_mises_hess_3_0< double* >(); return;}
            if( name == "von_mises_hess_3_1") { simplefunc = new von_mises_hess_3_1< double* >(); return;}
            if( name == "von_mises_hess_3_2") { simplefunc = new von_mises_hess_3_2< double* >(); return;}
            if( name == "von_mises_hess_3_3") { simplefunc = new von_mises_hess_3_3< double* >(); return;}
            if( name == "von_mises_hess_3_4") { simplefunc = new von_mises_hess_3_4< double* >(); return;}
            if( name == "von_mises_hess_3_5") { simplefunc = new von_mises_hess_3_5< double* >(); return;}
            if( name == "von_mises_hess_3_6") { simplefunc = new von_mises_hess_3_6< double* >(); return;}
            if( name == "von_mises_hess_3_7") { simplefunc = new von_mises_hess_3_7< double* >(); return;}
            if( name == "von_mises_hess_4_0") { simplefunc = new von_mises_hess_4_0< double* >(); return;}
            if( name == "von_mises_hess_4_1") { simplefunc = new von_mises_hess_4_1< double* >(); return;}
            if( name == "von_mises_hess_4_2") { simplefunc = new von_mises_hess_4_2< double* >(); return;}
            if( name == "von_mises_hess_4_3") { simplefunc = new von_mises_hess_4_3< double* >(); return;}
            if( name == "von_mises_hess_4_4") { simplefunc = new von_mises_hess_4_4< double* >(); return;}
            if( name == "von_mises_hess_4_5") { simplefunc = new von_mises_hess_4_5< double* >(); return;}
            if( name == "von_mises_hess_4_6") { simplefunc = new von_mises_hess_4_6< double* >(); return;}
            if( name == "von_mises_hess_4_7") { simplefunc = new von_mises_hess_4_7< double* >(); return;}
            if( name == "von_mises_hess_5_0") { simplefunc = new von_mises_hess_5_0< double* >(); return;}
            if( name == "von_mises_hess_5_1") { simplefunc = new von_mises_hess_5_1< double* >(); return;}
            if( name == "von_mises_hess_5_2") { simplefunc = new von_mises_hess_5_2< double* >(); return;}
            if( name == "von_mises_hess_5_3") { simplefunc = new von_mises_hess_5_3< double* >(); return;}
            if( name == "von_mises_hess_5_4") { simplefunc = new von_mises_hess_5_4< double* >(); return;}
            if( name == "von_mises_hess_5_5") { simplefunc = new von_mises_hess_5_5< double* >(); return;}
            if( name == "von_mises_hess_5_6") { simplefunc = new von_mises_hess_5_6< double* >(); return;}
            if( name == "von_mises_hess_5_7") { simplefunc = new von_mises_hess_5_7< double* >(); return;}
            if( name == "von_mises_hess_6_0") { simplefunc = new von_mises_hess_6_0< double* >(); return;}
            if( name == "von_mises_hess_6_1") { simplefunc = new von_mises_hess_6_1< double* >(); return;}
            if( name == "von_mises_hess_6_2") { simplefunc = new von_mises_hess_6_2< double* >(); return;}
            if( name == "von_mises_hess_6_3") { simplefunc = new von_mises_hess_6_3< double* >(); return;}
            if( name == "von_mises_hess_6_4") { simplefunc = new von_mises_hess_6_4< double* >(); return;}
            if( name == "von_mises_hess_6_5") { simplefunc = new von_mises_hess_6_5< double* >(); return;}
            if( name == "von_mises_hess_6_6") { simplefunc = new von_mises_hess_6_6< double* >(); return;}
            if( name == "von_mises_hess_6_7") { simplefunc = new von_mises_hess_6_7< double* >(); return;}
            if( name == "von_mises_hess_7_0") { simplefunc = new von_mises_hess_7_0< double* >(); return;}
            if( name == "von_mises_hess_7_1") { simplefunc = new von_mises_hess_7_1< double* >(); return;}
            if( name == "von_mises_hess_7_2") { simplefunc = new von_mises_hess_7_2< double* >(); return;}
            if( name == "von_mises_hess_7_3") { simplefunc = new von_mises_hess_7_3< double* >(); return;}
            if( name == "von_mises_hess_7_4") { simplefunc = new von_mises_hess_7_4< double* >(); return;}
            if( name == "von_mises_hess_7_5") { simplefunc = new von_mises_hess_7_5< double* >(); return;}
            if( name == "von_mises_hess_7_6") { simplefunc = new von_mises_hess_7_6< double* >(); return;}
            if( name == "von_mises_hess_7_7") { simplefunc = new von_mises_hess_7_7< double* >(); return;}
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
            throw std::runtime_error( "PSimpleBase< double*, double > " + name + " was not found in the PLibrary");
        }

        void PLibrary::checkout( std::string name, PFuncBase< std::vector<double>, double > *&func)
        {
            if( name == "neohook") { func = new neohook< std::vector<double> >(); return;}
            if( name == "stvenkir") { func = new stvenkir< std::vector<double> >(); return;}
            if( name == "linear_hardening") { func = new linear_hardening< std::vector<double> >(); return;}
            if( name == "Cu_hardening") { func = new Cu_hardening< std::vector<double> >(); return;}
            if( name == "von_mises") { func = new von_mises< std::vector<double> >(); return;}
            if( name == "quadlog") { func = new quadlog< std::vector<double> >(); return;}
            throw std::runtime_error( "PFuncBase< std::vector<double>, double > " + name + " was not found in the PLibrary");
        }
        void PLibrary::checkout( std::string name, PFuncBase< double*, double > *&func)
        {
            if( name == "neohook") { func = new neohook< double* >(); return;}
            if( name == "stvenkir") { func = new stvenkir< double* >(); return;}
            if( name == "linear_hardening") { func = new linear_hardening< double* >(); return;}
            if( name == "Cu_hardening") { func = new Cu_hardening< double* >(); return;}
            if( name == "von_mises") { func = new von_mises< double* >(); return;}
            if( name == "quadlog") { func = new quadlog< double* >(); return;}
            throw std::runtime_error( "PFuncBase< double*, double > " + name + " was not found in the PLibrary");
        }




}


#endif
