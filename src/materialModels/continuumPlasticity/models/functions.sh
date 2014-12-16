#!/bin/sh

fw --name von_mises -v beta1 beta2 beta3 tau_y q -d "First principle stress" "Second principle stress" "Third principle stress" "Yield stress" "Conjugate stress-like quantity" --sym "sqrt(pow(2./3.*beta1 - 1./3.*beta2 - 1./3.*beta3,2.) + pow(-1./3.*beta1 + 2./3.*beta2 - 1./3.*beta3,2.) + pow(-1./3.*beta1 - 1./3.*beta2 + 2./3.*beta3,2.)) - sqrt(2./3.)*(tau_y - q)" --grad --hess

fw --name quadlog -v lambda mu lambda1 lambda2 lambda3 K alpha -d "First Lame parameter" "Second Lame parameter" "First principle stretch" "Second principle stretch" "Third principle stretch" "Strain hardening coefficient" "Equivalent plastic strain" --sym "0.5*lambda*pow(log(lambda1) + log(lambda2) + log(lambda3),2.) + mu*(pow(log(lambda1),2.) + pow(log(lambda2),2.) + pow(log(lambda3),2.)) + 0.5*K*pow(alpha,2.)" --grad --hess

fw --name neohook -v lambda mu lambda1 lambda2 lambda3 K alpha -d "First Lame parameter" "Second Lame parameter" "First principle stretch" "Second principle stretch" "Third principle stretch" "Strain hardening coefficient" "Equivalent plastic strain" --sym "0.25*lambda*(pow(lambda1,2.)*pow(lambda2,2.)*pow(lambda3,2.) - 1.) - (0.5*lambda + mu)*log(lambda1*lambda2*lambda3) + 0.5*mu*(pow(lambda1,2.) + pow(lambda2,2.) + pow(lambda3,2.) - 3.) + 0.5*K*pow(alpha,2.)" --grad --hess

fw --name stvenkir -v lambda mu lambda1 lambda2 lambda3 K alpha -d "First Lame parameter" "Second Lame parameter" "First principle stretch" "Second principle stretch" "Third principle stretch" "Strain hardening coefficient" "Equivalent plastic strain" --sym "0.125*lambda*pow(pow(lambda1,2.) + pow(lambda2,2.) + pow(lambda3,2.) - 3.,2.) + 0.25*mu*(pow(lambda1,4.) + pow(lambda2,4.) + pow(lambda3,4.) - 2.*(pow(lambda1,2.) + pow(lambda2,2.) + pow(lambda3,2.)) + 3.) + 0.5*K*pow(alpha,2.)" --grad --hess

lw -h -d $PWD -v "std::vector<double>" "double*"

lw -c -d $PWD -v "std::vector<double>" "double*"
