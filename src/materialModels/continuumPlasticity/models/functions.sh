#!/bin/sh

fw --name "von_mises" -v "beta1" "beta2" "beta3" "tau_y" "q" "xi1" "xi2" "xi3" -d "First principle stress" "Second principle stress" "Third principle stress" "Yield stress" "Conjugate stress-like quantity" "First principle deviatoric back stress" "Second principle deviatoric back stress" "Third principle deviatoric back stress" --sym "sqrt(pow(2./3.*beta1 - 1./3.*beta2 - 1./3.*beta3 - xi1,2.) + pow(-1./3.*beta1 + 2./3.*beta2 - 1./3.*beta3 - xi2,2.) + pow(-1./3.*beta1 - 1./3.*beta2 + 2./3.*beta3 - xi3,2.)) - sqrt(2./3.)*(tau_y - q)" --grad --hess

fw --name "linear_hardening" -v "alpha" "K" -d "Equivalent plastic strain" "Hardening parameter" --sym "-K*alpha" --grad --hess

fw --name "Cu_hardening" -v "alpha" "K" -d "Equivalent plastic strain" "Hardening parameter" --sym "-1.1216965746072850e+09*pow(alpha,0.5)+3.2738656842295301e+08*alpha+3.2986102779887700e+08*pow(alpha,1./3.)" --grad --hess

fw --name "quadlog" -v "lambda" "mu" "lambda1" "lambda2" "lambda3" "K" "alpha" -d "First Lame parameter" "Second Lame parameter" "First principle stretch" "Second principle stretch" "Third principle stretch" "Strain hardening coefficient" "Equivalent plastic strain" --sym "0.5*lambda*pow(log(lambda1) + log(lambda2) + log(lambda3),2.) + mu*(pow(log(lambda1),2.) + pow(log(lambda2),2.) + pow(log(lambda3),2.)) + 0.5*K*pow(alpha,2.)" --grad --hess

fw --name "neohook" -v "lambda" "mu" "lambda1" "lambda2" "lambda3" "K" "alpha" -d "First Lame parameter" "Second Lame parameter" "First principle stretch" "Second principle stretch" "Third principle stretch" "Strain hardening coefficient" "Equivalent plastic strain" --sym "0.25*lambda*(pow(lambda1,2.)*pow(lambda2,2.)*pow(lambda3,2.) - 1.) - (0.5*lambda + mu)*log(lambda1*lambda2*lambda3) + 0.5*mu*(pow(lambda1,2.) + pow(lambda2,2.) + pow(lambda3,2.) - 3.) + 0.5*K*pow(alpha,2.)" --grad --hess

fw --name "stvenkir" -v "lambda" "mu" "lambda1" "lambda2" "lambda3" "K" "alpha" -d "First Lame parameter" "Second Lame parameter" "First principle stretch" "Second principle stretch" "Third principle stretch" "Strain hardening coefficient" "Equivalent plastic strain" --sym "0.125*lambda*pow(pow(lambda1,2.) + pow(lambda2,2.) + pow(lambda3,2.) - 3.,2.) + 0.25*mu*(pow(lambda1,4.) + pow(lambda2,4.) + pow(lambda3,4.) - 2.*(pow(lambda1,2.) + pow(lambda2,2.) + pow(lambda3,2.)) + 3.) + 0.5*K*pow(alpha,2.)" --grad --hess

lw -h -d $PWD -v "std::vector<double>"

lw -c -d $PWD -v "std::vector<double>"
