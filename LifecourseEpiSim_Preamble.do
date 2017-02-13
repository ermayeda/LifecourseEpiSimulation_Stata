* sim-001-preambling.do
* ALG ERM YG TF
* 20 Jun 2016
* -------------------

* paramcoefficients.do
global pexp = 0.4 // prevalence of graduating HS. reference: high education group
global seeed = 8675309
set seed $seeed
* new decision
* U1: genetic risk
* U2: pathology

* Effect of education to pathology (U2)
global h0 = .5

*iii. Variances and correlations
global s2z0 = 0.2     //variance of random cognitive intercept
global s2z1 = 0.005   //variance of random cognitive slope
global sz01 = -0.03   //covariance of random intercept and random slope
global s2e = 0.70     //variance of unexplained variation in Cij
global rho = 0.40     //correlation between noise terms for Cij
global s2d = 0.3     //variance of measurement error of Cij - .19

/*
scenarios 7/13/16
   1 unbiased situation
   2 u1 pathology affects survival. has less variability.
   3 u1 pathology education intx?
   4 u2 genetics on survival.
u2's variance is larger than u1.
*/


if `paramset'==1 {

   *parameters for Cij
   global b00 = 0        //cognitive intercept for unexposed
   global b01 = -.2     //effect of exposure on cognitive intercept
   global b02 = -0.05   //effect of u1 on cognitive intercept
   global b03 = -.2    //effect of U2 on cognitive intercept
   global b10 = -0.05    //cognitive slope for unexposed (high educ group)
   global b11 = -0.05    //direct effect of exposure on cognitive slope
   global b12 = -0.05   //effect of u1 on cognitive slope
   global b13 = -0.05   //effect of U2 on cognitive slope

   * New: model the survival to age 65
   *global g0 = ln(4) // overall odds of survival in the reference category (low educ group). 81%
   *global g0 = ln(4) - $pexp*$g1 // because we want to adjust g0 for the additional people who died due to low educ.
            * odds = p/(1-p)=.81/.19=4
            * .83*.98. 83% of newborns survive to age 65. US Life Tables 2006. 98% of newborns survive to age 20
   global g1 = -ln(2) // log OR of educ-->survival
   global g2 = -ln(1) // log OR of u1 on survival (pathology)
   global g3 = -ln(1) // log OR of u2 on survival (genetics)
   global g4 = -ln(1) // log OR of intx between educ and U1 on survival
   if $causalcondition==1 { // no bias
      di "Hi!"
      *global b11 = $b11 + $h0*$b12
      di -0.1+.3*-0.05
   }
   if $causalcondition==2 { // u1 affects survival
      global g2 = -ln(2) // log OR of u1 on survival (pathology)
   }
   else {
      global g2 = -ln(1) // log OR of u1 on survival (pathology)
   }
   if $causalcondition==3 { // it's like a silencing interaction when we keep g2==0. U2 only has an effect in the exposed.
      global g4 = -ln(2)   // log OR of intx between educ and U1 on survival
   }
   else { // it's like a silencing interaction when we keep g2==0. U2 only has an effect in the exposed.
      global g4 = -ln(1)   // log OR of intx between educ and U1 on survival
   }
   if $causalcondition==4 {
      global g3 = -ln(2) // log OR of u2 on survival (genetics)
   }
   else {
      global g3 = -ln(1) // log OR of u2 on survival (genetics)
   }
   if $causalcondition==5 {
      * to comment this out is like having a silencing intx.
      *global g2 = -ln(1.5) // log OR of u1 on survival (pathology)
      global g3 = -ln(2) // log OR of u2 on survival (genetics)
      global g4 = -ln(2)   // log OR of intx between educ and U1 on survival
   }
   else {
      global g3 = -ln(1) // log OR of u2 on survival (genetics)
      global g4 = -ln(1)   // log OR of intx between educ and U1 on survival
   }
}

if `paramset'==2 {

   *parameters for Cij
   global b00 = 0        //cognitive intercept for unexposed
   global b01 = -.2     //effect of exposure on cognitive intercept
   global b02 = -0.05   //effect of u1 on cognitive intercept
   global b03 = -.2    //effect of U2 on cognitive intercept
   global b10 = -0.05    //cognitive slope for unexposed (high educ group)
   global b11 = -0.1    //direct effect of exposure on cognitive slope
   global b12 = -0.1   //effect of u1 on cognitive slope
   global b13 = -0.1   //effect of U2 on cognitive slope

   * New: model the survival to age 65
   *global g0 = ln(4) // overall odds of survival in the reference category (low educ group). 81%
   *global g0 = ln(4) - $pexp*$g1 // because we want to adjust g0 for the additional people who died due to low educ.
            * odds = p/(1-p)=.81/.19=4
            * .83*.98. 83% of newborns survive to age 65. US Life Tables 2006. 98% of newborns survive to age 20
   global g1 = -ln(6) // log OR of educ-->survival
   global g2 = -ln(1) // log OR of u1 on survival (pathology)
   global g3 = -ln(1) // log OR of u2 on survival (genetics)
   global g4 = -ln(1)   // log OR of intx between educ and U1 on survival
   if $causalcondition==1 { // no bias
      di "Hi!"
      *global b11 = $b11 + $h0*$b12
      di -0.1+.3*-0.05
   }
   if $causalcondition==2 { // u1 affects survival
      global g2 = -ln(6) // log OR of u1 on survival (pathology)
   }
   else {
      global g2 = -ln(1) // log OR of u1 on survival (pathology)
   }
   if $causalcondition==3 { // it's like a silencing interaction when we keep g2==0. U2 only has an effect in the exposed.
      global g4 = -ln(6)   // log OR of intx between educ and U1 on survival
   }
   else { // it's like a silencing interaction when we keep g2==0. U2 only has an effect in the exposed.
      global g4 = -ln(1)   // log OR of intx between educ and U1 on survival
   }
   if $causalcondition==4 {
      global g3 = -ln(6) // log OR of u2 on survival (genetics)
   }
   else {
      global g3 = -ln(1) // log OR of u2 on survival (genetics)
   }
   if $causalcondition==5 {
      * to comment this out is like having a silencing intx.
      *global g2 = -ln(1.5) // log OR of u1 on survival (pathology)
      global g3 = -ln(6) // log OR of u2 on survival (genetics)
      global g4 = -ln(6)   // log OR of intx between educ and U1 on survival
   }
   else {
      global g3 = -ln(1) // log OR of u2 on survival (genetics)
      global g4 = -ln(1)   // log OR of intx between educ and U1 on survival
   }
}



* challenge: we want the same number of people living.
                  // bin   cont con  cont
*global g0 = ln(4) - $pexp*$g1 - $g2 - $pexp*$h0*$g3 - $g4


* New 6/9/16
*this is a do file to run the data generating do file with different values for the baseline mortality odds
* start with a value for the target that is definitely toolow,
* then it should notify you at the first increment where you are no longer
* too low
clear
local toolow =1
local target_p_surv65 =.8
*add lower bound guess here
local gg0 =1

forvalues i=1/100 {
   clear
   local toolow =1
   local target_p_surv65 =.8
   *add lower bound guess here
   local gg0 =1
   quietly forvalues x = 0(.01)5 {
      if `toolow'==1 {
         *set seed 8675309
         local gg0t =`gg0'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "ggo = `gg0t'"
         clear
         set obs 1500 //creates blank dataset with 1,500 observations
         gen exposure = runiform()< $pexp
         gen U2 = 0 + $h0 *exposure + rnormal()*.1
         gen U1 = rnormal()

         gen p_surv65= ///
            exp(`gg0t'+$g1*exposure + $g2*U1+ $g3*U2 + $g4*U1*exposure) ///
               / (1+exp(`gg0t'+$g1*exposure + $g2*U1+ $g3*U2 + $g4*U1*exposure) )
         label var p_surv65 "Probability of surviving to 65"
         gen surv65 = runiform() < p_surv65
         label var surv65 "Survived to 65"
         *label define surv65 0 "0 Died" 1 "1 Alive"
         sum surv65 , meanonly
         local mean =`r(mean)'
         if `r(mean)' > `target_p_surv65'+.001 {
            local toolow=0
            noisily di in red "I stopped at g0==`gg0t', surv65=`r(mean)'
         }
         else {
            noisily di in white "I did NOT stopp at g0==`gg0t', surv65=`r(mean)'
         }
      }
   }
   local g0`i' = `gg0t'
}
clear
set obs 100
gen g0=.
forvalues i=1/100 {
   replace g0=`g0`i'' in `i'
}
sum g0
global g0 = `r(mean)'












/*

$pexp = 0.6
1-$pexp=0.4
p0+p1=0.8

p0=exp($go)/(1+exp($g0))



/*
As g1 g2 g3 g4 go up, we want g0 to go down.
By how much?

gen p_surv65= ///
   exp(`g0'+`g1'*exposure + `g2'*U1+ `g3'*U2 + `g4'*U1*exposure) ///
      / (1+exp(`g0'+`g1'*exposure + `g2'*U1+ `g3'*U2 + `g4'*U1*exposure) )
label var p_surv65 "Probability of surviving to 65"

gen surv65 = runiform() < p_surv65
label var surv65 "Survived to 65"
tab surv65
cap label drop surv65
label define surv65 0 "0 Died" 1 "1 Alive"
label values surv65 surv65

*/

* http://www.stata.com/support/faqs/programming/system-of-nonlinear-equations/
discard
cap program drop nlfaq
program nlfaq
        syntax varlist(min=1 max=1) [if], at(name)
        tempname g0 p
        scalar `g0' = `at'[1, 1]
        scalar `p'  = `at'[1, 2]
        tempvar yh
        gen double `yh' = (1-$pexp)*exp(`g0')/(1+exp(`g0')) + $pexp*exp(`g0'+$g1)/(1+exp(`g0'+$g1)) -`p' in 1
        replace    `yh' = `p'+`g0'-1                                                                     in 2
        replace `varlist' = `yh'
end
clear
set obs 2
gen y=0
replace y=1 in 2
nl faq @ y , parameters( g0 b ) initials( g0 1 b 2 )


discard
cap program drop nlfaq
 program nlfaq
        syntax varlist(min=1 max=1) [if], at(name)
        tempname A B
        scalar `A' = `at'[1, 1]
        scalar `B' = `at'[1, 2]
        tempvar yh
        gen double `yh' = $pexp*exp(`A') + `B' - 2-1  in 1
        replace `yh' = `A'/`B' - log(`B')       in 2
        replace `varlist' = `yh'
 end
clear
set obs 2
generate y = 0
replace y = 1 in 1
nl faq @ y, parameters(A B ) initial(A 1 B 1 )
*          /A |    .848199          .        .       .            .           .
*          /B |   1.664563          .        .       .            .           .
di exp(.848199)+1.664563-2          -1
di .848199/1.664563 - log(1.664563)+1




* http://www.stata.com/support/faqs/programming/system-of-nonlinear-equations/
clear mata
mata:
void mysolver(todo,p,lnf,S,H)
   {
            g0  = p[1]
            lnf = (1-$pexp)*exp(g0)/(1+exp(g0))+$pexp*exp(g0+$g1)/(1+exp(g0+$g1))-0.8
   }
S=optimize_init()
optimize_init_evaluator(S, &mysolver())
optimize_init_evaluatortype(S, "v0")
optimize_init_params(S, (1))
optimize_init_which(S,  "min" )
optimize_init_tracelevel(S,"none")
optimize_init_conv_ptol(S, 1e-16)
optimize_init_conv_vtol(S, 1e-16)
p = optimize(S)
p
end

*odds = exp(odds) / (1+exp(odds))
*di 1- exp(1) / (1+exp(1))

di (1-$pexp)*exp(-50.56272752)/(1+exp(-50.56272752)) + $pexp*exp(-50.56272752+$g1)/(1+exp(-50.56272752+$g1)) -0.8

*      a=p[1]
*      lnf = (1-$pexp)*exp(`g0')/(1+exp(`g0')) + $pexp*exp(`g0'+$g1)/(1+exp(`g0'+$g1)) -`b'


* for condition1
(1-$pexp)*exp($g0)/(1+exp($g0)) + $pexp*exp($g0+$g1)/(1+exp($g0+$g1)) -0.8 = 0
Solve for $g0

* for condition2
(1-$pexp)*exp($g0)/(1+exp($g0)) + $pexp*exp($g0+$g1)/(1+exp($g0+$g1)) -0.8 = 0

*/

