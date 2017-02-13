
set more off , perm
clear

/*Create blank data set*/
set obs 1500 //creates blank dataset with 1,500 observations
gen id = _n  //creates id numbers


/*Step 1: Set parameters*/

*i. Specify time (in years) between cognitive assessments
local int_time = 1.5

*ii. Specify prevalence of binary exposure (exp)
local pexp = $pexp // 0.4 // prevalence of graduating HS. reference: high education group


*iii. Variances and correlations
local s2z0 =$s2z0 //variance of random cognitive intercept
local s2z1 =$s2z1   //variance of random cognitive slope
local sz01 =$sz01    //covariance of random intercept and random slope
local s2e = $s2e  //variance of unexplained variation in Cij
local rho = $rho  //correlation between noise terms for Cij
local s2d = $s2d  //variance of measurement error of Cij
*Notes for modifications:
*1. To eliminate random intercept and random slope terms, set
*sigma_square_zeta0(s2z0)=0, sigma_square_zeta1(s2z1)=0,
*and sigma_zeta_01(sz01)=0.
*2. To eliminate autoregressive covariance structure, set rho=0.
*3. To eliminate measurement error, set sigma_square_delta(s2d)=0.

*iv. Parameters for cognitive function (Cij)
* Commented out for when we run from the program.
local b00 = $b00 //cognitive intercept for unexposed
local b01 = $b01   //effect of exposure on cognitive intercept
local b02 = $b02    //effect of u1 on cognitive intercept
local b03 = $b03 //effect of U2 on cognitive intercept
local b10 = $b10 //cognitive slope for unexposed
local b11 = $b11    //effect of exposure on cognitive slope
local b12 = $b12   //effect of u1 on cognitive slope
local b13 = $b13    //effect of U2 on cognitive slope



* New: model the survival to age 65
local g0 =$g0 // ln(5) // overall odds of survival in the reference category (low educ group). 81%
         *   * odds = p/(1-p)=.81/.19=4
         *   * .83*.98. 83% of newborns survive to age 65. US Life Tables 2006. 98% of newborns survive to age 20
local g1 =$g1 // -ln(1.5) // log OR of educ-->survival
local g2 =$g2  //-ln(1.5) // log OR of u1 on survival (pathology)
local g3 =$g3 // -ln(1.3) // log OR of u2 on survival (genetics)
local g4 =$g4 //          // log OR of intx between educ and U1 on survival

* Effect of education to pathology (U1)
local h0 = $h0

/*Step 3: Generate exposure variable with prevalence pexp by generating a U[0,1)
distribution and setting exposure=0 if the random number<pexp, else exposure=1*/
gen exposure = runiform()<`pexp'
cap label drop exposure
label define exposure 0 "0 High education" 1 "1 Low education"
label values exposure exposure
tab exposure

/*Step 4: Generate U1, a continuous time-constant variable, U~N(0,1) pathology variable*/
* Effect of education to pathology, h0
gen U2 = 0 + `h0' *exposure + rnormal()*1
*gen U2 = rnormal()*1


/*Step 4x: Generate U2, a continuous time-constant variable, U~N(0,1) genetic polygenic risk score*/
*gen U1 = rnormal()
gen U1 = rnormal()

standardize U1, replace
standardize U2, replace

/*Step 5: Create visit times and ages for each cognitive assessment. The study
will include 7 cognitive assessment waves "int_time" years apart*/
*Generate visit times
gen time0 = 0
gen time1 = `int_time'
forvalues i=2/6 {
   local j=`i'-1
   gen time`i' = time`j' + `int_time'
}

/*Generate random terms for slope and intercept, zeta_0i (z0i) and zeta_1i (z1i),
where zeta_0i and zeta_1i covary*/
matrix c = (`s2z0', `sz01'\ `sz01', `s2z1')
drawnorm z0i z1i, cov(c)

/*Generate autoregressive noise term (unexplained variance in Cij) for each visit*/
gen sd_alpha = sqrt((1-`rho'*`rho')*`s2e')

gen epsilon0 = sqrt(`s2e')*rnormal()
forvalues i=1/6 {
   local j=`i'-1
   gen alpha`i' = sd_alpha*rnormal()
   gen epsilon`i' = `rho'*epsilon`j' + alpha`i'
}

/*Generate true cognitive function*/
gen true_cogfxn0 = `b00' +`b01'*exposure +`b02'*U1 +`b03'*U2 ///
   +(`b10' +`b11'*exposure+`b12'*U1 +`b13'*U2)*0 ///
   +z0i +z1i*0 +epsilon0

forvalues i=1/6 {
   gen true_cogfxn`i' = `b00' +`b01'*exposure+`b02'*U1 +`b03'*U2 ///
      +(`b10' +`b11'*exposure+`b12'*U1+`b13'*U2)*time`i' ///
      +z0i +z1i*time`i' +epsilon`i'
}

/*Generate measured cognitive function (true cognitive function measured with error)*/
*C_ij = C_ij + error_ij
*First, generate delta term (delta_ij = measurement error for Cij) for each visit
forvalues i=0/6 {
   gen delta`i' = sqrt(`s2d')*rnormal()
}

*Next, generate measured cognitive function as C_ij* = C_ij + delta_ij
forvalues i=0/6 {
   gen measured_cogfxn`i' = true_cogfxn`i' +delta`i'
}



* New 6/6/16 Generate survival to age 65
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




/******************************************************************************/
/***  End data generation                          ***/
/******************************************************************************/


/******************************************************************************/
/***  MODELS                                ***/
/******************************************************************************/

/******************************************/
*pull N
scalar N= _N
/******************************************/

/******************************************/
*pull int time
scalar int_time = `int_time'
/******************************************/

/******************************************/
*pull percentage of deaths
summarize exposure , meanonly
scalar p_low_educ = r(mean)
summarize surv65, meanonly
scalar p_surv65 = r(mean)
/******************************************/


/*Data generating model for cognitive function:
Cij = b00 + b01*exposurei + b02*age0i + b03*Ui + (b10 +
      + b11*exposurei + b12*age0i + b13*Ui)*timeij + z0i + z1i*timeij + epsilonij*/

/******************************************/
/*Conventional linear mixed effects model (linear mixed effects models with
random intercepts and slopes allowing for possible correlation of the random
intercept and slope and no additional within-person covariance structure)*/
reshape long measured_cogfxn time, i(id) j(year)
qui xtmixed measured_cogfxn i.exposure##c.time || id: time, mle cov(uns) , if surv65==1

* "true" model except U1 is conditional on education.
* just for sanity checking!
*xtmixed measured_cogfxn i.exposure##c.time c.time#c.U1 c.U1 c.time#c.U2 c.U2 || id: time, mle cov(uns)

*store fixed effect estimates as scalars;
matrix var=e(V)
*matrix list var
scalar mixed_b00_int = _b[_cons]
scalar mixed_b00_int_se = _se[_cons]
scalar mixed_b01_exposure = _b[1.exposure]
scalar mixed_b01_exposure_se = _se[1.exposure]
scalar mixed_b10_time = _b[time]
scalar mixed_b10_time_se = _se[time]
scalar mixed_b11_exptime = _b[1.exposure#c.time]
scalar mixed_b11_exptime_se = _se[1.exposure#c.time]
/******************************************/
*store random effect estimates as scalars;
matrix var=e(b)
_diparm lns1_1_2, f(exp(@)) d(exp(@))
scalar mixed_z0i_est = r(est)
scalar mixed_z0i_se = r(se)
_diparm lns1_1_1, f(exp(@)) d(exp(@))
scalar mixed_z1i_est = r(est)
scalar mixed_z1i_se = r(se)
_diparm atr1_1_1_2, f(tanh(@)) d(1-tanh(@)^2)
scalar mixed_rho_est = r(est)
scalar mixed_rho_se = r(se)
/******************************************/

