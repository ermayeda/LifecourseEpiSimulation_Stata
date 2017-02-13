
set more off , perm
clear

/*Create blank data set*/
set obs $N //creates blank dataset with $N observations
gen id = _n  //creates id numbers


/******************************************************************************/
/*Step 1: Set parameters*/

*i. Specify time (in years) between cognitive assessments
local int_time = 1.5

*ii. Specify prevalence of binary exposure (exp)
local pexp = $pexp // 0.4 // prevalence of <high school education; reference: >=high school education group

*iii. Variances and correlations
local s2z0 =$s2z0       //variance of random cognitive intercept
local s2z1 =$s2z1       //variance of random cognitive slope
local sz01 =$sz01       //covariance of random intercept and random slope
local s2e = $s2e        //variance of unexplained variation in Cij
local rho = $rho        //correlation between noise terms for Cij
local s2d = $s2d        //variance of measurement error of Cij
*Notes for modifications:
*1. To eliminate random intercept and random slope terms, set
*sigma_square_zeta0(s2z0)=0, sigma_square_zeta1(s2z1)=0,
*and sigma_zeta_01(sz01)=0.
*2. To eliminate autoregressive covariance structure, set rho=0.
*3. To eliminate measurement error, set sigma_square_delta(s2d)=0.

*iv. Parameters for cognitive function (Cij)
* Commented out for when we run from the program.
local b00 = $b00        //cognitive intercept for unexposed
local b01 = $b01        //effect of exposure on cognitive intercept
local b02 = $b02        //effect of U1 on cognitive intercept
local b03 = $b03        //effect of U2 on cognitive intercept
local b10 = $b10        //cognitive slope for unexposed
local b11 = $b11        //effect of exposure on cognitive slope
local b12 = $b12        //effect of U1 on cognitive slope
local b13 = $b13        //effect of U2 on cognitive slope

*v. Parameters for death prior to Age K 
local g0 =$g0           //log odds of death by age K for reference group
local g1 =$g1           //log OR for effect of exposure on death
local g2 =$g2           //log OR for effect of U1 on death 
local g3 =$g3           //log OR for effect of U2 on death 
local g4 =$g4           //log OR for interaction between exposure and U1 on death

*vi. Parameter for U2
local a0 = $a0          //effect of education on U2
/******************************************************************************/


/******************************************************************************/
/*Step 2: Generate exposure variable with prevalence pexp by generating a U[0,1)
distribution and setting exposure=0 if the random number<pexp, else exposure=1*/
gen exposure = runiform()<`pexp'
cap label drop exposure
label define exposure 0 "0 >=high school" 1 "1 <high school"
label values exposure exposure
tab exposure
/******************************************************************************/


/******************************************************************************/
/*Step 3: Generate U1, a continuous time-constant variable, U~N(0,1)*/
* Effect of U2 on exposure, a0
gen U2 = 0 + `a0' *exposure + rnormal()*1
/******************************************************************************/


/******************************************************************************/
/*Step 4: Generate U2, a continuous time-constant variable, U~N(0,1)*/
gen U1 = rnormal()

standardize U1, replace
standardize U2, replace
/******************************************************************************/


/*Step 5: Generate cognitive function values at each cognitive assessment. The study
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
/******************************************************************************/


/******************************************************************************/
/*Step 6: Generate death by age K
gen p_deathK= ///
   exp(`g0'+`g1'*exposure + `g2'*U1+ `g3'*U2 + `g4'*U1*exposure) ///
      / (1+exp(`g0'+`g1'*exposure + `g2'*U1+ `g3'*U2 + `g4'*U1*exposure) )
label var p_deathK "Probability of dying by age K"

gen deathK = runiform() < p_deathK
label var survK "Died by age K"
gen survK = 1-deathK
label var survK "Survived to age K"
tab survK
cap label drop survK
label define survK 0 "0 Died" 1 "1 Alive"
label values survK survK
/******************************************************************************/



/******************************************************************************/
/*** End data generation ***/
/******************************************************************************/


/******************************************************************************/
/*** Analyze data ***/
*The following steps include restricting data set to:
a) people alive at age K
b) random sample of n=$sample people alive at age K
/******************************************************************************/

/******************************************/
*pull int time
scalar int_time = `int_time'
/******************************************/

/******************************************/
*pull total N and n in each exposure group at age 20
sum exposure
scalar N=r(N)
scalar N_exposure1=r(mean)*r(N)
scalar N_exposure0=r(N)-r(mean)*r(N)

*pull proportion of cohort alive at age K
sum survK, meanonly
scalar p_survK = r(mean)
/******************************************/

/******************************************/
*before restricting sample to survivors to age K, estimate "true" data-generating model--sanity check
*xtmixed measured_cogfxn i.exposure##c.time c.time#c.U1 c.U1 c.time#c.U2 c.U2 || id: time, mle cov(uns)
/******************************************/

/******************************************/
*restrict data set to people alive at age K
drop if deathK==1
*pull total N and n in each exposure group alive at age K
sum exposure
scalar N_ageK=r(N)
scalar N_exposure1_ageK=r(mean)*r(N)
scalar N_exposure0_ageK=r(N)-r(mean)*r(N)
/******************************************/

/******************************************************************************/
/*restrict to random sample of n=$sample among those alive at age K
sample $sample, count
/******************************************************************************/

/******************************************/
*pull total N (should be n=$sample) and n in each exposure group alive at age K
*all following steps use this sample
summarize exposure
scalar N_sample=r(N)
scalar p_exposure1_sample = r(mean)
scalar N_exposure0_sample=r(N)-r(mean)*r(N)
/******************************************/

/******************************************/
*summarize mean U1 and U2 by exposure group
/******************************************/
*exposure=1 group
summarize U1 if exposure==1, meanonly 
scalar mean_U1_exposure1 = r(mean)
summarize U2 if exposure==1, meanonly
scalar mean_U2_exposure1 = r(mean)

*exposure=0 group
summarize U1 if exposure==0, meanonly 
scalar mean_U1_exposure0 = r(mean)
summarize U2 if exposure==0, meanonly
scalar mean_U2_exposure0 = r(mean)


/*Data generating model for cognitive function:
Cij = b00 + b01*exposurei + b02*age0i + b03*Ui + (b10 +
      + b11*exposurei + b12*age0i + b13*Ui)*timeij + z0i + z1i*timeij + epsilonij*/

/******************************************/
/*Linear mixed effects model with random intercepts and slopes allowing for possible 
correlation of the random intercept and slope and no additional within-person 
covariance structure)*/
reshape long measured_cogfxn time, i(id) j(year)
qui xtmixed measured_cogfxn i.exposure##c.time || id: time, mle cov(uns) , if survK==1

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



