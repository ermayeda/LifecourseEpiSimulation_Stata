set more off
set seed $seeed



global betalist "mixed_b00_int mixed_b10_time mixed_b01_exposure mixed_b11_exptime"
*local B =11 // number of simulations
*global B=4
*local N=1500


local simline ""
foreach x in p_low_educ p_surv65 N int_time { // dif_in_slope
   local simline "`simline' `x'=`x'"
}
foreach x in ///
   $betalist {
   local simline "`simline' `x'=`x' "
   local simline "`simline' `x'_se=`x'_se "
}

simulate `simline', ///
   reps($B) seed(8675309): do $analysis/sim-306-models-stata-sim

*number of individuals
local N = N
local int_time = int_time

*proportion died and avg f/u time
foreach x in p_low_educ p_surv65 { // dif_in_slope
   summarize `x', meanonly
   scalar mean_`x' = r(mean)
   scalar r_mean_`x' = round(mean_`x',0.01)
}
scalar pct_surv65 = 100*r_mean_p_surv65

*Mean
foreach b in $betalist {
   summarize `b', meanonly
   scalar mean_`b' = r(mean)
   local mean_`b' : di %5.3f r(mean)

   *SE
   qui summarize `b'_se, meanonly
   scalar mean_`b'_se = r(mean)
   local mean_`b'_se : di %5.3f r(mean)
   *RMSE
   *Step 1: Calculate empirical SE (SD of estimates of interest from all simulations);
   summarize `b'
   scalar empSE_`b' = r(sd)
   local empSE_`b' : di %5.3f r(sd)
   *Step 2: RMSE
   local foo= substr("`b'",7,3)
   local goo = ${`foo'}
   if "`b'"== "mixed_b11_exptime"  {
      local goo = $b11 + $h0*$b13
   }
   if "`b'"== "mixed_b01_exposure" {
      local goo = $b01 + $h0*$b03
   }
   scalar RMSE_`b' = sqrt((mean_`b'-`goo')*(mean_`b'-`goo')+empSE_`b'*empSE_`b')
   local RMSE_`b' : di %5.3f sqrt((mean_`b'-`goo')*(mean_`b'-`goo')+empSE_`b'*empSE_`b')
   di in red "`goo'"

   *95% CI coverage: in order to calc vars for coverage, need to calc LCI and UCI
   cap drop lb_`b'
   gen lb_`b' = `b'- 1.96*`b'_se
   cap drop ub_`b'
   gen ub_`b' = `b'+ 1.96*`b'_se
   cap drop coverage_`b'
   gen coverage_`b'= (ub_`b' > `goo' & lb_`b' < `goo')
   *Create summary variable for coverage for each parameter
   summarize coverage_`b', meanonly
   scalar Pcoverage_`b' = r(mean)
   local Pcoverage_`b' : di %5.3f r(mean)
   *Z test statistic and p-value for rejecting true value of each parameter
   scalar Z_`b' = ((mean_`b' - `goo')/(empSE_`b'/sqrt($B)))
   local Z_`b' = ((mean_`b' - `goo')/(empSE_`b'/sqrt($B)))
   local Z_`b' : di %5.3f `Z_`b''
   di in red "((mean_`b' - `goo' )/(empSE_`b'/sqrt($B)))"
   di in white "Z_`b':"
   scalar pvalue_`b'=2*(1-normal(abs(Z_`b')))
   local pvalue_`b': di %5.3f  2*(1-normal(abs(Z_`b')))
   * Percent Bias
   if ${`foo'} !=0 {
      scalar pb_`b' = 100*round(((mean_`b' - `goo' ) / `goo' ),.1)
      local pb_`b' : di %5.3f 100*round(((mean_`b' - `goo' ) / `goo' ),.1)
   }
   else {
      scalar pb_`b' = 0
      local pb_`b' = 0
   }
   scalar spb_`b' = 100*round(((mean_`b' - `goo' ) / empSE_`b' ),.1)
   local spb_`b' : di %5.3f 100*round(((mean_`b' - `goo' ) / empSE_`b' ),.1)

   pairplot lb_`b' ub_`b' , sort(lb_`b' ) title("95% Confidence coverage") ///
      xtitle(Replicate) ytitle("95% Confidence Interval, `b'") ///
      legend(off) yline( `goo' , lcolor(red) lwidth(thick) )
   graph export 95`b'_`causalcondition'_`paramset'.png , replace width(3000)
   graph box `b' , title("Boxplot of estimated effect of exposure on slope")
   graph export box`b'_`causalcondition'_`paramset'.png , replace width(3000)
}

*intercept (mixed_b00)
*time (mixed_b10)
*exposure(mixed_b01)
*exposure*time (mixed_b11)
local B=$B
local rp = $paramset - 1
if `rp'==0 {
   local rp ""
}
/***Mixed model results***/
putexcel ///
   A`rp'1=("STATA") ///
   ///Row headers
   A`rp'2=("Scenario `paramset'")       A`rp'3=("Intercept")           A`rp'4=("Time (years)")        A`rp'5=("Exposure")      A`rp'6=("Exposure*Time") ///
   ///Column headers
   B`rp'1=("avg est beta")        C`rp'1=("Percent Bias")       D`rp'1=("avg est SE")          E`rp'1=("empirical SE")                 F`rp'1=("RMSE")       G`rp'1=("95% CI coverage")     H`rp'1=("TS") I`rp'1=("p-value") J`rp'1=("Standardized bias")  ///
   ///
   B`rp'3=(mean_mixed_b00_int)          B`rp'4=(mean_mixed_b10_time)      B`rp'5=(mean_mixed_b01_exposure)        B`rp'6=(mean_mixed_b11_exptime) ///
   C`rp'3=(pb_mixed_b00_int)            C`rp'4=(pb_mixed_b10_time)        C`rp'5=(pb_mixed_b01_exposure)          C`rp'6=(pb_mixed_b11_exptime) ///
   D`rp'3=(mean_mixed_b00_int_se)       D`rp'4=(mean_mixed_b10_time_se)   D`rp'5=(mean_mixed_b01_exposure_se)     D`rp'6=(mean_mixed_b11_exptime_se) ///
   E`rp'3=(empSE_mixed_b00_int)         E`rp'4=(empSE_mixed_b10_time)     E`rp'5=(empSE_mixed_b01_exposure)       E`rp'6=(empSE_mixed_b11_exptime) ///
   F`rp'3=(RMSE_mixed_b00_int)          F`rp'4=(RMSE_mixed_b10_time)      F`rp'5=(RMSE_mixed_b01_exposure)        F`rp'6=(RMSE_mixed_b11_exptime) ///
   G`rp'3=(Pcoverage_mixed_b00_int)     G`rp'4=(Pcoverage_mixed_b10_time) G`rp'5=(Pcoverage_mixed_b01_exposure)   G`rp'6=(Pcoverage_mixed_b11_exptime) ///
   H`rp'3=(Z_mixed_b00_int)             H`rp'4=(Z_mixed_b10_time)         H`rp'5=(Z_mixed_b01_exposure)           H`rp'6=(Z_mixed_b11_exptime) ///
   I`rp'3=(pvalue_mixed_b00_int)        I`rp'4=(pvalue_mixed_b10_time)    I`rp'5=(pvalue_mixed_b01_exposure)      I`rp'6=(pvalue_mixed_b11_exptime) ///
   J`rp'3=(spb_mixed_b00_int)           J`rp'4=(spb_mixed_b10_time)       J`rp'5=(spb_mixed_b01_exposure)         J`rp'6=(spb_mixed_b11_exptime) ///
   ///B7=(mean_dif_in_slope) ///
   using c:\trash\Sim_N`N'_nreps_`B'_`causalcondition', sheet("MixedResults") modify

/**Pull death info**/
putexcel A1=("Scenario") B1=("Proportion survived to age 65") ///
   A2=("ScenarioA") B2=(mean_p_surv65)    ///
   using c:\trash\Sim_N`N'_nreps_`B'_`causalcondition', sheet("StataDeathInfo") modify

/**Pull # individuals (N) and # simulations (B) and time between visits**/
putexcel A1=("Scenario") B1=("N (# individuals)") C1=("B (# simulations)") D1=("time btwn visits") ///
   A2=("ScenarioA") B2=($N) C2=($B)  ///
   using c:\trash\Sim_N`N'_nreps_`B'_`causalcondition', sheet("StataNandBinfo") modify

* For the PDF report

tex \begin{landscape}
tex \newpage
tex \subsubsection{Stata result}
tex Table `++tablecounter'. Stata: Causal condition $causalcondition, parameter set $paramset. \\
tex \begin{longtable}{l rrrrr rrrr } \hline
tex Parameter & avg est beta & Percent Bias & avg est SE & empirical SE ///
   & RMSE & 95 CI coverage & TS & p-value & Standardized bias \\
tex \hline

tex Intercept        & `mean_mixed_b00_int'       & `pb_mixed_b00_int'         & `mean_mixed_b00_int_se'       & `empSE_mixed_b00_int'        & `RMSE_mixed_b00_int'         & `Pcoverage_mixed_b00_int'       & `Z_mixed_b00_int'       & `pvalue_mixed_b00_int'      & `spb_mixed_b00_int'       \\
tex Time (years)     & `mean_mixed_b10_time'      & `pb_mixed_b10_time'        & `mean_mixed_b10_time_se'      & `empSE_mixed_b10_time'       & `RMSE_mixed_b10_time'        & `Pcoverage_mixed_b10_time'      & `Z_mixed_b10_time'      & `pvalue_mixed_b10_time'     & `spb_mixed_b10_time'      \\
tex Exposure         & `mean_mixed_b01_exposure'  & `pb_mixed_b01_exposure'    & `mean_mixed_b01_exposure_se'  & `empSE_mixed_b01_exposure'   & `RMSE_mixed_b01_exposure'    & `Pcoverage_mixed_b01_exposure'  & `Z_mixed_b01_exposure'  & `pvalue_mixed_b01_exposure' & `spb_mixed_b01_exposure'  \\
tex Exposure x Time  & `mean_mixed_b11_exptime'   & `pb_mixed_b11_exptime'     & `mean_mixed_b11_exptime_se'   & `empSE_mixed_b11_exptime'    & `RMSE_mixed_b11_exptime'     & `Pcoverage_mixed_b11_exptime'   & `Z_mixed_b11_exptime'   & `pvalue_mixed_b11_exptime'  & `spb_mixed_b11_exptime'   \\


tex \hline
tex \end{longtable}
tex \end{landscape}

foreach b in $betalist {
   lstrfun btex , subinstr("`b'","_","\_",.)
   tex \subsubsection{Parameter `btex'}
   tex \begin{center}
   tex \includegraphics[width=.8\textwidth]{95`b'_`causalcondition'_`paramset'.png} \\
   tex \includegraphics[width=.8\textwidth]{box`b'_`causalcondition'_`paramset'.png} \\
   tex \end{center}
}



