* fh2016sim project
* sim-000-master.do
* ALG ERM YG TF
* 20 Jun 2016
* -------------------
* in line below, full UNC path would be best
*wftree , m(fh2016sim) p(TF ALG ERM YT)
*! junction "c:\work\fh2016sim" "C:\Users\aldengross\Dropbox\fh2016sim"
*cd c:\work\fh2016sim\posted\analysis
*automaster2 , m(fh2016sim) a(sim) u(ALG ERM YG TF)


wfenv fh2016sim , a(sim) alg
global N=100,000
local N=$N
global sample=2000
local sample=$sample
global B=2000
local B=$B
*global causalcondition = 5 // 1 2 3 4

include $analysis/sim-005-start-latex.do

foreach causalcondition in 1 2 3 4 {
   global causalcondition=`causalcondition'
   capture erase c:\trash\Sim_N`N'_nreps_`B'_`causalcondition'.xlsx
   tex section{Causal condition `causalcondition'}
   foreach paramset in 1 2 {
      global paramset=`paramset'
      tex \subsection{Parameter set `paramset'}

      include $analysis/sim-001-preambling.do // parameter coefficients
      *include $analysis/sim-105-call-source.do
      *include $analysis/sim-110-create-variables.do
      *include $analysis/sim-120-select-cases-for-analysis.do
      *include $analysis/sim-205-table1.do       // does stuff
      *include $analysis/sim-210-figure1-main-outcome-undadjusted.do
      include $analysis/sim-305-models-stata.do // does stuff

      **include $analysis/sim-307-models-mplus.do // does stuff
      *include $analysis/sim-310-table2.do
      *include $analysis/sim-320-figure2.do
      *include $analysis/sim-405-sensitivity-analyses.do
      *include $analysis/sim-505-stuff-in-text.do
   }
   copy c:\trash\Sim_N`N'_nreps_`B'_`causalcondition'.xlsx ///
      $text/Sim_N`N'_nreps_`B'_`causalcondition'.xlsx , replace
}

include $analysis/sim-990-close-latex.do

