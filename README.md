<h2>Hidden Markov regime inference model for error propagation mitigation in modular digital twins</h2>

This repo provides scripts to generate the figures in the manuscript. The article introduces a computational model for understanding error propagation in a modular digital-twin pipeline and for designing regime-aware mitigation strategies.

The repo contains the following items:

<ul>
  <li>Digital twin simulation + ARX surrogates across a six-module pipeline (actuator, plant, sensors, fusion, KPI).</li>

  <li>Residual construction and HMM regime inference, including posterior visualization and Viterbi decoding.</li>

  <li>  Deliberate fault-injection experiments with ground-truth regimes and Hungarian matching for state-to-regime mapping (regime identification accuracy figures).</li>

  <li> Model-order selection for the number of HMM states using BIC sweeps.</li>

  <li> Robustness/sensitivity analysis</li>
  <li> MCDA inspired mitigation action to consider costs, etc. </li>
</ul>


Please refer to <a href="https://rpubs.com/Annice/Hungarian">this Rmd notebook on Rpubs</a> for better understanding the Hungarian algorithm. 
