Code for the preprint _Preventing evolutionary rescue in cancer_, Srishti Patil, Armaan Ahmed, Yannick Viossat, Robert John Noble; bioRxiv 2023.11.22.568336; doi: https://doi.org/10.1101/2023.11.22.568336

## Simulation data
The folder ./simulation_data contains data from Gillepsie simulations used to generate extinction probability plots and heatmaps. Find further description in the subfolders.

## Simulation code
Code for the Gillespie simulations is in ./simulations

## Analytical model and plotting
The folder ./plotting contains scripts to calculate extinction probabilities using the analytical model and plot them along with simulation data:

1) ./plotting/EP_parameter_space.R and ./plotting/high_prob_regions.R generate the data and plots shown in Fig3A, Fig4B and FigA.4.
2) ./plotting/extinction_probability.R, ./plotting/plotting_Nq.R and ./plotting/plotting_extinction_probability.R generate the data and plots shown in Fig2, Fig3B,C, Fig4A, FigA.2, FigA.5 and FigA.6.
3) ./plotting/heatmaps_analytical.R and ./plotting/plotting_heatmaps.R generate the data and plots shown in FigA.7.
