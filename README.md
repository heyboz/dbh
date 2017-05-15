# Exploring Patterns of Cluster Size Distribution
### main.py
Clean data, define 3 simplest null models, define cluster size and find nearest adult tree of child tree, calculate "Neighborhood Crowding" index. Outputs are in "/Outputs/nearest/all/nc", which contain information of every adult tree of every species, including location, NC index, environment data. #Additionally, Paint contour, try dbh percolation...
### cs_all_species.py
count cluster size of every adult for every species, with different definitions which are counting number, the sum of children's basal areas, the sum of all trees' basal areas. Outputs are in "/Outputs/nearest/all/cs_scanner/". #Additionally, try climaxD to see how children distribute in clusters...
### weibull-gamma.R
fitting weibull, gamma and power-law distribution. Outputs are in "/Outputs/nearest/monoecious/fitting_parameter/". #Additionally, paint PDF...
### lm.R
linear reggression: y is cluster size, x are all factors it can be, including dbh, NC index, environmental and topological indices. Outputs are in "Outputs/nearest/monoecious/effect_size_numeric/" and "/effect_size_fig/".
### nc_effect.R
see the effect size of NC index for cluster size.
### passed-failed.R
see the difference of the fitting parameters between those species passed weibull & gamma fitting and those failed.
