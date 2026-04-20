
# Supplemental materials for our era-adjusted baseball statistics [website](https://eckeraadjustment.web.illinois.edu/)

Era adjustment is made via [Full House Modeling](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-19/issue-2/Comparing-baseball-players-across-eras-via-novel-Full-House-Modeling/10.1214/24-AOAS1992.short). This model computes era-adjusted statistics through a principled balancing of how well a player performed "vs their peers" and the size of the MLB talent pool.


## Version 3.0 (consistent with [website](https://eckeraadjustment.web.illinois.edu/))

We present an updated version of the website with:

- New [talent pool calculation](https://htmlpreview.github.io/?https://github.com/ecklab/era-adjustment-app-supplement/blob/main/writeups/MLB-talent-pool-v3.0.html). This estimate is similar to the [previous version](https://htmlpreview.github.io/?https://github.com/ecklab/era-adjustment-app-supplement/blob/main/writeups/MLBeligiblepop.html). The main change is that we change the reference group from aged 20-29 eligible males to aged 20-29 eligible white males. Results do not change too much.
- New stats: 2B, 3B, SLG, OPS, JAWS
- Data through 2025 season
- Better UI, comprehensive leaderboards with stat download, and player-specific linking 


## Version 2.1

We present an updated version of the website with:

- An updated [writeup](https://htmlpreview.github.io/?https://github.com/ecklab/era-adjustment-app-supplement/blob/main/writeups/era_adjusted_V2.1.html) that highlights some results and discussions based on our era-adjusted statistics.
- Data through 2024 season
- Better UI on mobile devices
- Added loess smoothing of home run counts for Dead Ball era players

## Version 2.0 (consistent with [paper](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-19/issue-2/Comparing-baseball-players-across-eras-via-novel-Full-House-Modeling/10.1214/24-AOAS1992.short))

We present a new version of the Full House Model with:- A new and better justified [estimation of the talent pool over time](https://htmlpreview.github.io/?https://github.com/ecklab/era-adjustment-app-supplement/blob/main/writeups/MLBeligiblepop.html)
- Version 2.0 of a [detailed writeup](https://htmlpreview.github.io/?https://github.com/ecklab/era-adjustment-app-supplement/blob/main/writeups/era_adjusted_V2_I.html) that highlights some results and discussions based on our era-adjusted statistics.
- A [technical report](https://htmlpreview.github.io/?https://github.com/ecklab/era-adjustment-app-supplement/blob/main/tech-reports/tech-report-2023.html) that goes through important calculations and additional analyses- Data through the 2023 season- Removed rotation adjustments
- Removed smoothing


## Version 1.0

The app contains several era-adjusted statistics for baseball players which are obtained from the Full House Model. This model computes era-adjusted statistics through a principled balancing of how well a player performed "vs their peers" and the size of the MLB eligible population. Under this model, great all-time statistics requires that an MLB player is both better than their peers and played during a time in which the MLB eligible population is large. In this way, the model constructs an even playing field that extends across eras.

For a fun read on the Full House Model see our three-part series: [article I](https://htmlpreview.github.io/?https://github.com/ecklab/era-adjustment-app-supplement/blob/main/writeups/article_I.html), [article II](https://htmlpreview.github.io/?https://github.com/ecklab/era-adjustment-app-supplement/blob/main/writeups/article_II.html), [article III](https://htmlpreview.github.io/?https://github.com/ecklab/era-adjustment-app-supplement/blob/main/writeups/article_III.html)

Our original estimation of the talent pool can be seen [here](https://htmlpreview.github.io/?https://github.com/ecklab/era-adjustment-app-supplement/blob/main/writeups/MLB_eligible_pop.html).



