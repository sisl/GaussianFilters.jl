# GaussianFilters

GaussianFilters implements methods to define and run **Kalman**, **Extended Kalman**, **Unscented Kalman**, and **Gaussian-Mixture Probability Hypothesis Density** Filters on simulated data. It also implements simulation functions for the Kalman-class filters.

The Kalman Filter and its nonlinear extensions (the Extended and Unscented KFs) are used commonly to localize state information from noisy observations \[1\]. [Wikipedia's Kalman Filter entry](https://en.wikipedia.org/wiki/Kalman_filter) is a useful resource for learning about the applications and advantages of the Kalman Filter, along with its nonlinear extensions.

The Gaussian Mixture Probability Hypothesis Density Filter \[2\] is the analogous recursive Gaussian filter used for estimating the "time-varying number of targets and their states from a sequence of observation sets in the presence of data association uncertainty, detection uncertainty, noise, and false alarms". The GM-PHDF excels at tracking multiple targets in a low signal-to-noise environment without the costliness of methods with explicit data association (e.g. multiple hypotheses tracking (MHT) and the joint probabilistic data association filter (JPDAF)).

GaussianFilters uses automatic forward differentiation in order to save users from having to manually calculate Jacobians for nonlinear prediction and measurement steps.

1. Thrun, S., Burgard, W., & Fox, D. (2005). Probabilistic robotics. MIT press.
2. [Vo, B. N., & Ma, W. K. (2006). The Gaussian mixture probability hypothesis density filter. IEEE Transactions on signal processing, 54(11), 4091-4104.](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=1710358)
