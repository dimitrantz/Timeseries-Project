# Timeseries-Project
Implementation of linear and non-linear analysis on eye movement datasets

## Project description:

The project deals with eye-movement signals measured in the Department of Laboratory Otolaryngology / Head & Neck Surgery, Haukeland University Hospital, Bergen, Norway (see related article on [link](https://www.sciencedirect.com/science/article/pii/S0010480997914415). The experiment requires the person to sitting comfortably in a chair and trying to focus on the center of a screen, which they constantly spend black-and-white vertical stripes with specific direction (right / left) and constant angular velocity (30 and 60% sec). This experiment is done to diagnose vertigo (imbalance in the ear maze). Each one Time series is about such an experiment for different times, speed and category of person (healthy / patient). In the experiment, the eye imperceptibly follows the stripes and returns, generating a continuous signal
formations, consisting of an upward slow voltage (as the eye tricks and follows the stripes), and a sharp drop (at the point where the eye returns). This movement of the eye is called an optomotor nystagmus (optokinetic nystagmus, OKN).
The two datasets (time series) are in the same conditions of direction and speed of the light stimulus but in two different individuals, which can be both healthy, both patients, or one healthy and the other patient. They may also haveand different length. All measurements were made by sampling at 100Hz. The files are given on his [website's course](http://users.auth.gr/dkugiu/Teach/TimeSeriesTHMMY/index.html). We can use the extremes.m function as explained below.
### Linear analysis:
1. For each of the two OKN time series, call the extremes.m function with appropriate values at
function parameters (mainly the parameter set by the local window). This function will you
give a new time series of local Across values ​​and the corresponding times they appear. Based
the local extreme price ranges and times, you can create the following time series:
local maximum values ​​(amplitude maxima, AMA), local minimum values ​​(amplitude minima, AMI),
the local minimum distance to the next local maximum, ie the difference of local maximum values
and minimum (amplitude min-max difference, AMD), lasting from local minimum to next local
time to maximum (TMA), the duration from the local maximum to the next local minimum (time to
minimum, TMI), the duration from local maximum to next local maximum (TBP).
You are asked to select one of these time series (AMA, AMI, AMD, TMI, TMA, TBP) to apply methods
of linear analysis. The same type of timer will select for both original OKN timings.
2. For each of the two extreme time series of the type you selected, perform the linear analysis,
basic steps should be: a) eliminating the voltage if you consider that the time series has a tendency, b)
calculation and graphs of the autocorrelation function and partial autocorrelation function; b) Investigation
suitable linear model (justification of model selection, calculation of some statistical error
adaptation, e.g. NRMSE), (c) predicting one and two steps ahead of its evaluation set
model predictability (calculation of a statistical forecast error, eg NRMSE).
3. Compare the two extreme time series based on the results of the linear analysis. Does any of them contain (stronger) correlations and can therefore be better adapted and / or predicted;
4. Select a second type of extremes and repeat steps 2 and 3. Compare your results.
### Non-linear analysis:
5. For each of the two OKN time series, consider that they are likely to derive from a causal potential
system (with noise?) and proceed to the non-linear analysis, the basic steps of which should be: a)
scatter diagrams in two and three dimensions, b) detection of suitable reconstruction parameters
(c) calculate the correlation dimension for differentdip dimensions and an appropriate lag parameter, and (d) exploring an appropriate local model based on nearest neighbors and predicting 1 to 10 steps ahead of an evaluation set
model predictability (calculation of a statistical forecast error, eg NRMSE).
6. Compare the two LCN time series based on the results of the nonlinear analysis. Looks like some of them
both contain (stronger) nonlinear dynamics and can be better adapted and / or predicted;
7. Comment on your conclusions from your linear analysis amd nonlinear analysis and also comment on which system each of the two time series can come from(stochastic / causal, linear / non-linear, low / high dimension, small / large
complexity), as well as if they appear to differ. Could you provide an answer if the time series come from people in the same category (healthy / sick) or not;

## Contributors:

Diamanti Maria -  mfdiamanti@ece.auth.gr

Ntzioni Dimitra - dntzioni@ece.auth.gr
