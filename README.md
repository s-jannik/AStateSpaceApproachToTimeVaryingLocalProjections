# A State Space Approach to Time-Varying Local Projections

This repository contains R code accompanying my paper titled A State Space Approach to Time-Varying Local Projections. The repository contains self-written implementations from scratch of the proposed time-varying parameter Local Projections method alongside scripts used to obtain the main results of the paper. 

To ensure correct implementation of the Kalman Filter and Kalman Filter Smoother recursions at the heart of the TVP-LP methodology, i replicated the figures on pages 16 and 22 of Durbin and Koopman (2012). The replications are reported in the bottom of this repository. 

## Overview of scripts and data files

### Scripts

| File name    | Description |
| :-------- | :------- |
|main.R| runs the simulation study and empirical study|
|kf_kfs.R| contains the Kalman Filter recursions, Kalman Filter Smoother recursions and implementations of state space models used to replicate figures from Durbin and Koopman (2012)|
|tvplp.R| contains the TVP-LP implementation, relying on the Kalman Filter and Kalman Filter Smoother implementations in kf_kfs.R|
|studies.R| contains function used for the simulation study and the empirical study|
|support_functions.R| contains various functions for data cleaning and manipulation as well as stationarity transformation functions|
|function_testing.R| contains code replicating figures from Durbin and Koopman (2012)|

### Data files 
| File name    | Description |
| :-------- | :------- |
|fred_qd.csv| contains data used for the empirical study|
|Nile.csv| contains data used for replicating figures from Durbin and Koopman (2012)|


## Replication figures
![replication_filter](https://github.com/s-jannik/AStateSpaceApproachToTimeVaryingLocalProjections/assets/143873944/1cc1ad61-6730-4c67-b48c-432ef9466a36)

![replication_smoother](https://github.com/s-jannik/AStateSpaceApproachToTimeVaryingLocalProjections/assets/143873944/34402f1f-549f-459d-a968-8a30a1fe0f02)


## References 
Koopman, S.J. and Durbin, J. (2012). Time Series Analysis by State Space Methods. Oxford Statistical Science Series. Oxford.
