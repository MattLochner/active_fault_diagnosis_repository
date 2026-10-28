# Active Fault Detection on a High-Pressure Leaching Process Simulation

This repository contains the code used in Matthew Lochner's masters project: *Active fault detection on a high-pressure leaching process simulation*.

The repository is split into 3 folders:

- **`scripts`** — contains the script files that are run to generate results related to the 3-tank problem and high-pressure leach case studies.
  Functions are imported from the source code.
  Functions related to the input design step are included in the scripts folder as they must be customized to each application, as
  control setpoints need to be manipulated inside the function.
- **`src`** (source code) — reusable functions are defined in the source code and exported through `full_simulation.jl`.
- **`results`** (source code) — the results for the MAPE analysis of state estimation performance and Monte Carlo analysis of input design solve time are reported in this   folder.

---

## Scripts

### `3_tank_problem_observability`
This script runs the observability analysis for the 3 tank system.

### `3_tank_problem_EKF.jl`
This script runs the Three tank problem ODEs with the extended Kalman filter integrated into the framework.
It was used in investigating the effectiveness of the state estimation on the 3 tank problem.

### `3_tank_problem_with_input_design.jl`
This script contains the full application of the AFD methodology to the 3-tank problem.
By commenting sections, the user can select between constant inputs, input design with input constraints and input design with input and state constraints.
The input design functions for applying the AFD methodology to the 3-tank problem are also included in this file.

### `BMR_EKF_Input_Design_Functions.jl`
This script contains the input design functions adapted for application to the high-pressure leach case study.
It is acknowledged that designing the input design functions as model agnostic would have been ideal. However, it is necessary to change
setpoints in the model inside the function. The input design functions therefore have to be customized to each model.

The following functions are defined here:

- `ExtendedStatePrediction(xhatmin1, Pposmin1, nx, ny, R, Q, params, tstart, tend, measurement_pred_function, model_ODE_function)`
  → returns `xhatneg, Pneg_2`
  - Propagates from current model estimate forward in time
- `LinearisedMeasurementPrediction(xhatneg, Pneg, params, R, measurement_pred_function, tend)`
  → returns `yhat, Py`
  - Calculates a measurement prediction where the mean is calculated by the nonlinear function and the mean is calculated
- `Input_Design_Function(u_vec, p)`
  → returns `s`
- `bhatt_dist(x, X, y, Y)`
  → returns `t1 + t2`
- `Calc_Obj_Function(beliefs, yhats, Pys, no_of_modes)`
  → returns `s`

### `BMR_EKF_plot_responses.jl`
This script contains the code used to plot the model's response to step changes in control setpoints and disturbances. Step changes are implemented
by commenting and commenting out lines of code.

### `BMR_EKF_test_input_design_time`
This script is the example script which was used to develop the `test_input_design_function`; the function is then submitted to the HPC in scripts
`Input_D_MC_1-100.jl`.

### `BMR_EKF_with_hourly_samples`
This script is the implementation of the EKF with the 36 online measurements and the 13 composition samples taken hourly.

### `BMR_EKF.jl`
This script is the implementation of the EKF with only the 36 online measurements.

### `BMR_EKF_with_input_design`
This script contains the full implementation of the AFD methodology to the high pressure leach.

### `BMR_gt_params`
This script contains the parameters for the ground-truth model.
This includes the autoregressive functions (ARFs) that generate LinearInterpolations for the stream 1, 2 and 18 temperatures and slurry solid fraction.
It exports `data`, which is a struct that contains all the model parameters.

### `BMR_kf_params`
This script contains the parameters for the process models used in the EKF.
As no model mismatch is assumed, there is no difference between the parameters in the process models and ground-truth models, except that the
parameters specified as ARFs are not included in `BMR_kf_params` as they are estimated as auxiliary states by the EKF.

### `BMR_initial_conditions`
Initial conditions obtained by running model to steady state.

### `BMR_observability`
Contains the script to run the observability analysis on the BMR model.

### `BMR_run_ODE.jl`
Runs the high-pressure leach ODEs.

### `Input_D_MC_xx.jl` / `Input_D_MC_xx.sh`
Contains the scripts used to run the Monte Carlo analysis of optimisation solve time on the HPC.

### `Week_Simulation_x.jl`
Contains the scripts used to run the week long simulation used to demonstrate how AFD can be applied in weekly 8 hour shifts.

---

## Source code (`src`)

### `3_tank_generate_measurements.jl`
Contains function which adds random noise to measurements for 3-tank process.
This generates the ground-truth measurements.

### `3_tank_measurement_prediction.jl`
Contains function which generates measurement predictions for 3-tank process. As no plant-model mismatch is assumed, this is identical to the `generate_measurements_3_tank_system` function, just without measurement noise.
The measurement prediction function is used as h(x) in the Kalman filter.

### `3_tank_problem_ODEs.jl`
Contains ODE function for ground-truth process and model process. As no plant-model mismatch is assumed, these models are identical.

### `BMR_gt_generate_measurements.jl`
This file contains 4 functions, which combined generate the ground-truth measurements at each timestep from the states propagated by the ODEs.

### `BMR_gt_ODEs.jl`
Contains the following 4 functions.
`if` statements are used so that the same `Calculate_gt_autoclave_odes` & `Calculate_model_autoclave_odes` function can be used
for normal operating conditions and both fault modes.

1. **`Calculate_gt_autoclave_odes`** — This function represents the ground-truth process.
2. **`Calculate_oxygen_solubility2`** — This function is a helper function used to calculate the oxygen solubility in the slurry.
3. **`Calculate_reaction_rates`** — This function is a helper function that calculates the reaction rates of selected reactions in the autoclave compartments.
4. **`Calculate_model_autoclave_odes`** — This function is the process model f(x) used in the EKF. The model is identical to the ground-truth process model except that
   the parameters modelled as autoregressive functions and fault mode magnitudes are modelled as additional states which are estimated by the EKF.

### `BMR_kf_generate_measurements`
Contains 2 functions that generate the measurement predictions used in the EKF and Input Design Functions.
`if` statements are used so that the same function can be used for all system modes.

- `BMR_kf_measurement_prediction(X, data, t)`
- `BMR_kf_measurement_prediction_with_samples(X, data, t)`

### `extended_kalman_filter.jl`
Contains the following 4 functions:

1. **`lyap!`**
   This function is the right-hand side of the continuous-time Lyapunov differential equation. It is written in the in-place form
   that DifferentialEquations.jl expects. When called by `ODEProblem` in the EKF, it propagates the covariance matrix of the state
   estimate through the linearised process model for 1 timestep.
2. **`check_covariance`**
   This function checks whether the covariance matrix is positive-definite.
   If any eigenvalues are negative, then the minimum (i.e. most negative) eigenvalue is identified
   and `1e-10 + min(λ)` is added to the diagonal.
   The symmetry of the matrix is enforced by `Psym = (P + P') / 2`.
   It is vital to the Kalman filter that the matrix is positive-definite.
3. **`ekf`**
   `ekf(xhatmin1, Pposmin1, nx, ny, R, Q, params, tstart, tend, measurement_pred_function, model_ODE_function, y)`
   → returns `yhat, Py, xhatpos, Ppos`
   This function uses the process model and measurement prediction function to generate
4. **`run_ODE_model`**
   `run_ODE_model(tStart, tEnd, x, p, model_ODE_function)`
   → returns `ODEsol2.u[2]`
   Simple function that runs an ODE given the start time, end time, initial state, parameters and function.
   The function returns the state at the end time.

### `full_simulation.jl`
Imports all functions defined in the source code and exports the names as references.

### `multiple_model_method.jl`
Contains the following functions:

1. **`multiple_model_method`**
  This function the updated system mode probabilities based on the measurement predictions generated by the ekf function, the prior probabilities and the measurements at the current timestep. 
2. **`flatten_prior`**
  This function flattens the prior distribution by a user selected factor. This is similar to a fading memory factor and means
  recent measurements are weighted more heavily by the multiple model method.
