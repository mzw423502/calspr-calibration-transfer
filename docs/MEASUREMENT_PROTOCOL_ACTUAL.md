# Actual measurement protocol used in CalSPR V2.1

## Scope

The executable study is an **angular-interrogation abstraction**, not a model of a named commercial instrument. It evaluates p-polarized Kretschmann reflectance at 633 nm from 60 to 80 degrees. The 201-bin inversion runs and the earlier 401-bin information scan are different numerical discretizations; they are not compared point-for-point as if they were the same experiment.

The angle bins may represent sequential angle samples or detector pixels in a simultaneous angular image. The present likelihood specifies only independent Poisson counts conditional on the means. It does not choose between those hardware implementations and does not model angular point-spread, read noise, saturation, or bin integration.

## Photon and background accounting

`total_incident_photons` is the total across **all exposures and all angular bins** in one virtual measurement unit. If an exposure receives fraction `f` of the budget and contains `J` bins, each bin receives `f N / J` incident photons. The mean reflected count in bin `j` is

`mu_j = (f N / J) [exp(log_gain) R(theta_j + angle_zero) + background_fraction]`.

Thus, `background_fraction = 0.001` is an additive background equal to 0.001 of the incident photons assigned to that bin. It scales with exposure time/budget in this abstraction. It is not a fraction of the reflected count. Changing 201 to 401 bins holds total photons fixed and changes the numerical sampling of a continuous response; the current model does not claim that this is equivalent to a physical change in pixel width or angular blur.

The target quantity is the sample refractive index `q` (or `Delta q` only in the historical same-device differential analysis). The known reference state is `q0 = 1.3330` by scenario definition. Its preparation uncertainty is not yet modeled; its photon cost is modeled explicitly.

## Transfer strategies and parameter sharing

| Strategy | Charged exposures | Known quantities | Fitted target | Shared fitted nuisance | Device-specific fitted nuisance | Fixed/omitted device difference |
|---|---|---|---|---|---|---|
| T0: rigid transfer | source-device known reference; target-device sample | reference index `q0`; nominal model family | target sample `q` | source thickness parameter under a rigid equality, common angle zero, common log gain | none | target-minus-source thickness fixed to zero |
| T1: explicit thickness difference | same as T0 | reference index `q0` | target sample `q` | source thickness, common angle zero, common log gain | target-minus-source thickness `delta` | target angle change fixed to zero in the stable-angle scenario |
| T2: target reference | source-device known reference; target-device known reference; target-device sample | both reference states have index `q0` | target sample `q` | source thickness, common angle zero, common log gain | target-minus-source thickness `delta` | target angle change fixed to zero in the stable-angle scenario |
| O: oracle | same exposures and allocation as T0 | true target-minus-source thickness is supplied to the estimator | target sample `q` | source thickness, common angle zero, common log gain | none | no thickness omission; this is an unattainable ideal comparison |
| T1-angle / T2-angle | corresponding T1/T2 exposures | reference index `q0` | target sample `q` | source thickness, common reader angle zero, common log gain | thickness difference and an added target-device angle offset | none of these two differences |

The shared log gain represents a controlled common multiplicative calibration of the abstract reader. It is a conditional assumption, not an empirical statement that separate devices or exposures have identical gain. Likewise, the stable-angle T0/T1/T2 scenario represents chips exchanged without an added angle-zero change. The separate angle-aware analysis represents reinstallation with a target-specific angular offset. Independent device gains and time-dependent drift remain later sensitivity analyses, not silently shared facts.

T2 is charged against the same total incident-photon budget as T0/T1. The exploratory allocations are `(source reference, target reference, target sample) = (0.20,0.40,0.40)`, `(0.33,0.33,0.34)`, and `(0.50,0.25,0.25)`. A target reference can constrain target-device nuisance parameters but reduces photons available to the target sample. Reuse of one target reference across samples would create a shared calibration error; the present unit contains one sample and does not treat such reuse as independent per-sample calibration.

## Historical regimes

The earlier absolute regime used one target sample and assigned the full photon budget to it. The earlier same-device differential regime used a known baseline and sample on one film, each receiving half the total photons; thickness, angle zero, and gain were shared. The earlier transfer regime used a source reference and target sample, each receiving half the total photons. All analyses use the joint likelihood of raw counts rather than treating a subtracted signal as independent Poisson data.

## Terminology correction

The historical CSV column `*_separation_fraction` stores `||r_q|| / ||u||`, the retained **derivative-norm ratio** after nuisance projection. The retained Fisher-information fraction is its square, `||r_q||^2 / ||u||^2`. For example, a stored value near 0.71 corresponds to about 0.50 retained information, not 71% retained information. Historical files are preserved; new outputs report both quantities under explicit names.

## Inference boundary

This protocol supports conditional computational statements about the stated angular-count model. It does not establish measured instrument drift, fabrication tolerances, manufacturing yield, achieved detection limit, or performance of all gold films.

