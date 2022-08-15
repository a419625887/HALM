# Horizon-assisted Lithologic Modeling

This is a python implementation of the Horizon-assisted Lithologic Modeling (HALM) method. The HALM method aims to contruct high-fidelity litholgic model given well log data dip variations of horizons (bedding surfaces). The HALM method involves horizon restoration, discretization, and interpolation.

## Example problem

This example problem is used to demonstrate procedures of the HALM method. In this problem, we try to model lithofacies within two horizons (bedding surfaces) given lithologic data in five synthetic well logs. The well logs contain bimodal lithology, shown in blue and yellow. The green portion is not constrained by the two horizons, and thus, is not used for modeling.
![alt text](https://github.com/a419625887/HALM/edit/main/Plots/Dip_domain_view.png)

First, the bedding surfaces and well logs are transformed into a non-dipping domain through a horizon restoration technique. The transfomred surfaces are placed horizontally and the transformed logs can be slightly inclined.

Then, an indicator natural neighbor interpolation is performed in the non-dipping domain given the well log data, which produces a litholgic modeling in the non-dipping domain.

Finally, the interpolated lithology is transformed back to the orginal domain.


