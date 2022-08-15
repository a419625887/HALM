# Horizon-assisted Lithologic Modeling

This is a python implementation of the Horizon-assisted Lithologic Modeling (HALM) method. The HALM method aims to contruct high-fidelity litholgic model given well log data and dip variations of horizons (bedding surfaces). The HALM method involves horizon restoration, discretization, and interpolation.

## Example problem

This example problem is used to demonstrate procedures of the HALM method. In this problem, we try to model lithofacies in the domain bounded by two horizons (bedding surfaces) using lithologic data in five synthetic well logs. The well logs contain bimodal lithology, shown in blue and yellow. The green portion is not constrained by the two horizons, and thus, is not used for modeling.
<img src="/Plots/Dip_domain_view.png" width="500">

First, the bedding surfaces and well logs are transformed into a non-dipping domain through a horizon restoration technique. The transfomred surfaces are placed horizontally and the transformed logs can be slightly inclined.
<img src="/Plots/Nondip_domain_view.png" width="500">

Then, an indicator natural neighbor interpolation is performed in the non-dipping domain given the well log data, which produces a litholgic modeling in the non-dipping domain.

<img src="/Plots/Nondip_lithology_view.png" width="500">

Finally, the interpolated lithology is transformed back to the orginal domain. As a result, the dip variations of modelled lithofacies follow the curvature of the horizons above and below.

<img src="/Plots/Dip_lithology_view.png" width="500">

## Requirements

The following requirements can be directly installed from PyPi:
* numpy
* math
* matplotlib
* scipy
* os
* joblib
* subprocess
* csv
* shapely
* pyvista
* metpy

## Data description

Input data for the above example problem is available in `Input` folder.

### Input data
* `top_surface_point.txt`: Points that define the top bedding surface. Columns are point ID, X, Y, Z coordinates
* `bot_surface_point.txt`: Points that define the bottom bedding surface. Columns are point ID, X, Y, Z coordinates
* `triangle_faces.txt`: IDs of points that define each triangle on the top or bottom bedding surface. Columns are number of point, point ID, point ID, point ID
* `domain_discretization.txt`: 2D grid that defines the model horizontal resolution. Columns are 2Dgrid ID, X, Y, Z of the top surface at the XY location, Z of the bottom surface at the XY location
* `Logs.csv`: Well log data consist of well names, well locations, well datum, and lithology at depths.

### Output data
* `restored_top_point_data.txt`: Points that define the transformed top bedding surface in the non-dipping domain. Columns are point ID, X, Y, Z coordinates
* `restored_bot_point_data.txt`: Points that define the transformed bottom bedding surface in the non-dipping domain. Columns are point ID, X, Y, Z coordinates
* `Interpolation_nondipping_results.txt`: 3D grid with lithofacies in the non-dipping domain. Columns are 2Dgrid ID, X, Y, Z, lithofacies(binary)
* `grid_lithology.txt`: 3d grid with lithofacies in the original domain. Columns are 2Dgrid ID, X, Y, Z, lithofacies(binary)
* `Mapped_logs.csv`: Well log data in terms of lithologic sequences in the non-dipping domain. Columns are well name, X, Y, Z, lithology
* `log_intersection.txt`: Intersection points between each well log and the bedding surfaces. For visualization only

## Usage

The usage of the code is quite simple. After preparing all input data and put them in `Input` folder, run the following command in a terminal:

python halm.py

The results can be then found in `Output` folder. Run the following command to visualize the results:

python visualization.py

This visualization tool is for a relatively small dataset (e.g. the example problem). Other softwares may be used for large datasets.
