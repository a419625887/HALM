# Horizon-assisted Lithologic Modeling

This is a python implementation of the Horizon-assisted Lithologic Modeling (HALM) method which is proposed by Yang et al. (2024). The HALM method aims to contruct high-fidelity litholgic model given well log data and structures of horizons (bedding surfaces). The HALM method involves horizon restoration, domain transformation, and lithologic modeling. The python codes in this repository run the HALM method in a serial mode.

## Demonstration

The HALM method needs horizon (surface) data and well log data as input. As shown in the figure below, the model domain is bounded by two horizons. The horizons are represented by triangular meshes. There are five well logs containing bimodal lithology, shown in blue and yellow. The green portion of logs is outside the domain, and thus, is not used for modeling.

<img src="/Demo/Dip_domain_view.png" width="500">

First, the horizons and well logs are transformed into a non-dipping domain through a flexural restoration technique. The transformed surfaces are placed horizontally and the transformed logs can be slightly inclined.
<img src="/Demo/Nondip_domain_view.png" width="500">

Then, an indicator interpolation approach is performed in the non-dipping domain given the well log data, which produces a lithologic model in the non-dipping domain.

<img src="/Demo/Nondip_lithology_view.png" width="500">

Finally, the non-dipping lithofacies is transformed back to the orginal domain. As a result, the dip variations of modeled lithofacies conform to the horizons above and below.

<img src="/Demo/Dip_lithology_view.png" width="500">

## Citation
Yang, S., Tsai, F.TC. & Yarbrough, L.D. Horizon-assisted lithologic modeling: understanding Mississippi Embayment and Coastal Lowlands aquifer systems in Louisiana and southwestern Mississippi, USA. Hydrogeol J (2024). https://doi.org/10.1007/s10040-024-02804-z

## Requirements

The following requirements can be directly installed from PyPi:
* numpy
* math
* matplotlib
* scipy
* os
* subprocess
* csv
* vtk
* pyvista
* metpy

## Examples

There are two examples included in this repository. In the `Example1`, lithofacies as results of uniform sedimentation is modeled in a folded structure. The `Example2` is to construct a lithologic model with more complicated deposition patterns in a syncline structure.

## Data description

Input data should be put in `Input` folder. Output of the code would be saved in `Output` folder.

### Input data
* `top_surface_point.txt`: Nodal points of the triangular mesh that represents the top bedding surface. Columns are X, Y, Z coordinates.
* `bot_surface_point.txt`: Nodal points of the triangular mesh that represents the bottom bedding surface. Columns are X, Y, Z coordinates. The X and Y in this file are identical to those in `top_surface_point.txt`, only Z values are different.
* `triangle_faces.txt`: IDs of points that define each triangle face on the top or bottom bedding surface. Columns are number of point (must be 3), 1st point ID, 2nd point ID, 3rd point ID.
* `domain_discretization.txt`: 2D grid that defines the horizontal resolution of lithologic model. Columns are 2Dgrid ID, X, Y, Z of the top surface at the XY location, Z of the bottom surface at the XY location.
* `Logs.csv`: Well log data consist of well names, well locations, well datum, and lithology at depths.

### Output data
* `Geo-model.dat`: The final lithologic model. Columns in the file are 2Dgrid ID, X, Y, Z, lithofacies (binary). This 2D grid is the same as that specified in `domain_discretization.txt`.
* `Geo-model_nondip.dat`: The lithologic model in the non-dipping domain. Columns are grid ID, X, Y, Z, lithofacies (binary).
* `Grid2D_nondip.txt`: Location of 2D grid center in the non-dipping domain. Columns are grid ID, X, Y.
* `Logs_nondip.csv`: Well log data in terms of lithologic sequences in the non-dipping domain. Columns are well name, X, Y, Z, lithology.
* `restored_top_point_data.txt`: Nodal points of the triangular mesh that represents the transformed top bedding surface in the non-dipping domain. Columns are point ID, X, Y, Z coordinates.
* `restored_bot_point_data.txt`: Nodal points of the triangular mesh that represents the transformed bottom bedding surface in the non-dipping domain. Columns are point ID, X, Y, Z coordinates.

## Usage

The HALM codes are in `Source` folder. The main script is `halm.py`. The main script calls the modules in `HALM_lib.py`.

Usage of the code is quite simple. After putting all python scripts from the `Source` into the same directory as the `Input` folder. run the following command in a terminal:
```
python halm.py
```
The results can be then found in `Output` folder. To visualize the model results, put all python scripts from the `Visual_tools` into the same directory as the `Output` folder, then run the following command:
```
python visualization.py
```
Figures are exported to `Plots` folder. This visualization tool is for a relatively small dataset (e.g. the example problems). Other softwares may be used to visualize large datasets.

## License

MIT
