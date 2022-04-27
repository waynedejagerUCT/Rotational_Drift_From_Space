# Rotational_Drift_From_Space

This repository contains the code used to produce the manuscript "Rotational drift in Antarctic sea ice: pronounced cyclonic features and differences between data products" submitted to the Copernicus peer reviewed journal The Cryosphere

de Jager, W. and Vichi, M.: Rotational drift in Antarctic sea ice: pronounced cyclonic features and differences between data products, The Cryosphere, 16, 925–940, https://doi.org/10.5194/tc-16-925-2022, 2022.




-------
The VorticityFeatureDetection_v001.py script is used to iterate through EUMETSAT OSI-405-c low resolution 48hr sea ice drift data and detect circular vorticity features in the sea ice. The script pools all detected features into yearly .csv files.
Feature parameters extracted include:
* Time0 of the displacement tracking window
* Time1 of the displacement tracking window
* Nan percentage (represented as a faction of 1)
* The x-coordinate of the centre of the feature, corresponding to the polar stereographic grid of the OSI-405-c product.
* The y-coordinate of the centre of the feature, corresponding to the polar stereographic grid of the OSI-405-c product.
* Mean Vorticity (in units "per second")
* Standard deviation of the vorticity (in units "per second")
* Mean uncertaity of the vorticity (in units "per second")
* Standard deviation of the uncertaity of vorticity (in units "per second")

Adjustable parameters include:
* Temporal range
* Spatial range (can work for both Arctic and Antarctic regions with suitable adjustments made to the boundary parameters)
* Sensor products (both single- and multi-sensor products)
* Radius of the feature
* Strictness of the valid vorticity threshold requirement applied to each featur (Note that a isnan_threshold value = 0.1 means that a maximum 10% of vorticity values can be NaN, i.e., a 90% validity threshold)
* Flag status according to OSI-405-c product index

All neccessary python packages are listed in the import list.
