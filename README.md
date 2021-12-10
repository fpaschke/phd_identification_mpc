General Information
===================
This project contains software, models and identification data which was developed during my phd thesis. It was developed and tested using MATLAB 2020b. It has th following content:
- Toolbox for system identification (folder: +idModels) 
- Examples (folder: examples)
- Modelica and compiled fmu binaries that are used as process models for mpc design. (folder: examples/Chap6)
- Some utility functions for plotting, formatting plots and so on (folder +util)

Dependencies
============
Software was written under MATLAB 2020b. The Toolbox and all examples require MATLAB with control systems and optimization toolbox (tested on 2020b). If optimization toolbox is not available OPTI toolbox can be used for identifcation as well (https://github.com/jonathancurrie/OPTI). The Modelica Models were created using Dymola 2020 and Buildings library (https://simulationresearch.lbl.gov/modelica/) and were compiled to FMU binaries under Windows. If you want to run the the examples in Chapter 6 you will need MATLAB, Simulink and a dymola license (The FMU was created with standard license of dymola (no XBM) and thus needs a valid license file.). Make sure to set a path variable DYMOLA_RUNTIME_LICENSE to valid path (e. g.*PATH_TO_LICENSEFILE*/dymola.lic) so that simulink can find the necessary license.

Usage
=====
Go to specific folder under examples and open the run_******.m file.