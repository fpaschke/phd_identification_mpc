General Information
===================

Both examples contain folders:
- ...\model -> modelica process model (*.mo), compiled FMU and resource files (weather: DEU_Berlin.103840_IWEC.mos and Temps of neighboring rooms RoomTempData.txt)
- ...\sim 	-> simulink Models, function to generate parameter structure for simulation (getParametes.m), scripts to run simulation and identification (run_sim and run_ident) and copy of above mentioned resource files
	- ...\sim\resource 	-> saved identified models
	- ...\sim\subfun\ctrl -> MPC and Reference Control algorithm 
	- ...\sim\subfun\ident-> functions that setup an model structure for identification estimate MPC models 
	- ...\sim\subfun\sig	-> functions to build Signals necessary for identification

A thorough description of the modelica process model can be found in the modelica files (*.mo)

Usage and Dependencies
======================

To run the examples you will need
- Windows environment. If you want to run on Linux you need to generate an FMU on a Linux machine from *.mo files.
- Matlab 2020b, or higher with Control System and Optimization Toolbox 
- Dymola (tested on dymola 2020 only). A Dymola license is also necessary for running the FMU within the simulink model.  

NOTE: Make sure to set Environment variable DYMOLA_RUNTIME_LICESE to corresponding value e.g. PATHTOLICENSEFILE\dymolalicense.LIC