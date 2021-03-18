18-03-2021 - Benjamin Dalton - dalton@zedat.fu-berlin.de

Simulation of active filament networks and analysis. Source code in C++, analysis scripts in MATLAB.

Code suitable for simulations similar to those presented in Dalton et. al, 2021.

Three distinct code versions included:

- 01_Active_Flow: Simulate the flow of cross-linked dynamics filaments in flow channel confinement
  driving by active motors interacting with rigid pinning field structure. Analysis scriots to 
  generate mpg-4 movies and calculate long-range velocity profiles
  
- 02_Active_MR: Simulate active microrheology of cross-linked filament networks. Analysis scripts 
  to generate mpg-4 movies and present trapped filament response trajectory
  
- 03_Branching_Nucleation: simulate the collision of two branched microtubule networks, interacting
  with active motors to generate sorter-polarity configureations. Generate mpg-4 movies.
  
Directories contain documentation: ReadMe.pdf

Run each code version with included Makefile:

- make clean
- make
- ./execute_dynamics
