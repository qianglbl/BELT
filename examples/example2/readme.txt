-------------------
backward tracking a linear accelerator
-------------------
Input files:
pts.in -- initial distribution information
eblt.in -- accelerator lattice input
-------------------
Output files:
-------------------
fort.2: 3 columns - distance, kin. energy, gamma , <dz>, rms dz, <delta gamma>, rms delta gamma
fort.100: initial current profile (3 columns) - bunch length (m), charge per cell, current (A)
fort.101: initial particle distribution (4 columns) - z coordinate (m), del_gamma, ptcl. weight, dE/E0
fort.200: final current profile (3 columns) - bunch length (m), charge per cell, current (A)
fort.201: final particle distribution (4 columns) - z coordinate (m), del_gamma, ptcl. weight, dE/E0
