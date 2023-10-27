Background
----------
The CMS has been a research and development area of the Coastal Inlets Research Program (CIRP) at the United States Army Corps of Engineers - Engineering Research and Development Center (USACE-ERDC), Coastal and Hydraulics Laboratory (CHL) since 2006. It was built from a group of numerical models that have been under development since 2002. Information on the CIRP and publications on the CMS can be found at the `CIRP Website <https://cirp.usace.army.mil>`_.

Key Features
------------
* Fully integrated system
	* CMS-Flow - Inline flow, sediment and salinity model
	* CMS-Wave - Spectral wave transformation model
* Finite Volume Method - Mass conservative
* Non-uniform Cartesian Grid - Easy to setup and run
* Quadtree Cartesian Grid - Flexible, efficient, and easier to generate than unstructured meshes
* Supports most common types of forcing and boundary conditions
* Robust numerical schemes for reliable, crash-free simulations
* Parallelization on desktop computers for fast computation
* User-friendly interface

System Components
-----------------

CMS-Flow
^^^^^^^^
CMS-Flow is a coupled hydrodynamic and sediment transport model capable of simulating depth-averaged circulation, salinity and sediment transport due to tides, wind and waves. The hydrodynamic model solves the conservative form of the shallow water equations and includes terms for the Coriolis force, wind stress, wave stress, bottom stress, vegetation flow drag, bottom and friction, and turbulent diffusion. 

There are three sediment transport models available in CMS: 

* a sediment mass balance model
* an equilibrium advection diffusion model
* a non-equilibrium advection-diffusion model. 

Salinity transport is simulated with the standard advection diffusion model and includes evaporation and precipitation. 

All equations are solved using the Finite Volume Method on a non-uniform or quadtree (telescoping) Cartesian grid. For additional information on CMS-Flow visit the `CMS-Flow Wiki Page <https://cirpwiki.info/wiki/CMS-Flow>`_.

CMS-Wave
^^^^^^^^
CMS-Wave is a spectral wave transformation model and solves the steady-state wave-action balance equation on a non-uniform Cartesian grid. It considers wind wave generation and growth, diffraction, reflection, dissipation due to bottom friction, whitecapping and breaking, wave-wave and wave-current interactions, wave runup, wave setup, and wave transmission through structures. For additional information information on CMS-Wave visit the `CMS-Wave Wiki Page <https://cirpwiki.info/wiki/CMS-Wave>`_.
