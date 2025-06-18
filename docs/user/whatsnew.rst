What's New
==========

v5.4.5 (18 June 2025)
-------------------------

Changes
^^^^^^^
* Bug fix - need to test for enabled sediment after reading in cards in prestart
* Enabled XMDF compilation on Linux. Updated workflow to compile CMS including XMDF libraries.
* Add canned test cases for XMDF that should work on Linux or Windows
* Moved third-party libraries from root folder to be under the 'external_libraries' folder.


v5.4.4.4 (19 May 2025)
-------------------------

Changes
^^^^^^^
* Minor changes to wave code to initialize a few variables.
* Added test to ensure initial mixing layer thickness is always less than thicknesses for bed layers.
* Added test for met station wind data import.
* Test fix for divide by zero in SHEARLUND routine.
* Added allocation for Observation and Nest points when using the newest format.
* (5.4.4.1) Bug fix for mixing layer thickness error.
* (5.4.4.2) Handles a '--non-interactive' flag as final argument to avoid the "Press any key to continue" message after errors.
* (5.4.4.3) Bug fix for array out of bounds when nesting is used.
* (5.4.4.4) Bug fix for mixing layer change in 5.4.4.1 when checking after reading in cards. 


v5.4.3 (18 March 2025)
-------------------------

Changes
^^^^^^^
* Optimized ADAK and ADBK dot product routines
* Found that time reference has been set to BCE instead of CE. Changing reference time to CE.


v5.4.2 (17 March 2025)
-------------------------

Changes
^^^^^^^
* Upgraded CMS-Wave variables to Dynamic allocation
* Implemented new extracted boundary condition for use with SMS 13.4+.
* Implemented partial CF compliance by adding HDF5 descriptive attributes for many solution types.
* Bug fix for deep water hard bottom values (deeper than 900m)

v5.4.0 (17 December 2024)
-------------------------

Changes
^^^^^^^
* Initial release for SMS 13.4
* Includes new implementations for CMS-Flow structures (Culverts, Weirs, Tide Gates)


v5.3.12.1 (10 December 2024)
----------------------------

Changes
^^^^^^^
* Bug fix for Wind Speed equal to 0.0 in some situations.


v5.3.12 (26 November 2024)
--------------------------

Changes
^^^^^^^
* Added new card to hot start a simulation, 'HOT_START_SIMULATION'.
* Added error message when either ORIGIN_X/Y card is missing.
* Fixed an issue where the CMS-Wave .wav file cannot be written to multiple times in some situations.
* Implemented new method for weirs, culverts, and tide gates in SMS 13.4+ and improved summary output for those structures.
* A few minor changes for compiling on Linux


v5.3.11 (19 August 2024)
------------------------

Changes
^^^^^^^
* Changed to alternate between recurring XMDF hot start files. ASCII files TBD.
* Added folder to start holding test cases for Continuous Integration implementation.


v5.3.10 (16 July 2024)
----------------------

Changes
^^^^^^^
* Added compile option for VS2022, Intel OneAPI 2024.x
* Fix for an erroneous Save Point warning that was always written.
* Improved status descrition update for Culverts in CMS.
* Added Error check for missing X/Y Origin cards in .cmcards file (SMS issue).
* Implemented n extra file to store Save Point cards when there are more than 10 defined.


v5.3.9 (16 April 2024)
----------------------

Changes
^^^^^^^
* Fix for Dredge module diagnostic file not found.
* Clarified some diagnostic output regarding simulation duration and computed residuals.
* Fix for steering variable (NHV_str to NH_str) possibly related to Flow/Wave steering issue.


v5.3.8 (08 January 2024) 
------------------------
 
Improvements
^^^^^^^^^^^^
* Update for internal path length for tidal database file locations.

Notice
^^^^^^
* Investigating an issue with internal mapping from UTM to Lat/Long in tidal database cases. 
  * Use State Plane if possible.


v5.3.7.1 (20 December 2023)
------------------------

Improvements
^^^^^^^^^^^^
* Began adding CF Compliant naming for HDF5 solution datasets (not yet implemented)
* Incremented CMS-Wave version number to 3.3. The code was modified in 2021 but the number remained the same.

Bug fixes
^^^^^^^^^
* Tidal boundary conditions where an offset was used. The offset was being applied twice which doubled the effect.
* Steering for ADCIRC tidal databases with waves.
* Some tidal database forcing issues.

Documentation
^^^^^^^^^^^^^
* Began adding new user documentation
* First version available as Open Source.


v5.3.6 (26 September 2023)
--------------------------

Improvements
^^^^^^^^^^^^
* Minor fixes and diagnostic output improvements for Rubble Mound Jetties with SMS 13.3+.

Bug fixes
^^^^^^^^^
*  Incorporate updates for explicit issues (Reed)


v5.3.5 (8 September 2023)
-------------------------

Notes
^^^^^
* Initial version of CMS released with SMS version 13.3.4 to public (out of beta).

Improvements
^^^^^^^^^^^^
* Improved summary of set up in diagnostic file and on screen.

Bug fixes
^^^^^^^^^
* Fix issues with Tidal Database operation with CMS. SMS 13.3.4+ will export correctly.
* Fixed issue with Grid Angle introduced with SMS 13.1.
* Minor GFortran fixes for compiling on some flavors of Linux Gnu compiler. 


v5.3.4 (16 May 2023)
--------------------

Improvements
^^^^^^^^^^^^
* Added a warning if user-specified boundary angle was too different than internally calculated angle for that boundary.

Bug fixes
^^^^^^^^^
* Fixed missing zero output time in files when increment was more than 100 hours.
* Fixed issue with infinity in certain situations within CMS-Wave GSR solver.


v5.3.3 (7 April 2023)
---------------------

Improvements
^^^^^^^^^^^^
* Split third-party code (spatial and datum transformations) out as separate linkable libraries for Open Source requirement.
* Implemented two new ADCIRC tidal databases (EC2015 and ENPAC2015).


v5.3.2 (11 August 2022)
-----------------------

Bug fixes
^^^^^^^^^
* Minor change to reading parameter file for Explicit scheme.
* Fixed array allocation error when choosing multi-sediment with D35,D50,D90.
* Fixed lookup-table issue in 'bs_init' routine.
* Fix in Tools (Option 4) when merging datasets written by SMS instead of CMS solution datasets.


v5.3.1 (7 July 2022)
--------------------

Improvements
^^^^^^^^^^^^
* Updated the array upper limit for several CMS-Wave variables.

Bug fixes
^^^^^^^^^
* Change to Weir Structure specification cards for integration into SMS 13.2+.
  * Add option to specify Cell IDs in the same manner as for Rubble Mound Jetties.


v5.3.0 (16 May 2022)
--------------------

Notes
^^^^^
* Initial version of CMS released with SMS version 13.2.

Improvements
^^^^^^^^^^^^
* Implemented C2Shore as a new option for Cross-shore sediment transport.
  * Requires CMS to be run with both waves and flow and will fail otherwise.


Previous changes 
----------------

See CMS-Releases on `CIRP Wiki <https://cirpwiki.info/wiki/CMS_Releases>`_.
