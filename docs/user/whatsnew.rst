What's New
==========

v5.4.8.1 (20 October 2026)
-------------------------

Changes
^^^^^^^
* Integrated changes from implicit developer

  * Implement sediment influx (Cohesive/Non-cohesive mixtures) on river boundary
  * A few OpenMP directive improvements
  * A few improvements for stability (relax(8) in flow.f90 amd other)
  * Eliminate oscillations in deep water when using Lund-CIRP formula

* Test each sediment percentile dataset when reading cards to ensure all non-negative values.
* Added function to return specific messages based on XMDF error codes.
* Fix a few issues reading cards with Explicit solution scheme
* Improve diagnostic output for in_xmdf_lib calls.
* Fix for floating overflow in LundCirp sediment transport code for currents.
* (5.4.8.1) Bug fix for Hardbottom checks when Rubble Mound cells specified and fix garbled output in status.

v5.4.7.3 (09 September 2026)
-------------------------

Changes
^^^^^^^
* Change fractional dataset names from "Datasets/Fraction_01 (1)" to "Datasets/Fraction_S1_L1" for updates in SMS for sediment mapping for readability. Retain previous format for backward compatibility
* Initialize some output logical variables related to sediment mapping.
* (5.4.7.1) Added functionality to write new naming convention for pertinent output datasets.
* (5.4.7.1) Added CF Compliance attributes to output datasets written bed layer thickness, bed layer fractions, and sediment percentiles
* (5.4.7.1) Now writing CF compliance attributes to each solution's root group instead of "PROPRERTIES" group.
* (5.4.7.1) Stop writing CF attributes for vector datasets.  Will implement on separated scalars at a later date.
* (5.4.7.1) Added new function to test a character to see if it is a numeric digit.
* (5.4.7.2) Corrected a few APPEND statements for compiler that is now more strict.
* (5.4.7.2) Added vs2026 solution
* (5.4.7.3) Added missing DEALLOCATE statements in wave code to potentially fix a 'insufficient virtual memory' issue.


v5.4.6 (03 July 2025)
-------------------------

Changes
^^^^^^^
* If hardbottom dataset seems inverted (elevations instead of depths), give warning and write out cells unless there are > 250.
* Fix potential rounding issue on hard bottom check if above Z.  Tel file keeps only 3 significant digits for Z but HB has more.

v5.4.5 (23 June 2025)
-------------------------

Changes
^^^^^^^
* Bug fix - need to test for enabled sediment after reading in cards in prestart
* Enabled XMDF compilation on Linux. Updated workflow to compile CMS including XMDF libraries.
* Add canned test cases for XMDF that should work on Linux or Windows
* Moved third-party libraries from root folder to be under the 'external_libraries' folder.
* Modified CI for CMS-Flow and added a few CMS-Wave tests (still to be added to CI)
* Add extra print line for Hot Start indicating the starting date/time.
* Bug fix in steering for printing out correct timestamp when it has only 8 digits instead of 12.
* Ignore certain Advanced cards if they are already handled. Print warning message.

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


Version 5.3.x
-------------

The What's New list for version 5.3.x is in the 5.3 documentation: `5.3 What's New <https://cms2d.readthedocs.io/en/cms_v5.3/user/whatsnew.html>`_ (also available from the version menu).


Previous changes 
----------------

See CMS-Releases on `CIRP Wiki <https://cirpwiki.info/wiki/CMS_Releases>`_.
