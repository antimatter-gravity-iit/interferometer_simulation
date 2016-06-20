# Change Log

Summary of changes made since Summer 2016. 

- Version indication is by date of last commit.

- Changes are summarized across sequences of commits, with the final commit from each sequence being identified by hash number.

## [2016-06-16](https://github.com/lnevesabrantes/interferometer_simulation/commit/9fede859df7ec4bcda48c86f521227b0ff8394de) (current)
Commit **9fede859df7ec4bcda48c86f521227b0ff8394de**

**One-line summary:** added info files, corrected math error in GSM simulation, and allowed choice of scale.

Files: **CHANGELOG.md**, **CONTRIBUTING.md**, **INSTALL**
- Added change log, coding conventions and installation instructions to the repository.

File: **README.md**
- Updated collaborator list and added reference to other info files.

File: **central_simulator.c**
- Commented a line that apparently forced usage of a regular scale in the simulation plot.

File: **Gratings.c**
- Corrected equations in lines 122 and 130 according to article in http://journals.aps.org/pra/abstract/10.1103/PhysRevA.78.013601 .

File: **Makefile**
- Changed standard ROOT path.

File: **PhaseShifts.c**
- Clarified indentation.

## [2016-06-14](https://github.com/lnevesabrantes/interferometer_simulation/commit/e518c780e3f1bf5c5165599272aecc56e7440afe)
Commit **e518c780e3f1bf5c5165599272aecc56e7440afe**

**One-line summary:** Makefile update; minor commentary and indentation edits.

File: **central_simulator.c**
- Updated collaborator list.

File: **Makefile**
- Added *-fext-numeric-literals* flag to Makefile targets. This fixed a compiler error on g++-5.
- Eliminated obsolete comments from a past implementation;
- Changed ROOT path and added other options for developer convenience;
- Added comment about path where ROOT is installed (must be altered according to the user's system);
- Created 3 new targets for debugging;
- Incorporated existing (but inactive) debug flags into new targets;
- Updated "clean" target;
- Updated target description.

File: **Misc.c**
- Added comment to ixgenerator function (copied from central_simulator.c). 

File: **PhaseShifts.c**
- Corrected mistake with commented variable initialization.

## [2016-06-07](https://github.com/lnevesabrantes/interferometer_simulation/commit/15041ff13e1f9e6d274aa64181e563c07ac30a0e) 
Commit **15041ff13e1f9e6d274aa64181e563c07ac30a0e**

**One-line summary:** extensive comment cleanup (with some additions); minor changes in control-flow structures.

Files: **central_simulator.c**, **BeamParams.c**, **BeamParams.h**, **Gratings.c**, **Gratings.h**, **Misc.c**, **Misc.h**, **PhaseShifts.c**
- Cleaned up most comments on files listed above. Removed obsolete implementations that survived as comments, stripped ancient personal commentary and obscure bits of unused code, and generally deleted anything that wasn't meant to be final;
- Block comments were also reformatted based on agreed-upon conventions.

File: **central_simulator.c**
- Significant reformatted the initial comment section. Updated collaborator list, moved various comments there and standardized the block comment structure;
- Added function descriptions by Isaac Gewarges to the code. These were previously on a separate spreadsheet;
- As suggested in the code itself, replaced *if* statement nested inside an *else* by standard *else-if* sequence;
- Re-standardized indent structure inside the main function.

File: **Misc.c**
- Rewrote the function *sinc* with a different *if-else* structure;
- Rewrote nested *if* sequence as a single 'if' block.   

File: **Misc.h**
- Removed commented variable initializations from parameter structure (actual initializations are on central_simulator.c).

## [2016-06-02](https://github.com/lnevesabrantes/interferometer_simulation/commit/2f25cd198c6595d801ff2745bce7926e9ff5da17) 
Commit **2f25cd198c6595d801ff2745bce7926e9ff5da17**

**One-line summary:** This is the code as it was before any modifications by the Summer 2016 team.
