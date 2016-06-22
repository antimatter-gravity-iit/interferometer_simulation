# Change Log

Summary of changes made since Summer 2016. 

- Version indication is by date of last commit.

- Changes are summarized across sequences of commits, with the final commit from each sequence being identified by hash number.

## [2016-06-21] (https://github.com/lnevesabrantes/interferometer_simulation/commit/663a437f54a43c88441e7b6c4a204903ba36610d) (current)
Commit **663a437f54a43c88441e7b6c4a204903ba36610d**

**One-line summary:** partial rewrite of function *gp2*, merge and correction of *w* and *el* functions and addition of timers.

Files: **BeamParams.c**, **BeamParams.h**, **central_simulator.c**, **Gratings.c**, **Gratings.h**, **Misc.h**
- Merged former functions *w* (that calculated the beamwidth) and *el* (that calculated the coherence width) into a single C function called *calculate_width*. 
The formulas for these parameters (found in [McMorran 2008](http://journals.aps.org/pra/abstract/10.1103/PhysRevA.78.013601))
are the same, except for a multiplicative factor that represents the initial value of the width being calculated.
In the previous versions of the code, both widths were always assigned the same value (i.e., **the beam was taken 
to be perfectly coherent**), because the coefficient in the two functions was identical.
The new function contains an extra parameter, which is used as the coefficient in the physical formula,
so it returns the correct value for each parameter when correctly called.
- Updated function declarations and calls, and added explanatory comments. Notice that the formal arguments of function *gp1* weren't changed, but a comment on their nomenclature was included.

Files: **central_simulator.c**, **Gratings.c**, **Misc.h**
- Renamed variables:
  - w0 → initial_beamwidth
  - el0 → initial_coherence_width
  - r0 → initial_radius_of_wavefront_curvature

File: **central_simulator.c**, **Gratings.c**
- A timer was added inside each of the functions *gp0*, *gp1* and *gp2* to compute the time spent in each function call.
The program prints both that information and a message indicating which iteration of the various *for* loops
it is currently in.

File: **central_simulator.c**
- Added copyright notice to central simulator source code;
- Removed line *int col = 2* (unused variable);
- Rewrote values in scientific notation.

File: **Gratings.c**, **Gratings.h**
- Changed *gp0*, *gp1* and *gp2* function types to void (they do not return anything).

File: **Gratings.c**, **PhaseShifts.c**
- Replaced definition of *pi* by standard *math.h* library *M_PI* constant.

File: **Gratings.c**
- The sequence of *if* statements on lines 125 - 129 was rewritten as an if-else sequence;
- Line 269: equation rewritten to account for the case when *d1* doesn't equal *d2*, folowing [McMorran 2008](http://journals.aps.org/pra/abstract/10.1103/PhysRevA.78.013601) equation 18e;
- Equations in *gp2* were rewritten: six new variables were created to make the code more readable and
to be more consistent with [McMorran 2008](http://journals.aps.org/pra/abstract/10.1103/PhysRevA.78.013601).
They are *argument_d*, *argument_f*, *argument_p*, and *argument_v*, representing the arguments of the
exponential functions 18b, c, d, and e, respectively, according to the numbering in the article.
The two other variables (*function_d* and *function_v*) hold the value of the exponentials
in 18b and 18e (lines 280 and 281). The former implementation
had typos and included some obscure terms (i.e. *G2_x*) taken from the original IGOR PRO code, the origin and meaning of
which is unclear. **All these equations now match with the article**.

Files: various
- Standardized indentation and block comments;
- Added variable explanations.
 
## [2016-06-16](https://github.com/lnevesabrantes/interferometer_simulation/commit/9fede859df7ec4bcda48c86f521227b0ff8394de)
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
