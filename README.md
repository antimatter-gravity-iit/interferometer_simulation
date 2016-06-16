# Interferometer simulation
Simulates phase shift of particles (due to various effects) moving through 
a 3-grating interferometer.

For installation instructions, refer to the INSTALL file. 
A summary of changes can be found in CHANGELOG.md. 


______________________
**Effects Accounted For:**
______________________
- Gravity

- VanDerWaals effect


______________________
**Authors:**
______________________
- Dr. Ben McMorran, Assistant Professor, University of Oregon

- Arthur Romero, IIT research student, summer 2015

- Adam Denchfield, IIT IPRO student, fall 2015

- Melanie Cornelius (nee Dooley), IIT IPRO and research student, fall 2014 
       through present

- Dr. Tom Roberts, IIT Research Professor of Physics

- Isaac Gewarges, IIT IPRO student, spring 2016

- Lucas Maia Rios, IIT research student, summer 2016

- Lucas Neves Abrantes, IIT research student, summer 2016

- Yuri Rossi Tonin, IIT research student, summer 2016

______________________
**How to Commit:**
______________________
- Place all source files in the source/ dir

- In general, keep functions pertaining to different effects in different files

- Name variables descriptively. If a variable pertains to a common physics term
       or symbol, place that symbol in a nearby comment

- Comment functions explaining the math in use 

- Add any new effect into this README


______________________
**General Git Etiquette:**
______________________
- **DO NOT commit executables** (.out files by default)

- DO NOT make dated folders in any directory except output/

- In general, avoid committing output files - the output directory is a 
       compromise between programming convention and needs of professors and 
       students in class.  Remember - you can always store your outputs in the 
       Google Drive or on iGroups.

- That said, DO NOT commit output files unless:

  - You have a good reason to track the output (IE, showing new functionality)

  - You have placed the output in a dated directory within the output/ dir

# [License](https://github.com/antimatter-gravity-iit/interferometer_simulation/blob/master/LICENSE)

   Copyright (C) 2016 Antimatter Gravity Interferometer Group, Illinois Institute of Technology (IIT). 

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
