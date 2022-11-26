# CLAS12Flow
Repo for representing MC::Lund collider events in a flowchart  
Uses Diagrams by @mingrammer https://github.com/mingrammer/diagrams
and hipopy by @mfmceneaney https://github.com/mfmceneaney/hipopy

Dependencies:
1. Diagrams
2. hipopy
3. shutil
4. numpy

## Use
1. Clone this repo
1. edit main.py
    1. Choose what particles you want in the endstate with userpid values
    1. Give the hipo file you want to use
1. run main.py
    1. This creates the copypythondiagram.py file
1. run FlowDev/copypythondiagram.py
    1. This creates the png of the diagram
1. The image of your diagram will now be in the repo directory
## Current Features
1. Reads MC::Lund bank from hipo files, extracts parent, daughter, pid, and index for each variable
1. Searches for specified hadron endstate (either exact counts or at least as many as specified)
1. Creates python file from template and adds definitions for particles in chosen event

1. Currently can create flowcharts for three types of collision events
    1. Quark-diquark result from collision
    1. Quark-Anti-quark-baryon result from collision
    1. Quark-diquark-meson result from collision
## Planned Features
1. More customization for event selection
1. Finishing python file building to create ready-to-run file that outputs flowchart png
1. Custom images for each particle
1. Adding kinematic data to flowchart for convenience
## Long-term Goals
1. Add support for data from other detectors (read in root ttrees and other files, not depend on clas12root)
1. Add module for REC::Particle bank;
