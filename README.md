# CLAS12Flow
Repo for representing MC::Lund collider events in a flowchart
Uses Diagrams by @mingrammer https://github.com/mingrammer/diagrams
## Current Features
1. Reads MC::Lund bank from hipo files, extracts parent, daughter, pid, and index for each variable
1. Searches for specified hadron endstate (either exact counts or at least as many as specified)
1. Creates python file from template and adds definitions for particles in chosen event
## Planned Features
1. More customization for event selection
1. Finishing python file building to create ready-to-run file that outputs flowchart png
1. Custom images for each particle
1. Adding kinematic data to flowchart for convenience
## Long-term Goals
1. Add support for data from other detectors (read in root ttrees and other files, not depend on clas12root)
1. Add module for REC::Particle bank;
