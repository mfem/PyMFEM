#!/bin/bash
# [-it] interactive mode
# [-e DISPLAY and -v /tmp/.X11-unix:/tmp/.X11-unix] view GLVis through docker
docker run -it -e DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix pymfem:dev
