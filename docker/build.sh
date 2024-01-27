#!/bin/bash
# [--no-cache] Will take significantly longer but can help with build issues
docker build -t pymfem:dev -f Dockerfile .
