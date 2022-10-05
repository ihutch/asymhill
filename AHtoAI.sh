#!/bin/bash
# Convert AsymHill.f90 to AsymIon.f90 by swapping the labelling suffixes.
sed -e 's/!de/!dq/g; s/!di/!dr/g; s/!dq/!di/g; s/!dr/!de/g '\
    AsymHill.f90 > AsymIon.f90
