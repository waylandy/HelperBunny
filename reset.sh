#!/bin/bash
for d in `find . -name "__pycache__"`; do rm -r $d; done
