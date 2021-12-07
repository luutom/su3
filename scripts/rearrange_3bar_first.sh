#!/bin/bash
for MOM in "0 0 0" "1 0 0" 
do
  python3 ./rearrange_3bar_first.py --input=scratch/3bar --NT=64 --output=3bar --mom="${MOM}" 
done
