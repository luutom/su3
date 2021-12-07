#!/bin/bash
for MOM in "0 0 0" "1 0 0" 
do
  python3 ./rearrange_3bar_third.py --input=scratch/3bar --NT=64 --output=3bar_test --mom="${MOM}" 
done
