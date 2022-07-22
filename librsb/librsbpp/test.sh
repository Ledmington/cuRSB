#!/bin/bash
set -e
set -x
#OMP_NUM_THREADS=2 ./rsbpp bayer02.mtx
#OMP_NUM_THREADS=2 ./rsbpp raefsky4.mtx
#OMP_NUM_THREADS=2 zcat raefsky3.mtx.gz | ./rsbpp -
#OMP_NUM_THREADS=2 ./rsbpp raefsky4.mtx
test $(./rsbpp Td,s G.mtx | grep Z-sort | wc -l ) = 54
test $(./rsbpp Td   G.mtx | grep Z-sort | wc -l ) = 27
test $(./rsbpp Td,z     G.mtx | grep Z-sort | wc -l ) = 54
test $(./rsbpp vTd,z G.mtx | grep Z-sort | wc -l ) = 54
test $(./rsbpp vTd,z  G.mtx | grep Z-sort | wc -l ) = 54
test $(./rsbpp vvvTd,z  G.mtx | grep Zorted | wc -l ) = 8
test $(./rsbpp vvTd,z  G.mtx | grep Z-sort | wc -l ) = 54
test $(./rsbpp vvTd,z  G.mtx | grep Range  | wc -l ) = 0
test $(./rsbpp vvvTd,z G.mtx | grep Range  | wc -l ) -gt 0
test $(./rsbpp vvvTd,z S.mtx | grep Range  | wc -l ) -eq 0
test $(./rsbpp vvvTd,z G.mtx | grep Range  | wc -l ) = 258
test $(OMP_NUM_THREADS=1 ./rsbpp m10M10I1r1,4,8sFv | grep spmm- | wc -l ) = 9
test $(OMP_NUM_THREADS=1 ./rsbpp    C1000m100M100I1r1,4,8sFv | grep spmm- | wc -l ) = 9
test $(OMP_NUM_THREADS=1 ./rsbpp    C1000m100M100I1r1sFvtN,T | grep spmm- | wc -l ) = 3
test $(OMP_NUM_THREADS=1 ./rsbpp    C1000m100M100I1r1vtN,TsF | grep spmm- | wc -l ) = 2
test $(OMP_NUM_THREADS=1 ./rsbpp    C1000m100M100I1r0vtN,TsF | grep spmm- | wc -l ) = 0
test $(OMP_NUM_THREADS=1 RSB_NUM_THREADS=1 ./rsbpp vvvC1000m100M100I1r1vtN,TorsF | grep Recursing | wc -l ) = 4
test $(OMP_NUM_THREADS=2 RSB_NUM_THREADS=2 ./rsbpp vvvC1000m100M100I1r1vtN,TorsF | grep Recursing | wc -l ) = 4
if ! test $(OMP_NUM_THREADS=1 RSB_NUM_THREADS=1 ./rsbpp vvvC1000m100M100I1r1vtN,ToRsF | grep Recursing | wc -l ) = 208 ; then echo "Note: testing oR is not yet deterministic"; fi
if ! test $(OMP_NUM_THREADS=2 RSB_NUM_THREADS=2 ./rsbpp vvvC1000m100M100I1r1vtN,ToRsF | grep Recursing | wc -l ) = 410 ; then echo "Note: testing oR is not yet deterministic"; fi
