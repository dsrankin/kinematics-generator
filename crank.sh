#!/bin/bash

for AMASS in {200..400..10}
do
  let "MAXCMASS = ${AMASS} - 100"
  echo "${MAXCMASS}"
  for CMASS in $(seq 20 10 ${MAXCMASS});
  do
    echo "Running for A mass ${AMASS} and C mass ${CMASS}..."
    cp test.cc.bak test.cc
    sed -i "s/    { 'A', 200. },/    { 'A', ${AMASS}. },/g" test.cc
    sed -i "s/    { 'C', 100. },/    { 'C', ${CMASS}. },/g" test.cc
    sed -i "s/test_Xq_Wqq_qqqq_M200/test_qX_qWY_qqqlv_X${AMASS}_Y${CMASS}/g" test.cc
    if ! [ -f ./test_qX_qWY_qqqlv_X${AMASS}_Y${CMASS}.root ]; then
      root -l -b -q test.cc+
    fi
  done
done
