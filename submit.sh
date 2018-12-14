#! /usr/bin/env bash

( ./compile_all.sh &> compile_all.out && ./driver.sh &> driver.out ) &
