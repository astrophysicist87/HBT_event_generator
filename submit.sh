#! /usr/bin/env bash

./compile_all.sh &> compile_all.out &

wait

./driver.sh &> driver.out &
