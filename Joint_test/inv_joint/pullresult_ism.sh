#!/bin/bash

rm result_inv/*
rsync -avz -e ssh 25cr05tm@ismxd:/home/25cr05tm/ActFEM/Joint_test/inv_joint/result_inv/ ./result_inv/