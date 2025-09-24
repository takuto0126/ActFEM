#!/bin/bash

rsync -avz -e ssh minami@10.35.22.56:/home/minami/ActFEM/Fwd_test/MTfwd_test_mpi/result_homo/ ./result_homo/