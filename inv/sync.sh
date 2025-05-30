#!/bin/bash

rsync -avz -e ssh --exclude result_ms0/* ./ minami@10.35.22.56:/home/minami/ActFEMv1.0/inv/
#rsync -avz -e ssh ../inv2014_15_base minami@10.35.22.56:/home/minami/ActFEMv1.0/inv2014_15_base/
#rsync -avz -e ssh ../mesh_aso_A04/ minami@10.35.22.56:/home/minami/ActFEMv1.0/mesh_aso_A04/
#rsync -avz -e ssh ../src/ minami@10.35.22.56:/home/minami/ActFEMv1.0/src/