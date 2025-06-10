#!/bin/bash

# for eic
#rsync -avz -e ssh tminami@eic.eri.u-tokyo.ac.jp:/home/tminami/volcano/FEM_edge/FEM_edge_Nakadake/inv_IUGG/result_m0/ result_m0/

#rsync -avz -e ssh tminami@eic.eri.u-tokyo.ac.jp:/home/tminami/volcano/FEM_edge/FEM_edge_Nakadake/inv_IUGG/result_ms0/ result_ms0/

#rsync -avz -e ssh tminami@eic.eri.u-tokyo.ac.jp:/home/tminami/volcano/FEM_edge/FEM_edge_Nakadake/inv_IUGG/result_fwd/ result_fwd/

rsync -avz -e "ssh -p 10122" ishibashi@130.54.189.251:/home/ishibashi/ActFEMv1.0/inv/result_ms0/ result_ms0/


