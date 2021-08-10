mpiexec -np 1 build/mcpartra RegressionTests/MonteCarlo3D_5Cycles.lua check_num_procs=false run_tape_basename=[["ZRunTape0"]] seed=0 > ZOut0.txt &
mpiexec -np 1 build/mcpartra RegressionTests/MonteCarlo3D_5Cycles.lua check_num_procs=false run_tape_basename=[["ZRunTape1"]] seed=1 > ZOut1.txt &
mpiexec -np 1 build/mcpartra RegressionTests/MonteCarlo3D_5Cycles.lua check_num_procs=false run_tape_basename=[["ZRunTape2"]] seed=2 > ZOut2.txt &
mpiexec -np 1 build/mcpartra RegressionTests/MonteCarlo3D_5Cycles.lua check_num_procs=false run_tape_basename=[["ZRunTape3"]] seed=3 > ZOut3.txt &