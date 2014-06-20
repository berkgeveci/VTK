add_test(NAME ${vtk-module}Python-TestParallelNumPy
  COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
    $<TARGET_FILE:pvtkpython> ${CMAKE_CURRENT_SOURCE_DIR}/TestParallelNumPy.py
    ${VTK_MPI_POSTFLAGS})