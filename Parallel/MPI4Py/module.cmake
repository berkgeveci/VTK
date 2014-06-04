vtk_module(vtkParallelMPI4Py
  GROUPS
    MPI
  TCL_NAME # Not actually useful for tcl, but there's no "only wrap for X" option.
    vtkParallelMPIFourPy
  DEPENDS
    vtkParallelMPI
    vtkPython
  COMPILE_DEPENDS
    vtkmpi4py
  )
