""" Test for various numpy_interface modules. Main goal is to test
parallel algorithms in vtk.numpy_interface.algorithms."""

import sys
try:
    import numpy
except ImportError:
    print "Numpy (http://numpy.scipy.org) not found.",
    print "This test requires numpy!"
    sys.exit(0)

import vtk
from vtk.test import Testing
import vtk.numpy_interface.dataset_adapter as dsa
import vtk.numpy_interface.algorithms as algs
from mpi4py import MPI

c = vtk.vtkMPIController()
#c.SetGlobalController(None)
rank = c.GetLocalProcessId()
size = c.GetNumberOfProcesses()

def PRINT(text, values):
    res = numpy.array(values, dtype=numpy.float64)
    MPI.COMM_WORLD.Allreduce([values, MPI.DOUBLE], [res, MPI.DOUBLE], MPI.SUM)
    assert numpy.abs(numpy.sum(res)) < 1E-5
    if rank == 0:
        print text, numpy.sum(res)

def testArrays(rtData, rtData2, grad, grad2, total_npts):
    " Test various parallel algorithms."
    if rank == 0:
        print '-----------------------'
    PRINT( "SUM ones:", algs.sum(rtData / rtData) - total_npts )

    PRINT( "SUM sin:", (algs.sum(algs.sin(rtData) + 1) - numpy.sum(numpy.sin(rtData2) + 1)) / numpy.sum(numpy.sin(rtData2) + 1) )

    PRINT( "rtData min:", algs.min(rtData) - numpy.min(rtData2) )
    PRINT( "rtData max:", algs.max(rtData) - numpy.max(rtData2) )
    PRINT( "rtData sum:", (algs.sum(rtData) - numpy.sum(rtData2)) / (2*numpy.sum(rtData2)) )
    PRINT( "rtData mean:", (algs.mean(rtData) - numpy.mean(rtData2)) / (2*numpy.mean(rtData2)) )
    PRINT( "rtData var:", (algs.var(rtData) - numpy.var(rtData2)) / numpy.var(rtData2) )
    PRINT( "rtData std:", (algs.std(rtData) - numpy.std(rtData2)) / numpy.std(rtData2) )

    PRINT( "grad min:", algs.min(grad) - numpy.min(grad2) )
    PRINT( "grad max:", algs.max(grad) - numpy.max(grad2) )
    PRINT( "grad min 0:", algs.min(grad, 0) - numpy.min(grad2, 0) )
    PRINT( "grad max 0:", algs.max(grad, 0) - numpy.max(grad2, 0) )
    PRINT( "grad min 1:", algs.sum(algs.min(grad, 1)) - numpy.sum(numpy.min(grad2, 1)) )
    PRINT( "grad max 1:", algs.sum(algs.max(grad, 1)) - numpy.sum(numpy.max(grad2, 1)) )
    PRINT( "grad sum 1:", algs.sum(algs.sum(grad, 1)) - numpy.sum(numpy.sum(grad2, 1)) )
    PRINT( "grad var:", (algs.var(grad) - numpy.var(grad2)) / numpy.var(grad2) )
    PRINT( "grad var 0:", (algs.var(grad, 0) - numpy.var(grad2, 0)) / numpy.var(grad2, 0) )

w = vtk.vtkRTAnalyticSource()
w.UpdateInformation()
# Update with ghost level because gradient needs it
# to be piece independent
w.SetUpdateExtent(rank, size, 1)
w.Update()

# The parallel arrays that we care about
ds = dsa.WrapDataObject(w.GetOutput())
rtData = ds.PointData['RTData']
grad = algs.gradient(rtData)
ds.PointData.append(grad, 'gradient')

# Crop the any ghost points out
org_ext = w.GetOutput().GetExtent()
ext = list(org_ext)
wext = w.GetOutputInformation(0).Get(vtk.vtkStreamingDemandDrivenPipeline.WHOLE_EXTENT())
for i in range(3):
    if ext[2*i] != wext[2*i]:
        ext[2*i] = ext[2*i] + 2
    if ext[2*i+1] != wext[2*i+1]:
        ext[2*i+1] = ext[2*i+1] - 1
if ext != list(org_ext):
    w.GetOutput().Crop(ext)

# Croppped arrays
rtData = ds.PointData['RTData']
grad = ds.PointData['gradient']

# The whole dataset so that we can compare
# against parallel algorithms.
w2 = vtk.vtkRTAnalyticSource()
w2.Update()

ds2 = dsa.WrapDataObject(w2.GetOutput())
rtData2 = ds2.PointData['RTData']
grad2 = algs.gradient(rtData2)

npts = numpy.array(numpy.int32(ds.GetNumberOfPoints()))
total_npts = numpy.array(npts)
MPI.COMM_WORLD.Allreduce([npts, MPI.INT], [total_npts, MPI.INT], MPI.SUM)

# Test simple distributed data.
testArrays(rtData, rtData2, grad, grad2, total_npts)

# Check that we can disable parallelism by using a dummy controller
# even when a global controller is set
assert algs.sum(rtData / rtData, controller=vtk.vtkDummyController()) != total_npts

# Test where arrays are NoneArray on one of the ranks.
if size > 1:
    if rank == 0:
        rtData3 = rtData2
        grad3 = grad2
    else:
        rtData3 = dsa.NoneArray
        grad3 = dsa.NoneArray

    testArrays(rtData3, rtData2, grad3, grad2, total_npts)

# Test composite arrays
rtData3 = dsa.VTKCompositeDataArray([rtData, dsa.NoneArray])
grad3 = dsa.VTKCompositeDataArray([dsa.NoneArray, grad])

testArrays(rtData3, rtData2, grad3, grad2, total_npts)

# Test where arrays are NoneArray on one of the ranks
# and composite on others.
if size > 1:
    if rank == 1:
        rtData3 = dsa.VTKCompositeDataArray([rtData2])
        grad3 = dsa.VTKCompositeDataArray([grad2])
    else:
        rtData3 = dsa.NoneArray
        grad3 = dsa.NoneArray

    testArrays(rtData3, rtData2, grad3, grad2, total_npts)

# Test composite arrays with multiple blocks.

# Split the local image to 2.
datasets = []
for i in range(2):
    image = vtk.vtkImageData()
    image.ShallowCopy(w.GetOutput())
    t = vtk.vtkExtentTranslator()
    wext = image.GetExtent()
    t.SetWholeExtent(wext)
    t.SetPiece(i)
    t.SetNumberOfPieces(2)
    t.PieceToExtent()
    ext = list(t.GetExtent())

    # Crop the any ghost points out
    for i in range(3):
        if ext[2*i] != wext[2*i]:
            ext[2*i] = ext[2*i] + 1
    if ext != list(org_ext):
        image.Crop(ext)

    datasets.append(dsa.WrapDataObject(image))

rtData3 = dsa.VTKCompositeDataArray([datasets[0].PointData['RTData'], datasets[1].PointData['RTData']])
grad3 = dsa.VTKCompositeDataArray([datasets[0].PointData['gradient'], datasets[1].PointData['gradient']])

testArrays(rtData3, rtData2, grad3, grad2, total_npts)
