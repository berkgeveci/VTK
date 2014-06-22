import dataset_adapter as dsa
import internal_algorithms as algs
import itertools
import numpy
try:
    from vtk.vtkParallelCore import vtkMultiProcessController
    from vtk.vtkParallelCore import vtkDummyController
    from vtk.vtkParallelMPI4Py import vtkMPI4PyCommunicator
except ImportError:
    vtkMultiProcessController = None
    vtkMPI4PyCommunicator = None

def apply_func2(func, array, args):
    if array is dsa.NoneArray:
        return []
    res = []
    for a in array.Arrays:
        res.append(func(a, *args))
    return res

def apply_func(func, array, args):
    return dsa.VTKCompositeDataArray(apply_func2(func, array, args), dataset = array.DataSet)

def make_ufunc(ufunc):
    def new_ufunc(array):
        if type(array) == dsa.VTKCompositeDataArray:
            res = []
            for a in array.Arrays:
                if a is dsa.NoneArray:
                    res.append(dsa.NoneArray)
                else:
                    res.append(ufunc(a))
            return dsa.VTKCompositeDataArray(res, dataset = array.DataSet)
        elif array is dsa.NoneArray:
            return dsa.NoneArray
        else:
            return ufunc(array)
    return new_ufunc

def make_dfunc(dfunc):
    def new_dfunc(array1, val2):
        if type(array1) == dsa.VTKCompositeDataArray and type(val2) == dsa.VTKCompositeDataArray:
            res = []
            for a1, a2 in itertools.izip(array1.Arrays, val2.Arrays):
                if a1 is dsa.NoneArray or a2 is dsa.NoneArray:
                    res.append(dsa.NoneArray)
                else:
                    res.append(dfunc(a1, a2))
            return dsa.VTKCompositeDataArray(res, dataset = array1.DataSet)
        elif type(array1) == dsa.VTKCompositeDataArray:
            res = []
            for a in array1.Arrays :
                if a is dsa.NoneArray:
                    res.append(dsa.NoneArray)
                else:
                    res.append(dfunc(a, val2))
            return dsa.VTKCompositeDataArray(res, dataset = array1.DataSet)
        elif array1 is dsa.NoneArray:
            return dsa.NoneArray
        else:
            return dfunc(array1, val2)
    return new_dfunc

def make_dsfunc(dsfunc):
    def new_dsfunc(array, ds=None):
        if type(array) == dsa.VTKCompositeDataArray:
            res = []
            for a in array.Arrays:
                if a is dsa.NoneArray:
                    res.append(dsa.NoneArray)
                else:
                    res.append(dsfunc(a, ds))
            return dsa.VTKCompositeDataArray(res, dataset = array.DataSet)
        elif array is dsa.NoneArray:
            return dsa.NoneArray
        else:
            return dsfunc(array, ds)
    return new_dsfunc

def make_dsfunc2(dsfunc):
    def new_dsfunc2(ds):
        if type(ds) == dsa.CompositeDataSet:
            res = []
            for dataset in ds:
                res.append(dsfunc(dataset))
            return dsa.VTKCompositeDataArray(res, dataset = dataset)
        else:
            return dsfunc(ds)
    return new_dsfunc2

def _lookup_mpi_type(ntype):
    from mpi4py import MPI
    if ntype == numpy.float64:
        return MPI.DOUBLE
    elif ntype == numpy.bool:
        return MPI.BOOL
    else:
        raise ValueError

def global_func(impl, array, axis, controller):
    if type(array) == dsa.VTKCompositeDataArray:
        if axis is None or axis == 0:
            res = impl.serial_composite(array, axis)
        else:
            res = apply_func(impl.op(), array, (axis,))
    else:
        res = impl.op()(array, axis)
        if res is not dsa.NoneArray:
            res = res.astype(numpy.float64)

    if axis is None or axis == 0:
        if controller is None and vtkMultiProcessController is not None:
            controller = vtkMultiProcessController.GetGlobalController()
        if controller and controller.IsA("vtkMPIController"):
            from mpi4py import MPI
            comm = vtkMPI4PyCommunicator.ConvertToPython(controller.GetCommunicator())

            ntuples = numpy.int32(0)
            if type(res) == dsa.VTKArray or type(res) == numpy.ndarray:
                shp = shape(res)
                if len(shp) == 0:
                    ntuples = numpy.int32(1)
                else:
                    ntuples = numpy.int32(res.shape[0])
            elif type(res) == numpy.float64:
                ntuples = numpy.int32(1)
            max_tuples = numpy.array(ntuples, dtype=numpy.int32)
            comm.Allreduce([ntuples, MPI.INT], [max_tuples, MPI.INT], MPI.MAX)

            if res is dsa.NoneArray:
                if max_tuples == 1:
                    # Weird trick to make the array look like a scalar
                    max_tuples = ()
                res = impl.default(max_tuples)

            res_recv = numpy.array(res)
            mpi_type = _lookup_mpi_type(res.dtype)
            comm.Allreduce([res, mpi_type], [res_recv, mpi_type], impl.mpi_op())
            res = res_recv

    return res

def sum(array, axis=None, controller=None):
    class SumImpl:
        def op(self):
            return algs.sum

        def mpi_op(self):
            from mpi4py import MPI
            return MPI.SUM

        def serial_composite(self, array, axis):
            res = None
            arrays = array.Arrays
            for a in arrays:
                if a is not dsa.NoneArray:
                    if res is None:
                        res = algs.sum(a, axis).astype(numpy.float64)
                    else:
                        res += algs.sum(a, axis)
            return res

        def default(self, max_tuples):
            return numpy.zeros(max_tuples, dtype=numpy.float64)

    return global_func(SumImpl(), array, axis, controller)

def max(array, axis=None, controller=None):
    class MaxImpl:
        def op(self):
            return algs.max

        def mpi_op(self):
            from mpi4py import MPI
            return MPI.MAX

        def serial_composite(self, array, axis):
            res = apply_func2(algs.max, array, (axis,))
            clean_list = []
            for a in res:
                if a is not dsa.NoneArray:
                    clean_list.append(a)
            if clean_list is []:
                return None
            return algs.max(clean_list, axis=0).astype(numpy.float64)

        def default(self, max_tuples):
            return numpy.ones(max_tuples, dtype=numpy.float64) * numpy.finfo(numpy.float64).min

    return global_func(MaxImpl(), array, axis, controller)

def min(array, axis=None, controller=None):
    class MinImpl:
        def op(self):
            return algs.min

        def mpi_op(self):
            from mpi4py import MPI
            return MPI.MIN

        def serial_composite(self, array, axis):
            res = apply_func2(algs.min, array, (axis,))
            clean_list = []
            for a in res:
                if a is not dsa.NoneArray:
                    clean_list.append(a)
            if clean_list is []:
                return None
            return algs.min(clean_list, axis=0).astype(numpy.float64)

        def default(self, max_tuples):
            return numpy.ones(max_tuples, dtype=numpy.float64) * numpy.finfo(numpy.float64).max

    return global_func(MinImpl(), array, axis, controller)

def max_per_block(array, axis=None, controller=None):
    t = type(array)
    if t == dsa.VTKArray or t == numpy.ndarray:
        return max(array, axis, controller)
    elif array is dsa.NoneArray:
        return dsa.NoneArray

    maxes = apply_func2(max, array, (axis, vtkDummyController()))

    if controller is None and vtkMultiProcessController is not None:
        controller = vtkMultiProcessController.GetGlobalController()
    if controller and controller.IsA("vtkMPIController"):
        from mpi4py import MPI
        comm = vtkMPI4PyCommunicator.ConvertToPython(controller.GetCommunicator())

        for res in maxes:
            if res is not dsa.NoneArray:
                break

        ntuples = numpy.int32(0)
        if type(res) == dsa.VTKArray or type(res) == numpy.ndarray:
            shp = shape(res)
            if len(shp) == 0:
                ntuples = numpy.int32(1)
            else:
                ntuples = numpy.int32(res.shape[0])
        elif type(res) == numpy.float64:
            ntuples = numpy.int32(1)
        max_tuples = numpy.array(ntuples, dtype=numpy.int32)
        comm.Allreduce([ntuples, MPI.INT], [max_tuples, MPI.INT], MPI.MAX)

        # Get all ids from dataset, including empty ones.
        it = array.DataSet.NewIterator()
        it.UnRegister(None)
        it.SetSkipEmptyNodes(False)
        ids = []
        max_id = 0
        while not it.IsDoneWithTraversal():
            _id = it.GetCurrentFlatIndex()
            max_id = numpy.max((max_id, _id))
            if it.GetCurrentDataObject() is not None:
                ids.append(_id)
            it.GoToNextItem()

        if max_id == 0:
            return dsa.VTKCompositeDataArray(maxes, dataset=array.DataSet)

        has_ids = numpy.zeros(max_id+1, dtype=numpy.int32)
        for _id in ids:
            has_ids[_id] = 1
        id_cout = numpy.array(has_ids)
        comm.Allreduce([has_ids, MPI.INT], [id_cout, MPI.INT], MPI.SUM)

        reduce_ids = []
        for _id in ids:
            if id_cout[_id] > 1:
                reduce_ids.append(_id)

        to_reduce = len(reduce_ids)
        if to_reduce == 0:
            return dsa.VTKCompositeDataArray(maxes, dataset=array.DataSet)

        lmaxes = numpy.empty(max_tuples*to_reduce)
        lmaxes.fill(numpy.finfo(numpy.float64).min)

        # Just get non-empty ids. Doing this again in case
        # the traversal above results in a different order.
        # We need the same order since we'll use izip below.
        it = array.DataSet.NewIterator()
        it.UnRegister(None)
        ids = []
        while not it.IsDoneWithTraversal():
            ids.append(it.GetCurrentFlatIndex())
            it.GoToNextItem()

        for _id, _max in itertools.izip(ids, maxes):
            try:
                loc = reduce_ids.index(_id)
                lmaxes[loc*max_tuples:(loc+1)*max_tuples] = _max
            except ValueError:
                pass

        rmaxes = numpy.array(lmaxes)
        comm.Allreduce([lmaxes, MPI.DOUBLE], [rmaxes, MPI.DOUBLE], MPI.MAX)

        for i in xrange(to_reduce):
            _id = reduce_ids[i]
            try:
                loc = ids.index(_id)
                a = maxes[loc]
                if len(a.shape) == 0:
                    maxes[loc] = dsa.VTKArray(rmaxes[i])
                else:
                    maxes[loc][:] = rmaxes[i*max_tuples:(i+1)*max_tuples]
            except ValueError:
                pass

    return dsa.VTKCompositeDataArray(maxes, dataset=array.DataSet)

def min_per_block(array, axis=None):
    l = apply_func2(min, array, (axis, vtkDummyController()))
    return dsa.VTKCompositeDataArray(l, dataset=array.DataSet)

def all(array, axis=None, controller=None):
    class MinImpl:
        def op(self):
            return algs.all

        def mpi_op(self):
            from mpi4py import MPI
            return MPI.LAND

        def serial_composite(self, array, axis):
            res = apply_func2(algs.all, array, (axis,))
            clean_list = []
            for a in res:
                if a is not dsa.NoneArray:
                    clean_list.append(a)
            if clean_list is []:
                return None
            return algs.all(clean_list, axis=0)

        def default(self, max_tuples):
            return numpy.ones(max_tuples, dtype=numpy.bool)

    return global_func(MinImpl(), array, axis, controller)

def _array_count(array, axis, controller):

    if array is dsa.NoneArray:
        size = numpy.int64(0)
    elif axis is None:
        size = numpy.int64(array.size)
    else:
        size = numpy.int64(shape(array)[0])

    if controller is None and vtkMultiProcessController is not None:
        controller = vtkMultiProcessController.GetGlobalController()

    if controller and controller.IsA("vtkMPIController"):
        from mpi4py import MPI
        comm = vtkMPI4PyCommunicator.ConvertToPython(controller.GetCommunicator())

        total_size = numpy.array(size, dtype=numpy.int64)
        comm.Allreduce([size, MPI.INT64_T], [total_size, MPI.INT64_T], MPI.SUM)
        size = total_size

    return size

def mean(array, axis=None, controller=None, size=None):

    if axis is None or axis == 0:
        if size is None:
            size = _array_count(array, axis, controller)
        return sum(array, axis) / size
    else:
        if type(array) == dsa.VTKCompositeDataArray:
            return apply_func(algs.mean, array, (axis,))
        else:
            return algs.mean(array, axis)

def var(array, axis=None, controller=None):

    if axis is None or axis == 0:
        size = _array_count(array, axis, controller)
        tmp = array - mean(array, axis, controller, size)
        return sum(tmp*tmp, axis, controller) / size
    else:
        if type(array) == dsa.VTKCompositeDataArray:
            return apply_func(algs.var, array, (axis,))
        else:
            return algs.var(array, axis)

def std(array, axis=None, controller=None):
    return sqrt(var(array, axis, controller))

def shape(array):
    if type(array) == dsa.VTKCompositeDataArray:
        shp = None
        for a in array.Arrays:
            if a is not dsa.NoneArray:
                if shp is None:
                    shp = list(a.shape)
                else:
                    tmp = a.shape
                    if (len(shp) != len(tmp)):
                        raise ValueError, "Expected arrays of same shape"
                    shp[0] += tmp[0]
                    for idx in range(1,len(tmp)):
                        if shp[idx] != tmp[idx]:
                            raise ValueError, "Expected arrays of same shape"
        return tuple(shp)
    elif array is dsa.NoneArray:
        return ()
    else:
        return numpy.shape(array)

def make_vector(arrayx, arrayy, arrayz=None):
    if type(arrayx) == dsa.VTKCompositeDataArray and type(arrayy) == dsa.VTKCompositeDataArray and (type(arrayz) == dsa.VTKCompositeDataArray or arrayz is None):
        res = []
        if arrayz is None:
            for ax, ay in itertools.izip(arrayx.Arrays, arrayy.Arrays):
                if ax is not dsa.NoneArray and ay is not dsa.NoneArray:
                    res.append(algs.make_vector(ax, ay))
                else:
                    res.append(dsa.NoneArray)
        else:
            for ax, ay, az in itertools.izip(arrayx.Arrays, arrayy.Arrays, arrayz.Arrays):
                if ax is not dsa.NoneArray and ay is not dsa.NoneArray and az is not dsa.NoneArray:
                    res.append(algs.make_vector(ax, ay, az))
                else:
                    res.append(dsa.NoneArray)
        return dsa.VTKCompositeDataArray(res, dataset = arrayx.DataSet)
    else:
        return algs.make_vector(arrayx, arrayy, arrayz)

sqrt = make_ufunc(numpy.sqrt)
exp = make_ufunc(numpy.exp)
floor = make_ufunc(numpy.floor)
ceil = make_ufunc(numpy.ceil)
round = make_ufunc(numpy.round)
sin = make_ufunc(numpy.sin)
cos = make_ufunc(numpy.cos)
tan = make_ufunc(numpy.tan)
arcsin = make_ufunc(numpy.arcsin)
arccos = make_ufunc(numpy.arccos)
arctan = make_ufunc(numpy.arctan)
sinh = make_ufunc(numpy.sinh)
cosh = make_ufunc(numpy.cosh)
tanh = make_ufunc(numpy.tanh)
arcsinh = make_ufunc(numpy.arcsinh)
arccosh = make_ufunc(numpy.arccosh)
arctanh = make_ufunc(numpy.arctanh)
where = make_ufunc(numpy.where)
expand_dims = make_dfunc(numpy.expand_dims)

abs = make_ufunc(algs.abs)
area = make_dsfunc2(algs.area)
aspect = make_dsfunc2(algs.aspect)
aspect_gamma = make_dsfunc2(algs.aspect_gamma)
condition= make_dsfunc2(algs.condition)
cross = make_dfunc(algs.cross)
curl = make_dsfunc(algs.curl)
divergence = make_dsfunc(algs.divergence)
det = make_ufunc(algs.det)
determinant = make_ufunc(algs.determinant)
diagonal = make_dsfunc2(algs.diagonal)
dot = make_dfunc(algs.dot)
eigenvalue = make_ufunc(algs.eigenvalue)
eigenvector = make_ufunc(algs.eigenvector)
gradient = make_dsfunc(algs.gradient)
inv = make_ufunc(algs.inv)
inverse = make_ufunc(algs.inverse)
jacobian = make_dsfunc2(algs.jacobian)
laplacian = make_dsfunc(algs.laplacian)
ln = make_ufunc(algs.ln)
log = make_ufunc(algs.log)
log10 = make_ufunc(algs.log10)
max_angle = make_dsfunc2(algs.max_angle)
mag = make_ufunc(algs.mag)
min_angle = make_dsfunc2(algs.min_angle)
norm = make_ufunc(algs.norm)
shear = make_dsfunc2(algs.shear)
skew = make_dsfunc2(algs.skew)
strain = make_dsfunc(algs.strain)
surface_normal = make_dsfunc2(algs.surface_normal)
trace = make_ufunc(algs.trace)
volume = make_dsfunc2(algs.volume)
vorticity = make_dsfunc(algs.vorticity)
vertex_normal = make_dsfunc2(algs.vertex_normal)
