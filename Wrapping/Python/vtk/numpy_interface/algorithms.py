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

def _apply_func2(func, array, args):
    """Apply a function to each member of a VTKCompositeDataArray.
    Returns a list of arrays.

    Note that this function is mainly for internal use by this module."""
    if array is dsa.NoneArray:
        return []
    res = []
    for a in array.Arrays:
        res.append(func(a, *args))
    return res

def _apply_func(func, array, args):
    """Apply a function to each member of a VTKCompositeDataArray.
    Returns a VTKCompositeDataArray of results.

    Note that this function is mainly for internal use by this module."""
    return dsa.VTKCompositeDataArray(_apply_func2(func, array, args), dataset = array.DataSet)

def _make_ufunc(ufunc):
    """ Given a ufunc, creates a closure that applies it to each member
    of a VTKCompositeDataArray.

    Note that this function is mainly for internal use by this module."""
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

def _make_dfunc(dfunc):
    """ Given a function that requires two arguments, creates a closure that
    applies it to each member of a VTKCompositeDataArray.

    Note that this function is mainly for internal use by this module."""
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

def _make_dsfunc(dsfunc):
    """ Given a function that requires two arguments (one array, one dataset),
    creates a closure that applies it to each member of a VTKCompositeDataArray.
    Note that this function is mainly for internal use by this module."""
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

def _make_dsfunc2(dsfunc):
    """ Given a function that requires a dataset, creates a closure that
    applies it to each member of a VTKCompositeDataArray.

    Note that this function is mainly for internal use by this module."""
    def new_dsfunc2(ds):
        if type(ds) == dsa.CompositeDataSet:
            res = []
            for dataset in ds:
                res.append(dsfunc(dataset))
            return dsa.VTKCompositeDataArray(res, dataset = ds)
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

def _reduce_dims(array, comm):
    from mpi4py import MPI
    dims = numpy.array([0, 0], dtype=numpy.int32)
    if array is not dsa.NoneArray:
        shp = shape(array)
        if len(shp) == 0:
            dims = numpy.array([1, 0], dtype=numpy.int32)
        elif len(shp) == 1:
            dims = numpy.array([shp[0], 0], dtype=numpy.int32)
        else:
            dims = numpy.array(shp, dtype=numpy.int32)
    max_dims = numpy.array(dims, dtype=numpy.int32)
    comm.Allreduce([dims, MPI.INT], [max_dims, MPI.INT], MPI.MAX)

    if max_dims[1] == 0:
        max_dims = numpy.array((max_dims[0],))
        size = max_dims[0]
    else:
        size = max_dims[0]*max_dims[1]

    if max_dims[0] == 1:
        max_dims = 1

    return (max_dims, size)

def _global_func(impl, array, axis, controller):
    if type(array) == dsa.VTKCompositeDataArray:
        if axis is None or axis == 0:
            res = impl.serial_composite(array, axis)
        else:
            res = _apply_func(impl.op(), array, (axis,))
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

            max_dims, size = _reduce_dims(res, comm)

            # All NoneArrays
            if size == 0:
                return dsa.NoneArray;

            if res is dsa.NoneArray:
                if max_dims is 1:
                    # Weird trick to make the array look like a scalar
                    max_dims = ()
                res = numpy.empty(max_dims)
                res.fill(impl.default())

            res_recv = numpy.array(res)
            mpi_type = _lookup_mpi_type(res.dtype)
            comm.Allreduce([res, mpi_type], [res_recv, mpi_type], impl.mpi_op())
            if array is dsa.NoneArray:
                return dsa.NoneArray
            res = res_recv

    return res

def sum(array, axis=None, controller=None):
    """Returns the sum of all values along a particular axis (dimension).
    Given an array of m tuples and n components:
    * Default is to return the sum of all values in an array.
    * axis=0: Sum values of all components and return a one tuple,
      n-component array.
    * axis=1: Sum values of all components of each tuple and return an
      m-tuple, 1-component array.

    When called in parallel, this function will sum across processes
    when a controller argument is passed or the global controller is
    defined. To disable parallel summing when running in parallel, pass
    a dummy controller as follows:

    sum(array, controller=vtk.vtkDummyController().
    """
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

        def default(self):
            return numpy.float64(0)

    return _global_func(SumImpl(), array, axis, controller)

def max(array, axis=None, controller=None):
    """Returns the max of all values along a particular axis (dimension).
    Given an array of m tuples and n components:
    * Default is to return the max of all values in an array.
    * axis=0: Return the max values of all components and return a
      one tuple, n-component array.
    * axis=1: Return the max values of all components of each tuple
      and return an m-tuple, 1-component array.

    When called in parallel, this function will compute the max across
    processes when a controller argument is passed or the global controller
    is defined. To disable parallel summing when running in parallel, pass a
    dummy controller as follows:

    max(array, controller=vtk.vtkDummyController().
    """
    class MaxImpl:
        def op(self):
            return algs.max

        def mpi_op(self):
            from mpi4py import MPI
            return MPI.MAX

        def serial_composite(self, array, axis):
            res = _apply_func2(algs.max, array, (axis,))
            clean_list = []
            for a in res:
                if a is not dsa.NoneArray:
                    clean_list.append(a)
            if clean_list is []:
                return None
            return algs.max(clean_list, axis=0).astype(numpy.float64)

        def default(self):
            return numpy.finfo(numpy.float64).min

    return _global_func(MaxImpl(), array, axis, controller)

def min(array, axis=None, controller=None):
    """Returns the min of all values along a particular axis (dimension).
    Given an array of m tuples and n components:
    * Default is to return the min of all values in an array.
    * axis=0: Return the min values of all components and return a one
      tuple, n-component array.
    * axis=1: Return the min values of all components of each tuple and
      return an m-tuple, 1-component array.

    When called in parallel, this function will compute the min across processes
    when a controller argument is passed or the global controller is defined.
    To disable parallel summing when running in parallel, pass a dummy controller as follows:

    min(array, controller=vtk.vtkDummyController().
    """
    class MinImpl:
        def op(self):
            return algs.min

        def mpi_op(self):
            from mpi4py import MPI
            return MPI.MIN

        def serial_composite(self, array, axis):
            res = _apply_func2(algs.min, array, (axis,))
            clean_list = []
            for a in res:
                if a is not dsa.NoneArray:
                    clean_list.append(a)
            if clean_list is []:
                return None
            return algs.min(clean_list, axis=0).astype(numpy.float64)

        def default(self):
            return numpy.finfo(numpy.float64).max

    return _global_func(MinImpl(), array, axis, controller)

def _global_per_block(impl, array, axis=None, controller=None):
    if axis > 0:
        return impl.op()(array, axis=axis, controller=controller)

    try:
        dataset = array.DataSet
    except AttributeError:
        dataset = None

    t = type(array)
    if t == dsa.VTKArray or t == numpy.ndarray:
        from vtk.vtkCommonDataModel import vtkMultiBlockDataSet
        array = dsa.VTKCompositeDataArray([array])
        ds = vtkMultiBlockDataSet()
        ds.SetBlock(0, dataset.VTKObject)
        dataset = ds

    results = _apply_func2(impl.op2(), array, (axis,))

    if controller is None and vtkMultiProcessController is not None:
        controller = vtkMultiProcessController.GetGlobalController()
    if controller and controller.IsA("vtkMPIController"):
        from mpi4py import MPI
        comm = vtkMPI4PyCommunicator.ConvertToPython(controller.GetCommunicator())

        # First determine the number of components to use
        # for reduction
        res = dsa.NoneArray
        for res in results:
            if res is not dsa.NoneArray:
                break

        max_dims, size = _reduce_dims(res, comm)

        # All NoneArrays
        if size == 0:
            return dsa.NoneArray;

        # Next determine the max id to use for reduction
        # operations

        # Get all ids from dataset, including empty ones.
        ids = []
        lmax_id = numpy.int64(0)
        if dataset is not None:
            it = dataset.NewIterator()
            it.UnRegister(None)
            it.SetSkipEmptyNodes(False)
            while not it.IsDoneWithTraversal():
                _id = it.GetCurrentFlatIndex()
                lmax_id = numpy.max((lmax_id, _id))
                if it.GetCurrentDataObject() is not None:
                    ids.append(_id)
                it.GoToNextItem()
        max_id = numpy.array(0, dtype=numpy.int64)
        comm.Allreduce([lmax_id, MPI.INT], [max_id, MPI.INT], MPI.MAX)

        has_ids = numpy.zeros(max_id+1, dtype=numpy.int32)
        for _id in ids:
            has_ids[_id] = 1
        id_count = numpy.array(has_ids)
        comm.Allreduce([has_ids, MPI.INT], [id_count, MPI.INT], MPI.SUM)

        if numpy.all(id_count <= 1):
            return dsa.VTKCompositeDataArray(results, dataset=dataset)

        # Now that we know which blocks are shared by more than
        # 1 rank. The ones that have a count of 2 or more.
        reduce_ids = []
        for _id in xrange(len(id_count)):
            if id_count[_id] > 1:
                reduce_ids.append(_id)

        to_reduce = len(reduce_ids)
        # If not block is shared, short circuit. No need to
        # communicate any more.
        if to_reduce == 0:
            return dsa.VTKCompositeDataArray(results, dataset=dataset)

        # Create the local array that will be used for
        # reduction. Set it to a value that won't effect
        # the reduction.
        lresults = numpy.empty(size*to_reduce)
        lresults.fill(impl.default())

        # Just get non-empty ids. Doing this again in case
        # the traversal above results in a different order.
        # We need the same order since we'll use izip below.
        if dataset is not None:
            it = dataset.NewIterator()
            it.UnRegister(None)
            ids = []
            while not it.IsDoneWithTraversal():
                ids.append(it.GetCurrentFlatIndex())
                it.GoToNextItem()

        # Fill the local array with available values.
        for _id, _res in itertools.izip(ids, results):
            success = True
            try:
                loc = reduce_ids.index(_id)
            except ValueError:
                success = False
            if success:
                if _res is not dsa.NoneArray:
                    lresults[loc*size:(loc+1)*size] = _res.flatten()

        # Now do the MPI reduction.
        rresults = numpy.array(lresults)
        comm.Allreduce([lresults, MPI.DOUBLE], [rresults, MPI.DOUBLE], impl.mpi_op())

        if array is dsa.NoneArray:
            return dsa.NoneArray

        # Fill in the reduced values.
        for i in xrange(to_reduce):
            _id = reduce_ids[i]
            success = True
            try:
                loc = ids.index(_id)
            except ValueError:
                success = False
            if success:
                if size == 1:
                    results[loc] = dsa.VTKArray(rresults[i])
                else:
                    results[loc] = rresults[i*size:(i+1)*size].reshape(max_dims)

    return dsa.VTKCompositeDataArray(results, dataset=dataset)

def max_per_block(array, axis=None, controller=None):
    """Returns the max of all values along a particular axis (dimension)
    for each block of a VTKCompositeDataArray.
    Given an array of m tuples and n components:
    * Default is to return the max of all values in an array.
    * axis=0: Return the max values of all components and return a one
      tuple, n-component array.
    * axis=1: Return the max values of all components of each tuple and return
      an m-tuple, 1-component array.

    When called in parallel, this function will compute the max across
    processes when a controller argument is passed or the global controller
    is defined. To disable parallel summing when running in parallel, pass a
    dummy controller as follows:

    max(array, controller=vtk.vtkDummyController().
    """
    class MaxPerBlockImpl:
        def op(self):
            return max

        def op2(self):
            return algs.max

        def mpi_op(self):
            from mpi4py import MPI
            return MPI.MAX

        def default(self):
            return numpy.finfo(numpy.float64).min

    return _global_per_block(MaxPerBlockImpl(), array, axis, controller)

def min_per_block(array, axis=None, controller=None):
    """Returns the min of all values along a particular axis (dimension)
    for each block of a VTKCompositeDataArray.
    Given an array of m tuples and n components:
    * Default is to return the min of all values in an array.
    * axis=0: Return the min values of all components and return a one
      tuple, n-component array.
    * axis=1: Return the min values of all components of each tuple and
      return an m-tuple, 1-component array.

    When called in parallel, this function will compute the min across
    processes when a controller argument is passed or the global controller
    is defined. To disable parallel summing when running in parallel, pass a
    dummy controller as follows:

    min(array, controller=vtk.vtkDummyController().
    """
    class MinPerBlockImpl:
        def op(self):
            return min

        def op2(self):
            return algs.min

        def mpi_op(self):
            from mpi4py import MPI
            return MPI.MIN

        def default(self):
            return numpy.finfo(numpy.float64).max

    return _global_per_block(MinPerBlockImpl(), array, axis, controller)

def all(array, axis=None, controller=None):
    """Returns True if all values of an array evaluate to True, returns
    False otherwise.
    This is useful to check if all values of an array match a certain
    condition such as:

    algorithms.all(array > 5)
    """
    class MinImpl:
        def op(self):
            return algs.all

        def mpi_op(self):
            from mpi4py import MPI
            return MPI.LAND

        def serial_composite(self, array, axis):
            res = _apply_func2(algs.all, array, (axis,))
            clean_list = []
            for a in res:
                if a is not dsa.NoneArray:
                    clean_list.append(a)
            if clean_list is []:
                return None
            return algs.all(clean_list, axis=0)

        def default(self, max_comps):
            return numpy.ones(max_comps, dtype=numpy.bool)

    return _global_func(MinImpl(), array, axis, controller)

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
    """Returns the mean of all values along a particular axis (dimension).
    Given an array of m tuples and n components:
    * Default is to return the mean of all values in an array.
    * axis=0: Return the mean values of all components and return a one
      tuple, n-component array.
    * axis=1: Return the mean values of all components of each tuple and
      return an m-tuple, 1-component array.

    When called in parallel, this function will compute the mean across
    processes when a controller argument is passed or the global controller
    is defined. To disable parallel summing when running in parallel, pass a
    dummy controller as follows:

    mean(array, controller=vtk.vtkDummyController().
    """

    if axis is None or axis == 0:
        if size is None:
            size = _array_count(array, axis, controller)
        return sum(array, axis) / size
    else:
        if type(array) == dsa.VTKCompositeDataArray:
            return _apply_func(algs.mean, array, (axis,))
        else:
            return algs.mean(array, axis)

def var(array, axis=None, controller=None):
    """Returns the variance of all values along a particular axis (dimension).
    Given an array of m tuples and n components:
    * Default is to return the variance of all values in an array.
    * axis=0: Return the variance values of all components and return a one
      tuple, n-component array.
    * axis=1: Return the variance values of all components of each tuple and
      return an m-tuple, 1-component array.

    When called in parallel, this function will compute the variance across
    processes when a controller argument is passed or the global controller
    is defined. To disable parallel summing when running in parallel, pass a
    dummy controller as follows:

    var(array, controller=vtk.vtkDummyController().
    """

    if axis is None or axis == 0:
        size = _array_count(array, axis, controller)
        tmp = array - mean(array, axis, controller, size)
        return sum(tmp*tmp, axis, controller) / size
    else:
        if type(array) == dsa.VTKCompositeDataArray:
            return _apply_func(algs.var, array, (axis,))
        else:
            return algs.var(array, axis)

def std(array, axis=None, controller=None):
    """Returns the standard deviation of all values along a particular
    axis (dimension).
    Given an array of m tuples and n components:
    * Default is to return the standard deviation of all values in an array.
    * axis=0: Return the standard deviation values of all components and
      return a one tuple, n-component array.
    * axis=1: Return the standard deviation values of all components of
      each tuple and return an m-tuple, 1-component array.

    When called in parallel, this function will compute the standard deviation
    across processes when a controller argument is passed or the global controller
    is defined. To disable parallel summing when running in parallel, pass a dummy
    controller as follows:

    std(array, controller=vtk.vtkDummyController().
    """
    return sqrt(var(array, axis, controller))

def shape(array):
    "Returns the shape (dimensions) of an array."
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
    """Given 2 or 3 scalar arrays, returns a vector array. If only
    2 scalars are provided, the third component will be set to 0."""
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

sqrt = _make_ufunc(numpy.sqrt)
sqrt.__doc__ = "Computes square root."

exp = _make_ufunc(numpy.exp)
exp.__doc__ = "The exponential function."

floor = _make_ufunc(numpy.floor)
floor.__doc__ = "Returns the floor of floating point values."

ceil = _make_ufunc(numpy.ceil)
ceil.__doc__ = "Returns the ceiling of floating point values."

round = _make_ufunc(numpy.round)
round.__doc__ = "Rounds floating points values to integers."

sin = _make_ufunc(numpy.sin)
sin.__doc__ = "Computes sine of values in radians."

cos = _make_ufunc(numpy.cos)
cos.__doc__ = "Computes cosine of values in radians."

tan = _make_ufunc(numpy.tan)
tan.__doc__ = "Computes tangent of values in radians."

arcsin = _make_ufunc(numpy.arcsin)
arcsin.__doc__ = "Computes inverse sine."

arccos = _make_ufunc(numpy.arccos)
arccos.__doc__ = "Computes inverse cosine."

arctan = _make_ufunc(numpy.arctan)
arctan.__doc__ = "Computes inverse tangent."

sinh = _make_ufunc(numpy.sinh)
sinh.__doc__ = "Computes hyperbolic sine."

cosh = _make_ufunc(numpy.cosh)
cosh.__doc__ = "Computes hyperbolic cosine."

tanh = _make_ufunc(numpy.tanh)
tanh.__doc__ = "Computes hyperbolic tangent."

arcsinh = _make_ufunc(numpy.arcsinh)
arcsinh.__doc__ = "Computes inverse hyperbolic sine."

arccosh = _make_ufunc(numpy.arccosh)
arccosh.__doc__ = "Computes inverse hyperbolic cosine."

arctanh = _make_ufunc(numpy.arctanh)
arctanh.__doc__ = "Computes inverse hyperbolic tangent."

where = _make_ufunc(numpy.where)
where.__doc__ = """Returns the location (indices) of an array where the given
expression is true. For scalars, it returns a single array of indices.
For vectors and matrices, it returns two arrays: first with tuple indices,
second with component indices. The output of this method can be used to
extract the values from the array also by using it as the index of the [] operator.

For example:

>>> algs.where(algs.array([1,2,3]) == 2)
(array([1]),)

>>> algs.where(algs.array([[1,2,3], [2,1,1]]) == 2)
(array([0, 1]), array([1, 0]))

>>> a = array([[1,2,3], [2,1,1]])
>>> indices = algs.where(a > 2)
>>> a[indices]
array([3])
"""

expand_dims = _make_dfunc(numpy.expand_dims)
expand_dims.__doc__ = """Insert a new dimension, corresponding to a given
position in the array shape. In VTK, this function's main use is to
enable an operator to work on a vector and a scalar field. For example,
say you want to devide each component of a vector by the magnitude of
that vector. You might try this:

>>> v
VTKArray([[ 1.,  1.,  1.],
       [ 1.,  1.,  1.],
       [ 1.,  1.,  1.],
       [ 1.,  1.,  1.],
       [ 1.,  1.,  1.]])
>>> algs.mag(v)
VTKArray([ 1.73205081,  1.73205081,  1.73205081,  1.73205081,  1.73205081])
>>> v / algs.mag(v)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
ValueError: operands could not be broadcast together with shapes (5,3) (5)

The division operator does not know how to map a scalar to a vector
due to a mismatch in dimensions. This can be solved by making the
scalar a vector of 1 component (increasing its dimension to 2) as follows:

>>> v / algs.expand_dims(algs.mag(v), 1)
VTKArray([[ 0.57735027,  0.57735027,  0.57735027],
       [ 0.57735027,  0.57735027,  0.57735027],
       [ 0.57735027,  0.57735027,  0.57735027],
       [ 0.57735027,  0.57735027,  0.57735027],
       [ 0.57735027,  0.57735027,  0.57735027]])"""

abs = _make_ufunc(algs.abs)
abs.__doc__ = "Returns the absolute values of an array of scalars/vectors/tensors."

area = _make_dsfunc2(algs.area)
area.__doc__ = "Returns the surface area of each 2D cell in a mesh."

aspect = _make_dsfunc2(algs.aspect)
aspect.__doc__ = "Returns the aspect ratio of each cell in a mesh. See Verdict documentation for details."

aspect_gamma = _make_dsfunc2(algs.aspect_gamma)
aspect_gamma.__doc__ = "Returns the aspect gamma of each cell in a mesh. This metric compares root-mean-square edge length to volume. See Verdict documentation for details."

condition = _make_dsfunc2(algs.condition)
condition.__doc__ = "Returns the condition number of each cell in a mesh. See Verdict documentation for details."

cross = _make_dfunc(algs.cross)
cross.__doc__ = "Return the cross product of two vectors."

curl = _make_dsfunc(algs.curl)
curl.__doc__ = "Returns the curl a vector field."

divergence = _make_dsfunc(algs.divergence)
divergence.__doc__ = "Returns the divergence of a vector field."

det = _make_ufunc(algs.det)
det.__doc__ = "Returns the determinant of 2D matrices."

determinant = _make_ufunc(algs.determinant)
determinant.__doc__ = "Returns the determinant of 2D matrices."

diagonal = _make_dsfunc2(algs.diagonal)
diagonal.__doc__ = "Returns the diagonal length of each cell in a dataset. See Verdict documentation for details"

dot = _make_dfunc(algs.dot)
dot.__doc__ = "Returns the dot product of two vectors."

eigenvalue = _make_ufunc(algs.eigenvalue)
eigenvalue.__doc__ = "Returns the eigenvalues of 3x3 matrices. Currently only works with symmetric matrices."

eigenvector = _make_ufunc(algs.eigenvector)
eigenvector.__doc__ = "Returns the eigenvectors of 3x3 matrices. Currently only works with symmetric matrices."

gradient = _make_dsfunc(algs.gradient)
gradient.__doc__ = "Returns the gradient of  scalars or vectors."

inv = _make_ufunc(algs.inv)
inv.__doc__ = "Returns the inverse of 3x3 matrices."

inverse = _make_ufunc(algs.inverse)
inverse.__doc__ = "Returns the inverse of 3x3 matrices."

jacobian = _make_dsfunc2(algs.jacobian)
jacobian.__doc__ = "Returns the Jacobian of a dataset."

laplacian = _make_dsfunc(algs.laplacian)
laplacian.__doc__ = "Returns the Laplacian of a scalar field."

ln = _make_ufunc(algs.ln)
ln.__doc__ = "Returns the natural logarithm of its input."

log = _make_ufunc(algs.log)
log.__doc__ = "Returns the natural logarithm of its input."

log10 = _make_ufunc(algs.log10)
log10.__doc__ = "Returns the base 10 logarithm of its input."

max_angle = _make_dsfunc2(algs.max_angle)
max_angle.__doc__ = "Returns the maximum angle of each cell in a dataset. See Verdict documentation for details"

mag = _make_ufunc(algs.mag)
mag.__doc__ = "Returns the magnitude of vectors."

min_angle = _make_dsfunc2(algs.min_angle)
min_angle.__doc__ = "Returns the minimum angle of each cell in a dataset."

norm = _make_ufunc(algs.norm)
norm.__doc__ = "Computes the normalized values of vectors."

shear = _make_dsfunc2(algs.shear)
shear.__doc__ = "Returns the shear of each cell in a dataset. See Verdict documentation for details."

skew = _make_dsfunc2(algs.skew)
skew.__doc__ = "Returns the skew of each cell in a dataset. See Verdict documentation for details."

strain = _make_dsfunc(algs.strain)
strain.__doc__ = "Given a deformation vector, this function computes the infinitesimal (Cauchy) strain tensor. It can also be used to compute strain rate if the input is velocity."

surface_normal = _make_dsfunc2(algs.surface_normal)
surface_normal.__doc__ = "Returns the surface normal of each cell in a dataset."

trace = _make_ufunc(algs.trace)
trace.__doc__ = "Returns the trace of square matrices."

volume = _make_dsfunc2(algs.volume)
volume.__doc__ = "Returns the volume of each cell in a dataset. Use sum to calculate total volume of a dataset."

vorticity = _make_dsfunc(algs.vorticity)
vorticity.__doc__ = "Given a velocity field, calculates vorticity."

vertex_normal = _make_dsfunc2(algs.vertex_normal)
vertex_normal.__doc__ = "Returns the normal at each vertex of a dataset, which is defined as the average of the cell normals of all cells containing that vertex."