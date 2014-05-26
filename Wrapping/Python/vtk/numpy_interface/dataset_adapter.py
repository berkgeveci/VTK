try:
    import numpy
except ImportError:
    raise RuntimeError("This module depends on the numpy module. Please make\
sure that it is installed properly.")

import itertools
import operator
from vtk.util import numpy_support

class ArrayAssociation :
    """Easy access to vtkDataObject.AttributeTypes"""
    POINT = 0
    CELL  = 1
    FIELD = 2
    ROW = 6

class VTKObjectWrapper(object):
    "Superclass for classes that wrap VTK objects with Python objects."
    def __init__(self, vtkobject):
        self.VTKObject = vtkobject

    def __getattr__(self, name):
        "Forwards unknown attribute requests to VTK object."
        return getattr(self.VTKObject, name)

def MakeObserver(numpy_array):
    "Internal function used to attach a numpy array to a vtk array"
    def Closure(caller, event):
        foo = numpy_array
    return Closure

def vtkDataArrayToVTKArray(array, dataset=None):
    "Given a vtkDataArray and a dataset owning it, returns a VTKArray."
    narray = numpy_support.vtk_to_numpy(array)

    # Make arrays of 9 components into matrices. Also transpose
    # as VTK store matrices in Fortran order
    shape = narray.shape
    if len(shape) == 2 and shape[1] == 9:
        narray = narray.reshape((shape[0], 3, 3)).transpose(0, 2, 1)

    return VTKArray(narray, array=array, dataset=dataset)

def numpyTovtkDataArray(array, name="numpy_array"):
    """Given a numpy array or a VTKArray and a name, returns a vtkDataArray.
    The resulting vtkDataArray will store a reference to the numpy array
    through a DeleteEvent observer: the numpy array is released only when
    the vtkDataArray is destroyed."""
    if not array.flags.contiguous:
        array = array.copy()
    vtkarray = numpy_support.numpy_to_vtk(array)
    vtkarray.SetName(name)
    # This makes the VTK array carry a reference to the numpy array.
    vtkarray.AddObserver('DeleteEvent', MakeObserver(array))
    return vtkarray

def make_tensor_array_contiguous(array):
    if array == None:
        return None
    if array.flags.contiguous:
        return array
    array = numpy.asarray(array)
    size = array.dtype.itemsize
    strides = array.strides
    if len(strides) == 3 and strides[1]/size == 1 and strides[2]/size == 3:
        return array.transpose(0, 2, 1)
    return array

class VTKArray(numpy.ndarray):
    """This is a sub-class of numpy matrix that stores a
    reference to a vtk array as well as the owning dataset.
    The numpy array and vtk array should point to the same
    memory location."""

    def __new__(cls, input_array, array=None, dataset=None):
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        obj = numpy.asarray(input_array).view(cls)
        obj.Association = ArrayAssociation.FIELD
        # add the new attributes to the created instance
        obj.VTKObject = array
        # if dataset:
        #     import weakref
        #     obj.DataSet = weakref.ref(dataset)
        obj.DataSet = dataset
        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self,obj):
        # Copy the VTK array only if the two share data
        slf = make_tensor_array_contiguous(self)
        obj2 = make_tensor_array_contiguous(obj)

        self.VTKObject = None
        try:
            # This line tells us that they are referring to the same buffer.
            # Much like two pointers referring to same memory location in C/C++.
            if buffer(slf) == buffer(obj2):
                self.VTKObject = getattr(obj, 'VTKObject', None)
        except TypeError:
            pass

        self.Association = getattr(obj, 'Association', None)
        self.DataSet = getattr(obj, 'DataSet', None)

    def __getattr__(self, name):
        "Forwards unknown attribute requests to VTK array."
        if not hasattr(self, "VTKObject") or not self.VTKObject:
            raise AttributeError("class has no attribute %s" % name)
        return getattr(self.VTKObject, name)

    def __mul__(self, other):
        return numpy.multiply(self, other)

    def __rmul__(self, other):
        return numpy.multiply(self, other)

    def __pow__(self, other):
        return numpy.power(self, other)


class VTKCompositeDataArray():

    def __metaclass__(name, parent, attr):
        """Simplify the implementation of the numeric/logical sequence API."""
        def add_numeric_op(attr_name, op):
            """Create an attribute named attr_name that calls
            _numeric_op(self, other, op)."""
            def closure(self, other):
                return VTKCompositeDataArray._numeric_op(self, other, op)
            closure.__name__ = attr_name
            attr[attr_name] = closure

        def add_reverse_numeric_op(attr_name, op):
            """Create an attribute named attr_name that calls
            _reverse_numeric_op(self, other, op)."""
            def closure(self, other):
                return VTKCompositeDataArray._reverse_numeric_op(self, other, op)
            closure.__name__ = attr_name
            attr[attr_name] = closure

        def add_default_reverse_numeric_op(op_name):
            """Adds '__r[op_name]__' attribute that uses operator.[op_name]"""
            add_reverse_numeric_op("__r%s__"%op_name, getattr(operator, op_name))

        def add_default_numeric_op(op_name):
            """Adds '__[op_name]__' attribute that uses operator.[op_name]"""
            add_numeric_op("__%s__"%op_name, getattr(operator, op_name))

        def add_default_numeric_ops(op_name):
            """Call both add_default_numeric_op and add_default_reverse_numeric_op."""
            add_default_numeric_op(op_name)
            add_default_reverse_numeric_op(op_name)

        add_default_numeric_ops("add")
        add_default_numeric_ops("sub")
        add_default_numeric_ops("mul")
        add_default_numeric_ops("div")
        add_default_numeric_ops("truediv")
        add_default_numeric_ops("floordiv")
        add_default_numeric_ops("mod")
        add_default_numeric_ops("pow")
        add_default_numeric_ops("lshift")
        add_default_numeric_ops("rshift")
        add_numeric_op("__and__", operator.and_)
        add_reverse_numeric_op("__rand__", operator.and_)
        add_default_numeric_ops("xor")
        add_numeric_op("__or__", operator.or_)
        add_reverse_numeric_op("__ror__", operator.or_)

        add_default_numeric_op("lt")
        add_default_numeric_op("le")
        add_default_numeric_op("eq")
        add_default_numeric_op("ne")
        add_default_numeric_op("ge")
        add_default_numeric_op("gt")
        return type(name, parent, attr)

    def GetSize(self):
        size = numpy.int64(0)
        for a in self.Arrays:
            try:
                size += a.size
            except AttributeError:
                pass
        return size

    size = property(GetSize)

    def __init__(self, arrays = []):
        self.Arrays = arrays

    def InitFromCompositeData(self, composite_data, array_name):
        self.Arrays = []
        for ds in composite_data:
            self.Arrays.append(ds.PointData[array_name])

    def __getitem__(self, index):
        if type(index) != tuple:
            index = (index,)
        res = []
        if type(index[0]) == VTKCompositeDataArray:
            for a, idx in itertools.izip(self.Arrays, index[0].Arrays):
                if a != None:
                    res.append(a.__getitem__((idx,)+index[1:]))
                else:
                    res.append(None)
        else:
            for a in self.Arrays:
                if a != None:
                    res.append(a.__getitem__(index))
                else:
                    res.append(None)
        return VTKCompositeDataArray(res)

    def _numeric_op(self, other, op):
        """Used to implement numpy-style numerical operations such as __add__,
        __mul__, etc."""
        res = []
        if type(other) == VTKCompositeDataArray:
            for a1, a2 in itertools.izip(self.Arrays, other.Arrays):
                if a1 != None and a2 != None:
                    res.append(op(a1,a2))
                else:
                    res.append(None)
        else:
            for a in self.Arrays:
                if a != None:
                    res.append(op(a, other))
                else:
                    res.append(None)
        return VTKCompositeDataArray(res)

    def _reverse_numeric_op(self, other, op):
        """Used to implement numpy-style numerical operations such as __add__,
        __mul__, etc."""
        res = []
        if type(other) == VTKCompositeDataArray:
            for a1, a2 in itertools.izip(self.Arrays, other.Arrays):
                if a1 != None and a2 != None:
                    res.append(op(a2,a1))
                else:
                    res.append(None)
        else:
            for a in self.Arrays:
                if a != None:
                    res.append(op(other, a))
                else:
                    res.append(None)
        return VTKCompositeDataArray(res)

    def __str__(self):
        return self.Arrays.__str__()

class DataSetAttributes(VTKObjectWrapper):
    """This is a python friendly wrapper of vtkDataSetAttributes. It
    returns VTKArrays. It also provides the dictionary interface."""

    def __init__(self, vtkobject, dataset, association):
        super(DataSetAttributes, self).__init__(vtkobject)
        # import weakref
        # self.DataSet = weakref.ref(dataset)
        self.DataSet = dataset
        self.Association = association

    def __getitem__(self, idx):
        """Implements the [] operator. Accepts an array name."""
        return self.GetArray(idx)

    def GetArray(self, idx):
        "Given an index or name, returns a VTKArray."
        vtkarray = self.VTKObject.GetArray(idx)
        if not vtkarray:
            vtkarray = self.VTKObject.GetAbstractArray(idx)
            if vtkarray:
                return vtkarray
            return None
        array = vtkDataArrayToVTKArray(vtkarray, self.DataSet)
        array.Association = self.Association
        return array

    def keys(self):
        """Returns the names of the arrays as a list."""
        kys = []
        narrays = self.VTKObject.GetNumberOfArrays()
        for i in range(narrays):
            name = self.VTKObject.GetArray(i).GetName()
            if name:
                kys.append(name)
        return kys

    def values(self):
        """Returns the arrays as a list."""
        vals = []
        narrays = self.VTKObject.GetNumberOfArrays()
        for i in range(narrays):
            a = self.VTKObject.GetArray(i)
            if a.GetName():
                vals.append(a)
        return vals

    def PassData(self, other):
        try:
            self.VTKObject.PassData(other)
        except TypeError:
            self.VTKObject.PassData(other.VTKObject)

    def append(self, narray, name):
        """Appends a new array to the dataset attributes."""

        if self.Association == ArrayAssociation.POINT:
            arrLength = self.DataSet.GetNumberOfPoints()
        elif self.Association == ArrayAssociation.CELL:
            arrLength = self.DataSet.GetNumberOfCells()

        # Fixup input array length:
        if not isinstance(narray, numpy.ndarray): # Scalar input
            narray = narray * numpy.ones((arrLength, 1))
        elif narray.shape[0] != arrLength: # Vector input
            components = reduce(operator.mul, narray.shape)
            narray = narray.flatten() * numpy.ones((arrLength, components))

        shape = narray.shape

        if len(shape) == 3:
            # Array of matrices. We need to make sure the order  in memory is right.
            # If column order (c order), transpose. VTK wants row order (fortran
            # order). The deep copy later will make sure that the array is contiguous.
            # If row order but not contiguous, transpose so that the deep copy below
            # does not happen.
            size = narray.dtype.itemsize
            if (narray.strides[1]/size == 3 and narray.strides[2]/size == 1) or \
                (narray.strides[1]/size == 1 and narray.strides[2]/size == 3 and \
                 not narray.flags.contiguous):
                narray  = narray.transpose(0, 2, 1)

        # If array is not contiguous, make a deep copy that is contiguous
        if not narray.flags.contiguous:
            narray = narray.copy()

        # Flatten array of matrices to array of vectors
        if len(shape) == 3:
            narray = narray.reshape(shape[0], shape[1]*shape[2])

        arr = numpyTovtkDataArray(narray, name)
        self.VTKObject.AddArray(arr)


class CompositeDataSetAttributes():
    """This is a python friendly wrapper for vtkDataSetAttributes for composite
    datsets. Since composite datasets themselves don't have attribute data, but
    the attribute data is associated with the leaf nodes in the composite
    dataset, this class simulates a DataSetAttributes interface by taking a
    union of DataSetAttributes assiciated with all leaf nodes."""

    def __init__(self, dataset, association):
        # import weakref
        # self.DataSet = weakref.ref(dataset)
        self.DataSet = dataset
        self.Association = association
        self.ArrayNames = []

        # build the set of arrays available in the composite dataset. Since
        # composite datasets can have partial arrays, we need to iterate over
        # all non-null blocks in the dataset.
        self.__determine_arraynames()

    def __determine_arraynames(self):
        array_set = set()
        array_list = []
        for dataset in self.DataSet:
            dsa = dataset.GetAttributes(self.Association)
            for array_name in dsa.keys():
                if array_name not in array_set:
                    array_set.add(array_name)
                    array_list.append(array_name)
        self.ArrayNames = array_list

    def keys(self):
        """Returns the names of the arrays as a list."""
        return self.ArrayNames

    def __getitem__(self, idx):
        """Implements the [] operator. Accepts an array name."""
        return self.GetArray(idx)

    def append(self, narray, name):
        """Appends a new array to the composite dataset attributes."""
        for ds, array in itertools.izip(self.DataSet, narray.Arrays):
            if array != None:
                if self.Association == ArrayAssociation.POINT:
                    ds.PointData.append(array, name)
                elif self.Association == ArrayAssociation.CELL:
                    ds.CellData.append(array, name)
                elif self.Association == ArrayAssociation.FIELD:
                    ds.FieldData.append(array, name)
                elif self.Association == ArrayAssociation.ROW:
                    ds.RowData.append(array, name)

    def GetArray(self, idx):
        """Given an index or name, returns a VTKCompositeArray."""
        if type(idx) == int:
            arrayname = self.ArrayNames[idx]
        else:
            arrayname = idx
        if arrayname not in self.ArrayNames:
            return None
        array = VTKCompositeDataArray()
        array.InitFromCompositeData(self.DataSet, arrayname)
        return array

    def PassData(self, other):
        """Emulate PassData for composite datasets."""
        for this,that in zip(self.DataSet, other.DataSet):
            for assoc in [ArrayAssociation.POINT, ArrayAssociation.CELL]:
                this.GetAttributes(assoc).PassData(that.GetAttributes(assoc))

class CompositeDataIterator(object):
    """Wrapper for a vtkCompositeDataIterator class to satisfy
       the python iterator protocol.
       """

    def __init__(self, cds):
        self.Iterator = cds.NewIterator()
        if self.Iterator:
            self.Iterator.UnRegister(None)
            self.Iterator.GoToFirstItem()

    def __iter__(self):
        return self

    def next(self):
        if not self.Iterator:
            raise StopIteration

        if self.Iterator.IsDoneWithTraversal():
            raise StopIteration
        retVal = self.Iterator.GetCurrentDataObject()
        self.Iterator.GoToNextItem()
        return WrapDataObject(retVal)

    def __getattr__(self, name):
        """Returns attributes from the vtkCompositeDataIterator."""
        return getattr(self.Iterator, name)

class MultiCompositeDataIterator(CompositeDataIterator):
    def __init__(self, cds):
        CompositeDataIterator.__init__(self, cds[0])
        self.Datasets = cds

    def next(self):
        if not self.Iterator:
            raise StopIteration

        if self.Iterator.IsDoneWithTraversal():
            raise StopIteration
        retVal = []
        retVal.append(WrapDataObject(self.Iterator.GetCurrentDataObject()))
        if len(self.Datasets) > 1:
            for cd in self.Datasets[1:]:
                retVal.append(WrapDataObject(cd.GetDataSet(self.Iterator)))
        self.Iterator.GoToNextItem()
        return retVal

class DataObject(VTKObjectWrapper):

    def GetAttributes(self, type):
        """Returns the attributes specified by the type as a DataSetAttributes
         instance."""
        return DataSetAttributes(self.VTKObject.GetAttributes(type), self, type)

    def GetFieldData(self):
        "Returns the field data as a DataSetAttributes instance."
        return DataSetAttributes(self.VTKObject.GetFieldData(), self, ArrayAssociation.FIELD)

    FieldData = property(GetFieldData, None, None, "This property returns \
        the field data of a data object.")

class Table(DataObject):
    def GetRowData(self):
        "Returns the row data as a DataSetAttributes instance."
        return self.GetAttributes(ArrayAssociation.ROW)

    RowData = property(GetRowData, None, None, "This property returns \
        the row data of the table.")

class CompositeDataSet(DataObject):
    def __iter__(self):
        "Creates an iterator for the contained datasets."
        return CompositeDataIterator(self)

    def GetNumberOfElements(self, assoc):
        result = 0
        for dataset in self:
            result += dataset.GetNumberOfElements(assoc)
        return int(result)

    def GetNumberOfPoints(self):
        return self.GetNumberOfElements(ArrayAssociation.POINT)

    def GetNumberOfCells(self):
        return self.GetNumberOfElements(ArrayAssociation.CELL)

    def GetAttributes(self, type):
        """Returns the attributes specified by the type as a
        CompositeDataSetAttributes instance."""
        return CompositeDataSetAttributes(self, type)

    def GetPointData(self):
        "Returns the point data as a DataSetAttributes instance."
        return self.GetAttributes(ArrayAssociation.POINT)

    def GetCellData(self):
        "Returns the cell data as a DataSetAttributes instance."
        return self.GetAttributes(ArrayAssociation.CELL)

    PointData = property(GetPointData, None, None, "This property returns \
        the point data of the dataset.")
    CellData = property(GetCellData, None, None, "This property returns \
        the cell data of a dataset.")

class DataSet(DataObject):
    """This is a python friendly wrapper of a vtkDataSet that defines
    a few useful properties."""

    def GetPointData(self):
        "Returns the point data as a DataSetAttributes instance."
        return self.GetAttributes(ArrayAssociation.POINT)

    def GetCellData(self):
        "Returns the cell data as a DataSetAttributes instance."
        return self.GetAttributes(ArrayAssociation.CELL)

    PointData = property(GetPointData, None, None, "This property returns \
        the point data of the dataset.")
    CellData = property(GetCellData, None, None, "This property returns \
        the cell data of a dataset.")

class PointSet(DataSet):
    def GetPoints(self):
        """Returns the points as a VTKArray instance. Returns None if the
        dataset has implicit points."""
        if not self.VTKObject.GetPoints():
            return None
        return vtkDataArrayToVTKArray(
            self.VTKObject.GetPoints().GetData(), self)

    Points = property(GetPoints, None, None, "This property returns the \
        point coordinates of dataset.")

class PolyData(PointSet):
    def GetPolygons(self):
        """Returns the points as a VTKArray instance. Returns None if the
        dataset has implicit points."""
        if not self.VTKObject.GetPolys():
            return None
        return vtkDataArrayToVTKArray(
            self.VTKObject.GetPolys().GetData(), self)

    Polygons = property(GetPolygons, None, None, "This property returns the \
        connectivity of polygons.")

def WrapDataObject(ds):
    if ds.IsA("vtkPolyData"):
        return PolyData(ds)
    elif ds.IsA("vtkPointSet"):
        return PointSet(ds)
    elif ds.IsA("vtkDataSet"):
        return DataSet(ds)
    elif ds.IsA("vtkCompositeDataSet"):
        return CompositeDataSet(ds)
    elif ds.IsA("vtkTable"):
        return Table(ds)
