import dataset_adapter as dsa
import internal_algorithms as algs
import itertools
import numpy

def apply_func2(func, array, args):
    res = []
    for a in array.Arrays:
        if a != None:
            res.append(func(a, *args))
        else:
            res.append(None)
    return res

def apply_func(func, array, args):
    return dsa.VTKCompositeDataArray(apply_func2(func, array, args))

def make_ufunc2(ufunc):
    def new_ufunc2(array):
        if type(array) == dsa.VTKCompositeDataArray:
            res = []
            for a in array.Arrays:
                if a != None:
                    res.append(ufunc(a))
                else:
                    res.append(None)
            return res
        else:
            return ufunc(array)
    return new_ufunc2

def make_ufunc(ufunc):
    def new_ufunc(array):
        if type(array) == dsa.VTKCompositeDataArray:
            res = []
            for a in array.Arrays:
                if a != None:
                    res.append(ufunc(a))
                else:
                    res.append(None)
            return dsa.VTKCompositeDataArray(res)
        else:
            return ufunc(array)
    return new_ufunc

def make_dfunc(dfunc):
    def new_dfunc(array1, val2):
        if type(array1) == dsa.VTKCompositeDataArray and type(val2) == dsa.VTKCompositeDataArray:
            res = []
            for a1, a2 in itertools.izip(array1.Arrays, val2.Arrays):
                if a1 != None and a2 != None:
                    res.append(dfunc(a1, a2))
                else:
                    res.append(None)
            return dsa.VTKCompositeDataArray(res)
        elif type(array1) == dsa.VTKCompositeDataArray:
            res = []
            for a in array1.Arrays :
                if a != None:
                    res.append(dfunc(a, val2))
                else:
                    res.append(None)
            return dsa.VTKCompositeDataArray(res)
        else:
            return dfunc(array1, val2)
    return new_dfunc

def make_dsfunc(dsfunc):
    def new_dsfunc(array, ds=None):
        if type(array) == dsa.VTKCompositeDataArray:
            res = []
            for a in array.Arrays:
                if a != None:
                    res.append(dsfunc(a, ds))
                else:
                    res.append(None)
            return dsa.VTKCompositeDataArray(res)
        else:
            return dsfunc(array, ds)
    return new_dsfunc

def make_dsfunc2(dsfunc):
    def new_dsfunc2(ds):
        if type(ds) == dsa.CompositeDataSet:
            res = []
            for dataset in ds:
                res.append(dsfunc(dataset))
            return dsa.VTKCompositeDataArray(res)
        else:
            return dsfunc(ds)
    return new_dsfunc2

def sum(array, axis=None):
    if type(array) == dsa.VTKCompositeDataArray:
        if axis is None or axis == 0:
            res = None
            arrays = array.Arrays
            for a in arrays:
                if a != None:
                    if res == None:
                        res = numpy.sum(a, axis)
                    else:
                        res += numpy.sum(a, axis)
            return res
        else:
            return apply_func(numpy.sum, array, (axis,))
    else:
        return numpy.sum(array, axis)

def max(array, axis=None):
    if type(array) == dsa.VTKCompositeDataArray:
        l = apply_func2(numpy.max, array, (axis,))
        if axis is None or axis == 0:
            # Reduce over the list
            return numpy.max(l, axis=0)
        else:
            return l
    else:
        return numpy.max(array)

def min(array, axis=None):
    if type(array) == dsa.VTKCompositeDataArray:
        l = apply_func2(numpy.min, array, (axis,))
        if axis is None or axis == 0:
            # Reduce over the list
            return numpy.min(l, axis=0)
        else:
            return l
    else:
        return numpy.max(array)

def mean(array, axis=None):
    if type(array) == dsa.VTKCompositeDataArray:
        if axis == None or axis == 0:
            return sum(array, axis) / array.size
        else:
            return apply_func(numpy.mean, array, (axis,))
    else:
        return numpy.mean(array)

def var(array, axis=None):
    if type(array) == dsa.VTKCompositeDataArray:
        if axis is None:
            tmp = array - mean(array)
            return sum(tmp*tmp) / array.size
        elif axis == 0:
            mn = mean(array, axis=0)
            tmp = array - mn
            return sum(tmp*tmp, axis=0) / shape(array)[0]
        else:
            return apply_func(numpy.var, array, (axis,))
    else:
        return numpy.var(array)

def std(array, axis=None):
    if type(array) == dsa.VTKCompositeDataArray:
        return sqrt(var(array,axis))
    else:
        return numpy.std(array,axis)

def shape(array):
    if type(array) == dsa.VTKCompositeDataArray:
        shp = None
        for a in array.Arrays:
            if a != None:
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
        return shp
    else:
        return numpy.shape(array)

def make_vector(arrayx, arrayy, arrayz=None):
    res = []
    if arrayz is None:
        for ax, ay in itertools.izip(arrayx.Arrays, arrayy.Arrays):
            if ax != None and ay != None:
                res.append(algs.make_vector(ax, ay))
            else:
                res.append(None)
    else:
        for ax, ay, az in itertools.izip(arrayx.Arrays, arrayy.Arrays, arrayz.Arrays):
            if ax != None and ay != None and az != None:
                res.append(algs.make_vector(ax, ay, az))
            else:
                res.append(None)
    return dsa.VTKCompositeDataArray(res)

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
