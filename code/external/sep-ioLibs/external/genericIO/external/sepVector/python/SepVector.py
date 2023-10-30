import pyHypercube
import pySepVector
import Hypercube
import numpy
import pyVector
from sys import version_info


dtypeToSepVecType={
 "float32":"dataFloat",
 "int32":"dataInt",
 "float64":"dataDouble",
 "uint8":"dataByte",
 "complex64":"dataComplex",
 "complex128":"dataComplexDouble"

 }

class vector(pyVector.vectorIC):
    """Generic sepVector class"""

    def __init__(self):
        """Initialize a vector object"""
        if "fromCpp" in self.kw:
            self.cppMode = self.kw["fromCpp"]
        elif "fromHyper" in self.kw:
            self.cppMode = getCppSepVector(self.kw["fromHyper"], self.kw)
        super().__init__(self.getNdArray())

    def getCpp(self):
        """Return the associated cpp object"""
        return self.cppMode

    def clip(self, bclip, eclip):
        """Clip dataset
                bclip - Minimum value to clip to
                eclip - Maximum value to clip to"""
        self.cppMode.clip(bclip, eclip)

    def cent(self, pct, jsamp=1):
        """Calculate the percentile of a dataset
                pct - Percentile of the dataset
                jsamp - Sub-sampling of dartaset"""
        return self.cppMode.cent(pct, jsamp)

    def getStorageType(self):
        """Return type of storage"""
        return self.storage

    def zero(self):
        """Function to zero out a vector"""
        self.cppMode.set(0.)
        return self

    def set(self, val):
        """Function to a vector to a value"""
        self.cppMode.set(val)
        return self

    def getHyper(self):
        """Return the hypercube associated with the vector"""
        return Hypercube.hypercube(hypercube=self.cppMode.getHyper())

    def getNdArray(self):
        """Return a numpy version of the array (same memory"""
        return numpy.array(self.cppMode, copy=False)

    def adjustHyper(self,hyper):
        """Adjust the hypercube associated with vector. Does not reallocate. Must be same dims/size"""
        hyperOld=self.getHyper()
        if hyperOld.getN123() != hyper.getN123():
            raise Exception("Trying to reset with a different sized hyper")
        if len(hyperOld.axes) != len(hyper.axes):
            raise Exception("Trying to reset with a different sized hyper")
        self.cppMode.setHyper(hyper.cppMode)


    def windowInternal(self, **kw):
        """Window a vector return another vector (of the same dimension
            specify min1..min6, max1...max6, f1...f6, j1...j6, n1...n6, or
            any of these by lists.
            Can not specify n and min or max """
        axes = self.getHyper().axes
        ndim = len(axes)
        nw = []
        fw = []
        jw = []
        for i in range(1, ndim + 1):
            nset = False
            fset = False
            jset = False
            if "n%d" % i in kw:
                nset = True
                n = int(kw["n%d" % i])
            if "f%d" % i in kw:
                fset = True
                f = int(kw["f%d" % i])
            if "j%d" % i in kw:
                jset = True
                j = int(kw["j%d" % i])
            biSet = False
            eiSet = False
            if "min%d" % i in kw:
                bi = int(float(kw["min" %
                                  i] - axes[i - 1].o) / axes[i - 1].d + .5)
                biSet = True
            if "max%d" % i in kw:
                ei = int(float(kw["max" %
                                  i] - axes[i - 1].o) / axes[i - 1].d + .5)
                eiSet = True
            if fset:
                if axes[i - 1].n <= f:
                    raise Exception(
                        "invalid f%d=%d parameter n%d of data=%d" %
                        (i, f, i, axes[
                            i - 1].n))
            if nset:
                if axes[i - 1].n < n:
                    raise Exception(
                        "invalid n%d=%d parameter n%d of data=%d" %
                        (i, n, i, axes[
                            i - 1].n))
            if jset and j <= 0:
                raise Exception("invalid j%d=%d " % (i, j))
            if not jset:
                j = 1
            if not nset:
                if not fset:
                    if not biSet:
                        f = 0
                    elif bi < 0 or bi >= axes[i - 1].n:
                        raise Exception("Invalid min%d=%f" %
                                        (i, kw["min%d" % i]))
                    else:
                        f = bi
                if eiSet:
                    if ei <= f or ei >= axes[i - 1].n:
                        raise Exception("Invalid max%d=%f" %
                                        (i, kw["max%d" % i]))
                    else:
                        n = (ei - f - 1) / j + 1
                else:
                    n = (axes[i - 1].n - f - 1) / j + 1

                if not biSet and not eiSet and not jset and not fset:
                    n = axes[i - 1].n
                    j = 1
                    f = 0
            elif not fset:
                if not biSet:
                    f = 0
                elif bi < 0 or bi >= axes[i - 1].n:
                    raise Exception("Invalid min%d=%f" % (i, kw["min%d" % i]))
                else:
                    f = fi
            if axes[i - 1].n < (1 + f + j * (n - 1)):
                raise Exception("Invalid window parameter")
            nw.append(int(n))
            fw.append(int(f))
            jw.append(int(j))
        return self.cppMode.window(nw, fw, jw)

    def checkSame(self, vec2):
        """Function to check if two vectors belong to the same vector space"""
        return self.getCpp().checkSame(vec2.getCpp())


class floatVector(vector):
    """Generic float vector class"""

    def __init__(self, **kw):
        self.kw = kw
        super().__init__()
        self.storage = "dataFloat"

    def norm(self, N=2):
        """Function to compute vector N-norm"""
        self.cppMode.norm(N)
        return self
    def scale(self, sc):
        """Function to scale a vector"""
        self.cppMode.scale(sc)
        return self

    def __repr__(self):
        """Override print method"""
        return "f loatVector\n%s"%str(self.getHyper())

    def calcHisto(self, nelem, mn, mx):
        """Calculate histogram
           mn - Minimum for cell
           mx - Maximum for cell
           nelem - Return number of elements in histogram

           @return Histogram """
        histo = getSepVector(ns=[nelem], storage="dataInt")
        self.cppMode.calcHisto(histo.getCpp(), mn, mx)
        return histo

    def rand(self):
        """Function to fill with random numbers"""
        self.cppMode.rand()
        return self

    def clone(self):
        """Function to clone (deep copy) a vector"""
        return floatVector(fromCpp=self.cppMode.clone())

    def clipVector(self, low, high):
        """Clip vector element by element vec=min(high,max(low,vec))"""
        self.cppMode.clipVector(low.cppMode, high.cppMode)
        return self

    def cloneSpace(self):
        """Funtion tor return the space of a vector"""
        return floatVector(fromCpp=self.cppMode.cloneSpace())

    def scaleAdd(self, vec2, sc1=1., sc2=1.):
        """self = self * sc1 + sc2 * vec2"""
        self.cppMode.scaleAdd(vec2.cppMode, sc1, sc2)
        return self

    def copy(self, vec2):
        """Copy vec2 into self"""
        self.scaleAdd(vec2, 0., 1.)
        return self

    def dot(self, vec2):
        """Compute dot product of two vectors"""
        return self.cppMode.dot(vec2.cppMode)

    def multiply(self, vec2):
        """self = vec2 * self"""
        self.cppMode.mult(vec2.cppMode)
        return self

    def norm(self, nrm=2):
        """Return the norm of a vector"""
        return self.cppMode.norm(nrm)

    def window(self, **kw):
        """Window a vector return another vector (of the same dimension
            specify min1..min6, max1...max6, f1...f6, j1...j6, n1...n6, or
            any of these by lists.
            Can not specify n and min or max """
        return floatVector(fromCpp=self.windowInternal(**kw))


class doubleVector(vector):
    """Generic double vector class"""

    def __init__(self, **kw):
        self.kw = kw
        super().__init__()
        self.storage = "dataDouble"

    def norm(self, N=2):
        """Function to compute vector N - norm"""
        self.cppMode.norm(N)
        return self

    def scale(self, sc):
        """Function to scale a vector"""
        self.cppMode.scale(sc)
        return self
    
    def __repr__(self):
        """Override print method"""
        return "doubleVector\n%s"%str(self.getHyper())
    
    def rand(self):
        """Function to fill with random numbers"""
        self.cppMode.rand()
        return self

    def clipVector(self, low, high):
        """Clip vector element by element vec=min(high,max(low,vec))"""
        self.cppMode.clipVector(low.cppMode, high.cppMode)
        return self

    def calcHisto(self, nelem, mn, mx):
        """Calculate histogram
           mn - Minimum for cell
           mx - Maximum for cell
           nelem - Return number of elements in histogram

           @return Histogram """
        histo = getSepVector(ns=[nelem], storage="dataInt")
        self.cppMode.calcHisto(histo.getCpp(), mn, mx)
        return histo

    def clone(self):
        """Function to clone (deep copy) a vector"""
        return doubleVector(fromCpp=self.cppMode.clone())

    def copy(self, vec2):
        """Copy vec2 into self"""
        self.scaleAdd(vec2, 0., 1.)
        return self

    def cloneSpace(self):
        """Funtion tor return the space of a vector"""
        return doubleVector(fromCpp=self.cppMode.cloneSpace())

    def scaleAdd(self, vec2, sc1=1., sc2=1.):
        """self = self * sc1 + sc2 * vec2"""
        self.cppMode.scaleAdd(vec2.cppMode, sc1, sc2)
        return self

    def dot(self, vec2):
        """Compute dot product of two vectors"""
        return self.cppMode.dot(vec2.cppMode)

    def multiply(self, vec2):
        """self = vec2 * self"""
        self.cppMode.mult(vec2.cppMode)
        return self

    def norm(self, nrm=2):
        """Return the norm of a vector"""
        return self.cppMode.norm(nrm)

    def window(self, **kw):
        """Window a vector return another vector (of the same dimension
            specify min1..min6, max1...max6, f1...f6, j1...j6, n1...n6, or
            any of these by lists.
            Can not specify n and min or max """
        return doubleVector(fromCpp=self.windowInternal(**kw))


class intVector(vector):
    """Generic int vector class"""

    def __init__(self, **kw):
        self.kw = kw
        super().__init__()
        self.storage = "dataInt"

    def __repr__(self):
        """Override print method"""
        return "intVector\n%s"%str(self.getHyper())


class complexVector(vector):
    """Generic complex vector class"""

    def __init__(self, **kw):
        self.kw = kw
        super().__init__()
        self.storage = "dataComplex"

    def cloneSpace(self):
        """Funtion tor return the space of a vector"""
        return complexVector(fromCpp=self.cppMode.cloneSpace())

    def norm(self, N=2):
        """Function to compute vector N-norm"""
        return self.cppMode.norm(N)

    def zero(self):
        """Function to zero out a vector"""
        self.cppMode.set(0.)
        return self

    def multiply(self, vec2):
        """self = vec2 * self"""
        self.cppMode.mult(vec2.cppMode)
        return self

    def rand(self):
        """Fill with random numbers"""
        self.cppMode.rand()
        return self

    def copy(self, vec2):
        """Copy vec2 into self"""
        self.scaleAdd(vec2, 0., 1.)
        return self
    def __repr__(self):
        """Override print method"""
        return "complexVector\n%s"%str(self.getHyper())

    def clone(self):
        """clone a vector"""
        return complexVector(fromCpp=self.cppMode.clone())

    def clipVector(self, low, high):
        """Clip vector element by element vec=min(high,max(low,vec))"""
        self.cppMode.clipVector(low.cppMode, high.cppMode)
        return self

    def dot(self, vec2):
        """Compute dot product of two vectors"""
        return self.cppMode.dot(vec2.cppMode)

    def scale(self, sc):
        """Function to scale a vector"""
        self.cppMode.scale(sc)
        return self

    def scaleAdd(self, vec2, sc1=1., sc2=1.):
        """self = self * sc1 + sc2 * vec2"""
        self.cppMode.scaleAdd(vec2.cppMode, sc1, sc2)
        return self

    def window(self, **kw):
        """Window a vector return another vector (of the same dimension
            specify min1..min6, max1...max6, f1...f6, j1...j6, n1...n6, or
            any of these by lists.
            Can not specify n and min or max """
        return complexVector(fromCpp=self.windowInternal(**kw))

class complexDoubleVector(vector):
    """Generic complex vector class"""

    def __init__(self, **kw):
        self.kw = kw
        super().__init__()
        self.storage = "dataComplexDouble"

    def cloneSpace(self):
        """Funtion tor return the space of a vector"""
        return complexDoubleVector(fromCpp=self.cppMode.cloneSpace())

    def norm(self, N=2):
        """Function to compute vector N-norm"""
        return self.cppMode.norm(N)
    def __repr__(self):
        """Override print method"""
        return "complexDoubleVector\n%s"%str(self.getHyper())
    def zero(self):
        """Function to zero out a vector"""
        self.cppMode.set(0.)
        return self

    def multiply(self, vec2):
        """self = vec2 * self"""
        self.cppMode.mult(vec2.cppMode)
        return self

    def rand(self):
        """Fill with random numbers"""
        self.cppMode.rand()
        return self

    def copy(self, vec2):
        """Copy vec2 into self"""
        self.scaleAdd(vec2, 0., 1.)
        return self

    def clone(self):
        """clone a vector"""
        return complexDoubleVector(fromCpp=self.cppMode.clone())

    def clipVector(self, low, high):
        """Clip vector element by element vec=min(high,max(low,vec))"""
        self.cppMode.clipVector(low.cppMode, high.cppMode)
        return self

    def dot(self, vec2):
        """Compute dot product of two vectors"""
        return self.cppMode.dot(vec2.cppMode)

    def scale(self, sc):
        """Function to scale a vector"""
        self.cppMode.scale(sc)
        return self

    def scaleAdd(self, vec2, sc1=1., sc2=1.):
        """self = self * sc1 + sc2 * vec2"""
        self.cppMode.scaleAdd(vec2.cppMode, sc1, sc2)
        return self

    def window(self, **kw):
        """Window a vector return another vector (of the same dimension
            specify min1..min6, max1...max6, f1...f6, j1...j6, n1...n6, or
            any of these by lists.
            Can not specify n and min or max """
        return complexDoubleVector(fromCpp=self.windowInternal(**kw))


class byteVector(vector):
    """Generic byte vector class"""

    def __init__(self, **kw):
        self.kw = kw
        super().__init__()
        self.storage = "dataByte"

    def calcHisto(self, nelem, mn, mx):
        """Calculate histogram
           mn - Minimum for cell
           mx - Maximum for cell
           nelem - Return number of elements in histogram

           @return Histogram """
        histo = getSepVector(ns=[nelem], storage="dataInt")
        self.cppMode.calcHisto(histo.getCpp(), mn, mx)
        return histo
    def __repr__(self):
        """Override print method"""
        return "byteVector\n%s"%str(self.getHyper())

def getSepVector(*args, **keys):
    """ Get a sepVector object
            Option 1 (supply hypercube):
                    hyper, kw args
            Option 2 (build hypercube):
                    ns = [] - list of sizes
                    os = [] - list of origins
                    ds = [] - list os sampling
                    labels = [] list of labels
                    axes = [] list of axes
            Option 3 (get 2 - D slice)
                    vector - vector to grab from
                    iax1, iax2 - axes to grab
                    rev1, rev1 - whether or not to reverse axes
                    beg, end - beg and end position for all axes(lists)
            storage = StorageType(dataFloat[default], dataComplex,dataComplexDouble,
                        dataDouble,dataInt,dataByte)

            Option 4 (numpy)
                Provide hyper, ns, os, or ds,label,s axes


    """
    myt = "dataFloat"
    haveHyper=False
    haveNumpy=False
    array=None
    if len(args) == 1:
        if isinstance(args[0],Hypercube.hypercube):
            haveHyper=True
            hyper = args[0]
        elif isinstance(args[0],numpy.ndarray):
            haveNumpy=True
            array=args[0]
            if "hyper" in keys:
                hyper=keys["hyper"]
            elif "axes" in keys or "ns" in keys:
                hyper = Hypercube.hypercube(**keys)
            else:
                nt=list(array.shape)
                ns=[]
                for i in range(len(nt)):
                    ns.append(nt[len(nt)-1-i])
                hyper =Hypercube.hypercube(ns=ns)
        else:
            raise Exception("First argument must by a hypercube or numpy array")
    elif len(args) == 0:
        if "axes" in keys or "ns" in keys:
            hyper = Hypercube.hypercube(**keys)
        elif "vector" in keys:
            if "iax1" in keys and "iax2" in keys and "rev1" in keys and "rev2" in keys and "ipos" in keys and "beg" in keys and "end" in keys:
                if isinstance(keys["vector"], floatVector):
                    return floatVector(fromCpp=pySepVector.float2DReg(keys["vector"].getCpp(), keys["iax1"], keys["rev1"],
                                                                      keys["iax2"], keys["rev2"], keys["ipos"], keys["beg"], keys["end"]))
                elif isinstance(keys["vector"], doubleVector):
                    return doubleVector(fromCpp=pySepVector.double2DReg(keys["vector"].cppMode, keys["iax1"], keys["rev1"],
                                                                        keys["iax2"], keys["rev2"], keys["ipos"], keys["beg"], keys["end"]))
                elif isinstance(keys["vector"], intVector):
                    return intVector(fromCpp=pySepVector.int2DReg(keys["vector"].cppMode, keys["iax1"], keys["rev1"],
                                                                  keys["iax2"], keys["rev2"], keys["ipos"], keys["beg"], keys["end"]))
                elif isinstance(keys["vector"], byteVector):
                    return byteVector(fromCpp=pySepVector.byte2DReg(keys["vector"].cppMode, keys["iax1"], keys["rev1"],
                                                                    keys["iax2"], keys["rev2"], keys["ipos"], keys["beg"], keys["end"]))
            else:
                raise Exception("Must supply iax1,iax2,rev1,rev2,ipos,beg,end")
        else:
            raise Exception("Must supply Hypercube,vector  or ns/axes")
    else:
        raise Exception(
            "Only understand 0 or 1 (hypercube) non-keyword arguments")



    if haveNumpy:
        if not str(array.dtype) in dtypeToSepVecType:
            raise Exception("Numpy array not supported",str(array.dtype))
        myt=dtypeToSepVecType[str(array.dtype)]
    else:
        myt="dataFloat"
        if "storage" in keys:
            myt = keys["storage"]

    if myt == "dataFloat":
        x = getFloatVector(hyper)
        y= floatVector(fromCpp=x)
    elif myt == "dataComplexDouble":
        x = getComplexDoubleVector(hyper)
        y=complexDoubleVector(fromCpp=x)
    elif myt == "dataComplex":
        x = getComplexVector(hyper)
        y=complexVector(fromCpp=x)
    elif myt == "dataDouble":
        x = getDoubleVector(hyper)
        y= doubleVector(fromCpp=x)
    elif myt == "dataInt":
        x = getIntVector(hyper)
        y= intVector(fromCpp=x)
    elif myt == "dataByte":
        x = getByteVector(hyper)
        y= byteVector(fromCpp=x)
    else:
        raise Exception("Unknown type %s" % myt)
    if haveNumpy:
        numpy.copyto(y.getNdArray(),array)
    return y

def getCppSepVector(hyper, **keys):
    h = hyper.getCpp()
    myt = "dataFloat"
    if "storage" in keys:
        myt = keys["storage"]
    if myt == "datFloat":
        x = getFloatVector(hyper)
    elif myt == "dataComplex":
        x = getComplexVector(hyper)
    elif myt == "dataComplexDouble":
        x = getComplexDoubleVector(hyper)
    elif myt == "dataDouble":
        x = getDoubleVector(hyper)
    elif myt == "dataInt":
        x = getIntVector(hyper)
    elif myt == "dataByte":
        x = getByteVector(hyper)
    else:
        raise Exception("Unknown type" % myt)


def getFloatVector(hyper):
    h = hyper.getCpp()
    if len(hyper.axes) == 1:
        return pySepVector.float1DReg(h)
    elif len(hyper.axes) == 2:
        return pySepVector.float2DReg(h)
    elif len(hyper.axes) == 3:
        return pySepVector.float3DReg(h)
    elif len(hyper.axes) == 4:
        return pySepVector.float4DReg(h)
    elif len(hyper.axes) == 5:
        return pySepVector.float5DReg(h)
    elif len(hyper.axes) == 6:
        return pySepVector.float6DReg(h)


def getComplexVector(hyper):
    h = hyper.getCpp()

    if len(hyper.axes) == 1:
        return pySepVector.complex1DReg(h)
    elif len(hyper.axes) == 2:
        return pySepVector.complex2DReg(h)
    elif len(hyper.axes) == 3:
        return pySepVector.complex3DReg(h)
    elif len(hyper.axes) == 4:
        return pySepVector.complex4DReg(h)
    elif len(hyper.axes) == 5:
        return pySepVector.complex5DReg(h)
    elif len(hyper.axes) == 6:
        return pySepVector.complex6DReg(h)
def getComplexDoubleVector(hyper):
    h = hyper.getCpp()

    if len(hyper.axes) == 1:
        return pySepVector.complexDouble1DReg(h)
    elif len(hyper.axes) == 2:
        return pySepVector.complexDouble2DReg(h)
    elif len(hyper.axes) == 3:
        return pySepVector.complexDouble3DReg(h)
    elif len(hyper.axes) == 4:
        return pySepVector.complexDouble4DReg(h)
    elif len(hyper.axes) == 5:
        return pySepVector.complexDouble5DReg(h)
    elif len(hyper.axes) == 6:
        return pySepVector.complexDouble6DReg(h)



def getByteVector(hyper):
    h = hyper.getCpp()

    if len(hyper.axes) == 1:
        return pySepVector.byte1DReg(h)
    elif len(hyper.axes) == 2:
        return pySepVector.byte2DReg(h)
    elif len(hyper.axes) == 3:
        return pySepVector.byte3DReg(h)
    elif len(hyper.axes) == 4:
        return pySepVector.byte4DReg(h)
    elif len(hyper.axes) == 5:
        return pySepVector.byte5DReg(h)
    elif len(hyper.axes) == 6:
        return pySepVector.byte6DReg(h)


def getDoubleVector(hyper):
    h = hyper.getCpp()
    if len(hyper.axes) == 1:
        return pySepVector.double1DReg(h)
    elif len(hyper.axes) == 2:
        return pySepVector.double2DReg(h)
    elif len(hyper.axes) == 3:
        return pySepVector.double3DReg(h)
    elif len(hyper.axes) == 4:
        return pySepVector.double4DReg(h)
    elif len(hyper.axes) == 5:
        return pySepVector.double5DReg(h)
    elif len(hyper.axes) == 6:
        return pySepVector.double6DReg(h)


def getIntVector(hyper):
    h = hyper.getCpp()
    if len(hyper.axes) == 1:
        return pySepVector.int1DReg(h)
    elif len(hyper.axes) == 2:
        return pySepVector.int2DReg(h)
    elif len(hyper.axes) == 3:
        return pySepVector.int3DReg(h)
    elif len(hyper.axes) == 4:
        return pySepVector.int4DReg(h)
    elif len(hyper.axes) == 5:
        return pySepVector.int5DReg(h)
    elif len(hyper.axes) == 6:
        return pySepVector.int6DReg(h)


class rectFilter1D(floatVector):

    def __init__(self, n, f, pef=False):
        """Initialize a rectFilter1D
                n - size of box(list 1 - D)
                f - location of the 1 (list 1 - D) (set to - 1)
                pef - Whether or not this is a PEF"""
        self.cppMode = pySepVector.rectFilter1D(n, f, pef)


class rectFilter2D(floatVector):

    def __init__(self, *arg, **kw):
        """Initialize a rectFilter2D
            Option 1:
                n - size of box(list 2 - D)
                f - location of the 1 (list 2 - D) (set to - 1)
                pef- Whether or not the filter is a PEF
            Option 2:
                clone=True/False True-with memory, False without"""
        if len(arg) == 3:
            self.cppMode = pySepVector.rectFilter2D(arg[0], arg[1], arg[2])
        elif len(arg) == 2:
            self.cppMode = pySepVector.rectFilter2D(arg[0], arg[1], False)
        elif "clone" in kw:
            if "space" in kw:
                if kw["space"]:
                    self.cppMode = kw["clone"].cloneSpace()
                else:
                    self.cppMode = kw["clone"].clone()
            else:
                self.cppMode = kw["clone"].clone()
        else:
            raise Exception("Unknown way to create rectFilter2D")

    def clone(self):
        """Function to clone (deep copy) a vector"""
        return rectFilter2D(clone=self.cppMode)

    def cloneSpace(self):
        """Funtion tor return the space of a vector"""
        return rectFilter2D(clone=self.cppMode, space=False)


def readColTextFile(file):
    f = open(file)
    lines = f.readlines()
    array = []
    for line in lines:
        parts = line.split()
        vec1 = []
        for x in parts:
            try:
                y = float(x)
                vec1.append(y)
            except:
                y = x
        array.append(vec1)
    array2 = []
    for i in range(len(array[0])):
        lst = []
        for j in range(len(array)):
            lst.append(array[j][i])
        array2.append(lst)
    return array2
