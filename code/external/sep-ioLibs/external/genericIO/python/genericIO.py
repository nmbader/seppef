import pyGenericIO
import SepVector
import Hypercube
import re
storageConvert = {"dataFloat": pyGenericIO.dataType.dataFloat,
                  "dataInt": pyGenericIO.dataType.dataInt,
                  "dataByte": pyGenericIO.dataType.dataByte,
                  "dataComplex": pyGenericIO.dataType.dataComplex,
                  "dataShort": pyGenericIO.dataType.dataShort,

                  "dataDouble": pyGenericIO.dataType.dataDouble,
                  "dataComplexDouble": pyGenericIO.dataType.dataComplexDouble
                  }
usageConvert = {"usageIn": pyGenericIO.usage_code.usageIn,
                "usageOut": pyGenericIO.usage_code.usageOut,
                "usageInOut": pyGenericIO.usage_code.usageInOut,
                "usageScr": pyGenericIO.usage_code.usageScr}

ioModes = pyGenericIO.ioModes([""])


class pythonParams:

    def __init__(self, inPars):
        """Create a parameter object that reads from a python dictionary"""
        self.pars = {}
        self.inPars = {}
        self._hyper=None
        if isinstance(inPars, list):
            prse = re.compile("(.+)=(.+)")
            for l in inPars:
                m = prse.search(l)
                if m:
                    self.inPars[m.group(1)] = m.group(2)
        elif isinstance(inPars, dict):
            self.inPars = inPars
        else:
            raise Exception("Must initialize with either list or dictionary")
        for k, v in self.inPars.items():
            found = False
            if isinstance(v, int) or isinstance(
                    v, float) or isinstance(v, str):
                vout = str(v)
                found = True
            elif isinstance(v, list):
                d = None
                found = True
                for l in v:
                    if not d:
                        d = l
                    else:
                        d += "," + l
            if found:
                self.pars[k] = vout
        self.cppMode = pyGenericIO.pythonParams(self.pars)
        ioModes.changeParamObj(self.cppMode)
        defaultIO.setParamObj(self.cppMode)

    def getCpp(self):
        return self.cppMode


class argParseParams:

    def __init__(self, inpars):
        """Create a parameter object that uses info from argsparse module"""
        self.pars = {}
        self.pars["junk"] = "jUnK"
        for k, v in inpars.items():
            if v:
                if isinstance(v, str):
                    self.pars[k] = v
            found = False
            if isinstance(v, int) or isinstance(
                    v, float) or isinstance(v, str):
                vout = str(v)
                found = True
            elif isinstance(v, list):
                d = None
                found = True
                for l in v:
                    if not d:
                        d = l
                    else:
                        d += "," + l
            if found:
                self.pars[k] = vout
        self.cppMode = pyGenericIO.pythonParams(self.pars)
        ioModes.changeParamObj(self.cppMode)
        defaultIO.setParamObj(self.cppMode)


class regFile:

    def __init__(self, ioM, tag, **kw):
        """Get regular file python object
                Requiered:
                    iom  - IO mode
                        tag  - Tag for file
                Optional:
                        fromHyper - Hypercube describing file (also must specify storage)
                        storage  - float,complex,byte,double, or int (defaults to float)
                        fromVector - Vector
                        usage - Defaults to IN for from file OUT for everything else
                        ndims - Minimum dimensions for the file
        """
        self.tag = tag
        self.usage = None
        ndimMax = -1
        self._type=ioM.getType()
        if "ndims" in kw:
            ndimMax = kw["ndims"]
        if "usage" in kw:
            if kw["usage"] == "INPUT":
                self.usage = "UsageIn"
            elif kw["usage"] == "output":
                self.usage = "usageOut"
            elif kw["usage"] == "inout":
                self.usage = "usageInOut"
            else:
                raise Exception(
                    "Only understand input,output, and inout for usage")
        self.storage = "dataFloat"
        if "storage" in kw:
            if kw["storage"] == "dataFloat":
                self.storage = "dataFloat"
            elif kw["storage"] == "dataInt":
                self.storage = "dataInt"
            elif kw["storage"] == "dataByte":
                self.storage = "dataByte"
            elif kw["storage"] == "dataDouble":
                self.storage = "dataDouble"
            elif kw["storage"] == "dataComplex":
                self.storage = "dataComplex"
            elif kw["storage"] == "dataComplexDouble":
                self.storage = "dataComplexDouble"
            else:
                raise Exception(
                    "Only understand dataFloat,dataInt,dataDouble,dataByte, dataComplexDouble,and dataComplex for storage not=%s"%storage)

        if "fromHyper" in kw:
            if not isinstance(kw["fromHyper"], Hypercube.hypercube):
                raise Exception(
                    "Must pass a hypercube object when creating a file from a hypercube")
            if not self.usage:
                self.usage = "usageOut"
            elif self.usage == "usageIn":
                raise Exception(
                    "Can not have usageIn when creating from Hypercube")
            if isinstance(ioM,io):
                self.cppMode = ioM.cppMode.getRegFile(
                    self.tag, usageConvert[
                        self.usage], ndimMax)
            else:
                self.cppMode = ioM.getRegFile(
                    self.tag, usageConvert[
                        self.usage], ndimMax)
            self.cppMode.setHyper(kw["fromHyper"].getCpp())
            self.cppMode.setDataType(storageConvert[self.storage])
            self.cppMode.writeDescription()
        elif "fromVector" in kw:
            if not isinstance(kw["fromVector"], SepVector.vector):
                raise Exception(
                    "When creating a file from a vector must be inherited from SepVector.vector")
            self.storage = kw["fromVector"].getStorageType()
            if not self.usage:
                self.usage = "usageOut"
            elif self.usage == "usageIn":
                raise Exception(
                    "Can not have usageIn when creating from Hypercube")
            
            self.cppMode = ioM.getRegFile(
                self.tag, usageConvert[
                    self.usage], ndimMax)
            self.cppMode.setHyper(kw["fromVector"].getCpp().getHyper())
            self.cppMode.setDataType(storageConvert[self.storage])
            self.cppMode.writeDescription()
        else:  # Assuming from file
            if not self.usage:
                self.usage = "usageIn"
            elif self.usage == "usageOut":
                raise Exception(
                    "Can not specify usageOut when creating from a file")
            self.cppMode = ioM.getRegFile(
                self.tag, usageConvert[
                    self.usage], ndimMax)
            self.cppMode.readDescription(ndimMax)
            inv_map = {v: k for k, v in storageConvert.items()}
            self.storage = inv_map[self.cppMode.getDataType()]


    def __repr__(self):
        """Print information about file"""
        x="Regular file type=%s \t "%self._type
        if self.storage=="dataFloat":
            x+="Data type=float\n"
        elif  self.storage=="datInt":
            x+="Data type=integer\n"
        elif  self.storage=="dataByte":
            x+="Data type=byte\n"
        elif  self.storage=="dataDouble":
            x+="Data type=double\n"
        elif  self.storage=="dataComplex":
            x+="Data type=complex\n"
        elif  self.storage=="dataComplexDouble":
            x+="Data type=complexDouble\n"
        else:
            x+="Data type=UKNOWN\n"
        x+="Binary=%s\n"%self.cppMode.getBinary()
        x+=str(self.getHyper())
        return x
    def remove(self):
        """Remove the given dataset"""
        self.cppMode.remove()
        
         
    def getEsize(self):
        """Return element size"""
        if  self.storage=="dataByte":
           return 1
        elif  self.storage=="dataDouble" or self.storage=="dataComplex":
            return 8
        elif  self.storage=="dataComplexDouble":
            return 16
        else:
            return 4
    def getStorageType(self):
        """Return storage type"""
        return self.storage

    def getInt(self, tag, *arg):
        """Get integer from a given IO"""
        if(len(arg) == 1):
            return self.cppMode.getInt(tag, arg[0])
        return self.cppMode.getInt(tag)

    def getFloat(self, tag, *arg):
        """Get a float from a given IO"""
        if(len(arg) == 1):
            return self.cppMode.getFloat(tag, arg[0])
        return self.cppMode.getFloat(tag)

    def getBool(self, tag, *arg):
        """Get a boolean from a given IO"""
        if(len(arg) == 1):
            return self.cppMode.getBool(tag, arg[0])
        return self.cppMode.getBool(tag)

    def getString(self, tag, *arg):
        """Get a string from a given IO"""
        if(len(arg) == 1):
            return self.cppMode.getString(tag, arg[0])
        return self.cppMode.getString(tag)

    def getInts(self, tag, *arg):
        """Get integers from a given IO"""
        if(len(arg) == 1):
            return self.cppMode.getInts(tag, arg[0])
        return self.cppMode.getInts(tag)

    def getInts(self, tag, *arg):
        """Get integers from a given IO"""
        if(len(arg) == 1):
            return self.cppMode.getInts(tag, arg[0])
        return self.cppMode.getInts(tag)

    def getFloats(self, tag, *arg):
        """Get floats from a given IO"""
        if(len(arg) == 1):
            return self.cppMode.getFloats(tag, arg[0])
        return self.cppMode.getFloats(tag)

    def readWindow(self, vec, **kw):
        """Read a window of a file into the vector
                Required:
                  vec - sepVector
                Optional:
                  n,f,j - Standard windowing parameters"""
        if not kw:
            raise Exception("Must supply windowing parameters")
        nw, fw, jw = self.getWindowParam(**kw)
        if not isinstance(vec, SepVector.vector):
            raise Exception("vec must be deriverd SepVector.vector")
        n123 = 1
        for n in nw:
            n123 = n123 * n
        if n123 != vec.getHyper().getN123():
            raise Exception(
                "Data size of vector not the same size as window parameters window=%s vec=%s"%(str(nw),str(vec)))
        elif self.storage == "dataFloat":
            if not isinstance(vec, SepVector.floatVector):
                raise Exception("File type is float vector type is not")
            self.cppMode.readFloatWindow(nw, fw, jw, vec.cppMode)
        elif self.storage == "dataInt":
            if not isinstance(vec, SepVector.intVector):
                raise Exception("File type is int vector type is not")
            self.cppMode.readIntWindow(nw, fw, jw, vec.cppMode)
        elif self.storage == "dataComplex":
            if not isinstance(vec, SepVector.complexVector):
                raise Exception("File type is complex vector type is not")
            self.cppMode.readComplexWindow(nw, fw, jw, vec.cppMode)
        elif self.storage == "dataComplexDouble":
            if not isinstance(vec, SepVector.complexDoubleVector):
                raise Exception("File type is complex double vector type is not")
            self.cppMode.readComplexDoubleWindow(nw, fw, jw, vec.cppMode)
        elif self.storage == "dataByte":
            if not isinstance(vec, SepVector.byteVector):
                raise Exception("File type is nyte vector type is not")
            self.cppMode.readByteWindow(nw, fw, jw, vec.cppMode)
        elif self.storage == "dataDouble":
            if not isinstance(vec, SepVector.doubleVector):
                raise Exception("File type is double vector type is not")
            self.cppMode.readDoubleWindow(nw, fw, jw, vec.cppMode)
        else:
            print("Unknown or unhandled storage type "%self.storage)
    def writeWindow(self, vec, **kw):
        """Write  a window of a file into the vector
                Required:
                  vec - sepVector
                Optional:
                  n,f,j - Standard windowing parameters"""
        nw, fw, jw = self.getWindowParam(**kw)
        if not isinstance(vec, SepVector.vector):
            raise Exception("vec must be deriverd SepVector.vector")
        n123 = 1
        for n in nw:
            n123 = n123 * n
        if n123 != vec.getHyper().getN123():
            raise Exception(
                "Data size of vector not the same size as window parameters")
        elif self.storage == "dataFloat":
            if not isinstance(vec, SepVector.floatVector):
                raise Exception("File type is float vector type is not")
            self.cppMode.writeFloatWindow(nw, fw, jw, vec.cppMode)
        elif self.storage == "dataInt":
            if not isinstance(vec, SepVector.intVector):
                raise Exception("File type is int vector type is not")
            self.cppMode.writeIntWindow(nw, fw, jw, vec.cppMode)
        elif self.storage == "dataComplex":
            if not isinstance(vec, SepVector.complexVector):
                raise Exception("File type is complex vector type is not")
            self.cppMode.writeComplexWindow(nw, fw, jw, vec.cppMode)
        elif self.storage == "dataComplexDouble":
            if not isinstance(vec, SepVector.complexDoubleVector):
                raise Exception("File type is complexDouble vector type is not")
            print("Write complex double window")
            self.cppMode.writeComplexDoubleWindow(nw, fw, jw, vec.cppMode)        
        elif self.storage == "dataByte":
            if not isinstance(vec, SepVector.byteVector):
                raise Exception("File type is nyte vector type is not")
            self.cppMode.writeByteWindow(nw, fw, jw, vec.cppMode)
        elif self.storage == "dataDouble":
            if not isinstance(vec, SepVector.doubleVector):
                raise Exception("File type is double vector type is not")
            self.cppMode.writeDoubleWindow(nw, fw, jw, vec.cppMode)
        else:
            print("Unknown or unhandled storage type "%self.storage)
    def getWindowParam(self, **kw):
        """Return window parameters
                must supply n,f,or j"""
        axes = self.getHyper().axes
        ndim = len(axes)
        for a in ["f", "j", "n"]:
            if a in kw:
                if not isinstance(kw[a], list):
                    raise Exception(
                        "Expecting %s to be a list the dimensions of your file" %
                        a)
                if len(kw[a]) != ndim:
                    raise Exception(
                        "Expecting %s to be a list the dimensions of your file" %
                        a)
        if "j" in kw:
            js = kw["j"]
        else:
            js = [1] * ndim
        if "f" in kw:
            fs = kw["f"]
            for i in range(len(fs)):
                if fs[i] >= axes[i].n:
                    raise Exception(
                        "Invalid f parameter f(%d)>=ndata(%d) for axis %d" %
                        (fs[i], axes[i].n, i + 1))

        else:
            fs = [0] * ndim
        if "n" in kw:
            ns = kw["n"]
            for i in range(len(fs)):
                if ns[i] > axes[i].n:
                    raise Exception(
                        "Invalid n parameter n(%d) > ndata(%d) for axes %d" %
                        (ns[i], axes[i].n, i + 1))
        else:
            ns = []
            for i in range(ndim):
                ns.append(int((axes[i].n - 1 - fs[i]) / js[i] + 1))
        for i in range(ndim):
            if axes[i].n < (1 + fs[i] + js[i] * (ns[i] - 1)):
                raise Exception(
                    "Invalid window parameter (outside axis range) f=%d j=%d n=%d iax=%d ndata=%d" %
                    (fs[i], js[i], ns[i], i + 1, axes[i].n))
        return ns, fs, js

    def close(self):
        """Close file"""
        self.cppMode.close()

    def putInt(self, tag, arg):
        """Put a float into description"""
        self.cppMode.putInt(tag, arg)

    def putFloat(self, tag, arg):
        """Put a float into description"""
        self.cppMode.putFloat(tag, arg)

    def putString(self, tag, arg):
        """Put a string into description"""
        self.cppMode.putString(tag, arg)

    def putBool(self, tag, arg):
        """Put a boolean into description"""
        self.cppMode.putBool(tag, arg)

    def putFloats(self, tag, arg):
        """Put floats into description"""
        self.cppMode.putFloats(tag, arg)

    def putInts(self, tag, arg):
        """Put ints into description"""
        self.cppMode.putInts(tag, arg)

    def getCpp(self):
        """Return cpp object"""
        return self.cppMode

    def getHyper(self):
        """Return hypercube for the file"""
        return Hypercube.hypercube(hypercube=self.cppMode.getHyper())


class AppendFile:
    """Class for append files"""

    def __init__(self, io, tag, vec, maxLength, flush):
        """io -  IO
                tag - filename
                vec - Vector
                maxLength - Maximum number of frames
                flush - How often to flush a file"""
        self.shape = vec.cloneSpace()
        self.vecs = []
        self.flush = flush
        self.hyper = vec.getHyper()
        self.hyper.addAxis(Hypercube.axis(n=maxLength))
        self.nmax = maxLength
        if isinstance(vec, SepVector.floatVector):
            storage = "float"
        elif isinstance(vec, SepVector.doubleVector):
            storage = "double"
        elif isinstance(vec, SepVector.intVector):
            storage = "int"
        elif isinstance(vec, SepVector.complexVector):
            storage = "complex"
        elif isinstance(vec, SepVector.byteVector):
            storage = "byte"
        self.file = regFile(io, tag, fromHyper=self.hyper, storage=storage)
        self.icount = 0

    def addVector(self, vec):
        """Adds a vector to this of vectors being held"""
        self.vecs.append(vec.clone())
        self.icount += 1
        if len(self.vecs) == self.flush:
            return True
        return False

    def flushVectors(self):
        """Return the list of vectors, zero out"""
        vs = []
        for v in self.vecs:
            vs.append(v)
        self.vecs = []
        return vs

    def finish(self, iextra=0):
        """Fix the correct number of frames in a file and update description"""
        self.hyper.axes[len(self.hyper.axes) - 1].n = self.icount + iextra
        self.hyper.buildCpp()
        self.file.cppMode.setHyper(self.hyper.cppMode)
        self.file.cppMode.writeDescription()
        self.nmax = self.icount + iextra


class io:

    def __init__(self, *arg, **kw):
        """Initialize IO"""
        if len(arg) > 0:
            self.cppMode = ioModes.getIO(arg[0])
        else:
            self.cppMode = ioModes.getDefaultIO()
        self.param = self.cppMode.getParamObj()
        if "params" in kw:
            if isinstance(kw["params"], dict):
                self.param.addParams(kw["params"])
            elif isinstance(kw["params"], list):
                self.param.addParams(kw["params"])
        self.appendFiles = {}
        self._files = {}

    def getType(self):
        """Return type"""
        return self.cppMode.getType()

    def getInt(self, tag, *arg):
        """Get integer from a given IO"""
        if(len(arg) == 1):
            return self.param.getInt(tag, arg[0])
        return self.param.getInt(tag)

    def getFloat(self, tag, *arg):
        """Get a float from a given IO"""
        if(len(arg) == 1):
            return self.param.getFloat(tag, arg[0])
        return self.param.getFloat(tag)

    def getBool(self, tag, *arg):
        """Get a boolean from a given IO"""
        if(len(arg) == 1):
            return self.param.getBool(tag, arg[0])
        return self.param.getBool(tag)

    def getString(self, tag, *arg):
        """Get a string from a given IO"""
        if(len(arg) == 1):
            return self.param.getString(tag, arg[0])
        return self.param.getString(tag)

    def getInts(self, tag, *arg):
        """Get integers from a given IO"""
        if(len(arg) == 1):
            return self.param.getInts(tag, arg[0])
        return self.param.getInts(tag)

    def getInts(self, tag, *arg):
        """Get integers from a given IO"""
        if(len(arg) == 1):
            return self.param.getInts(tag, arg[0])
        return self.param.getInts(tag)

    def getFloats(self, tag, *arg):
        """Get floats from a given IO"""
        if(len(arg) == 1):
            return self.param.getFloats(tag, arg[0])
        return self.param.getFloats(tag)

    def getRegFile(self, tag, **kw):
        """Get regular file python object
                Requiered:
                        tag  - Tag for file
                Optional:
                        fromHyper - Hypercube describing file (also must specify storage)
                        storage  - float,complex,byte,double, or int (defaults to float)
                        fromVector - Vector
                        ndims   - Minimum dimension of file
        """
        x = regFile(self.cppMode, tag, **kw)
        self._files[tag] = x
        return x

    def getFile(self, tag):
        """Return file assumed it has been read through io

            tag - Tag associated with file
        """
        if tag not in self._files:
            raise Exception("Requested tag not loaded into this IO")
        return self._files[tag]

    def getVector(self, tag, **kw):
        """Get vector from a file and read its contents
           Optional
             ndims - Force the hypercube to at least ndim axes"""
        file = self.getRegFile(tag, **kw)
        self._files[tag] = file
        hyper = file.getHyper()
        nw = file.getHyper().getNs()
        fw = [0] * len(nw)
        jw = [1] * len(nw)

        vec = SepVector.getSepVector(hyper, storage=file.storage)
        if file.storage == "dataFloat":
            file.getCpp().readFloatWindow(nw, fw, jw, vec.getCpp())
        elif file.storage == "dataComplex":
            file.getCpp().readComplexWindow(nw, fw, jw, vec.getCpp())
        elif file.storage == "dataComplexDouble":
            file.getCpp().readComplexDoubleWindow(nw, fw, jw, vec.getCpp())
        elif file.storage == "dataByte":
            file.getCpp().readByteWindow(nw, fw, jw, vec.getCpp())
        elif file.storage == "dataInt":
            file.getCpp().readIntWindow(nw, fw, jw, vec.getCpp())
        elif file.storage == "dataDouble":
            file.getCpp().readDoubleWindow(nw, fw, jw, vec.getCpp())
        file.close()
        return vec

    def writeVector(self, tag, vec):
        """Write entire sepVector to disk
           tag - File to write to
           vec - Vector to write"""
        file = regFile(self.cppMode, tag, fromVector=vec)
        nw = file.getHyper().getNs()
        fw = [0] * len(nw)
        jw = [1] * len(nw)
        if file.storage == "dataFloat":
            file.getCpp().writeFloatWindow(nw, fw, jw, vec.getCpp())
        elif file.storage == "dataComplex":
            file.getCpp().writeComplexWindow(nw, fw, jw, vec.getCpp())
        elif file.storage == "dataComplexDouble":
            file.getCpp().writeComplexDoubleWindow(nw, fw, jw, vec.getCpp())
        elif file.storage == "dataByte":
            file.getCpp().writeByteWindow(nw, fw, jw, vec.getCpp())
        elif file.storage == "dataInt":
            file.getCpp().writeIntWindow(nw, fw, jw, vec.getCpp())
        elif file.storage == "dataDouble":
            file.getCpp().writeDoubleWindow(nw, fw, jw, vec.getCpp())

    def writeVectors(self, file, vecs, ifirst):
        """Write a collection of vectors to a file
                file - regFile
                vecs - Vectors
                ifirst - Position in file for first vector"""
        nw = file.getHyper().getNs()
        fw = [0] * len(nw)
        jw = [1] * len(nw)
        nw[len(nw) - 1] = 1
        iouter = ifirst
        for vec in vecs:
            fw[len(fw) - 1] = iouter
            iouter += 1
            if file.storage == "dataFloat":
                file.getCpp().writeFloatWindow(nw, fw, jw, vec.getCpp())
            elif file.storage == "dataComplex":
                file.getCpp().writeComplexWindow(nw, fw, jw, vec.getCpp())
            elif file.storage == "dataComplexDouble":
                file.getCpp().writeComplexDoubleWindow(nw, fw, jw, vec.getCpp())
            elif file.storage == "dataByte":
                file.getCpp().writeByteWindow(nw, fw, jw, vec.getCpp())
            elif file.storage == "dataInt":
                file.getCpp().writeIntWindow(nw, fw, jw, vec.getCpp())
            elif file.storage == "dataDouble":
                file.getCpp().writeDoubleWindow(nw, fw, jw, vec.getCpp())

    def appendVector(self, tag, vec, maxLength=1000, flush=1):
        """Write entire sepVector to disk
           tag - File to write to
           vec - Vector to write
           maxLength - Maximum number of appended frames
           flush - How often to flush the files"""
        if tag not in self.appendFiles:
            self.appendFiles[tag] = AppendFile(
                self.cppMode, tag, vec, maxLength, flush)
        if self.appendFiles[tag].addVector(vec):
            vs = self.appendFiles[tag].flushVectors()
            if self.appendFiles[tag].icount > self.appendFiles[tag].nmax:
                self.appendFiles[tag].finish()
            self.writeVectors(
                self.appendFiles[tag].file,
                vs, self.appendFiles[tag].icount - len(vs))

    def closeAppendFile(self, tag):
        """Close an append file and fix the description to the right number of frames"""
        if tag not in self.appendFiles:
            raise Exception("No record of appended file")
        vs = self.appendFiles[tag].flushVectors()
        self.writeVectors(
            self.appendFiles[tag].file,
            vs, self.appendFiles[tag].icount - len(vs))
        self.appendFiles[tag].finish(0)

    def setParamObj(self, par):
        """Set parameter object"""
        self.param = par
defaultIO = io()
