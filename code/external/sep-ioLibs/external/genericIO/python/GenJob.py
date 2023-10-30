import Hypercube
import SepVector
import copy

esizeFromType={"dataFloat":4,"dataDouble":8,"dataInt":4,"dataComplex":8,"dataShort":2,"dataComplexDouble":16,"dataByte":1} 


class regSpace:
    """Base class for regular sampled datasets"""
    def __init__(self,process,mem,minOutputDim,inputJob=None,inputType="dataFloat",outputType="dataFloat", hasInput=True, hasOutput=True):
        """Initialization of the regspace job
            mem - Memory in megabytes needed by the job
            minOutputDim - Minimum output dimensions
            inputJob - Input job 
            inputType - Input data type
            outputType - Output data Type
            process - Function to process data
            inOnly - Only has input
            outOnly - Only has output

        """

        self.process=process
        self._mem=mem
        self._minDim=minOutputDim
        self._inputType=inputType
        self._outputType=outputType
        self._ein=esizeFromType[inputType]
        self._eout=esizeFromType[outputType]
        self._hyperOut=None
        self._hyperIn=None
        self._hasInput=hasInput
        self._hasOutput=hasOutput
        if inputJob:
            if not self.hasInput:
                raise Exception("Job specified as out only and inputJob passed in")
            if not isinstance(regSpace,inputJob):
                raise Exception("Expecting genericJob.regSpace instance to be passed in")
        self._inputJob=inputJob
        self._nw=[]
        self._fw=[]
        self._jw=[]
        self._inputFile=None
        self._outputFile=None
        self._inputBuffer=None
        self._outputBuffer=None
        self._ioBufferIn=None
        self._ioBufferOut=None

    
    def allocateBuffer(self,hyperOut,iwind):
        """Allocate buffer needed for this window
            hyperOut - Output hypercube
            iwind    - Block number
        """
        hyperIn=self._hyperIn.subCube(self._nw[iwind],self._fw[iwind],self._jw[iwind])
        if self._hasOutput:
            self._outputBuffer=self.reallocBuffer(self._outputBuffer,hyperOut,self._outputType)
        if not self._inputJob:
            if self._hasInput:
                self._inputBuffer=self.reallocBuffer(self._inputBuffer,hyperIn,self._inputType)
        else:
            self._inputJob.allocateBuffer(hyperIn,iwind)
    
    def allocateIOBufferOut(self,hyperOut,iwind):
        """Allocate buffer needed for this window
            hyperOut - Output hypercube
            iwind    - Block number

            @returnIOBuffer Out, file object
        """
        self._ioBufferOut=self.reallocBuffer(self._ioBufferOut,hyperOut,self._outputType)
        return self._ioBufferOut,self._outputFile
    
    def allocateIOBufferIn(self,hyperOut,iwind):
        """Allocate buffer needed for this window
            hyperOut - Output hypercube
            iwind    - Block number

            @returns IOBufferIn,File object
        """
        hyperIn=self._hyperIn.subCube(self._nw[iwind],self._fw[iwind],self._jw[iwind])
        if not self._inputJob:
            self._ioBufferIn=self.reallocBuffer(self._ioBufferIn,hyperIn,self._inputType)
            return self._ioBufferIn,self._inputFile,self._nw[iwind],self._fw[iwind],self._jw[iwind]
        else:
            return self._inputJob.allocateIOBufferIn(hyperIn,iwind),self._nw[iwind],self._fw[iwind],self._jw[iwind]



    def reallocBuffer(self,buf,hyper,typ):
        alloc=False
        if not buf:
            alloc=True 
        else:
            nc=buf.getHyper().getNs()
            nn=hyper.getNs()
            for i in range(len(nn)):
                if nn[i] != nc[i]:
                    alloc=True 
        if alloc:
            x= SepVector.getSepVector(hyper,storage=typ)
            return x
        else:
            buf.adjustHyper(hyper)
            return buf

    def deallocateBuffers(self):
        """Deallocate buffers"""
        self._outputBuffer=None
        self._inputBuffer=None
        self._ioBufferIn=None
        self._ioBufferOut=None
        if self._inputJob:
            self._inputJob.deallocateBuffers()

    def swapIOBufferPtrsOut(self):
        """SWap buffer pointers to allow for IO overlap"""
        tmp=self._ioBufferOut
        self._ioBufferOut=self._outputBuffer
        self._outputBuffer=tmp
        return self._ioBufferOut
    
    def swapIObufferPtrsIn(self):
        """SWap buffer pointers to allow for IO overlap"""
        tmp=self._ioBufferIn
        self._ioBufferIn=self._inputBuffer
        self._inputBuffer=tmp


    def processBuffer(self ):
        """Function that actuall does the work"""
        self.process(self._inputBuffer,self._outputBuffer)




    def setOutputFile(self,outputFile):
        """Set the output file"""
        self._outputFile=outputFile
    
    def setInputFile(self,inputFile):
        """Set input file"""
        self._inputFile=inputFile

    def checkLogic(self,first=True):
        if first:
            if self._hasOutput:
                if not self._outputFile:
                    raise Exception("Output file must exist at the end of the chain")
                if self._outputFile.getStorageType()!=self._outputType:
                    raise Exception("Output file type and and outputType set must match")
        if not self._inputJob:
            if self._hasInput:
                if not self._inputFile:
                    raise Exception("Input file must exist at the end of the chain")
                if self._inputFile.getStorageType()!=self._inputType:
                    raise Exception("Input file type and and outputType set must match") 
        else:
            if self._inType != self._inputJob.outType:
                raise Exception("Output/Input type mismatch in pipe")
            self._inputJob.checkLogic(False)

    def returnBaseMemory(self):
        """Return base memory needed by job"""
        x=0
        if self._inputJob:
            x= self._inputJob.returnBaseMemory()
        return x+self._mem

    def calcInputWindow(self,nw,fw,jw):
        """Calculate input window size from output window size

           nw,fw,jw - Standard window parameters"""
        self._nw.append(copy.deepcopy(nw))
        self._fw.append(copy.deepcopy(fw))
        self._jw.append(copy.deepcopy(jw))
        if self._inputJob:
            self._inputJob.calcInputWindow(nw,fw,jw)



    def setCompleteHyperOut(self,hyper):
        """Set output hypercube"""
        self._hyperOut=hyper
        self._hyperIn=self.createHyperIn(self._hyperOut)
        if self._inputJob:
            self._inputJob.setCompleteHyperOut(self._hyperIn)
        

    def getHasInput(self):
        """Return True if job does have input"""
        return self._hasInput

    def getHasOutput(self):
        """Return True if job has an output"""
        return self._hasOutput 

    def getCompleteHyperOut(self):
        """Return hypercube out"""
        if not self._hyperOut:
            raise Exception("Hyperout has not been set")
        return self._hyperOut

    def getCompleteHyperIn(self):
        """Return hypercube in"""
        if not self._hyperIn:
            raise Exception("Hyperin has not been set")
        return self._hyperIn
        
    def createHyperOut(self,hyperIn):
        """Return hypercube out given hypercube in.
        Defaults to same hypercube"""
        return self._hyperIn

    def createHyperIn(self,hyperOut):
        """Return hypercube in given hypercube out.
            Defaults to same hypercube"""
        return hyperOut
    def returnMinDim(self,minDimOut):
        """
            Return the minimum dimension needed in output space.
            Asumes input and ouput same size. Override if not true.   
        """
        dimOut=minDimOut
        if self._inputJob:
            dimOut=max(self._inputJob.returnMinDim(),dimOut)
        return dimOut
            

    def returnSize(self,ns):
        """Return the input dimension given output dimension.
        
            ns - Output size

            Assumes same size

            @return inDim, nbytes
        """
        nbytes=self._eout
        for n in ns:
            nbytes*=n

        if self._inputJob:
           nbytes+=self._inputJob.returnSize(ns)
        else:
            nb=self._ein
            for n in ns:
                nb=nb*n
            nbytes+=nb
        return nbytes

