#!/usr/bin/env python3
import argparse
import genericIO
import GenJob
import GenSplit
import numpy as np
from numba import jit,prange

class ddJob(GenJob.regSpace):
    def __init__(self,inputType,outputType,real):
        """Intialize object

            inputType - Input type
            outputType - Output type
            real - Convert real (rather than imaginary) of complex number to/from 
        """
        #super().__init__(self.convertBuf,0,0,inputType=inputType,outputType=outputType)
        self._real=real

    
    def convertBuf(self,ina,outa):
        """Convert a buffer from one type to another

        ina - Input vector
        outa - Output vector
        """

        n123=ina.getHyper().getN123()
        inN=np.reshape(ina.getNdArray(),(n123,))
        outN=np.reshape(outa.getNdArray(),(n123,))

        if self._inputType=="dataComplex" or self._inputType=="dataComplexDouble":
            if self._outputType=="dataShort":
                complex2Short(inN,outN,self._real)
            elif self._outputType=="dataInt":
                complex2Int(inN,outN,self._real)
            elif self._outputType=="dataFloat" or self._outputType=="dataDouble":
                complex2Real(inN,outN,self._real)
            elif self._outputType=="dataComplex" or self._outputType=="dataComplexDouble":
                complex2Complex(inN,outN)
            elif self._outType=="dataByte":
                complex2Byte(inN,outN,self._real)
            else:
                print("Uknown conversion %s"%self._outputType)
        else:
            if self._outputType=="dataShort":
                real2Short(inN,outN)
            elif self._outputType=="dataInt":
                real2Int(inN,outN)
            elif self._outputType=="dataFloat" or self._outputType=="dataDouble":
                real2Real(inN,outN)
            elif self._outputType=="dataComplex" or self._dataType=="dataComplexDouble":
                real2Complex(inN,outN,self_real)
            elif self._outputType=="dataByte":
                real2Byte(inN,outN)
            else:
                print("Unknown conversion %s to %s"%(self._inputType,self._outputType))

@jit(nopython=True, parallel=True)
def complex2Short(inA,outA,realFlag):
    if realFlag:
        for i in prange(outA.shape[0]):
            outA[i]=min(255,max(0, int(.5+real(inA[i]))))
    else:
        for i in prange(outA.shape[0]):
            outA[i]=min(255,max(0, int(.5+imag(inA[i]))))
            
@jit(nopython=True, parallel=True)
def complex2Short(inA,outA,realFlag):
    if realFlag:
        for i in prange(outA.shape[0]):
            outA[i]=min(32567,max(-32567, int(.5+real(inA[i]))))
    else:
        for i in prange(outA.shape[0]):
            outA[i]=min(32567,max(-32567, int(.5+imag(inA[i]))))

@jit(nopython=True, parallel=True)
def complex2Int(inA,outA,realFlag):
    if realFlag:
        for i in prange(outA.shape[0]):
            outA[i]=min(2147483647,max(-2147483647, int(.5+real(inA[i]))))
    else:
        for i in prange(outA.shape[0]):
            outA[i]=min(2147483647,max(-2147483647, int(.5+imag(inA[i]))))


@jit(nopython=True, parallel=True)
def complex2Real(inA,outA,realFlag):
    if realFlag:
        for i in prange(outA.shape[0]):
            outA[i]=real(inA[i])
    else:
        for i in prange(outA.shape[0]):
            outA[i]=imag(inA[i])


@jit(nopython=True, parallel=True)
def complex2Complex(inA,outA):
    for i in prange(outA.shape[0]):
        outA[i]=inA[i]

@jit(nopython=True, parallel=True)
def real2Byte(inA,outA):
    for i in prange(outA.shape[0]):
        outA[i]=min(255,max(0, int(.5+inA[i])))

@jit(nopython=True, parallel=True)
def real2Short(inA,outA):
    for i in prange(outA.shape[0]):
        outA[i]=min(32567,max(-32567, int(.5+inA[i])))

@jit(nopython=True, parallel=True)
def real2Int(inA,outA):
    for i in prange(outA.shape[0]):
        outA[i]=min(2147483647,max(-2147483647, int(.5+inA[i])))

@jit(nopython=True, parallel=True)
def real2Real(inA,outA):
    for i in prange(outA.shape[0]):
        outA[i]=inA[i]


@jit(nopython=True, parallel=True)
def real2Complex(inA,outA,realFlag):
    if realFlag:
        for i in prange(outA.shape[0]):
            outA[i]=inA[i]
    else:
        for i in prange(outA.shape[0]):
            outA[i]=complex(0,inA[i])



        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert file types')
    parser.add_argument('input', metavar='Input', type=str,
                        help='Input file')
    parser.add_argument('output', metavar='Output', type=str,
                        help='Output file')                   
    parser.add_argument("--ioIn", type=str,choices=[@GEN_IO_TYPES@], help='IO type. Defaults to defaultIO')
    parser.add_argument("--ioOut", type=str,choices=[@GEN_IO_TYPES@], help='IO type. Defaults to defaultIO')
    parser.add_argument("--storage",type=str,choices=["dataByte","dataInt","dataFloat","dataComplex","dataShort",
    "dataComplexDouble","dataDouble"],default="dataFloat")
    parser.add_argument("--memory",type=float,help="Memory in terms of GB",default=.5)
    parser.add_argument("--real",type=bool, help="Convert float to real portion of complex")
    parser.add_argument("--print_pct",type=float,help="Print progress every X pct (above 100 means no printing)",default=101)

    args = parser.parse_args()

    ioIn=genericIO.defaultIO
    ioOut=ioIn

    if args.ioIn:
        ioIn=genericIO.io(args.ioIn)

    if args.ioOut:
        ioOut=genericIO.io(args.ioOut)

    inFile=ioIn.getRegFile(args.input)
    outFile=genericIO.regFile(ioOut,args.output,storage=args.storage,fromHyper=inFile.getHyper())
    job=ddJob(inFile.getStorageType(),outFile.getStorageType(),args.real)
    job.setOutputFile(outFile)
    job.setCompleteHyperOut(outFile.getHyper())
    job.setInputFile(inFile)
    split=GenSplit.serialRegSpace(job, args.memory)
    split.loop(args.print_pct)









