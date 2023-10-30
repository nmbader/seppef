#!/usr/bin/env python3
import argparse
import genericIO
import GenJob
import GenSplit
import numpy as np
from numba import jit,prange

class cpJob(GenJob.regSpace):
    def __init__(self,fileType):
        """Intialize object
        fileType - Type of vector
        
        """
        super().__init__(self.cpBuf,0,0,inputType=fileType,outputType=fileType)
 
    
    def cpBuf(self,ina,outa):
        """Copy a vector

        ina - Input vector
        outa - Output vector
        """
        outa.scaleAdd(ina,0.,1.)
       
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Move file')
    parser.add_argument('input', metavar='Input', type=str,
                        help='Input file')
    parser.add_argument('output', metavar='Output', type=str,
                        help='Output file')                   
    parser.add_argument("--ioIn", type=str,choices=[@GEN_IO_TYPES@], help='IO type. Defaults to defaultIO')
    parser.add_argument("--ioOut", type=str,choices=[@GEN_IO_TYPES@], help='IO type. Defaults to defaultIO')
    parser.add_argument("--memory",type=float,help="Memory in terms of GB",default=.2)
    parser.add_argument("--print_pct",type=float,help="Print progress every X pct (above 100 means no printing)",default=101)

    args = parser.parse_args()

    ioIn=genericIO.defaultIO
    ioOut=ioIn

    if args.ioIn:
        ioIn=genericIO.io(args.ioIn)

    if args.ioOut:
        ioOut=genericIO.io(args.ioOut)

    inFile=ioIn.getRegFile(args.input)
    outFile=genericIO.regFile(ioOut,args.output,storage=inFile.getStorageType(),fromHyper=inFile.getHyper())
    job=cpJob(inFile.getStorageType())
    job.setOutputFile(outFile)
    job.setCompleteHyperOut(outFile.getHyper())
    job.setInputFile(inFile)
    split=GenSplit.serialRegSpace(job, args.memory)
    split.loop(args.print_pct)

    inFile.remove()









