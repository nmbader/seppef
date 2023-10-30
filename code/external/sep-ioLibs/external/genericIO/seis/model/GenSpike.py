#!/usr/bin/env python3
import argparse
import genericIO
import GenJob
import GenSplit
import Hypercube
import math
import numpy as np
from numba import jit,prange

class spikeJob(GenJob.regSpace):
    def __init__(self,outputType,events):
        """Intialize object

            outputType - Output type
            events - List of events
        """
        super().__init__(self.spikeBuf,0,0,outputType=outputType,hasInput=False)
        self._events=events
    
    def spikeBuf(self,ina,outa):
        """Convert a buffer from one type to another

        ina - Input vector
        outa - Output vector
        """
    
        f=[0]*6
        n=[1]*6

        axesVec=outa.getHyper().axes
        axesOut=self.getCompleteHyperOut().axes
        for i in range(len(axesOut)):
            f[i]=int(round((axesVec[i].o-axesOut[i].o)/axesOut[i].d))
            n[i]=axesVec[i].n 
        outa.scale(0.)
        outN=np.reshape(outa.getNdArray(),(n[5],n[4],n[3],n[2],n[1],n[0]))
        for ev in self._events:
            fill(outN,ev._mag,f[0],f[1],f[2],f[3],f[4],f[5],ev._k1,ev._k2,ev._k3,ev._k4,ev._k5,ev._k6)

@jit(nopython=True)
def fill(ar,mag,f1,f2,f3,f4,f5,f6,k1,k2,k3,k4,k5,k6):
    for i6 in k6:
        j6=i6-f6
        if j6 >=0 and j6 < ar.shape[0]:
            for i5 in k5:
                j5=i5-f5
                if j5 >=0 and j5 <ar.shape[1]:
                    for i4 in k4:
                        j4=i4-f4
                        if j4 >=0 and j4 < ar.shape[2]:
                            for i3 in k3:
                                j3=i3-f3
                                if j3>=0 and j3 < ar.shape[3]:
                                    for i2 in k2:
                                        j2=i2-f2
                                        if j2>=0 and j2<ar.shape[4]:
                                            for i1 in k1:
                                                j1=i1-f1
                                                if j1 >=0 and j1<ar.shape[5]:
                                                    ar[j6,j5,j4,j3,j2,j1]=mag
    

def hyperFromArgs(args):
    """Create a hypercube from arguments"""
    axs=[]
    axs.append(Hypercube.axis(n=args.n1,o=args.o1,d=args.d1,label=args.label1))
    axs.append(Hypercube.axis(n=args.n2,o=args.o2,d=args.d2,label=args.label2))
    axs.append(Hypercube.axis(n=args.n3,o=args.o3,d=args.d3,label=args.label3))
    axs.append(Hypercube.axis(n=args.n4,o=args.o4,d=args.d4,label=args.label4))
    axs.append(Hypercube.axis(n=args.n5,o=args.o5,d=args.d5,label=args.label5))
    axs.append(Hypercube.axis(n=args.n6,o=args.o6,d=args.d6,label=args.label6))

    return Hypercube.hypercube(axes=axs)

class event:
    """Create an event for the spike"""
    def __init__(self,ievent,hyper,storage,k1,k2,k3,k4,k5,k6,mag):
        """Event 
            hyper - Hypercube of output
            storage - Storage type
            ievent - Event number
            k1   - Location first axis (fortran)
            k2   - Location second axis (fortran)
            k3   - Location third axis (fortran)
            k4   - Location fourth axis (fortran)
            k5   - Location fifth axis (fortran)
            k6   - Location sixth axis (fortran)
  
            mag  - Magnitude for event"""
        if ievent>= len(mag):
            raise Exception("Event number larger than the magnitude list")


        if storage=="dataInt" or storage=="dataInt" or storage=="dataShort":
            try:
                self._mag=int(mag[ievent])
            except:
                raise Exception("Trouble converting to int:",self._mag[ievent])
            if storage=="dataByte":
                if self._mag < 0 or self._mag >255:
                    raise Exception("For byte expecting values between 0 and 255")
            if storage=="dataShort":
                if self._mag < -32767 or self._mag >32767:
                    raise Exception("For short expecting values between -32767 and 32767")
            if storage=="dataInt":
                if self._mag < -2147483648 or self._mag >2147483648:
                    raise Exception("For short expecting values between -2147483648 and 2147483648")
        elif storage=="dataFloat" or storage=="dataDouble":
            try:
                self._mag=float(mag[ievent])
            except:
                raise Exception("Trouble converting to float:",self._mag[ievent])
        elif storage=="dataComplex" or storage=="dataComplexDouble":
            try:
                self._mag=complex(mag[ievent])
            except:
                raise Exception("Trouble converting to complex:",self._mag[ievent])   
        else:
            raise Exception("Unknown data type ",storage)             

        if ievent >= len(k1):
            self._k1=np.asarray(range(hyper.axes[0].n),dtype=np.int32)
        else:
            k1[ievent]=int(k1[ievent])
            if k1[ievent]<=0 or k1[ievent] > hyper.axes[0].n:
                raise Exception("Illegal value of k1=%d must be between 1 and %d"%(k1[ievent],hyper.axes[0].n))
            self._k1=np.asarray([k1[ievent]-1],dtype=np.int32)

        if ievent >= len(k2):
            self._k2=np.asarray(range(hyper.axes[1].n),dtype=np.int32)
        else:
            k2[ievent]=int(k2[ievent])
            if k2[ievent]<=0 or k2[ievent] > hyper.axes[1].n:
                raise Exception("Illegal value of k2=%d must be between 1 and %d"%(k2[ievent],hyper.axes[1].n))
            self._k2=np.asarray([k2[ievent]-1],dtype=np.int32)

        if ievent >= len(k3):
            self._k3=np.asarray(range(hyper.axes[2].n),dtype=np.int32)
        else:
            k3[ievent]=int(k3[ievent])
            if k3[ievent]<=0 or k3[ievent] > hyper.axes[2].n:
                raise Exception("Illegal value of k3=%d must be between 1 and %d"%(k3[ievent],hyper.axes[2].n))
            self._k3=np.asarray([k3[ievent]-1],dtype=np.int32)

        if ievent >= len(k4):
            self._k4=np.asarray(range(hyper.axes[3].n),dtype=np.int32)
        else:
            k4[ievent]=int(k4[ievent])
            if k4[ievent]<=0 or k4[ievent] > hyper.axes[3].n:
                raise Exception("Illegal value of k4=%d must be between 1 and %d"%(k4[ievent],hyper.axes[3].n))
            self._k4=np.asarray([k4[ievent]-1],dtype=np.int32)
        if ievent >= len(k5):
            self._k5=np.asarray(range(hyper.axes[4].n),dtype=np.int32)
        else:
            k5[ievent]=int(k5[ievent])

            if k5[ievent]<=0 or k5[ievent] > hyper.axes[4].n:
                raise Exception("Illegal value of k5=%d must be between 1 and %d"%(k5[ievent],hyper.axes[4].n))
            self._k5=np.asarray([k5[ievent]-1],dtype=np.int32)

        if ievent >= len(k6):
            self._k6=np.asarray(range(hyper.axes[5].n),dtype=np.int32)
        else:
            k6[ievent]=int(k6[ievent])
            if k6[ievent]<=0 or k6[ievent] > hyper.axes[5].n:
                raise Exception("Illegal value of k6=%d must be between 1 and %d"%(k6[ievent],hyper.axes[5].n))
            self._k6=np.asarray([k6[ievent]-1],dtype=np.int32)




        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make delta functions and impulsive plane waves ')
    parser.add_argument('output', metavar='Output', type=str,
                        help='Output file')                   
    parser.add_argument("--io", type=str,choices=[@GEN_IO_TYPES@], help='IO type. Defaults to defaultIO')
    parser.add_argument("--storage",type=str,choices=["dataByte","dataInt","dataFloat","dataComplex","dataShort",
    "dataComplexDouble","dataDouble"],default="dataFloat")
    parser.add_argument("--memory",type=float,help="Memory in terms of GB",default=.2)
    parser.add_argument("--print_pct",type=float,help="Print progress every X pct (above 100 means no printing)",default=101)
    parser.add_argument("--k1",type=int,help="Fortran index for delta function (first axis), if absent, all axis values, set to value",nargs="+",default=[])
    parser.add_argument("--k2",type=int,help="Fortran index for delta function (second axis), if absent, all axis values, set to value",nargs="+",default=[])
    parser.add_argument("--k3",type=int,help="Fortran index for delta function (third axis), if absent, all axis values, set to value",nargs="+",default=[])
    parser.add_argument("--k4",type=int,help="Fortran index for delta function (fourth axis), if absent, all axis values, set to value",nargs="+",default=[])
    parser.add_argument("--k5",type=int,help="Fortran index for delta function (fift axis), if absent, all axis values, set to value",nargs="+",default=[])
    parser.add_argument("--k6",type=int,help="Fortran index for delta function (six axis), if absent, all axis values, set to value",nargs="+",default=[])
    parser.add_argument("--n1",type=int,help="Number of elements first axis",default=1)
    parser.add_argument("--n2",type=int,help="Number of elements second axis",default=1)
    parser.add_argument("--n3",type=int,help="Number of elements third axis",default=1)
    parser.add_argument("--n4",type=int,help="Number of elements fourth axis",default=1)
    parser.add_argument("--n5",type=int,help="Number of elements fifth axis",default=1)
    parser.add_argument("--n6",type=int,help="Number of elements sixth axis",default=1)
    parser.add_argument("--o1",type=float,help="First location first axis",default=0.)
    parser.add_argument("--o2",type=float,help="First location second axis",default=0.)
    parser.add_argument("--o3",type=float,help="First location third axis",default=0.)
    parser.add_argument("--o4",type=float,help="First location fourth axis",default=0.)
    parser.add_argument("--o5",type=float,help="First location fifth axis",default=0.)
    parser.add_argument("--o6",type=float,help="First location sixth axis",default=0.)
    parser.add_argument("--d1",type=float,help="Sampling first axis",default=1.)
    parser.add_argument("--d2",type=float,help="Sampling second axis",default=1.)
    parser.add_argument("--d3",type=float,help="Sampling third axis",default=1.)
    parser.add_argument("--d4",type=float,help="Sampling fourth axis",default=1.)
    parser.add_argument("--d5",type=float,help="Sampling fifth axis",default=1.)
    parser.add_argument("--d6",type=float,help="Sampling sixth axis",default=1.)
    parser.add_argument("--label1",type=str,help="Label first axis",default="")
    parser.add_argument("--label2",type=str,help="Label second axis",default="")
    parser.add_argument("--label3",type=str,help="Label third axis",default="")
    parser.add_argument("--label4",type=str,help="Label fourth axis",default="")
    parser.add_argument("--label5",type=str,help="Label fifth axis",default="")
    parser.add_argument("--label6",type=str,help="Label sixth axis",default="")
    parser.add_argument("--magnitude",required=True,type=list,help="Magnitude of spikes")

    args = parser.parse_args()
    ioOut=genericIO.defaultIO

    if args.io:
        ioOut=genericIO.io(args.io)

    hyper=hyperFromArgs(args)
    events=[]
    for i in range(len(args.magnitude)):
        events.append(event(i,hyper,args.storage,args.k1,args.k2,args.k3,args.k4,args.k5,args.k6,args.magnitude))

    outFile=genericIO.regFile(ioOut,args.output,storage=args.storage,fromHyper=hyper)
    job=spikeJob(outFile.getStorageType(),events)
    job.setOutputFile(outFile)
    job.setCompleteHyperOut(outFile.getHyper())
    split=GenSplit.serialRegSpace(job, args.memory)
    split.loop(args.print_pct)









