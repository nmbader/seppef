import pyHypercube


class axis:

    def __init__(self, **kw):
        """Axis
                defaults to n=1, o=0., d=1."""
        self.n = 1
        self.o = 0.
        self.d = 1.
        self.label = ""
        self.unit = ""
        if "n" in kw:
            self.n = kw["n"]
        if "o" in kw:
            self.o = kw["o"]
        if "d" in kw:
            self.d = kw["d"]
        if "label" in kw:
            self.label = kw["label"]
        if "unit" in kw:
            self.unit = kw["unit"]
        if "axis" in kw:
            self.n = kw["axis"].n
            self.o = kw["axis"].o
            self.d = kw["axis"].d
            self.label = kw["axis"].label
            self.unit = kw["axis"].unit

    def __repr__(self):
        """Define print method for class"""
        if self.unit!="":
            return "n=%d\to=%f\td=%f\tlabel=%s\tunit=%s"%(self.n,self.o,self.d,self.label,self.unit)
        elif self.label!="":
            return "n=%d\to=%f\td=%f\tlabel=%s"%(self.n,self.o,self.d,self.label)
        else:
            return "n=%d\to=%f\td=%f"%(self.n,self.o,self.d)


    def getCpp(self):
        """Return a c++ version of the python representation"""
        return pyHypercube.axis(self.n, self.o, self.d, self.label, self.unit)


class hypercube:

    def __init__(self, **kw):
        """initialize with
          - axes=[] A list of Hypercube.axis
          - ns=[] An list of integers (optionally lists of os,ds,labels)
          - hypercube From another hypercube"""
        isSet = False
        self.axes = []
        if "axes" in kw:
            for ax in kw["axes"]:
                self.axes.append(axis(n=ax.n, o=ax.o, d=ax.d, label=ax.label))
                isSet = True
        elif "ns" in kw:
            for n in kw["ns"]:
                self.axes.append(axis(n=n))
                isSet = True
        if "os" in kw:
            for i in range(len(kw["os"]) - len(self.axes)):
                self.axes.append(axis(n=1))
            for i in range(len(kw["os"])):
                self.axes[i].o = kw["os"][i]
        if "ds" in kw:
            for i in range(len(kw["ds"]) - len(self.axes)):
                self.axes.append(axis(n=1))
            for i in range(len(kw["ds"])):
                self.axes[i].d = kw["ds"][i]
        if "labels" in kw:
            for i in range(len(kw["labels"]) - len(self.axes)):
                self.axes.append(axis(n=1))
            for i in range(len(kw["labels"])):
                self.axes[i].label = kw["labels"][i]
        if "units" in kw:
            for i in range(len(kw["units"]) - len(self.axes)):
                self.axes.append(axis(n=1))
            for i in range(len(kw["units"])):
                self.axes[i].unit = kw["units"][i]
        if "hypercube" in kw:
            for i in range(kw["hypercube"].getNdim()):
                a = axis(axis=kw["hypercube"].getAxis(i + 1))
                self.axes.append(a)
        self.buildCpp()

    def clone(self):
        """Clone hypercube"""
        x=hypercube(hypercube=self.getCpp().clone())
        return x
    
    def subCube(self,nw,fw,jw):
        """Return a sub-cube"""
        axes=[]
        for i in range(len(self.axes)):
            axes.append(axis(n=nw[i],o=self.axes[i].o+self.axes[i].d*fw[i],d=self.axes[i].d*jw[i],label=self.axes[i].label,unit=self.axes[i].unit))
        return hypercube(axes=axes)

    def getNdim(self):
        """Return the number of dimensions"""
        return len(self.axes)

    def getAxis(self, i):
        """Return an axis"""
        return self.axes[i - 1]

    def getN123(self):
        """Get the number of elements"""
        n123 = 1
        for ax in self.axes:
            n123 = n123 * ax.n
        return n123

    def getNs(self):
        """Get a list of the sizes of the axes"""
        ns = []
        for a in self.axes:
            ns.append(a.n)
        return ns

    def getCpp(self):
        """Get the c++ representation"""
        return self.cppMode

    
    def checkSame(self,hyper):
        """CHeck to see if hypercube is the same space"""
        return self.cppMode.checkSame(self.cppMode)

    def __repr__(self):
        """Define print method for hypercube class"""
        x=""
        for i in range(len(self.axes)):
            x+="Axis %d: %s\n"%(i+1,str(self.axes[i]))
        return x

    def addAxis(self, axis):
        """Add an axis to the hypercube"""
        self.axes.append(axis)
        self.cppMode.addAxis(axis.getCpp())

    def buildCpp(self):
        """Internal function to build c++ version"""
        ax1 = self.axes[0].getCpp()
        if len(self.axes) > 1:
            ax2 = self.axes[1].getCpp()
            if len(self.axes) > 2:
                ax3 = self.axes[2].getCpp()
                if len(self.axes) > 3:
                    ax4 = self.axes[3].getCpp()
                    if len(self.axes) > 4:
                        ax5 = self.axes[4].getCpp()
                        if len(self.axes) > 5:
                            ax6 = self.axes[5].getCpp()
                            if len(self.axes) >6:
                                ax7 = self.axes[6].getCpp()
                                self.cppMode=pyHypercube.hypercube(
                                    ax1,ax2,ax3,ax4,ax5,ax6,ax7)
                            else:
                                self.cppMode = pyHypercube.hypercube(
                                    ax1, ax2, ax3, ax4, ax5, ax6)
                        else:
                            self.cppMode = pyHypercube.hypercube(
                                ax1, ax2, ax3, ax4, ax5)
                    else:
                        self.cppMode = pyHypercube.hypercube(
                            ax1, ax2, ax3, ax4)
                else:
                    self.cppMode = pyHypercube.hypercube(ax1, ax2, ax3)
            else:
                self.cppMode = pyHypercube.hypercube(ax1, ax2)
        else:
            self.cppMode = pyHypercube.hypercube(ax1)
