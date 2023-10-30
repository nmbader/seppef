import Hypercube
import pytest


def checkAxis(ax, n, o, d, label, unit):
    assert(n == ax.n)
    assert(abs(o - ax.o) < .0001)
    assert(abs(d - ax.d) < .0001)
    assert(label == ax.label)
    assert(unit == ax.unit)


def testAxisN():
    a = Hypercube.axis(n=3)
    checkAxis(a, 3, 0., 1., "", "")


def testAxisNO():
    a = Hypercube.axis(n=3, o=1.)
    checkAxis(a, 3, 1., 1., "", "")


def testAxisNOD():
    a = Hypercube.axis(n=3, o=1., d=2.)
    checkAxis(a, 3, 1., 2., "", "")


def testAxisLabelUnit():
    a = Hypercube.axis(n=3, o=1., d=2., label="blah", unit="dsd")
    checkAxis(a, 3, 1., 2., "blah", "dsd")


def testAxisAxis():
    a = Hypercube.axis(n=3, o=1., d=2., label="blah", unit="dsd")
    b = Hypercube.axis(axis=a)
    checkAxis(b, 3, 1., 2., "blah", "dsd")


def compare(h, ns, os, ds, labels, units):
    for i in range(len(h.axes)):
        assert(ns[i] == h.axes[i].n)
        assert(abs(os[i] - h.axes[i].o) < .001)
        assert(abs(ds[i] - h.axes[i].d) < .001)
        assert(labels[i] == h.axes[i].label)
        assert(units[i] == h.axes[i].unit)


def testGetCpp():
    a = Hypercube.axis(n=3, o=1., d=2., label="blah", unit="dsd")
    b = a.getCpp()
    assert(b.n == 3)
    assert(b.o == 1.)
    assert(b.d == 2.)
    assert(b.label == "blah")
    assert(b.unit == "dsd")


nk = [2, 3, 4, 1, 2]
ok = [1., 2., 3., 4., 5.]
dk = [.2, .3, .4, .4, .5]
labelk = ["one", "two", "three", "four", "five"]
unitk = ["xx", "yy", "zz", "aa", "bb"]


def testHyperN():
    m = Hypercube.hypercube(ns=nk)
    compare(
        m, nk, [
            0., 0., 0., 0., 0.], [
            1., 1., 1., 1., 1.], [
                "", "", "", "", ""], [
                    "", "", "", "", ""])


def testHyperAll():
    m = Hypercube.hypercube(ns=nk, os=ok, ds=dk, labels=labelk, units=unitk)
    compare(m, nk, ok, dk, labelk, unitk)


def testHyperHyper():
    m = Hypercube.hypercube(ns=nk, os=ok, ds=dk, labels=labelk, units=unitk)
    x = Hypercube.hypercube(hypercube=m)
    compare(x, nk, ok, dk, labelk, unitk)


def testNdim():
    m = Hypercube.hypercube(ns=nk, os=ok, ds=dk, labels=labelk, units=unitk)
    assert(m.getNdim() == 5)


def testGetAxis():
    m = Hypercube.hypercube(ns=nk, os=ok, ds=dk, labels=labelk, units=unitk)
    ax = m.getAxis(2)
    checkAxis(ax, 3, 2., .3, "two", "yy")


def testGetNs():
    m = Hypercube.hypercube(ns=nk, os=ok, ds=dk, labels=labelk, units=unitk)
    nt = m.getNs()
    for i in range(len(nt)):
        assert(nt[i] == nk[i])


def testGetCpp2():
    m = Hypercube.hypercube(ns=nk, os=ok, ds=dk, labels=labelk, units=unitk)
    k = m.getCpp()
    for i in range(5):
        ax = k.getAxis(i + 1)
        checkAxis(ax, nk[i], ok[i], dk[i], labelk[i], unitk[i])
