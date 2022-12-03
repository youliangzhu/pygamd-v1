#!/usr/bin/python
import MagicTools as MT

system = MT.GALAMOSTTopology(inpMagiC='01.magic.inp', xmol='03.peo10-10.i001.start.xmol',eps=1.0)
pot = MT.ReadPot('03.peo10-10.i010.pot', quiet=True)
pot.CutTail(RcutNB=15.0)
MT.PotsExport2GALAMOST(pot, npoints=1500, Rmaxtable=1.5, noplot=True, sigma=0.55)

MT.dcd2xtc("trj.dcd","topology.xml")
