#!/usr/bin/python
import MagicTools as MT 
rdf_repro_GALAMOST = MT.ReadRDF('peo10-10.repro.galamost.rdf')
rdf_ref = MT.ReadRDF('peo10-10.rdf',quiet=True)
rdf_ref.SetPlotProperty('linewidth',3)
MT.MultPlot([ rdf_ref, rdf_repro_GALAMOST])
