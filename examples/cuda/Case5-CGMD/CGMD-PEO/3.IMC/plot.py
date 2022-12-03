#!/usr/bin/python 
import MagicTools as MT
MT.Deviation(['01.magic.out'])
rdfs_imc = MT.ReadMagiC(['01.magic.out'], DFType='RDF')
rdf_ref = MT.ReadRDF('peo10-10.rdf', quiet=True)
rdf_ref_imc = MT.ReadMagiC('01.magic.out', DFType='RDFref', quiet=True)
rdf_ref.SetPlotProperty('linewidth', 3)
for rdf_ in rdfs_imc[0:10]:
    rdf_.SetPlotProperty('linestyle', '--')
rdfs_imc_filtered = [rdfs_imc[i] for i in [0,2,4,6,8,9]]
MT.MultPlot(rdfs_imc_filtered + [rdf_ref], hardcopy=False)
#pot_corr_imc = MT.ReadMagiC(['02.magic.out','03.magic.out'], DFType='PotCorr', iters = (1,3,5,7,9), quiet=True)
# make potentials from 02.magic.out dashed
#for pot_ in pot_corr_imc[0:5]:
#    pot_.SetPlotProperty('linestyle', '--') 

#MT.MultPlot(pot_corr_imc)
