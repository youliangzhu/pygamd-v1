#!/user/bin/python 
import MagicTools as MT 
# Simple read and plot, no additional parameters provided
#rdf_ref = MT.ReadRDF('peo10-10.rdf', quiet=True)
#MT.MultPlot(rdf_ref)

# More advanced plot with custom templates for legend and title strings, and larger fontsizes 
# See the HowToPlot tutorial for more examples
rdf_ref = MT.ReadRDF('peo10-10.rdf', quiet=True, Name='Reference RDF')
MT.MultPlot(rdf_ref, title_template='DFsetName', legend_template='Type.Name', title_fontsize=20, legend_fontsize=20)
system = MT.System()
MT_PEO = MT.MolType('peo10-10.CG.mcm', system)
system.Box = [100, 100, 100] # some large enough box
system.WriteLAMMPSData('LAMMPS.data')

