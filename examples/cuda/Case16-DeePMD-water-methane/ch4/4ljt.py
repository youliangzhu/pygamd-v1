import pygamd
import sys

"""
    4ljt.py - Script for running molecular dynamics simulations of methane using PyGAMD's 

    Usage:
        Please ensure to provide valid model file and POSCAR file paths.
        How to run: python 4ljt.py < MACHINE LEARNING FORCE FIELD MODEL file path > < POSCAR file path >

    Parameters:
        info: Information object containing simulation details.
        model_file: Path to the machine learning force field model file.
        poscar_file: Path to the POSCAR file containing atomic positions.
"""

model_file = sys.argv[1]
poscar_file = sys.argv[2]

xml = pygamd.s_napshot.read("UNK.xml")

app = pygamd.application.dynamics(info=xml, dt=0.001)
fn = pygamd.force.dpk(info=xml, model_file=model_file, poscar_file=poscar_file)
app.add(fn)

inn = pygamd.integration.nve(info=xml, group='all')
app.add(inn)

dd = pygamd.dump.data(info=xml, group='all', file='data.log', period=1)
app.add(dd)

dm = pygamd.dump.xml(info=xml, group='all', file='dpk.xml', period=1)
app.add(dm)

app.run(5000)
