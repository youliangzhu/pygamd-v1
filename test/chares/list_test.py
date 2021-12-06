# 导入 pygamd 包
import pygamd

mst=pygamd.snapshot.read("example.mst")

# nl=pygamd.plists.clist.clist(info=mst, rcut=1.5)
nl=pygamd.plists.nlist.nlist(info=mst, rcut=1.5, rbuff=0.1)
nl.calculate()
