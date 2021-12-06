# 导入 pygamd 包
import pygamd

mst=pygamd.snapshot.read("example.mst")

# nl=pygamd.plists.clist.clist(info=mst, rcut=1.5)
gg=pygamd.chare.particle_set(info=mst, group="all")
print(gg.member)
