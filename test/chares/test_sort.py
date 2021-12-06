import pygamd
	
mst = pygamd.snapshot.read("lj.mst")

fn = pygamd.tinker.sort(info=mst, period=1)
fn.data.calculate(0)

