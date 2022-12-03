neighbor_listSR = galamost.NeighborList(all_info, rcut, rbuffer)
neighbor_listSR.addExclusionsFromBonds()
neighbor_listSR.addExclusionsFromAngles()
