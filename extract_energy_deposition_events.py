# Description: This script extracts the energy deposition events from a PTRAC file.
import mcnptools


# define a class to store the energy deposition event
class EnergyDepositionEvent:
    def __init__(self, NPS, ZAID, X, Y, Z, ReactionType, DepositedEnergy, TimeStamp, Weight):
        self.nps = NPS
        self.src_nuclide_ZAID = ZAID
        self.x = X
        self.y = Y
        self.z = Z
        self.reaction_type = ReactionType
        self.deposited_energy = DepositedEnergy
        self.time_stamp = TimeStamp
        self.weight = Weight

    def __str__(self):
        return f"{self.nps} {self.src_nuclide_ZAID} {self.x} {self.y} {self.z} {self.reaction_type} {self.deposited_energy} {self.time_stamp} {self.weight}"


# save the energy deposition events to a file
def save_to_file(photon_energy_deposition_events, file_name):
    with open(file_name, 'w') as file:
        # write the header
        file.write("NPS ZAID X Y Z ReactionType DepositedEnergy TimeStamp Weight\n")
        for ede in photon_energy_deposition_events:
            file.write(str(ede) + "\n")


if __name__ == "__main__":
    # explicitly open the file as an ascii ptrac
    p = mcnptools.Ptrac("ptrac", mcnptools.Ptrac.ASC_PTRAC)

    DETECTOR_CELL_NUMBER = 1
    photon_energy_deposition_events = []
    tmp = []
    prev_energy = 0.0
    ZAID = None
    hists = p.ReadHistories(10000)
    while len(hists) > 0:
        # loop over all histories
        for h in hists:
            # loop over all events in the history
            for e in range(h.GetNumEvents()):
                event = h.GetEvent(e)
                curr_energy = event.Get(mcnptools.Ptrac.ENERGY)
                if event.Type() == mcnptools.Ptrac.BNK: # bank event, new photon
                    ZAID = int(event.Get(mcnptools.Ptrac.ZAID))
                    tmp.clear()
                elif event.Type() == mcnptools.Ptrac.COL: # collision event
                    if event.Get(mcnptools.Ptrac.CELL) == DETECTOR_CELL_NUMBER: # interaction in detector
                        # special case: if reaction type is -4 (pair production), set the photon energy to 0
                        if int(event.Get(mcnptools.Ptrac.RXN)) == -4:  # pair production
                            curr_energy = 0.0
                        ede = EnergyDepositionEvent(h.GetNPS().NPS(),
                                                    ZAID,
                                                    event.Get(mcnptools.Ptrac.X),
                                                    event.Get(mcnptools.Ptrac.Y),
                                                    event.Get(mcnptools.Ptrac.Z),
                                                    int(event.Get(mcnptools.Ptrac.RXN)),
                                                    prev_energy - curr_energy,
                                                    event.Get(mcnptools.Ptrac.TIME),
                                                    event.Get(mcnptools.Ptrac.WEIGHT))
                        tmp.append(ede)
                elif event.Type() == mcnptools.Ptrac.TER: # termination event
                    if len(tmp) == 1:
                        photon_energy_deposition_events.append(tmp[0])
                    elif len(tmp) > 1:
                        # there are multiple energy deposition events for the same photon
                        # merge the energy deposition events
                        total_deposited_energy = sum(ed.deposited_energy for ed in tmp)
                        # compute energy-weighted average position
                        x = sum(ed.x * ed.deposited_energy for ed in tmp) / total_deposited_energy
                        y = sum(ed.y * ed.deposited_energy for ed in tmp) / total_deposited_energy
                        z = sum(ed.z * ed.deposited_energy for ed in tmp) / total_deposited_energy
                        ede = EnergyDepositionEvent(tmp[0].nps,
                                                    tmp[0].src_nuclide_ZAID,
                                                    x, y, z,
                                                    tmp[0].reaction_type,
                                                    total_deposited_energy,
                                                    tmp[0].time_stamp,
                                                    tmp[0].weight)
                        photon_energy_deposition_events.append(ede)
                else:
                    raise RuntimeError("Unknown event type: " + str(event.Type()))
                prev_energy = curr_energy
        hists = p.ReadHistories(10000)

    # sort by time stamp
    # photon_energy_deposition_events.sort(key=lambda ede: ede.time_stamp)
    save_to_file(photon_energy_deposition_events, "photon_energy_deposition_events.txt")
