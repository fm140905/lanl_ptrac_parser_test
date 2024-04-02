#include "mcnptools/McnpTools.hpp"
#include <iomanip>
#include <iostream>
#include <vector>


class EnergyDepositionEvent
{
private:
    /* data */
public:
    int src_nuclide_ZAID;
    int nps;
    double x, y, z;
    int reaction_type;
    double deposited_energy;
    double time_stamp;
    double weight;
    EnergyDepositionEvent(const int NPS, const int ZAID,
                            const double X, const double Y, const double Z, 
                            const int ReactionType,
                            const double DepositedEnergy, 
                            const double TimeStamp,
                            const double Weight)
     : nps(NPS), src_nuclide_ZAID(ZAID),
       x(X), y(Y), z(Z), 
       reaction_type(ReactionType), 
       deposited_energy(DepositedEnergy),
       time_stamp(TimeStamp),
       weight(Weight) {}
    
    friend std::ostream& operator<< (std::ostream& stream, const EnergyDepositionEvent& ede) 
    {
        stream << ede.nps << " " 
               << ede.src_nuclide_ZAID << " " 
               << ede.x << " " 
               << ede.y << " " 
               << ede.z << " " 
               << ede.reaction_type << " " 
               << ede.deposited_energy << " " 
               << ede.time_stamp << " " 
               << ede.weight;
        return stream;
    }
};

int save_to_file(const std::vector<EnergyDepositionEvent>& photon_energy_deposition_events, const std::string& file_name)
{
    std::ofstream file(file_name);
    if (file.is_open())
    {
        // write the header
        file << "NPS ZAID X Y Z ReactionType DepositedEnergy TimeStamp Weight\n";
        for (auto ede : photon_energy_deposition_events)
        {
            file << ede << "\n";
        }
        file.close();
        return 0;
    }
    else
    {
        return 1;
    }
}


int main() {

    std::cout << std::scientific << std::setprecision(5);

    // explicitly open the file as a ascii ptrac
    mcnptools::Ptrac p("/media/ming/Elements/lanl_ptrac_parser_test/cpp/data/ptrac", mcnptools::Ptrac::ASC_PTRAC);

    // initialize counter
    unsigned int cnt = 0;

    // // read histories in batches of 10000
    // std::vector<mcnptools::PtracHistory> hists = p.ReadHistories(100);
    // std::cout << "frist nps: " << hists.at(0).GetNPS().NPS() << "\n";
    // std::cout << "Number of histories: " << hists.size() << "\n";
    // while (hists.size() > 0) {
    //     // loop over all histories
    //     for (unsigned int h = 0; h < hists.size(); h++) {
    //         // loop over all events in the history
    //         for (unsigned int e = 0; e < hists.at(h).GetNumEvents(); e++) {
    //             mcnptools::PtracEvent event = hists.at(h).GetEvent(e);
    //             if (event.Type() == mcnptools::Ptrac::BNK &&
    //                 event.BankType() == mcnptools::Ptrac::BNK_N_XG &&
    //                 event.Get(mcnptools::Ptrac::ZAID) == 6000) {
    //                 cnt += 1;
    //                 std::cout << "NPS: " << hists.at(h).GetNPS().NPS() << "\n";
    //                 std::cout << std::setw(13) << cnt << std::setw(13)
    //                           << event.Get(mcnptools::Ptrac::ENERGY) << std::endl;
    //             }
    //         }
    //     }
    //     hists = p.ReadHistories(10000);
    // }

    const int DETECTOR_CELL_NUMBER = 1;
    std::vector<mcnptools::PtracHistory> hists = p.ReadHistories(100);
    std::vector<EnergyDepositionEvent> photon_energy_deposition_events;
    double prev_energy, curr_energy;
    int ZAID;
    std::vector<EnergyDepositionEvent> tmp={};
    while (hists.size() > 0) 
    {
        // loop over all histories
        for (unsigned int h = 0; h < hists.size(); h++) 
        {
            // loop over all events in the history
            for (unsigned int e = 0; e < hists.at(h).GetNumEvents(); e++) {
                mcnptools::PtracEvent event = hists.at(h).GetEvent(e);
                curr_energy = event.Get(mcnptools::Ptrac::ENERGY);
                if (event.Type() == mcnptools::Ptrac::BNK) // bank
                {
                    // new photon
                    ZAID = event.Get(mcnptools::Ptrac::ZAID);
                    tmp.clear();
                }
                else if (event.Type() == mcnptools::Ptrac::COL) // collion
                {
                    if (event.Get(mcnptools::Ptrac::CELL) == DETECTOR_CELL_NUMBER) // interaction in detector
                    {
                        // special case: if reaction type is -4 (pair production)
                        // set current energy to 0
                        if (event.Get(mcnptools::Ptrac::RXN) == -4)
                        {
                            curr_energy = 0.0;
                        }
                        EnergyDepositionEvent ede(hists.at(h).GetNPS().NPS(),
                                                  ZAID,
                                                  event.Get(mcnptools::Ptrac::X),
                                                  event.Get(mcnptools::Ptrac::Y),
                                                  event.Get(mcnptools::Ptrac::Z),
                                                  event.Get(mcnptools::Ptrac::RXN),
                                                  prev_energy-curr_energy,
                                                  event.Get(mcnptools::Ptrac::TIME),
                                                  event.Get(mcnptools::Ptrac::WEIGHT));
                        tmp.push_back(ede);
                    }
                }
                else if (event.Type() == mcnptools::Ptrac::TER)// termination event
                {
                    if (tmp.size() == 1)
                    {
                        photon_energy_deposition_events.push_back(tmp.at(0));
                    }
                    else if (tmp.size() > 1)
                    {
                        // merge the energy deposition events
                        double total_deposited_energy = 0.0;
                        for (auto ede : tmp)
                        {
                            total_deposited_energy += ede.deposited_energy;
                        }
                        // compute energy-weighted average position
                        double x = 0.0, y = 0.0, z = 0.0;
                        for (auto ede : tmp)
                        {
                            x += ede.x * ede.deposited_energy;
                            y += ede.y * ede.deposited_energy;
                            z += ede.z * ede.deposited_energy;
                        }
                        x /= total_deposited_energy;
                        y /= total_deposited_energy;
                        z /= total_deposited_energy;
                        EnergyDepositionEvent ede(tmp.at(0).nps, 
                                                  tmp.at(0).src_nuclide_ZAID,
                                                  x, y, z, 
                                                  tmp.at(0).reaction_type, 
                                                  total_deposited_energy, 
                                                  tmp.at(0).time_stamp, 
                                                  tmp.at(0).weight);
                        photon_energy_deposition_events.push_back(ede);
                    }
                }
                else
                {
                    // throw std::runtime_error("Unknown event type");
                    throw std::runtime_error("Unknown event type: " + std::to_string(event.Type()));
                }
                prev_energy = curr_energy;
            }
        }
        hists = p.ReadHistories(10000);
    }
    // // sort the energy deposition events by time stamp
    // std::sort(photon_energy_deposition_events.begin(), photon_energy_deposition_events.end(), 
    //           [](const EnergyDepositionEvent& a, const EnergyDepositionEvent& b) -> bool
    //             {
    //                 return a.time_stamp < b.time_stamp;
    //             }
    // );
    save_to_file(photon_energy_deposition_events, "photon_energy_deposition_events.txt");

    return 0;
}
