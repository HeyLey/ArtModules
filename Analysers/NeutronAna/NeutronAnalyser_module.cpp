
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Simulation/Particle.h"
#include "TTree.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include <vector>
#include <math.h>

namespace leyla {
    class NeutronAnalyzer;
}

class leyla::NeutronAnalyzer : public art::EDAnalyzer {

public:

    explicit NeutronAnalyzer(fhicl::ParameterSet const & p);

    NeutronAnalyzer(NeutronAnalyzer const &) = delete;
    NeutronAnalyzer(NeutronAnalyzer &&) = delete;
    NeutronAnalyzer & operator = (NeutronAnalyzer const &) = delete;
    NeutronAnalyzer & operator = (NeutronAnalyzer &&) = delete;

    void analyze(art::Event const & e) override;
    void beginJob();

private:

    int pdg;
    TTree* NeutronAnaTree;
    TLorentzVector *start_pos;
    TLorentzVector *end_pos;
    // TLorentzVector *gamma_pos;
    // TLorentzVector *gamma_start_mom;
    // TLorentzVector *gamma_end_mom;
    std::uint32_t event_id;
    // double energy;
    // double sum_energy;
    int n_ne;
    double time;
    double len;
};

void leyla::NeutronAnalyzer::beginJob() {
    start_pos = new TLorentzVector();
    end_pos = new TLorentzVector();
    // gamma_pos = new TLorentzVector();
    // gamma_start_mom = new TLorentzVector();
    // gamma_end_mom = new TLorentzVector();
    art::ServiceHandle<art::TFileService> tfs;
    NeutronAnaTree = tfs->make<TTree>("NeutronAnaTree", "Info_from_event");
    NeutronAnaTree->Branch("event_id", (& event_id), "event_id/I");
    NeutronAnaTree->Branch("len", (& len), ("len/D"));
    NeutronAnaTree->Branch("time", (& time), ("time/D"));
    NeutronAnaTree->Branch("start_pos", (& start_pos));
    NeutronAnaTree->Branch("end_pos", (& end_pos));
    // NeutronAnaTree->Branch("gamma_pos", (& gamma_pos));
    // NeutronAnaTree->Branch("gamma_start_mom", (& gamma_start_mom));
    // NeutronAnaTree->Branch("gamma_end_mom", (& gamma_end_mom));
    NeutronAnaTree->Branch("n_ne", (& n_ne), "n_ne/I");
    // NeutronAnaTree->Branch("energy", (& energy), "energy/D");
    // NeutronAnaTree->Branch("sum_energy", (& sum_energy), "sum_energy/D");
}

leyla::NeutronAnalyzer::NeutronAnalyzer(fhicl::ParameterSet const & p) : EDAnalyzer(p) {
    pdg = p.get<int>("pdg");

}


void leyla::NeutronAnalyzer::analyze(art::Event const & e) {
    auto particles = e.getValidHandle<std::vector<sim::Particle>>("geantgen:");
    // sum_energy = 0;
    n_ne = 0;
    event_id = e.event();
    std::cout <<"Event number: "<<event_id<<std::endl;
    for (const auto& p: *particles) {
        pdg = p.PdgCode();
        if (pdg == 2112 && p.Process() == "Primary") {

            *start_pos = p.Position(0);
            *end_pos = p.EndPosition();
            len = sqrt( pow((p.EndX() - p.Vx(0)), 2) + pow((p.EndY() - p.Vy(0)), 2) + pow((p.EndZ() - p.Vz(0)), 2) );
            time = p.EndT() - p.T(0);
            n_ne += 1;
        }

        /* if (pdg == 22) {

            *gamma_pos = p.Position(0);
            *gamma_start_mom = p.Momentum(0);
            *gamma_end_mom = p.EndMomentum();
            energy = p.E();
            sum_energy += energy;
            n_ph += 1;
        //    NeutronAnaTree -> Fill();
          }*/
        NeutronAnaTree -> Fill();
    }
    // std::cout<<"Summary energy of all photons: "<<sum_energy<<"\n"<<std::endl;
}

namespace leyla {
      DEFINE_ART_MODULE(leyla::NeutronAnalyzer);
}