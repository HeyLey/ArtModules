
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
    class MuonAnalyzer;
}

class leyla::MuonAnalyzer : public art::EDAnalyzer {

public:

    explicit MuonAnalyzer(fhicl::ParameterSet const & p);
    // The destructor generated by the compiler is fine for classes
    // without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    MuonAnalyzer(MuonAnalyzer const &) = delete;
    MuonAnalyzer(MuonAnalyzer &&) = delete;
    MuonAnalyzer & operator = (MuonAnalyzer const &) = delete;
    MuonAnalyzer & operator = (MuonAnalyzer &&) = delete;

    // Required functions.
    void analyze(art::Event const & e) override;
    void beginJob();

private:

    int pdg;
    TTree* MuonAnaTree;
    TLorentzVector *start_pos;
    TLorentzVector *end_pos;
    TLorentzVector *start_mom;
    TLorentzVector *end_mom;
    std::uint32_t event_id;
    double energy;
    double len;
};

void leyla::MuonAnalyzer::beginJob() {
    start_pos = new TLorentzVector();
    end_pos = new TLorentzVector();
    start_mom = new TLorentzVector();
    end_mom = new TLorentzVector();
    art::ServiceHandle<art::TFileService> tfs;
    MuonAnaTree = tfs->make<TTree>("MuonAnalyzerTree", "Info_from_event");
    MuonAnaTree->Branch("event_id", (& event_id), "event_id/I");
    MuonAnaTree->Branch("start_pos", (& start_pos));
    MuonAnaTree->Branch("end_pos", (& end_pos));
    MuonAnaTree->Branch("start_mom", (& start_mom));
    MuonAnaTree->Branch("end_mom", (& end_mom));
    MuonAnaTree->Branch("len", (& len), "len/D");
    MuonAnaTree->Branch("energy", (& energy), "energy/D");

}

leyla::MuonAnalyzer::MuonAnalyzer(fhicl::ParameterSet const & p) : EDAnalyzer(p) {
    pdg = p.get<int>("pdg");
    //std::cout<<"pdg="<<pdg<<std::endl;
}


void leyla::MuonAnalyzer::analyze(art::Event const & e) {
    auto particles = e.getValidHandle<std::vector<sim::Particle>>("geantgen:");
    // double summary_energy = 0;
    event_id = e.event();
    std::cout <<"Event number: "<<event_id<<std::endl;
    for (const auto& p: *particles) {
        pdg = p.PdgCode();
        //lifetime = p.EndT() - p.T(0);

        if (pdg == 13 && p.Mother() == 0) {
            *start_pos = p.Position(0);
            *end_pos = p.EndPosition();
            *start_mom = p.Momentum(0);
            len = sqrt(pow((p.EndX()-p.Vx(0)), 2) + pow((p.EndY()-p.Vy(0)), 2) + pow((p.EndZ()-p.Vz(0)), 2));
            *end_mom = p.EndMomentum();
            energy = p.E();
            //  summary_energy += energy;
            MuonAnaTree -> Fill();
            std::cout << "Energy: " << energy << " Length: " <<len<<"\n"<<std::endl;
        }


        // std::cout<<"PDG of daughter element with process 'nCapture': "<<pdg<<"\nStart energy: "<<p.E(0)<<"\n"<<std::endl;

        // std::cout<<"PDG of daughter element with process 'nCapture': "<<pdg<<"\nStart energy: "<<p.E(0)<<"\n"<<std::endl;
    }
    //std::cout<<"Summary energy of all photons: "<<summary_energy<<"\n"<<std::endl;
}

DEFINE_ART_MODULE(leyla::MuonAnalyzer)