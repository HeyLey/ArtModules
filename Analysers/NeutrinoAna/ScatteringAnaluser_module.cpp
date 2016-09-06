#include "SimulationBase/GTruth.h"
#include "art/Framework/Core/FindOne.h"
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
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCNeutrino.h"
#include "TTree.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include <vector>
#include <math.h>
#include "art/Framework/Core/FindMany.h"
#include <string.h>

namespace leyla {
    class ScatteringAnalyzer;
}

class leyla::ScatteringAnalyzer : public art::EDAnalyzer {

public:

    explicit ScatteringAnalyzer(fhicl::ParameterSet const & p);

    ScatteringAnalyzer(ScatteringAnalyzer const &) = delete;
    ScatteringAnalyzer(ScatteringAnalyzer &&) = delete;
    ScatteringAnalyzer & operator = (ScatteringAnalyzer const &) = delete;
    ScatteringAnalyzer & operator = (ScatteringAnalyzer &&) = delete;

    void analyze(art::Event const & e) override;
    void beginJob();

private:

    TTree* ScatteringTree;
    std::vector<sim::Particle> pions;
    std::vector<sim::Particle> muons;
    std::vector<sim::Particle> neus;
    int n_muons;
    int n_pions;
    int n_neus;
    double sum_e_m;
    double sum_e_p;
    double sum_e_n;
    std::vector<double> e_muons;
    std::vector<double> e_pions;
    std::vector<double> px_pions;
    std::vector<double> py_pions;
    std::vector<double> pz_pions;
    std::vector<double> l_pions;
    std::vector<double> xf_pions;
std::vector<double> z_pions;
std::vector<double> e_neus;
int N_pdg;
std::uint32_t event_id;
double neutrino_energy;
double momentum4;
double energy_transfer;
double recoil_mass;
double cm_energy2;
double x;
double y;
int inter_type;
int xC;
double xSec;
double e_lep;
double e_jet;
double theta;
};

void leyla::ScatteringAnalyzer::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    ScatteringTree = tfs->make<TTree>("ScatteringTree", "Info_from_event");
    ScatteringTree->Branch("event_id", (& event_id), "event_id/I");
    ScatteringTree->Branch("N_pdg", (& N_pdg), "N_pdg/I");
    ScatteringTree->Branch("inter_type", (& inter_type), "inter_type/I");
    ScatteringTree->Branch("xC", (& xC), "xC/I");
    ScatteringTree->Branch("neutrino_energy",(& neutrino_energy), "neutrino_energy/D");
    ScatteringTree->Branch("n_pions", (& n_pions ), "n_pions/I");
ScatteringTree->Branch("n_muons", (& n_muons), "n_muons/I");
ScatteringTree->Branch("e_pions", (& e_pions));
ScatteringTree->Branch("px_pions", (& px_pions));
ScatteringTree->Branch("py_pions", (& py_pions));
ScatteringTree->Branch("pz_pions", (& pz_pions));
ScatteringTree->Branch("l_pions", (& l_pions));
ScatteringTree->Branch("xf_pions", (& xf_pions));
ScatteringTree->Branch("z_pions", (& z_pions));
ScatteringTree->Branch("e_muons", (& e_muons), "e_muons");
ScatteringTree->Branch("n_neus", (& n_neus), "n_neus");
ScatteringTree->Branch("e_neus", (& e_neus), "e_neus");
ScatteringTree->Branch("E_pi", (& sum_e_p), "E_pi/D");
ScatteringTree->Branch("E_mu", (& sum_e_m), "E_mu/D");
ScatteringTree->Branch("E_neu", (& sum_e_n), "E_neu/D");
ScatteringTree->Branch("E_lep", (& e_lep), "E_lep/D");
ScatteringTree->Branch("xB", (& x), "xB/D");
ScatteringTree->Branch("yB", (& y), "yB/D");
ScatteringTree->Branch("QSqr", (& momentum4), "QSqr/D");
ScatteringTree->Branch("WSqr", (& recoil_mass), "WSqr/D");
ScatteringTree->Branch("xSec", (& xSec), "xSec/D");
ScatteringTree->Branch("S",(& cm_energy2), "S/D");
ScatteringTree->Branch("energy_transfer", (& energy_transfer), "energy_transfer/D");
ScatteringTree->Branch("e_jet", (& e_jet), "e_jet/D");
ScatteringTree->Branch("theta", (& theta), "theta/D");
}

leyla::ScatteringAnalyzer::ScatteringAnalyzer(fhicl::ParameterSet const & p) : EDAnalyzer(p) {
//  pdg = p.get<int>("pdg");
}

void leyla::ScatteringAnalyzer::analyze(art::Event const & e) {

    event_id = e.event();
    std::cout <<"******************Event number:****************** "<<event_id<<std::endl;

    auto tracks = e.getValidHandle<std::vector<simb::MCTruth>>("generator");
    std::vector<simb::MCTruth> const &tracks_vector = *tracks;

    std::string trackCreatorModule = "geantgen";
    art::FindMany<sim::Particle> fmh(tracks, e, trackCreatorModule);
    auto gtr = e.getValidHandle<std::vector<simb::GTruth>>("generator");
std::vector<simb::GTruth> const &gtr_vector = *gtr;

for (size_t t = 0; t < tracks_vector.size(); ++t) {
pions.clear();
muons.clear();
neus.clear();
e_pions.clear();
e_muons.clear();
e_neus.clear();
sum_e_m = 0;
sum_e_p = 0;
sum_e_n = 0;

px_pions.clear();
py_pions.clear();
pz_pions.clear();
l_pions.clear();
xf_pions.clear();
z_pions.clear();

std::cout << "MCTruth: " << (t + 1) << std::endl;
const std::vector<const sim::Particle*> hits = fmh.at(t);

simb::MCTruth mct = tracks_vector[t];

simb::GTruth g = gtr_vector[t];

simb::MCNeutrino n = mct.GetNeutrino();

xSec = g.fXsec;
xC = n.CCNC();
N_pdg = n.Nu().PdgCode();
inter_type = n.InteractionType() ;
x = n.X();
y = n.Y();
neutrino_energy = n.Nu().E(0);
e_lep = (1 - y) * neutrino_energy;
energy_transfer = neutrino_energy - e_lep;
momentum4 = n.QSqr();
recoil_mass = pow(n.W(), 2);
double m = momentum4 / (2 * x * energy_transfer);
cm_energy2 = pow(m, 2) + momentum4 / (x * y);
e_jet = energy_transfer + m;
theta = n.Theta();

std::cout << "Act_type: " << inter_type << std::endl;
std::cout << "xC: " << xC << std::endl;
std::cout << "Neu_pdg: " << N_pdg << std::endl;
std::cout << "Neu_E: " << neutrino_energy << std::endl;
std::cout << "xSec: " << xSec << std::endl;

for (const sim::Particle* p: hits) {
if (p->Process() == "Primary") {

if (p->PdgCode() == abs(211)) {
pions.push_back(sim::Particle(*p));
}
if (p->PdgCode() == abs(13)) {
muons.push_back(sim::Particle(*p));
}
if (p->PdgCode() == abs(14)) {
neus.push_back(sim::Particle(*p));
}
}
}

n_neus = neus.size();
n_pions = pions.size();
n_muons = muons.size();
std::cout << "n_pi: " << n_pions << std::endl;
std::cout << "n_mu: " << n_muons << std::endl;
std::cout << "n_ne: " << n_neus << std::endl;

for (auto& pi: pions) {
std::cout << "Pi_pdg: " << pi.PdgCode() << std::endl;
std::cout << "Pi_E0: " << pi.E(0) << std::endl;
sum_e_p += pi.E(0);
e_pions.push_back(pi.E(0));
px_pions.push_back(pi.Px(0));
py_pions.push_back(pi.Py(0));
pz_pions.push_back(pi.Pz(0));
double l = sqrt( pow(pi.EndX()-pi.Vx(0), 2) + pow(pi.EndY()-pi.Vy(0), 2) + pow(pi.EndZ()-pi.Vz(0), 2) );
l_pions.push_back(l);
//xf_pions.push_back( /*...*/);
z_pions.push_back(pi.E(0)/e_jet);
}

for (auto& mu: muons) {
std::cout << "Mu_pdg: " << mu.PdgCode() << std::endl;
std::cout << "Mu_E0: " << mu.E(0) << std::endl;
sum_e_m += mu.E(0);
e_muons.push_back(mu.E(0));
}

for (auto& ne: neus) {
std::cout << "Neu_pdg: " << ne.PdgCode() << std::endl;
std::cout << "Neu_E0: " << ne.E(0) << std::endl;
sum_e_n += ne.E(0);
e_neus.push_back(ne.E(0));
}
//double lepton_energy = sum_p_enrg + sum_m_enrg;

ScatteringTree -> Fill();
std::cout << "" << std::endl;
}
std::cout << "" << std::endl;
}


DEFINE_ART_MODULE(leyla::ScatteringAnalyzer)