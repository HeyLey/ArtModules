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
#include "TLorentzVector.h"
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCParticle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "Geometry/Geometry.h"
#include "EventGeneratorBase/evgenbase.h"
#include "SummaryData/RunData.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandBreitWigner.h"
#include "TTree.h"
#include "TH1.h"
#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <vector>
#include "TDatabasePDG.h"
#include "RecoBase/Prong.h"
#include "RecoBase/RecoHit.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "SimulationBase/MCTrajectory.h"

namespace mytest {
    class MyTestAna;
}

class mytest::MyTestAna : public art::EDAnalyzer {

public:

explicit MyTestAna(fhicl::ParameterSet const & p);

void analyze(art::Event const & e) override;
void beginJob() override ;

private:

const long pdg1 = 310L;
const long pdg2 = 211L;
double E[100];
TLorentzVector IP1[100];
TLorentzVector K;
TLorentzVector XN1,XN2,XK1,XK2;
int npart=0;
double XBegin[100],YBegin[100],ZBegin[100],XEnd[100],YEnd[100],ZEnd[100] ;
double Px[100];
double Py[100];
double Pz[100];
double L[100];
double CosThMom[100],PhiMom[100];
double CosThTrack[100],PhiTrack[100];
int Mother[100];
long Pdg[100];
TTree* MyTestAna_Tree;
TH1D  *H;
double Minv;
unsigned int y;
//  long pdg;
};



mytest::MyTestAna::MyTestAna(fhicl::ParameterSet const & p) : EDAnalyzer(p) {
//   pdg = p.get<long>("pdg");
}



void mytest::MyTestAna::beginJob() {
    npart = 0;
    art::ServiceHandle <art::TFileService> tfs;
    MyTestAna_Tree = tfs->make<TTree>("MyTestAna_Tree", "Information from Generator");
    MyTestAna_Tree->Branch("npart", &npart, "npart/I");
    MyTestAna_Tree->Branch("Minv", &Minv, "Minv/D");

    MyTestAna_Tree->Branch("Px", Px, "Px[100]/D");
    MyTestAna_Tree->Branch("Py", Py, "Py[100]/D");
    MyTestAna_Tree->Branch("Pz", Pz, "Pz[100]/D");
    MyTestAna_Tree->Branch("XBegin", XBegin, "XBegin[100]/D");
    MyTestAna_Tree->Branch("YBegin", YBegin, "YBegin[100]/D");
    MyTestAna_Tree->Branch("ZBegin", ZBegin, "ZBegin[100]/D");
    MyTestAna_Tree->Branch("XEnd", XEnd, "XEnd[100]/D");
    MyTestAna_Tree->Branch("YEnd", YEnd, "YEnd[100]/D");
    MyTestAna_Tree->Branch("ZEnd", ZEnd, "ZEnd[100]/D");
    MyTestAna_Tree->Branch("L", L, "L[100]/D");
    MyTestAna_Tree->Branch("CosThMom", CosThMom, "CosThMom[100]/D");
    MyTestAna_Tree->Branch("PhiMom", PhiMom, "PhiMom[100]/D");
    MyTestAna_Tree->Branch("CosThTrack", CosThTrack, "CosThTrack[100]/D");
    MyTestAna_Tree->Branch("Pdg", Pdg, "Pdg[100]/L");
    MyTestAna_Tree->Branch("Mother", Mother, "Mother[100]/I");
    MyTestAna_Tree->Branch("PhiTrack", PhiTrack, "PhiTrack[100]/D");
    H = tfs->make<TH1D>("Minv", ";Invariant Mass ", 100, 0, 5);
}

void mytest::MyTestAna::analyze(art::Event const & e) {
    std::cout<<npart<<std::endl;

    auto particles = e.getValidHandle<std::vector<sim::Particle>>("geantgen:");

    for(const auto& p: *particles) {
        y=p.NumberTrajectoryPoints();
        XBegin[npart]=p.Position().X();
        YBegin[npart]=p.Position().Y();
        ZBegin[npart]=p.Position().Z();
        XEnd[npart]=p.Position(y-1).X();
        YEnd[npart]=p.Position(y-1).Y();
        ZEnd[npart]=p.Position(y-1).Z();
        Pdg[npart]=p.PdgCode();
        Mother[npart]=p.Mother();
        E[npart]=p.Momentum().E();
        Px[npart]=p.Momentum().Px();
        Py[npart]=p.Momentum().Py();
        Pz[npart]=p.Momentum().Pz();
        L[npart]=sqrt(TMath::Power((XEnd[npart]-XBegin[npart]),2)+TMath::Power((YEnd[npart]-YBegin[npart]),2)+TMath::Power((ZEnd[npart]-ZBegin[npart]),2));
        CosThMom[npart]=Pz[npart]/sqrt(Px[npart]*Px[npart]+Py[npart]*Py[npart]+Pz[npart]*Pz[npart]);
        PhiMom[npart]=TMath::ATan2(Px[npart],Py[npart]);
        CosThTrack[npart]=(ZEnd[npart]-ZBegin[npart])/L[npart];
        PhiTrack[npart]=TMath::ATan2(XEnd[npart]-XBegin[npart],YEnd[npart]-YBegin[npart]);

        if ( p.PdgCode() == pdg1 && p.Mother() == 0) {
            IP1[0]=p.Momentum();
        }
        if(p.PdgCode() == pdg2 && p.Mother() == 0) {
            IP1[1]=p.Momentum();
        }

        npart++;
    }

    K = IP1[0] + IP1[1];
    H -> Fill(K.M());
    Minv = K.M();
    std::cout << npart << std::endl;
    MyTestAna_Tree -> Fill();
    std::cout << npart << std::endl;
    npart = 0;

}

namespace mytest {
    DEFINE_ART_MODULE(mytest::MyTestAna);
}