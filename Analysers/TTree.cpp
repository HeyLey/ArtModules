#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <vector>

#include "TDatabasePDG.h"
#include "TTree.h"


#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "Geometry/Geometry.h"
#include "SimulationBase/MCParticle.h"
#include "SimulationBase/MCTruth.h"
#include "EventGeneratorBase/evgenbase.h"
#include "SummaryData/RunData.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandBreitWigner.h"

namespace evgen {
    class SingleParticle;

    class SingleGen : public art::EDProducer {
    public:
        explicit SingleGen(fhicl::ParameterSet const& pset);
        virtual ~SingleGen();

        void produce(art::Event& evt);
        void reconfigure(const fhicl::ParameterSet& pset);
        void beginRun(art::Run& run);

    private:

        void Sample(simb::MCTruth &truth);

        double getMomentum(const double& pkf, const double& m) const;

        static const int    kUNIF = 0;
        static const int    kGAUS = 1;
        std::vector<int>    fPDG;
        int                 fPMeaning;
        std::vector<double> fP0;
        std::vector<double> fSigmaP;
        std::vector<int>    fPDist;

        std::vector<double> fX0;
        std::vector<double> fY0;
        std::vector<double> fZ0;
        std::vector<double> fT0;
        std::vector<double> fSigmaX;
        std::vector<double> fSigmaY;
        std::vector<double> fSigmaZ;
        std::vector<double> fSigmaT;
        std::vector<int>    fPosDist;
        std::vector<double> fCosZ0;

        std::vector<double> fPhiXY0;
        std::vector<double> fSigmaCosZ;
        std::vector<double> fSigmaPhiXY;
        std::vector<int>    fAngleDist;

        int npart;
        TLorentzVector Ptot;
        TTree* SingleGen_Tree;
        double E[10];
        double Px[10];
        double Py[10];
        double Pz[10];
        double M[10];
        int pdg[10];
        double Minv;

    };
};



namespace evgen{


    SingleGen::SingleGen(fhicl::ParameterSet const& pset):
            fPMeaning(0)
    {
        this->reconfigure(pset);


        int seed = pset.get< unsigned int >("Seed", evgb::GetRandomNumberSeed());

        createEngine( seed );

        produces< std::vector<simb::MCTruth> >();
        produces< sumdata::RunData, art::InRun >();
    }

    SingleGen::~SingleGen() { }


    void SingleGen::reconfigure(const fhicl::ParameterSet& pset)
    {
        fPDG        = pset.get< std::vector<int>    >("PDG");


        fPMeaning = pset.get< int >("PMeaning", 0);
        if(fPMeaning == 1) {
            std::cout<<"SingleGen: Using Kinetic energy for the meaning of P0, SigmaP and PDist\n";
        }
        else if(fPMeaning == 2) {
            std::cout<<"SingleGen: Using Total energy for the meaning of P0, SigmaP and PDist\n";
        }

        fP0         = pset.get< std::vector<double> >("P0");
        fSigmaP     = pset.get< std::vector<double> >("SigmaP");
        fPDist      = pset.get< std::vector<int>    >("PDist");
        fX0         = pset.get< std::vector<double> >("X0");
        fY0         = pset.get< std::vector<double> >("Y0");
        fZ0         = pset.get< std::vector<double> >("Z0");
        fT0         = pset.get< std::vector<double> >("T0");
        fSigmaX     = pset.get< std::vector<double> >("SigmaX");
        fSigmaY     = pset.get< std::vector<double> >("SigmaY");
        fSigmaZ     = pset.get< std::vector<double> >("SigmaZ");
        fSigmaT     = pset.get< std::vector<double> >("SigmaT");
        fPosDist    = pset.get< std::vector<int>    >("PosDist");
        fCosZ0      = pset.get< std::vector<double> >("CosZ0");
        fPhiXY0     = pset.get< std::vector<double> >("PhiXY0");
        fSigmaPhiXY = pset.get< std::vector<double> >("SigmaPhiXY");
        fSigmaCosZ  = pset.get< std::vector<double> >("SigmaCosZ");
        fAngleDist  = pset.get< std::vector<int>    >("AngleDist");
    }

    void SingleGen::beginRun(art::Run& run)
    {


        art::ServiceHandle<geo::Geometry> geo;

        std::unique_ptr<sumdata::RunData> runcol(new sumdata::RunData(geo->DetId(),
                                                                      geo->FileBaseName(),
                                                                      geo->ExtractGDML()));

        run.put(std::move(runcol));


        art::ServiceHandle<art::TFileService> tfs;
        SingleGen_Tree = tfs->make<TTree>("SingleGen_Tree","Information from Generator");
        SingleGen_Tree -> Branch("npart", &npart, "npart/I");
        SingleGen_Tree -> Branch("E",   E,   "E[npart]/D");
        SingleGen_Tree -> Branch("Px",  Px,  "Px[npart]/D");
        SingleGen_Tree -> Branch("Py",  Py,  "Py[npart]/D");
        SingleGen_Tree -> Branch("Pz",  Pz,  "Pz[npart]/D");
        SingleGen_Tree -> Branch("M",   M,   "M[npart]/D");
        SingleGen_Tree -> Branch("pdg", pdg, "pdg[npart]/I");
        SingleGen_Tree -> Branch("Minv", &Minv, "Minv/D");

        return;
    }


    void SingleGen::produce(art::Event& evt)
    {

        std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);

        simb::MCTruth truth;
        truth.SetOrigin(simb::kSingleParticle);
        Sample(truth);

        truthcol->push_back(truth);

        evt.put(std::move(truthcol));

        return;
    }

    void SingleGen::Sample(simb::MCTruth& mct)
    {

        art::ServiceHandle<art::RandomNumberGenerator> rng;
        CLHEP::HepRandomEngine &engine = rng->getEngine();
        CLHEP::RandFlat   flat(engine);
        CLHEP::RandGaussQ gauss(engine);
        CLHEP::RandBreitWigner bw(engine);


        npart=0;
        Ptot.SetPxPyPzE(0,0,0,0);
        for (unsigned int i=0; i<fPDG.size(); ++i) {
            double pkf = 0.0;
            if (fPDist[i] == kGAUS) {
                pkf = gauss.fire(fP0[i], fSigmaP[i]);
            }
            else {
                pkf = flat.fire(fP0[i]-fSigmaP[i], fP0[i]+fSigmaP[i]);
            }

            const TDatabasePDG* pdgt = TDatabasePDG::Instance();
            const TParticlePDG* pdgp = pdgt->GetParticle(fPDG[i]);

            double    m = 0.0;
            if (pdgp) m = pdgp->Mass();
            double   w=0.0;
            if(pdgp) w=pdgp->Width();

            const double p = getMomentum(pkf, m);



            double x[4];
            if (fPosDist[i] == kGAUS) {
                x[0] = gauss.fire(fX0[i], fSigmaX[i]);
                x[1] = gauss.fire(fY0[i], fSigmaY[i]);
                x[2] = gauss.fire(fZ0[i], fSigmaZ[i]);
                x[3] = gauss.fire(fT0[i], fSigmaT[i]);
            }
            else {
                x[0] = flat.fire(fX0[i]-fSigmaX[i], fX0[i]+fSigmaX[i]);
                x[1] = flat.fire(fY0[i]-fSigmaY[i], fY0[i]+fSigmaY[i]);
                x[2] = flat.fire(fZ0[i]-fSigmaZ[i], fZ0[i]+fSigmaZ[i]);
                x[3] = flat.fire(fT0[i]-fSigmaT[i], fT0[i]+fSigmaT[i]);
            }

            const TLorentzVector pos(x[0], x[1], x[2], x[3]);

            double ms=bw.fire(m,w);



            double cosz, phi;
            unsigned int itry;
            for (itry=0; itry<1000000; ++itry) {
                if (fAngleDist[i] == kGAUS) {
                    cosz = gauss.fire(fCosZ0[i],  fSigmaCosZ[i]);
                    phi  = gauss.fire(fPhiXY0[i], fSigmaPhiXY[i]);
                }
                else {
                    cosz = flat.fire(fCosZ0[i]-fSigmaCosZ[i],
                                     fCosZ0[i]+fSigmaCosZ[i]);

                    phi  = flat.fire(fPhiXY0[i]-fSigmaPhiXY[i],
                                     fPhiXY0[i]+fSigmaPhiXY[i]);
                }
                if (cosz>=-1.0 && cosz<=1.0) break;
            }
            if (cosz<-1.0 || cosz>1.0) {
                mf::LogError("SingleGen") << __FILE__ << ":" << __LINE__
                                          << " Failed to find an acceptable cos(theta_z)"
                                          << " after many tries.\n"
                                          << " Please adjust CosZ0 and SigmaCosZ0 in your"
                                          << " SingleGen.fcl file and rerun";
                abort();
            }

            const double sinz    = sqrt(1.0-cosz*cosz);
            const double sinphi  = sin(phi*M_PI/180.0);
            const double cosphi  = cos(phi*M_PI/180.0);


            M[npart] =m;

            const int    trackid = -1*(i+1);

            const std::string primary("primary");
            simb::MCParticle part(trackid, fPDG[i], primary);

            const TLorentzVector pvec(cosphi*sinz*p,
                                      sinphi*sinz*p,
                                      cosz*p,
                                      sqrt(p*p+ms*ms));

            part.AddTrajectoryPoint(pos, pvec);
            mct.Add(part);

            E[npart]   = pvec.E();
            Px[npart]  = pvec.Px();
            Py[npart]  = pvec.Py();
            Pz[npart]  = pvec.Pz();
            pdg[npart] = fPDG[i];

            npart++;
            Ptot += pvec;

        }

        Minv = Ptot.M();

        SingleGen_Tree -> Fill();

    }


    double SingleGen::getMomentum(const double& pkf, const double& m) const{

        if     (fPMeaning == 0) return pkf;

        double total_energy = 0.0;
        if     (fPMeaning == 1){
            total_energy = pkf + m;
        }
        else if(fPMeaning == 2){
            total_energy = pkf;
        }
        else{
            total_energy = sqrt(pkf*pkf + m*m);
        }

        return sqrt(total_energy*total_energy - m*m);
    }


    DEFINE_ART_MODULE(SingleGen);

}
