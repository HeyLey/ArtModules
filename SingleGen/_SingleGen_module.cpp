#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <vector>

// ROOT includes
#include "TDatabasePDG.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// NOvASoft includes
#include "Geometry/Geometry.h"
#include "SimulationBase/MCParticle.h"
#include "SimulationBase/MCTruth.h"
#include "EventGeneratorBase/evgenbase.h"
#include "SummaryData/RunData.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"


namespace evgen {
    class SingleParticle;

    /// A module to check the results from the Monte Carlo generator
    class SingleGen : public art::EDProducer {
    public:
        explicit SingleGen(fhicl::ParameterSet const& pset);
        virtual ~SingleGen();

        void produce(art::Event& evt);
        void reconfigure(const fhicl::ParameterSet& pset);
        void beginRun(art::Run& run);

    private:

        void Sample(simb::MCTruth &truth);

        /// Depending on configuration of fPMeaning
        /// calculate the momentum of the particle
        double getMomentum(const double& pkf, const double& m) const;

        static const int    kUNIF = 0;
        static const int    kGAUS = 1;
        std::vector<int>    fPDG;        ///< PDG code of particles to generate
        int                 fPMeaning;   ///< Meaning of P0, fSigmaP and fPDist.
        ///< By default (fP0Meaning=0), P0 and sigmaP is momentum
        ///< If fPMeaning=1, then P0 and sigmaP is kinetic energy in GeV
        ///< If fPMeaning=2, then P0and sigmaP is total energy in GeV
        std::vector<double> fP0;         ///< Central momentum (GeV/c) to generate
        std::vector<double> fSigmaP;     ///< Variation in momenta (GeV/c)
        std::vector<int>    fPDist;      ///< How to distribute momenta (0=uniform, 1=gaussian)
        std::vector<double> fX0;         ///< Central x position (cm)
        std::vector<double> fY0;         ///< Central y position (cm)
        std::vector<double> fZ0;         ///< Central z position (cm)
        std::vector<double> fT0;         ///< Central t position (ns)
        std::vector<double> fSigmaX;     ///< Variation in x position (cm)
        std::vector<double> fSigmaY;     ///< Variation in y position (cm)
        std::vector<double> fSigmaZ;     ///< Variation in z position (cm)
        std::vector<double> fSigmaT;     ///< Variation in t position (ns)
        std::vector<int>    fPosDist;    ///< How to distribute xyz (0=uniform, 1=gaussian)
        std::vector<double> fCosZ0;      ///< Cosine of central angle wrt z-axis
        std::vector<double> fPhiXY0;     ///< Central angle in the x-y plane (degrees)
        std::vector<double> fSigmaCosZ;  ///< Size of variation of cosz
        std::vector<double> fSigmaPhiXY; ///< Size of variation in phixy (degrees)
        std::vector<int>    fAngleDist;  ///< How to distribute angles (0=uniform, 1=gaussian)

    };
};

////////////////////////////////////////////////////////////////////////

namespace evgen{

    //____________________________________________________________________________
    SingleGen::SingleGen(fhicl::ParameterSet const& pset):
            fPMeaning(0)
    {
        this->reconfigure(pset);

        // get the random number seed, use a random default if not specified
        // in the configuration file.
        int seed = pset.get< unsigned int >("Seed", evgb::GetRandomNumberSeed());

        createEngine( seed );
        produces< std::vector<simb::MCTruth> >();
        produces< sumdata::RunData, art::InRun >();
    }

    //____________________________________________________________________________
    SingleGen::~SingleGen() { }

    //____________________________________________________________________________
    void SingleGen::reconfigure(const fhicl::ParameterSet& pset)
    {
        fPDG        = pset.get< std::vector<int>    >("PDG");

        // Default to momentum
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

    //____________________________________________________________________________
    void SingleGen::beginRun(art::Run& run)
    {

        // grab the geometry object to see what geometry we are using
        art::ServiceHandle<geo::Geometry> geo;

        std::unique_ptr<sumdata::RunData> runcol(new sumdata::RunData(geo->DetId(),
                                                                      geo->FileBaseName(),
                                                                      geo->ExtractGDML()));

        run.put(std::move(runcol));

        return;
    }

    //____________________________________________________________________________
    void SingleGen::produce(art::Event& evt)
    {

        ///unique_ptr allows ownership to be transferred to the art::Event after the ut statement
        std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);
        simb::MCTruth truth;
        truth.SetOrigin(simb::kSingleParticle);
        Sample(truth);

        //     std::cout << "put mctruth into the vector" << std::endl;
        truthcol->push_back(truth);

        //     std::cout << "add vector to the event " << truthcol->size() << std::endl;
        evt.put(std::move(truthcol));

        return;
    }

    //____________________________________________________________________________
    void SingleGen::Sample(simb::MCTruth& mct)
    {
        // get the random number generator service and make some CLHEP generators
        art::ServiceHandle<art::RandomNumberGenerator> rng;
        CLHEP::HepRandomEngine &engine = rng->getEngine();
        CLHEP::RandFlat   flat(engine);
        CLHEP::RandGaussQ gauss(engine);

        //
        // every event will have one of each particle species in the fPDG array
        //
        for (unsigned int i=0; i<fPDG.size(); ++i) {
            // std::cout << "picking a " << fPDG[i] << std::endl
            //    << fP0[i] << " " << fSigmaP[i] << " (" << fX0[i]
            //    << "/" << fSigmaX[i] << ", " << fY0[i]
            //    << "/" << fSigmaY[i] << ", " << fZ0[i]
            //    << "/" << fSigmaZ[i] << ") " << fTheta0XZ[i]
            //    << "/" << fSigmaThetaXZ[i] << " " << fTheta0YZ[i]
            //    << "/" << fSigmaThetaYZ[i] << std::endl;
            // std::cout << fPDist[i] << " " << fPosDist[i] << " " << fAngleDist[i] << std::endl;

            // Momentum, kinetic energy or full energy
            double pkf = 0.0;
            if (fPDist[i] == kGAUS) {
                pkf = gauss.fire(fP0[i], fSigmaP[i]);
            }
            else {
                pkf = flat.fire(fP0[i]-fSigmaP[i], fP0[i]+fSigmaP[i]);
            }

            // std::cout << "get the mass" << std::endl;
            const TDatabasePDG* pdgt = TDatabasePDG::Instance();
            const TParticlePDG* pdgp = pdgt->GetParticle(fPDG[i]);
            // Mass in GeV
            double    m = 0.0;
            if (pdgp) m = pdgp->Mass();

            // Momentum in GeV/c
            const double p = getMomentum(pkf, m);

            //std::cout<<"Momentum = "<<p<<"\n";

            // Choose position
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
            // std::cout << "set the position" << std::endl;
            const TLorentzVector pos(x[0], x[1], x[2], x[3]);

            // Choose angles

            double       cosz, phi;
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
            //
            // set track id to -i as these are all primary particles and have
            // id <= 0
            //
            const int    trackid = -1*(i+1);

            const std::string primary("primary");
            simb::MCParticle part(trackid, fPDG[i], primary);

            const TLorentzVector pvec(cosphi*sinz*p,
                                      sinphi*sinz*p,
                                      cosz*p,
                                      sqrt(p*p+m*m));

            part.AddTrajectoryPoint(pos, pvec);

            mct.Add(part);
        } //end loop over particles
    }

//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
    DEFINE_ART_MODULE(SingleGen);

} //end namespace




