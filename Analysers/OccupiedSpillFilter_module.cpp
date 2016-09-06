#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// NOvA includes
#include "RawData/RawDigit.h"
#include "SimulationBase/MCTruth.h"
#include "Simulation/FLSHitList.h"

namespace evgen
{
    class OccupiedSpillFilter: public art::EDFilter
    {
    public:
        explicit OccupiedSpillFilter(const fhicl::ParameterSet& pset);
        ~OccupiedSpillFilter();

        bool filter(art::Event& evt);

        void reconfigure(const fhicl::ParameterSet& pset);

    private:
        std::string fGenLabel;    ///< Default "generator"
        std::string fG4GenLabel;  ///< Default "geantgen"
        std::string fDAQLabel;    ///< Default "daq"

        // Switches to configure checking of different objects
        bool fCheckG4Gen;
        bool fCheckDAQ;

    };
}


////////////////////////////////////////////////////////////////////////

namespace evgen
{
    //......................................................................
    OccupiedSpillFilter::OccupiedSpillFilter(const fhicl::ParameterSet& pset) :
            fCheckG4Gen(0),
            fCheckDAQ  (0)
    {
        reconfigure(pset);
    }

    //......................................................................
    OccupiedSpillFilter::~OccupiedSpillFilter()
    {
    }

    //......................................................................
    void OccupiedSpillFilter::reconfigure(const fhicl::ParameterSet& pset)
    {
        fGenLabel    = pset.get< std::string >("GenLabel"  );
        fG4GenLabel  = pset.get< std::string >("G4GenLabel");
        fDAQLabel    = pset.get< std::string >("DAQLabel"  );

        fCheckG4Gen  = pset.get< bool        >("CheckG4Gen");
        fCheckDAQ    = pset.get< bool        >("CheckDAQ"  );
    }
    //......................................................................
    bool OccupiedSpillFilter::filter(art::Event& evt)
    {

        // get the MCTruth list
        art::Handle<std::vector<simb::MCTruth> > mcTruths;
        evt.getByLabel(fGenLabel, mcTruths);
        // Filter if no truth info
        if(mcTruths.failedToGet()) {
            mf::LogError("OccupiedSpillFilter") << "No MCTruth object. \n";
            return false;
        }
        // Filter if truth info is empty
        if(mcTruths->empty()) {
            mf::LogError("OccupiedSpillFilter") << "Empty MCTruth object. \n";
            return false;
        }

        const simb::MCTruth mcTruth = (*mcTruths)[0];


        if(fCheckG4Gen){
            // get the FLSHit list
            art::Handle< std::vector<sim::FLSHitList> > FLScol;
            evt.getByLabel(fG4GenLabel, FLScol);
            // Filter if no FLShit info
            if (FLScol.failedToGet()) {
                mf::LogError("OccupiedSpillFilter") << "No FLSHit object. \n";
                return false;
            }
            // Filter if FLShit info is empty
            if(FLScol->empty()) {
                mf::LogError("OccupiedSpillFilter") << "Empty FLSHit object. \n";
                return false;
            }
        }

        if(fCheckDAQ){
            // get the RawDigit list
            art::Handle< std::vector<rawdata::RawDigit> > digitcol;
            evt.getByLabel(fDAQLabel, digitcol);
            // Filter if no rawdata info
            if(digitcol.failedToGet()) {
                mf::LogError("OccupiedSpillFilter") << "No RawDigit object. \n";
                return false;
            }
        }

        // Pass if there are Truth particles of course
        // after all subsequent checks
        return (mcTruth.NParticles() != 0);
    }
} // namespace

////////////////////////////////////////////////////////////////////////
namespace evgen
{
    DEFINE_ART_MODULE(OccupiedSpillFilter);
}