
#include <cassert>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <unistd.h>

// ROOT includes
#include "TStopwatch.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


// NOvA includes
#include "EventGeneratorBase/GENIE/GENIEHelper.h"
#include "Geometry/Geometry.h"
#include "SimulationBase/GTruth.h"
#include "SimulationBase/MCFlux.h"
#include "SimulationBase/MCTruth.h"
#include "SummaryData/POTSum.h"
#include "SummaryData/SpillData.h"
#include "SummaryData/RunData.h"
#include "Utilities/AssociationUtil.h"

///Monte Carlo event generation
namespace evgen {

    /// A module to check the results from the Monte Carlo generator
    class GENIEGen : public art::EDProducer {

    public:

        explicit GENIEGen(fhicl::ParameterSet const &pset);
        virtual ~GENIEGen();

        void produce(art::Event& evt);
        void beginJob();
        void beginRun(art::Run &run);
        void endSubRun(art::SubRun &sr);

    private:

        evgb::GENIEHelper  *fGENIEHelp;        ///< GENIEHelper object
        int                 fPassEmptySpills;  ///< whether or not to allow for events with no interactions to be put into the event
        TStopwatch          fStopwatch;        ///< keep track of how long it takes to run the job
        int                 fSpillCounter;
        double              fPOTPerSpill;
        double              fEventsPerSpill;
        double              fTotalExposure;

        double              fTotalPOTLimit;
    };
};

namespace evgen {

    //___________________________________________________________________________
    GENIEGen::GENIEGen(fhicl::ParameterSet const& pset)
            : fGENIEHelp       (0)
            , fPassEmptySpills (pset.get< bool >("PassEmptySpills"))
            , fSpillCounter    (0)
            , fPOTPerSpill     (pset.get< double >("POTPerSpill",    5.0e13))
            , fEventsPerSpill  (pset.get< double >("EventsPerSpill", 0))
            , fTotalExposure   (0)
            , fTotalPOTLimit   (pset.get< double >("TotalPOTLimit"))
    {
        fStopwatch.Start();

        produces< std::vector<simb::MCTruth> >();
        produces< std::vector<simb::MCFlux>  >();
        produces< std::vector<simb::GTruth>  >();
        produces< sumdata::SpillData >();
        produces< sumdata::POTSum, art::InSubRun  >();
        produces< sumdata::RunData, art::InRun    >();
        // Associate every truth with the flux it came from
        produces< art::Assns<simb::MCTruth, simb::MCFlux> >();
        produces< art::Assns<simb::MCTruth, simb::GTruth> >();

        art::ServiceHandle<geo::Geometry> geo;
        fGENIEHelp = new evgb::GENIEHelper(pset,
                                           geo->ROOTGeoManager(),
                                           geo->ROOTFile(),
                                           geo->TotalMass(pset.get< std::string>("TopVolume").c_str()));
    }

    //___________________________________________________________________________
    GENIEGen::~GENIEGen()
    {
        fStopwatch.Stop();
        mf::LogInfo("GENIEGen") << "real time to produce file: "
                                << fStopwatch.RealTime();
        delete fGENIEHelp; // clean up, and let dtor do its thing
    }

    //___________________________________________________________________________
    void GENIEGen::beginJob()
    {
    }

    //___________________________________________________________________________
    void GENIEGen::beginRun(art::Run& run)
    {
        // grab the geometry object to see what geometry we are using
        art::ServiceHandle<geo::Geometry> geo;

        std::unique_ptr<sumdata::RunData>
                runcol(new sumdata::RunData(geo->DetId(),
                                            geo->FileBaseName(),
                                            geo->ExtractGDML()));

        run.put(std::move(runcol));

        // initialize the GENIEHelper here rather than in beginJob to
        // avoid problems with the Geometry reloading at a run boundary.
        // If we ever make more than one run in a single job we will have
        // to re-evaluate
        fGENIEHelp->Initialize();
        fTotalExposure = 0.0;

        return;
    } //___________________________________________________________________________
    void GENIEGen::endSubRun(art::SubRun &sr)
    {
        std::unique_ptr< sumdata::POTSum > p(new sumdata::POTSum);

        // p->totpot     = fGENIEHelp->TotalExposure();
        // p->totgoodpot = fGENIEHelp->TotalExposure();
        p->totpot     = fTotalExposure;
        p->totgoodpot = fTotalExposure;
        p->totspills  = fSpillCounter;
        p->goodspills = fSpillCounter;
        p->Print(std::cout);

        sr.put(std::move(p));
    }

    //___________________________________________________________________________
    void GENIEGen::produce(art::Event& evt)
    {
        // A temporary value is needed to store the spill exposure that GENIEHelper uses.  GENIEGen
        // needs to remember the number after GENIEHelper has reset it to zero for the purposes of
        // updating fTotalExposure.
        double SpillExpTemp = 0.0;

        std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);
        std::unique_ptr< std::vector<simb::MCFlux>  > fluxcol (new std::vector<simb::MCFlux >);
        std::unique_ptr< std::vector<simb::GTruth>  > gtruthcol (new std::vector<simb::GTruth >);
        std::unique_ptr< art::Assns<simb::MCTruth, simb::GTruth> > tgtassn(new art::Assns<simb::MCTruth, simb::GTruth>);
        std::unique_ptr< art::Assns<simb::MCTruth, simb::MCFlux> > assns(new art::Assns<simb::MCTruth, simb::MCFlux>);

        bool passedPOTLimit = false;

        // keep throwing spills until we get one with something in it
        while ( truthcol->size() < 1 ) {

            // POT limit. Maybe we passed the limit in a previous event, so now we
            // have to keep skipping until the end of the file. Or maybe we're in the
            // process of breaking out of the inner loop, and need to break again.
            if(fTotalPOTLimit > 0 &&
               fTotalExposure + fGENIEHelp->SpillExposure() >= fTotalPOTLimit){
                passedPOTLimit = true;
                break;
            }

            while ( ! fGENIEHelp->Stop() ) {

                simb::MCTruth truth;
                simb::MCFlux  flux;
                simb::GTruth  gTruth;

                // GENIEHelper returns a false in the sample method if
                // either no neutrino was generated, or the interaction
                // occurred beyond the detector's z extent - ie something we
                // would never see anyway.
                if ( fGENIEHelp->Sample(truth, flux, gTruth ) ) {

                    // POT limit
                    if(fTotalPOTLimit > 0 &&
                       fTotalExposure + fGENIEHelp->SpillExposure() > fTotalPOTLimit){
                        // Top the total exposure up to the right value
                        SpillExpTemp = fTotalPOTLimit - fTotalExposure;
                        break;
                    }

                    // When running in "POT per Spill" mode, this will prevent the last neutrino (that was produced after
                    // the POT limit had been passed) from being put into the event.
                    //
                    // This is only a temporary fix to correct the problem of GENIEHelper always making one more neutrino
                    // than it should.
                    if((fGENIEHelp->SpillExposure() <= fPOTPerSpill && fEventsPerSpill == 0) || fEventsPerSpill > 0) {
                        truthcol ->push_back(truth);
                        gtruthcol->push_back(gTruth);
                        fluxcol  ->push_back(flux);

                        util::CreateAssn(*this, evt, *truthcol, *fluxcol, *assns,
                                         fluxcol->size()-1, fluxcol->size());

                        util::CreateAssn(*this, evt, *truthcol, *gtruthcol, *tgtassn,
                                         gtruthcol->size()-1, gtruthcol->size());
                        SpillExpTemp = fGENIEHelp->SpillExposure();
                    }

                } // end if genie was able to make an event

            } // end event generation loop

            // check to see if we are to pass empty spills
            if ( truthcol->size() < 1 && fPassEmptySpills ) {
                LOG_DEBUG("GENIEGen")
                        << "no events made for this spill "
                        << " continue on until we get a spill with events";
                break;
            }

        } // end loop while no interactions are made



        // put the collections in the event
        evt.put(std::move(truthcol));
        evt.put(std::move(fluxcol));
        evt.put(std::move(gtruthcol));
        evt.put(std::move(assns));
        evt.put(std::move(tgtassn));

        ++fSpillCounter;

        // If fEventsPerSpill > 0, then we should update by what we actually used.
        // Otherwise, we are in "POT per Spill" mode so we should update by the cutoff
        // value, fPOTPerSpill.

        double ThisSpillPoT = fPOTPerSpill;
        if(fEventsPerSpill > 0 || passedPOTLimit) ThisSpillPoT=SpillExpTemp;
        fTotalExposure += ThisSpillPoT;

        std::unique_ptr<sumdata::SpillData> sd(new sumdata::SpillData);
        sd->spillpot = ThisSpillPoT;
        sd->goodbeam = 1;

        // Fill in some typical values to hopefully keep good beam cuts happy
        // without having to special-case for MC.
        sd->deltaspilltimensec = 0;
        sd->hornI = -199.5;
        sd->posx = 0.9;
        sd->posy = 1.3;
        sd->widthx = 1.2;
        sd->widthy = 1.2;
        //based on data
/*
    sd->bposx[0]=-0.8;
    sd->bposx[1]=-0.2;
    sd->bposx[2]=-0.04;
    sd->bposx[3]=0.05;
    sd->bposx[4]=-0.008;
    sd->bposx[5]=-0.007;

    for(int i=0;i<6;i++){
      sd->bposy[i]=0.5;
      //assuming no slip-stacking
      sd->intx[i]=ThisSpillPoT/6;
      sd->inty[i]=ThisSpillPoT/6;
    }
 */
        evt.put(std::move(sd));
    }

}// namespace

namespace evgen { DEFINE_ART_MODULE(GENIEGen); }