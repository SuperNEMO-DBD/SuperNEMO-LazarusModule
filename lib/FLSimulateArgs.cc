// Ourselves
#include "FLSimulateArgs.h"

// Third party:
// - Boost:
#include <boost/algorithm/string.hpp>
// - Bayeux:
#include <bayeux/datatools/kernel.h>
#include <bayeux/datatools/urn.h>
#include <bayeux/datatools/urn_query_service.h>
#include <bayeux/datatools/utils.h>
#include <bayeux/mygsl/random_utils.h>

// This Project:
#include "FLSimulateCommandLine.h"
#include "FLSimulateErrors.h"
#include "falaise/config/property_reader.h"

namespace FLSimulate {

void do_postprocess(FLSimulateArgs& flSimParameters);

// static
FLSimulateArgs FLSimulateArgs::makeDefault() {
  FLSimulateArgs params;

  // Application specific parameters:
  params.logLevel = datatools::logger::PRIO_ERROR;
  params.userProfile = "normal";
  params.numberOfEvents = 1;
  params.doSimulation = true;
  params.doDigitization = false;
  // Identification of the experimental setup:
  params.experimentalSetupUrn = "";

  // Identification of the simulation setup:
  params.simulationSetupUrn = "urn:snemo:demonstrator:simulation:2.1";
  // Simulation manager internal parameters:
  params.simulationManagerParams.set_defaults();
  params.simulationManagerParams.interactive = false;
  params.simulationManagerParams.logging = "error";
  params.simulationManagerParams.manager_config_filename = "";
  // Seeding is auto (from system) unless explicit file supplied
  params.simulationManagerParams.input_prng_seeds_file = "";
  params.simulationManagerParams.vg_seed =
      mygsl::random_utils::SEED_AUTO;  // PRNG for the vertex generator
  params.simulationManagerParams.eg_seed =
      mygsl::random_utils::SEED_AUTO;  // PRNG for the primary event generator
  params.simulationManagerParams.shpf_seed =
      mygsl::random_utils::SEED_AUTO;  // PRNG for the back end MC hit
                                       // processors
  params.simulationManagerParams.mgr_seed =
      mygsl::random_utils::SEED_AUTO;  // PRNG for the Geant4 engine itself
  params.simulationManagerParams.output_profiles_activation_rule = "";

  // Variants support:
  params.variantConfigUrn = "";
  params.variantProfileUrn = "";
  params.variantSubsystemParams.config_filename = "";

  // Service support:
  params.servicesSubsystemConfigUrn = "urn:snemo:demonstrator:geometry:4.0";
  params.servicesSubsystemConfig = "";

  // I/O control:
  params.outputMetadataFile = "";
  params.embeddedMetadata = false;
  params.outputFile = "";

  return params;
}

void FLSimulateArgs::print(std::ostream& out_) const {
  static const std::string tag("|-- ");
  static const std::string last_tag("`-- ");
  out_ << "FLSimulate setup parameters: " << std::endl;
  out_ << tag << "logLevel                   = "
       << datatools::logger::get_priority_label(this->logLevel) << std::endl;
  out_ << tag << "userProfile                = " << userProfile << std::endl;
  out_ << tag << "numberOfEvents             = " << numberOfEvents << std::endl;
  out_ << tag << "doSimulation               = " << std::boolalpha
       << doSimulation << std::endl;
  out_ << tag << "doDigitization             = " << std::boolalpha
       << doDigitization << std::endl;
  out_ << tag << "experimentalSetupUrn       = " << experimentalSetupUrn
       << std::endl;
  out_ << tag << "simulationSetupUrn         = " << simulationSetupUrn
       << std::endl;
  out_ << tag << "simulationSetupConfig      = "
       << simulationManagerParams.manager_config_filename << std::endl;
  out_ << tag << "digitizationSetupUrn       = " << digitizationSetupUrn
       << std::endl;
  out_ << tag << "variantConfigUrn           = " << variantConfigUrn
       << std::endl;
  out_ << tag << "variantProfileUrn          = " << variantProfileUrn
       << std::endl;
  out_ << tag << "variantSubsystemParams     = "
       << variantSubsystemParams.config_filename << std::endl;
  out_ << tag << "servicesSubsystemConfigUrn = " << servicesSubsystemConfigUrn
       << std::endl;
  out_ << tag << "servicesSubsystemConfig    = " << servicesSubsystemConfig
       << std::endl;
  out_ << tag << "outputMetadataFile         = " << outputMetadataFile
       << std::endl;
  out_ << tag << "embeddedMetadata           = " << std::boolalpha
       << embeddedMetadata << std::endl;
  out_ << last_tag << "outputFile                 = " << outputFile
       << std::endl;
  return;
}

}  // namespace FLSimulate
