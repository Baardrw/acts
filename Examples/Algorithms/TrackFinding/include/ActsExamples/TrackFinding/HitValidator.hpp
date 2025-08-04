#pragma once 

#include "Acts/Definitions/Algebra.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include <sys/types.h>


namespace ActsExamples {


struct HitWrite {
  Acts::Vector3 surfacePosition;
  Acts::Vector3 hitPosition;
  Acts::Vector3 pointOfClosestApproach;
  Acts::Vector3 hitDirection;
  uint64_t eventID;
  uint64_t volumeID;
  uint64_t time;  // Time of the hit, if available
};

class HitValidator final : public IAlgorithm {
 public:
  struct Config {
    // Empty for now
    std::string simHits;
    std::string particlesSimulated;
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
  };

  const Config& config() const { return m_cfg; }

  HitValidator(Config cfg, Acts::Logging::Level lvl);
  ProcessCode execute(const AlgorithmContext& ctx) const final;

 private:
   Config m_cfg;

  ReadDataHandle<SimHitContainer> m_simHits{this, "simHits"};
  ReadDataHandle<SimParticleContainer> m_outputParticles{this, "simOutputParticles"};
};

}  // namespace ActsExamples
