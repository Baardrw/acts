#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>

#include <sys/types.h>

namespace ActsExamples {

struct MDTWrite {
  Acts::Vector3 surfacePosition;
  Acts::Vector3 hitPosition;
  Acts::Vector3 pointOfClosestApproach;

  // General hit information for validation
  Acts::Vector3 hitDirection;
  uint64_t eventID;
  uint64_t volumeID;
  uint64_t time;
};

struct RPCWrite {
  Acts::Vector3 stripPosition;
  Acts::Vector3 stripDirection;
  Acts::Vector3 stripNormal;

  // General hit information for validation
  Acts::Vector3 hitPosition;
  Acts::Vector3 pointOfClosestApproach;
  uint64_t eventID;
  uint64_t volumeID;
  uint64_t time;
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

  void processMDT(const SimHit& hit, const Acts::Surface* surface,
                  const AlgorithmContext& ctx,
                  const Acts::Transform3& toReferenceFrame,
                  const Acts::SquareMatrix3& toReferenceFrameLinear,
                  std::vector<MDTWrite>& mdtWrites) const;

  void processRPC(const SimHit& hit, const Acts::Surface* surface,
                              const AlgorithmContext& ctx,
                              const Acts::Transform3& toReferenceFrame,
                              const Acts::SquareMatrix3& toReferenceFrameLinear,
                              std::vector<RPCWrite>& rpcWrites) const;
  
 private:
  Config m_cfg;

  ReadDataHandle<SimHitContainer> m_simHits{this, "simHits"};
  ReadDataHandle<SimParticleContainer> m_outputParticles{this,
                                                         "simOutputParticles"};
};

}  // namespace ActsExamples
