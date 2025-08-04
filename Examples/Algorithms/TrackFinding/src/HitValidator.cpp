#include "ActsExamples/TrackFinding/HitValidator.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>

using namespace Acts;
namespace ActsExamples {

HitValidator::HitValidator(HitValidator::Config cfg, Acts::Logging::Level lvl)
    : IAlgorithm("HitValidator", lvl), m_cfg(std::move(cfg)) {
  m_simHits.initialize(m_cfg.simHits);
  m_outputParticles.initialize(m_cfg.particlesSimulated);
}

ProcessCode HitValidator::execute(const AlgorithmContext& ctx) const {
  const SimHitContainer& simHits = m_simHits(ctx);
  const SimParticleContainer& particles = m_outputParticles(ctx);

  std::vector<HitWrite> hitWrites;
  hitWrites.reserve(simHits.size());

  bool first = true;
  uint64_t chosenVolume = 0;
  uint64_t chosenParticleId = 0;
  const Surface* referenceSurface = nullptr;
  const Surface* localOriginSurface = nullptr;

  for (const SimHit& hit : simHits) {
    const Surface* surface =
        m_cfg.trackingGeometry->findSurface(hit.geometryId());

    if (first) {
      first = false;
      chosenVolume = hit.geometryId().volume();
      chosenParticleId = hit.particleId().value();
      referenceSurface = surface;
      std::cout << "Chosen volume: " << chosenVolume
                << ", Chosen particle ID: " << chosenParticleId << std::endl;
    }
    if (chosenParticleId != hit.particleId().value()) {
      continue;  // Skip hits not in the chosen volume
    }

    // Transform to local coordinates
    const Vector3 mdtCenter = Vector3(0, 0, 0);
    const Vector3 hitSamplePos =
        surface->transform(ctx.geoContext).inverse() * hit.position();
    const Vector3 hitDirection =
        surface->transform(ctx.geoContext).linear().inverse() * hit.direction();
    const Vector3 k = -hitSamplePos;

    const double lambda =
        (k.dot(hitDirection) -
         k.dot(Vector3::UnitZ()) * Vector3::UnitZ().dot(hitDirection)) /
        (1 - Vector3::UnitZ().dot(hitDirection) *
                 Vector3::UnitZ().dot(hitDirection));

    const Vector3 pointOfClosestApproach = hitSamplePos + hitDirection * lambda;

    Transform3 toReferenceFrame =
        referenceSurface->transform(ctx.geoContext).inverse() *
        surface->transform(ctx.geoContext);
    Acts::SquareMatrix3 toReferenceFrameLinear = toReferenceFrame.linear();

    ACTS_VERBOSE("HitValidator: " << "Surface: " << surface->geometryId()
                                  << ", Position: " << hit.position()
                                  << ", Direction: " << hit.direction()
                                  << ", Particle ID: " << hit.particleId());

    HitWrite hitWrite;
    hitWrite.surfacePosition = toReferenceFrame * mdtCenter;
    hitWrite.hitPosition = toReferenceFrame * hitSamplePos;
    hitWrite.pointOfClosestApproach = toReferenceFrame * pointOfClosestApproach;
    hitWrite.hitDirection = toReferenceFrameLinear * hitDirection;
    hitWrite.eventID = ctx.eventNumber;
    hitWrite.volumeID = hit.geometryId().volume();
    hitWrite.time = hit.time();

    // Swap x and z coordinates to match the expected format of the math
    hitWrite.surfacePosition =
        Vector3(hitWrite.surfacePosition.z(), hitWrite.surfacePosition.y(),
                hitWrite.surfacePosition.x());
    hitWrite.hitPosition =
        Vector3(hitWrite.hitPosition.z(), hitWrite.hitPosition.y(),
                hitWrite.hitPosition.x());
    hitWrite.pointOfClosestApproach =
        Vector3(hitWrite.pointOfClosestApproach.z(),
                hitWrite.pointOfClosestApproach.y(),
                hitWrite.pointOfClosestApproach.x());
    hitWrite.hitDirection =
        Vector3(hitWrite.hitDirection.z(), hitWrite.hitDirection.y(),
                hitWrite.hitDirection.x());

    hitWrites.push_back(hitWrite);
  }

  // Sort by hitTime to ensure correct ordering in output
  std::sort(
      hitWrites.begin(), hitWrites.end(),
      [](const HitWrite& a, const HitWrite& b) { return a.time < b.time; });

  bool fileExists = std::filesystem::exists("hit_validation.csv");
  std::ofstream outFile("hit_validation.csv", std::ios::app);
  if (!fileExists) {
    outFile << "\nSurfacePositionX,SurfacePositionY,SurfacePositionZ,"
               "HitPositionX,HitPositionY,HitPositionZ,"
               "PointOfClosestApproachX,PointOfClosestApproachY,"
               "PointOfClosestApproachZ,HitDirectionX,HitDirectionY,"
               "HitDirectionZ,eventID,volumeID\n";
  }

  for (const auto& hitWrite : hitWrites) {
    outFile << hitWrite.surfacePosition.x() << ","
            << hitWrite.surfacePosition.y() << ","
            << hitWrite.surfacePosition.z() << "," << hitWrite.hitPosition.x()
            << "," << hitWrite.hitPosition.y() << ","
            << hitWrite.hitPosition.z() << ","
            << hitWrite.pointOfClosestApproach.x() << ","
            << hitWrite.pointOfClosestApproach.y() << ","
            << hitWrite.pointOfClosestApproach.z() << ","
            << hitWrite.hitDirection.x() << "," << hitWrite.hitDirection.y()
            << "," << hitWrite.hitDirection.z() << "," << hitWrite.eventID
            << "," << hitWrite.volumeID << "\n";
  }

  outFile.close();

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
