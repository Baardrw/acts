#include "ActsExamples/TrackFinding/HitValidator.hpp"

#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Utilities/GridAccessHelpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

using namespace Acts;
using SurfaceType = Acts::RegularSurface::SurfaceType;

namespace ActsExamples {

HitValidator::HitValidator(HitValidator::Config cfg, Acts::Logging::Level lvl)
    : IAlgorithm("HitValidator", lvl), m_cfg(std::move(cfg)) {
  m_simHits.initialize(m_cfg.simHits);
  m_outputParticles.initialize(m_cfg.particlesSimulated);
}

void HitValidator::processMDT(const SimHit& hit, const Surface* surface,
                              const AlgorithmContext& ctx,
                              const Transform3& toReferenceFrame,
                              const Acts::SquareMatrix3& toReferenceFrameLinear,
                              std::vector<MDTWrite>& mdtWrites) const {
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
      (1 -
       Vector3::UnitZ().dot(hitDirection) * Vector3::UnitZ().dot(hitDirection));

  const Vector3 pointOfClosestApproach = hitSamplePos + hitDirection * lambda;

  ACTS_VERBOSE("HitValidator: " << "Surface: " << surface->geometryId()
                                << ", Position: " << hit.position()
                                << ", Direction: " << hit.direction()
                                << ", Particle ID: " << hit.particleId());

  MDTWrite mdtWrite;
  mdtWrite.surfacePosition = toReferenceFrame * mdtCenter;
  mdtWrite.hitPosition = toReferenceFrame * hitSamplePos;
  mdtWrite.pointOfClosestApproach = toReferenceFrame * pointOfClosestApproach;
  mdtWrite.hitDirection = toReferenceFrameLinear * hitDirection;
  mdtWrite.eventID = ctx.eventNumber;
  mdtWrite.volumeID = hit.geometryId().volume();
  mdtWrite.time = hit.time();

  // Swap x and z coordinates to match the expected format of the math
  mdtWrite.surfacePosition =
      Vector3(mdtWrite.surfacePosition.z(), mdtWrite.surfacePosition.y(),
              mdtWrite.surfacePosition.x());
  mdtWrite.hitPosition =
      Vector3(mdtWrite.hitPosition.z(), mdtWrite.hitPosition.y(),
              mdtWrite.hitPosition.x());
  mdtWrite.pointOfClosestApproach = Vector3(
      mdtWrite.pointOfClosestApproach.z(), mdtWrite.pointOfClosestApproach.y(),
      mdtWrite.pointOfClosestApproach.x());
  mdtWrite.hitDirection =
      Vector3(mdtWrite.hitDirection.z(), mdtWrite.hitDirection.y(),
              mdtWrite.hitDirection.x());

  mdtWrites.push_back(mdtWrite);

  ACTS_VERBOSE("HitValidator: " << "Surface: " << surface->geometryId()
                                << ", MDT_Center: " << mdtCenter << "\n"
                                << ", MDT_Center: "
                                << mdtWrite.surfacePosition);
}

void HitValidator::processRPC(const SimHit& hit, const Surface* surface,
                              const AlgorithmContext& ctx,
                              const Transform3& toReferenceFrame,
                              const Acts::SquareMatrix3& toReferenceFrameLinear,
                              std::vector<RPCWrite>& rpcWrites) const {
  
  Vector3 stripDirection = Vector3::UnitZ(); // Unit Z and Y are along the beamline
  Vector3 stripNormal = Vector3::UnitY();
  
  const Vector3 rpcCenter = surface->transform(ctx.geoContext).inverse() *
                            surface->center(ctx.geoContext);
  const Vector3 hitSamplePos =
      surface->transform(ctx.geoContext).inverse() * hit.position();
  const Vector3 hitDirection =
      surface->transform(ctx.geoContext).linear().inverse() * hit.direction();

  // Calculate the point of closest approach as the point along the
  // directionVector that intersects the plane Since plane normal is unit Z, the
  // t that gives the line the same Z coordinate as the RPC center is the
  // solution

  const float t = (-hitSamplePos.z()) / hitDirection.z();
  const Vector3 pointOfClosestApproach = hitSamplePos + hitDirection * t;

  RPCWrite rpcWrite;
  rpcWrite.stripPosition = toReferenceFrame * rpcCenter;
  rpcWrite.stripDirection =  toReferenceFrameLinear * stripDirection;
  rpcWrite.stripNormal =  toReferenceFrameLinear * stripNormal;
  rpcWrite.hitPosition = toReferenceFrame * hitSamplePos;
  rpcWrite.pointOfClosestApproach = toReferenceFrame * pointOfClosestApproach;
  rpcWrite.eventID = ctx.eventNumber;
  rpcWrite.volumeID = hit.geometryId().volume();
  rpcWrite.time = hit.time();

  // Swap x and z coordinates to match the expected format of the math
  rpcWrite.stripPosition =
      Vector3(rpcWrite.stripPosition.z(), rpcWrite.stripPosition.y(),
              rpcWrite.stripPosition.x());
  rpcWrite.stripDirection =
      Vector3(rpcWrite.stripDirection.z(), rpcWrite.stripDirection.y(),
              rpcWrite.stripDirection.x());
  rpcWrite.stripNormal =
      Vector3(rpcWrite.stripNormal.z(), rpcWrite.stripNormal.y(),
              rpcWrite.stripNormal.x());
  rpcWrite.hitPosition =
      Vector3(rpcWrite.hitPosition.z(), rpcWrite.hitPosition.y(),
              rpcWrite.hitPosition.x());
  rpcWrite.pointOfClosestApproach = Vector3(
      rpcWrite.pointOfClosestApproach.z(), rpcWrite.pointOfClosestApproach.y(),
      rpcWrite.pointOfClosestApproach.x());
  rpcWrites.push_back(rpcWrite);
}

ProcessCode HitValidator::execute(const AlgorithmContext& ctx) const {
  const SimHitContainer& simHits = m_simHits(ctx);
  // const SimParticleContainer& particles = m_outputParticles(ctx);

  std::vector<MDTWrite> mdtWrites;
  std::vector<RPCWrite> rpcWrites;
  std::vector<SimHit>
      rpcBacklog;  // rpcs must be processed after the MDTs such that they can
                   // be expressed in the same reference frame

  bool first = true;
  uint64_t chosenParticleId = 0;
  uint64_t chosenVolume;

  const Surface* referenceSurface = nullptr;

  for (const SimHit& hit : simHits) {
    const Surface* surface =
        m_cfg.trackingGeometry->findSurface(hit.geometryId());

    if (first && surface->type() == SurfaceType::Straw) {  // Must be a MDT hit
      first = false;
      chosenVolume = hit.geometryId().volume();
      chosenParticleId = hit.particleId().value();
      referenceSurface = surface;

      // Filter out rpc hits not for the chosen ID
      std::vector<SimHit> filteredHits;
      std::copy_if(rpcBacklog.begin(), rpcBacklog.end(),
                   std::back_inserter(filteredHits),
                   [chosenParticleId](const SimHit& h) {
                     return h.particleId().value() == chosenParticleId;
                   });
      rpcBacklog = filteredHits;

    } else if (first && surface->type() == SurfaceType::Plane) {
      // Add all first to backlog and filter out later
      rpcBacklog.push_back(hit);
      continue;  // Process RPCs after MDTs
    }
    if (chosenParticleId != hit.particleId().value()) {
      continue;  // Skip hits not in the chosen volume
    }

    if (surface->type() == SurfaceType::Plane) {
      rpcBacklog.push_back(hit);
      continue;  // Process RPCs after MDTs
    }

    Transform3 toReferenceFrame =
        referenceSurface->transform(ctx.geoContext).inverse() *
        surface->transform(ctx.geoContext);
    Acts::SquareMatrix3 toReferenceFrameLinear = toReferenceFrame.linear();
    processMDT(hit, surface, ctx, toReferenceFrame, toReferenceFrameLinear,
               mdtWrites);
  }

  // Process RPCs
  for (const SimHit& backlogHit : rpcBacklog) {
    // Process the backlog RPC hits with the chosen reference surface
    const Surface* surface =
        m_cfg.trackingGeometry->findSurface(backlogHit.geometryId());

    Transform3 toReferenceFrame =
        referenceSurface->transform(ctx.geoContext).inverse() *
        surface->transform(ctx.geoContext);

    Acts::SquareMatrix3 toReferenceFrameLinear = toReferenceFrame.linear();

    HitValidator::processRPC(backlogHit, surface, ctx, toReferenceFrame,
                             toReferenceFrameLinear, rpcWrites);
  }

  // Sort by hitTime to ensure correct ordering in output
  std::sort(
      mdtWrites.begin(), mdtWrites.end(),
      [](const MDTWrite& a, const MDTWrite& b) { return a.time < b.time; });

  std::sort(
      rpcWrites.begin(), rpcWrites.end(),
      [](const RPCWrite& a, const RPCWrite& b) { return a.time < b.time; });

  bool fileExists = std::filesystem::exists("mdt_hit_validation.csv");
  std::ofstream outFile("mdt_hit_validation.csv", std::ios::app);
  if (!fileExists) {
    outFile << "\nSurfacePositionX,SurfacePositionY,SurfacePositionZ,"
               "HitPositionX,HitPositionY,HitPositionZ,"
               "PointOfClosestApproachX,PointOfClosestApproachY,"
               "PointOfClosestApproachZ,HitDirectionX,HitDirectionY,"
               "HitDirectionZ,eventID,volumeID\n";
  }

  for (const auto& mdtWrite : mdtWrites) {
    outFile << mdtWrite.surfacePosition.x() << ","
            << mdtWrite.surfacePosition.y() << ","
            << mdtWrite.surfacePosition.z() << "," << mdtWrite.hitPosition.x()
            << "," << mdtWrite.hitPosition.y() << ","
            << mdtWrite.hitPosition.z() << ","
            << mdtWrite.pointOfClosestApproach.x() << ","
            << mdtWrite.pointOfClosestApproach.y() << ","
            << mdtWrite.pointOfClosestApproach.z() << ","
            << mdtWrite.hitDirection.x() << "," << mdtWrite.hitDirection.y()
            << "," << mdtWrite.hitDirection.z() << "," << mdtWrite.eventID
            << "," << mdtWrite.volumeID << "\n";
  }

  outFile.close();

  fileExists = std::filesystem::exists("rpc_hit_validation.csv");
  std::ofstream rpcOutFile("rpc_hit_validation.csv", std::ios::app);
  if (!fileExists) {
    rpcOutFile << "\nStripPositionX,StripPositionY,StripPositionZ,"
                  "StripDirectionX,StripDirectionY,StripDirectionZ,"
                  "StripNormalX,StripNormalY,StripNormalZ,"
                  "HitPositionX,HitPositionY,HitPositionZ,"
                  "PointOfClosestApproachX,PointOfClosestApproachY,"
                  "PointOfClosestApproachZ,eventID,volumeID\n";
  }

  for (const auto& rpcWrite : rpcWrites) {
    rpcOutFile << rpcWrite.stripPosition.x() << ","
               << rpcWrite.stripPosition.y() << ","
               << rpcWrite.stripPosition.z() << ","
               << rpcWrite.stripDirection.x() << ","
               << rpcWrite.stripDirection.y() << ","
               << rpcWrite.stripDirection.z() << "," << rpcWrite.stripNormal.x()
               << "," << rpcWrite.stripNormal.y() << ","
               << rpcWrite.stripNormal.z() << "," << rpcWrite.hitPosition.x()
               << "," << rpcWrite.hitPosition.y() << ","
               << rpcWrite.hitPosition.z() << ","
               << rpcWrite.pointOfClosestApproach.x() << ","
               << rpcWrite.pointOfClosestApproach.y() << ","
               << rpcWrite.pointOfClosestApproach.z() << "," << rpcWrite.eventID
               << "," << rpcWrite.volumeID << "\n";
  }

  rpcOutFile.close();
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
