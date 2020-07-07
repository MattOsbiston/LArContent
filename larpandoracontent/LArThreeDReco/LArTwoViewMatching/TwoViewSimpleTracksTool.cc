/**
 *  @file   larpandoracontent/LArThreeDReco/LArShowerMatching/TwoViewSimpleTracksTool.cc
 *
 *  @brief  Implementation of the clear showers tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArThreeDReco/TwoViewMatching/TwoViewSimpleTracksTool.h"

using namespace pandora;

namespace lar_content
{

TwoViewSimpleTracksTool::TwoViewSimpleTracksTool() :
    m_minMatchedFraction(0.2f),
    m_minMatchingScore(0.98f),
    m_minMatchedSamplingPoints(40),
    m_minXOverlapFraction(0.5f)
    //m_minMatchedSamplingPointRatio(2)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewSimpleTracksTool::Run(ThreeViewShowersAlgorithm *const pAlgorithm, MatrixType &overlapMatrix)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    ProtoParticleVector protoParticleVector;
    this->FindBestTrack(overlapMatrix, protoParticleVector);

    const bool particlesMade(pAlgorithm->CreateThreeDParticles(protoParticleVector));
    return particlesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewSimpleTracksTool::FindBestTrack(const MatrixType &overlapMatrix, ProtoParticleVector &protoParticleVector) const
{
    ClusterVector sortedKeyClusters;
    overlapMatrix.GetSortedKeyClusters(sortedKeyClusters);

    for (const Cluster *const pKeyCluster : sortedKeyClusters)
    {
        if (!pKeyCluster->IsAvailable())
            continue;

        unsigned int nU(0), nV(0);
        MatrixType::ElementList elementList;
        overlapMatrix.GetConnectedElements(pKeyCluster, true, elementList, nU, nV);

        if (elementList.empty())
            continue;

        MatrixType::Element bestElement(elementList.back());

        if (!bestElement.GetOverlapResult().IsInitialized())
            continue;

        if ((NULL == bestElement.GetCluster1()) || (NULL == bestElement.GetCluster2()))
            continue;

        ProtoParticle protoParticle;
        protoParticle.m_clusterList.push_back(bestElement.GetCluster1());
        protoParticle.m_clusterList.push_back(bestElement.GetCluster2());
        protoParticleVector.push_back(protoParticle);

        return;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewSimpleTracksTool::PassesElementCuts(MatrixType::ElementList::const_iterator eIter) const
{
    if (eIter->GetOverlapResult().GetLocallyMatchedFraction() < m_minMatchedFraction)
        return false;

    if (eIter->GetOverlapResult().GetMatchingScore() < m_minMatchingScore)
        return false;

    if (eIter->GetOverlapResult().GetNMatchedReUpsampledSamplingPoints() < m_minMatchedSamplingPoints)
        return false;

    const TwoViewXOverlap &xOverlap(eIter->GetOverlapResult().GetTwoViewXOverlap());

    if ((xOverlap.GetXSpan0() > std::numeric_limits<float>::epsilon()) && (xOverlap.GetTwoViewXOverlapSpan() / xOverlap.GetXSpan0() > m_minXOverlapFraction) &&
        (xOverlap.GetXSpan1() > std::numeric_limits<float>::epsilon()) && (xOverlap.GetTwoViewXOverlapSpan() / xOverlap.GetXSpan1() > m_minXOverlapFraction))
    {
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewSimpleTracksTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedFraction", m_minMatchedFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchingScore", m_minMatchingScore));


    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPoints", m_minMatchedSamplingPoints));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlapFraction", m_minXOverlapFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPointRatio", m_minMatchedSamplingPointRatio));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
