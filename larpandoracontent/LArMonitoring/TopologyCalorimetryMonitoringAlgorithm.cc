/**
 * @file   larpandoracontent/LArMonitoring/TopologyCalorimetryMonitoringAlgorithm.cc
 *
 * @brief  Implementation file for the topology and calorimetry monitoring algorithm.
 *
 * $Log: 
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Helpers/MCParticleHelper.h"
#include "Objects/MCParticle.h"

#include "larpandoracontent/LArMonitoring/TopologyCalorimetryMonitoringAlgorithm.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"

using namespace pandora;

namespace lar_content
{

TopologyCalorimetryMonitoringAlgorithm::TopologyCalorimetryMonitoringAlgorithm() : m_treename{"mc"}
{
}
//------------------------------------------------------------------------------------------//

StatusCode TopologyCalorimetryMonitoringAlgorithm::Run()
{
    for (const std::string &clusterListName : m_inputClusterListNames)
    {
        try
        {
            const ClusterList *pClusterList = nullptr;
            PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

            if (!pClusterList || pClusterList->empty())
            {
                if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                    std::cout << "ClusterMergingAlgorithm: unable to find cluster list " << clusterListName << std::endl;

                continue;
            }

             this->TopologyCalorimetryMonitoring(pClusterList);

        }
        catch (StatusCodeException &statusCodeException)
        {
            throw statusCodeException;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//--------------------------------------------------------------------------------------------//
void TopologyCalorimetryMonitoringAlgorithm::TopologyCalorimetryMonitoring(const pandora::ClusterList *const pClusterList) const 
{
    for(const Cluster *const pCluster : *pClusterList)
    {
        CaloHitList clusterCaloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterCaloHitList);

        CartesianVector centroid(0.f, 0.f, 0.f);
        LArPcaHelper::EigenVectors eigenVecs;
        LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
        LArPcaHelper::RunPca(clusterCaloHitList, centroid, eigenValues, eigenVecs);

        if(!eigenVecs.empty())
        {
            const float pcaX{eigenVecs.front().GetX()};
            const float pcaY{eigenVecs.front().GetY()};
            const float pcaZ{eigenVecs.front().GetZ()};

            for(const CartesianVector &axis : eigenVecs)
            {
                std::cout << axis  << std::endl;
            }

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "PcaX", pcaX));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "PcaY", pcaY));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "PcaZ", pcaZ));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
        }
    }

    return;
}

//-------------------------------------------------------------------------------------------------

StatusCode TopologyCalorimetryMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MCParticleListName", m_mcParticleListName));
 
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TreeName", m_treename));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
