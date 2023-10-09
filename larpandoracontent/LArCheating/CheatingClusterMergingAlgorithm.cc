/**
 * @file   larpandoracontent/LArCheating/CheatingClusterMergingAlgorithm.cc
 *
 * @brief  Implementation file for the cheating cluster merging algorithm.
 *
 * $Log: 
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Helpers/MCParticleHelper.h"
#include "Objects/MCParticle.h"

#include <cmath>

#include "larpandoracontent/LArCheating/CheatingClusterMergingAlgorithm.h"
#include "larpandoracontent/LArTrackShowerId/ShowerGrowingAlgorithm.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"

using namespace pandora;

namespace lar_content
{

CheatingClusterMergingAlgorithm::CheatingClusterMergingAlgorithm() :
    m_maxClusterFraction(0.25f),
    m_minNCaloHits(1),
    m_writeTree{true}
{
    if (m_maxClusterFraction < 0)
        std::cout << "WARN: m_maxClusterFraction can not be below 0!" << std::endl;
}
//------------------------------------------------------------------------------------------//

StatusCode CheatingClusterMergingAlgorithm::Run()
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

            this->CheatedClusterMerging(pClusterList, clusterListName);
            std::cout << "---------------------------------------------------------------------------------" << std::endl;   
     
        }
        
        catch (StatusCodeException &statusCodeException)
        {
            throw statusCodeException;
        }
        
        if (m_writeTree)
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "RECREATE"))
         }

    }

    return STATUS_CODE_SUCCESS;
}

//--------------------------------------------------------------------------------------------//
const MCParticle* CheatingClusterMergingAlgorithm::GetMCForCluster(const Cluster *const cluster, std::map<const Cluster*,
    const MCParticle*> &clusterToMCMap) const
{
    const MCParticle* clusterMC = nullptr;

    if (clusterToMCMap.count(cluster) > 0)
    {
        clusterMC = clusterToMCMap.at(cluster);
    }
    else
    {
        try
        {
            clusterMC = MCParticleHelper::GetMainMCParticle(cluster);
            clusterToMCMap[cluster] = clusterMC;
        }
        catch (StatusCodeException e)
        {
            std::cout << "Failed to get MC particle for cluster of " << cluster->GetOrderedCaloHitList().size()
                      << " : " << e.ToString() << std::endl;
        }
    }

    return clusterMC;
}

//---------------------------------------------------------------------------------------------//

bool CheatingClusterMergingAlgorithm::IsValidToUse(const Cluster *const cluster, std::map<const Cluster*, bool> &clusterIsUsed) const
{

    if (!cluster->IsAvailable())
        return false;

    if (cluster->GetNCaloHits() < m_minNCaloHits)
        return false;

    if (clusterIsUsed.count(cluster) > 0)
        return false;

    return true;
}

//-----------------------------------------------------------------------------------------------//

void CheatingClusterMergingAlgorithm::CheatedClusterMerging(const pandora::ClusterList *const pClusterList, const std::string &listName) const
{
    std::map<const Cluster*, const MCParticle*> clusterToMCParticleMap;

    std::map<const Cluster*, bool> clusterIsUsed;
    std::map<const Cluster*, ClusterVector> clustersToMerge;
    
    std::cout << "Cluster List Name: " << listName << "   //   Clusters in List: " << pClusterList->size() <<  std::endl;
     
    for (auto it = pClusterList->begin(); it != pClusterList->end(); ++it)
    {
        const Cluster *cluster(*it);      
 
        std::cout << "  --- Hits in Cluster: " <<  cluster->GetNCaloHits() << "  // Cluster Length: " << LArClusterHelper::GetLength(cluster) << std::endl;     
        
        CartesianVector firstInnerCoordinate(0.f,0.f,0.f);
        CartesianVector firstOuterCoordinate(0.f,0.f,0.f);
         
        LArClusterHelper::GetExtremalCoordinates(cluster, firstInnerCoordinate, firstOuterCoordinate);   
       
        std::cout << "**1**  Inner Extremal Coordinate: " << firstInnerCoordinate << std::endl;        
        std::cout << "**1**  Outer Extremal Coordinate: " << firstOuterCoordinate << std::endl; 
        std::cout << std::endl;
  
        if (!this->IsValidToUse(cluster, clusterIsUsed))
            continue;

        //const MCParticle *clusterMC(this->GetMCForCluster(cluster, clusterToMCParticleMap));
        const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(cluster));
 
        std::cout << "~~~TEST: " << pMCParticle->GetParticleId() << std::endl;

        for (auto it2 = std::next(it); it2 != pClusterList->end(); ++it2) {
            const Cluster *const otherCluster(*it2);
            
            CartesianVector secondInnerCoordinate(0.f,0.f,0.f);
            CartesianVector secondOuterCoordinate(0.f,0.f,0.f);

            LArClusterHelper::GetExtremalCoordinates(otherCluster, secondInnerCoordinate, secondOuterCoordinate);

            std::cout << "   **2** Inner Extremal Coordinate: " << secondInnerCoordinate << std::endl;
            std::cout << "   **2** Outer Extremal Coordinate: " << secondOuterCoordinate << std::endl;
            std::cout << std::endl;

            const float diff1{Distance(firstInnerCoordinate,secondInnerCoordinate)};
            const float diff2{Distance(firstInnerCoordinate,secondOuterCoordinate)};
            const float diff3{Distance(firstOuterCoordinate,secondInnerCoordinate)};
            const float diff4{Distance(firstOuterCoordinate,secondOuterCoordinate)};
            const float distance{std::min({diff1, diff2, diff3, diff4})}; 
 
            const double ang1{Angle(firstInnerCoordinate,firstOuterCoordinate)};
            const double ang2{Angle(secondInnerCoordinate,secondOuterCoordinate)};
            const double angle{std::abs(ang1 - ang2)};
            
            std::cout << "   ****  Distance between clusters: " << distance << std::endl;
            std::cout << "   ****  Angle between clusters: " << angle << " Ang1: " << ang1 << " Ang2: " << ang2 << std::endl;

            if (!this->IsValidToUse(otherCluster, clusterIsUsed))
                continue;

            //const MCParticle *otherClusterMC(this->GetMCForCluster(otherCluster, clusterToMCParticleMap));
            const MCParticle *const pMCOtherParticle(MCParticleHelper::GetMainMCParticle(otherCluster));
            const int Pdg1{pMCParticle->GetParticleId()};
            const int Pdg2{pMCOtherParticle->GetParticleId()};
            bool sigBac{false};          
   
            std::cout << "   ****  MC Cluster 1:  " << Pdg1 << "  MC Cluster 2: " << Pdg2 << std::endl; 
  
            if (distance < 1000 && angle < 1)
            {
                sigBac = true;
                std::cout << "              !!!   These clusters should be merged !!!  " << sigBac << std::endl;
                clusterIsUsed[cluster] = true;
                clusterIsUsed[otherCluster] = true;
                clustersToMerge[cluster].push_back(otherCluster);
            }
           
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Distance", distance));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Angle", angle));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "FirstClusterMC", Pdg1));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "SecondClusterMC", Pdg2));
           // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "SignalvBackground", sigBac));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
        }
    }

 // What to do with clusters that are never merged?
 //     - Have a fallback check for second main MC, if over X%?
 
 for (auto clusterToMergePair : clustersToMerge)
    {
        const Cluster *currentCluster = clusterToMergePair.first;
        const auto clusters = clusterToMergePair.second;

        for (auto clusterToMerge : clusters)
        {
            if (! clusterToMerge->IsAvailable())
                continue;

            try
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, currentCluster, clusterToMerge, listName, listName));
            } catch (StatusCodeException) {}
        }
    }

    return;
}

//-------------------------------------------------------------------------------------------------

float CheatingClusterMergingAlgorithm::Distance(const CartesianVector vector1, const CartesianVector vector2) const
{
    const float dx{vector1.GetX() - vector2.GetX()};
    const float dz{vector1.GetZ() - vector2.GetZ()};
    float distance{std::sqrt(dx * dx + dz * dz)};
 
    return distance;
}

//-------------------------------------------------------------------------------------------------

double CheatingClusterMergingAlgorithm::Angle(const CartesianVector vector1, const CartesianVector vector2) const
{
    const double dx{vector1.GetX() - vector2.GetX()};
    const double dz{vector1.GetZ() - vector2.GetZ()};
    double angle{tan(dx / dz)};

    return angle;
}
 
//----------------------------------------------------------------------------------------------------

StatusCode CheatingClusterMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MCParticleListName", m_mcParticleListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxClusterFraction", m_maxClusterFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinNCaloHits", m_minNCaloHits));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TreeName", m_treeName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FileName", m_fileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteTree", m_writeTree));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
