/**
 *  @file   larpandoracontent/LArMonitoring/TopologyCalorimetryMonitoring.h
 *
 *  @brief  Header file for the Topology and Calorimetry Monitoring Algorithm 
 *
 *  $Log: $
 */

#ifndef LAR_TOPOLOGY_CALORIMETRY_MONITORING_ALGORITHM_H
#define LAR_TOPOLOGY_CALORIMETRY_MONITORING__ALGORITHM_H 1

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  Topology and Calorimetry Monitoring class
 */
class TopologyCalorimetryMonitoringAlgorithm : public pandora::Algorithm
{
public:
/**
 *  @brief  Default constructor
 */
    TopologyCalorimetryMonitoringAlgorithm();
private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Topology and Calorimetry Monitoring. Uses Calohits and pca to fit direction of clusters in 2 views and extracts calorimetric information corresponding to axis.
     *
     *  @param  pClusterList the list of clusters
     *  @param  listName the name of the current cluster list (ClustersU, ClustersV, ClustersW normally)
     */
    void TopologyCalorimetryMonitoring(const pandora::ClusterList *const pClusterList) const;

     /**
      *  @brief  Get the MC particle for a given cluster, caching to a map.
      *
      *  @param  cluster         The current cluster to lookup
      *  @param  clusterToMCMap  Map from Cluster to MC to cache results.
      */
    const pandora::MCParticle* GetMCForCluster(const pandora::Cluster *const cluster, std::map<const pandora::Cluster*,
        const pandora::MCParticle*> &clusterToMCMap) const;

    pandora::StringVector  m_inputClusterListNames; ///< The names of the input cluster lists.
    std::string            m_mcParticleListName;    ///< Input MC particle list name.
    std::string            m_treename;              ///< Input name for tree 

};

} // namespace lar_content
#endif // #ifndef LAR_TOPOLOGY_CALORIMETRY_MONITORING_ALGORITHM_H
