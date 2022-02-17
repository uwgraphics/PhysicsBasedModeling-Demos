#include <iostream>

#include "FiniteElementMesh.h"

#include <Eigen/Dense>
#include <random>

template<class T>
struct LatticeMesh : public FiniteElementMesh<T>
{
    using Base = FiniteElementMesh<T>;

    // from AnimatedMesh
    using Base::m_meshElements;
    using Base::m_particleX;
    using Base::initializeUSD;
    using Base::initializeTopology;
    using Base::initializeParticles;
    using Vector2 = typename Base::Vector2;

    // from FiniteElementMesh
    using Base::initializeUndeformedConfiguration;
    using Base::m_stepEndTime;

    std::array<int, 2> m_cellSize; // dimensions in grid cells
    T m_gridDX;

    const int m_pinchRadius;

    std::vector<Vector2> m_particleUndeformedX;
    std::vector<int> m_leftHandleIndices;
    std::vector<int> m_rightHandleIndices;
    Vector2 m_leftHandleDisplacement;
    Vector2 m_rightHandleDisplacement;
    
    LatticeMesh()
        :Base(1.e2, 1., 4., .05), m_pinchRadius(1)
    {
        m_leftHandleDisplacement  = Vector2(-1.5, 0.);
        m_rightHandleDisplacement = Vector2( 1.5, 0.);
    }

    void initialize()
    {
        initializeUSD("MembraneStretchPD.usda");

        // Create a Cartesian lattice topology
        for(int cell_i = 0; cell_i < m_cellSize[0]; cell_i++)
        for(int cell_j = 0; cell_j < m_cellSize[1]; cell_j++){
            m_meshElements.emplace_back(
                std::array<int, 3>{
                    gridToParticleID(cell_i  , cell_j  ), 
                    gridToParticleID(cell_i+1, cell_j  ),
                    gridToParticleID(cell_i+1, cell_j+1)
                }
            );
            m_meshElements.emplace_back(
                std::array<int, 3>{
                    gridToParticleID(cell_i  , cell_j  ), 
                    gridToParticleID(cell_i+1, cell_j+1),
                    gridToParticleID(cell_i  , cell_j+1)
                }
            );
        }
        initializeTopology();

        // Also initialize the associated particles
        for(int node_i = 0; node_i <= m_cellSize[0]; node_i++)
        for(int node_j = 0; node_j <= m_cellSize[1]; node_j++)
            m_particleX.emplace_back(m_gridDX * (T)node_i, m_gridDX * (T)node_j);
        initializeParticles();

        // Check particle indexing in mesh
        for(const auto& element: m_meshElements)
            for(const auto vertex: element)
                if(vertex < 0 || vertex >= m_particleX.size())
                    throw std::logic_error("mismatch between mesh vertex and particle array");

        // Initialize rest shape matrices and particle mass
        initializeUndeformedConfiguration();

        // Also record rest shape
        m_particleUndeformedX = m_particleX;

        // Identify particles on left and right handles
        for(int node_j = 0; node_j <= m_cellSize[1]; node_j++){
            m_leftHandleIndices.push_back(gridToParticleID(0, node_j));
            m_rightHandleIndices.push_back(gridToParticleID(m_cellSize[0], node_j));
        }
    }

    void applyPerturbation()
    {
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> dis(-1.0, +1.0);

        for (auto& X : m_particleX)
            for (int d = 0; d < 2; d++)
                X[d] += dis(gen);
    }

    void clearConstrainedParticles(std::vector<Vector2>& x) override
    { 
        for(const auto v: m_leftHandleIndices)
            x[v] = Vector2::Zero();
        for(const auto v: m_rightHandleIndices)
            x[v] = Vector2::Zero();       
    }

    void setBoundaryConditions() override
    { 
        T effectiveTime = std::min<T>(m_stepEndTime, 1.0);
        
        for(const auto v: m_leftHandleIndices){
            m_particleX[v] = m_particleUndeformedX[v] + effectiveTime * m_leftHandleDisplacement;
        }
        for(const auto v: m_rightHandleIndices){
            m_particleX[v] = m_particleUndeformedX[v] + effectiveTime * m_rightHandleDisplacement;
        }
    }

private:
    inline int gridToParticleID(const int i, const int j) const { return i * (m_cellSize[1]+1) + j; }
};

int main(int argc, char *argv[])
{
    LatticeMesh<float> simulationMesh;
    simulationMesh.m_cellSize = { 40, 40 };
    simulationMesh.m_gridDX = 0.05;
    simulationMesh.m_nFrames = 23;
    simulationMesh.m_subSteps = 1;
    simulationMesh.m_frameDt = 0.05;

    // Initialize the simulation example
    simulationMesh.initialize();

    // Output the initial shape of the mesh
    simulationMesh.writeFrame(0);

    // Perform the animation, output results at each frame
    for(int frame = 1; frame <= simulationMesh.m_nFrames; frame++){
        simulationMesh.simulateFrame(frame);
        if (frame%8 == 0)
            simulationMesh.applyPerturbation();
        simulationMesh.writeFrame(frame);
    }

    // Write the entire timeline to USD
    simulationMesh.writeUSD();

    return 0;
}

