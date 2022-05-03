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
    using Base::m_particleV;
    using Base::m_particleMass;
    using Base::initializeUndeformedConfiguration;
    using Base::m_stepEndTime;

    std::array<int, 2> m_cellSize; // dimensions in grid cells
    T m_gridDX;

    const int m_pinchRadius;

    std::vector<Vector2> m_particleUndeformedX;
    std::vector<int> m_leftHandleIndices;
    std::vector<int> m_rightHandleIndices;
    Vector2 m_leftHandleVelocity;
    Vector2 m_rightHandleVelocity;
    
    LatticeMesh()
        :Base(1.e0, 5., 4., .05), m_pinchRadius(1)
    {
        // m_leftHandleVelocity  = Vector2(-.2, 0.);
        // m_rightHandleVelocity = Vector2( .2, 0.);
        m_leftHandleVelocity  = Vector2::Zero();
        m_rightHandleVelocity = Vector2::Zero();
    }

    void initialize()
    {
        initializeUSD("GroundCollision1.usda");

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

        // Also resize the velocities to match
        m_particleV.resize(m_particleX.size(), Vector2::Zero());

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

    void initializeDeformation()
    {
        // No need to apply any deformation; this example is driven by moving handles
    }

    virtual void addExternalForce(std::vector<Vector2>& f)
    {
        Vector2 aGravity( 0., -9.81 );

        for (int p = 0; p < m_particleX.size(); p++)
            f[p] += m_particleMass[p] * aGravity;
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
            m_particleX[v] = m_particleUndeformedX[v] + effectiveTime * m_leftHandleVelocity;
            m_particleV[v] = m_leftHandleVelocity;
        }
        for(const auto v: m_rightHandleIndices){
            m_particleX[v] = m_particleUndeformedX[v] + effectiveTime * m_rightHandleVelocity;
            m_particleV[v] = m_rightHandleVelocity;
        }
    }

private:
    inline int gridToParticleID(const int i, const int j) const { return i * (m_cellSize[1]+1) + j; }
};

int main(int argc, char *argv[])
{
    LatticeMesh<float> simulationMesh;
    simulationMesh.m_cellSize = { 20, 20 };
    simulationMesh.m_gridDX = 0.05;
    simulationMesh.m_nFrames = 30;
    simulationMesh.m_subSteps = 1;
    simulationMesh.m_frameDt = 0.1;

    // Initialize the simulation example
    simulationMesh.initialize();
    simulationMesh.initializeDeformation();

    // Output the initial shape of the mesh
    simulationMesh.writeFrame(0);

    // Perform the animation, output results at each frame
    for(int frame = 1; frame <= simulationMesh.m_nFrames; frame++){
        simulationMesh.simulateFrame(frame);
        simulationMesh.writeFrame(frame);
    }

    // Write the entire timeline to USD
    simulationMesh.writeUSD();

    return 0;
}

