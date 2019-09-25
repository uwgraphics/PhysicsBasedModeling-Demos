#include "FiniteElementMesh.h"

#include <Eigen/Dense>

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
    using Base::initializeUndeformedConfiguration;
    
    std::array<int, 2> m_cellSize; // dimensions in grid cells
    T m_gridDX;

    const int m_pinchRadius;
    
    LatticeMesh()
        :Base(1.e3, 1., 1., .1), m_pinchRadius(1)
    {}

    void initialize()
    {
        initializeUSD("MembranePinch6.usda");

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
    }

    void initializeDeformation()
    {
        // Pinch membrane as initial configuration
        for(int node_i = m_cellSize[0]/2-m_pinchRadius; node_i <= m_cellSize[0]/2+m_pinchRadius; node_i++)
        for(int node_j = m_cellSize[1]/2-m_pinchRadius; node_j <= m_cellSize[1]/2+m_pinchRadius; node_j++)
            m_particleX[gridToParticleID(node_i,node_j)] = Vector2(
                m_gridDX * (T)node_i + m_gridDX * 0.1 * (T) (m_cellSize[0]+m_cellSize[1]),
                m_gridDX * (T)node_j + m_gridDX * 0.1 * (T) (m_cellSize[0]+m_cellSize[1])
            );

        // Relax every interior node
        for(int iteration = 0; iteration < 1000; iteration++)
            for(int node_i = 1; node_i < m_cellSize[0]; node_i++)
            for(int node_j = 1; node_j < m_cellSize[1]; node_j++){

                if( std::abs(node_i - m_cellSize[0]/2) <= m_pinchRadius &&
                    std::abs(node_j - m_cellSize[1]/2) <= m_pinchRadius ) continue;

                int pCenter = gridToParticleID(node_i  ,node_j  );
                int pPlusX  = gridToParticleID(node_i+1,node_j  );
                int pMinusX = gridToParticleID(node_i-1,node_j  );
                int pPlusY  = gridToParticleID(node_i  ,node_j+1);
                int pMinusY = gridToParticleID(node_i  ,node_j-1);
                
                m_particleX[pCenter] = .25 * ( m_particleX[pPlusX] + m_particleX[pMinusX] + m_particleX[pPlusY] + m_particleX[pMinusY]);
            }
    }

private:
    inline int gridToParticleID(const int i, const int j) const { return i * (m_cellSize[1]+1) + j; }
};

int main(int argc, char *argv[])
{
    LatticeMesh<float> simulationMesh;
    simulationMesh.m_cellSize = { 40, 40 };
    simulationMesh.m_gridDX = 0.025;
    simulationMesh.m_nFrames = 100;
    simulationMesh.m_subSteps = 1;
    simulationMesh.m_frameDt = 0.5;

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

