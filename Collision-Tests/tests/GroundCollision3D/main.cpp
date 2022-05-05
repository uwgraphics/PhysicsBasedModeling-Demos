#include <iostream>

#include "FiniteElementMesh.h"

#include <Eigen/Dense>
#include <map>

template<class T>
struct LatticeMesh : public FiniteElementMesh<T>
{
    using Base = FiniteElementMesh<T>;

    // from AnimatedTetrahedonMesh
    using Base::m_meshElements;
    using Base::m_meshBoundaryElements;
    using Base::m_particleX;
    using Base::initializeUSD;
    using Base::initializeTopology;
    using Base::initializeParticles;
    using Vector3 = typename Base::Vector3;

    // from FiniteElementMesh
    using Base::m_surfaceParticles;
    using Base::m_particleV;
    using Base::m_particleMass;
    using Base::initializeUndeformedConfiguration;
    using Base::m_stepEndTime;

    std::array<int, 3> m_cellSize; // dimensions in grid cells
    T m_gridDX;

    std::vector<std::array<int, 3>> m_activeCells; // Marks the "active" cells in the lattice
    std::map<std::array<int, 3>, int> m_activeNodes; // Maps the "active" nodes to their particle index

    std::vector<Vector3> m_particleUndeformedX;
    std::vector<int> m_leftHandleIndices;
    std::vector<int> m_rightHandleIndices;
    Vector3 m_leftHandleVelocity;
    Vector3 m_rightHandleVelocity;
    
    LatticeMesh()
        :Base(1.e0, 5., 4., .05)
    {
        // m_leftHandleVelocity  = Vector3(-.2, 0., 0.);
        // m_rightHandleVelocity = Vector3( .2, 0., 0.);
        m_leftHandleVelocity  = Vector3::Zero();
        m_rightHandleVelocity = Vector3::Zero();
    }

    void initialize()
    {
        initializeUSD("GroundCollision3D.usda");

        // Activate cells within a sphere of radius m_radius (in cells)

        for(int cell_i = 0; cell_i < m_cellSize[0]; cell_i++)
        for(int cell_j = 0; cell_j < m_cellSize[1]; cell_j++)
        for(int cell_k = 0; cell_k < m_cellSize[2]; cell_k++)
                m_activeCells.push_back(std::array<int, 3>{cell_i, cell_j, cell_k});
        std::cout << "Created a model including " << m_activeCells.size() << " lattice cells" <<std::endl;

        // Create (uniquely numbered) particles at the node corners of active cells

        for(const auto& cell: m_activeCells){
            std::array<int, 3> node;
            for(node[0] = cell[0]; node[0] <= cell[0]+1; node[0]++)
            for(node[1] = cell[1]; node[1] <= cell[1]+1; node[1]++)
            for(node[2] = cell[2]; node[2] <= cell[2]+1; node[2]++){
                auto search = m_activeNodes.find(node);
                if(search == m_activeNodes.end()){ // Particle not yet created at this lattice node location -> make one
                    m_activeNodes.insert({node, m_particleX.size()});
                    m_particleX.emplace_back(m_gridDX * T(node[0]), m_gridDX * T(node[1]), m_gridDX * T(node[2]));
                }
            }
        }
        std::cout << "Model contains " << m_particleX.size() << " particles" << std::endl;

        // Make tetrahedra out of all active cells (6 tetrahedra per cell)

        for(const auto& cell: m_activeCells){
            int vertexIndices[2][2][2];
            for(int i = 0; i <= 1; i++)
            for(int j = 0; j <= 1; j++)
            for(int k = 0; k <= 1; k++){
                std::array<int, 3> node{cell[0] + i, cell[1] + j, cell[2] + k};
                auto search = m_activeNodes.find(node);
                if(search != m_activeNodes.end())
                    vertexIndices[i][j][k] = search->second;
                else
                    throw std::logic_error("particle at cell vertex not found");
            }

            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][0][0], vertexIndices[1][1][0], vertexIndices[1][1][1]});
            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][0][0], vertexIndices[1][1][1], vertexIndices[1][0][1]});
            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][0][1], vertexIndices[1][1][1], vertexIndices[0][0][1]});
            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][1][1], vertexIndices[0][1][1], vertexIndices[0][0][1]});
            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][1][1], vertexIndices[0][1][0], vertexIndices[0][1][1]});
            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][1][0], vertexIndices[0][1][0], vertexIndices[1][1][1]});
        }
        
        // Perform the USD-specific initialization of topology & particles
        // (this will also create a boundary *surface* to visualuze

        initializeTopology();
        initializeParticles();

        // Record surface particles (for collisions)
        std::vector<bool> particleOnSurface(m_particleX.size());
        for(auto element: m_meshBoundaryElements) {
            for (auto p: element)
                if (!particleOnSurface[p]){
                    m_surfaceParticles.push_back(p);
                    particleOnSurface[p] = true;
                }
        }

        // Check particle indexing in mesh

        for(const auto& element: m_meshElements)
            for(const auto vertex: element)
                if(vertex < 0 || vertex >= m_particleX.size())
                    throw std::logic_error("mismatch between mesh vertex and particle array");

        // Also resize the velocities to match
        m_particleV.resize(m_particleX.size(), Vector3::Zero());

        // Initialize rest shape matrices and particle mass
        initializeUndeformedConfiguration();

        // Also record rest shape
        m_particleUndeformedX = m_particleX;

        // Identify particles on left and right handles
#if 0
        for(int node_j = 0; node_j <= m_cellSize[1]; node_j++)
            for(int node_k = 0; node_k <= m_cellSize[2]; node_k++) {
                {
                    std::array<int, 3> node = {0, node_j, node_k}; // left side
                    auto search = m_activeNodes.find(node);
                    
                    if(search == m_activeNodes.end())
                        throw std::logic_error("expected to find particle at left side");
                    else
                        m_leftHandleIndices.push_back(search->second);
                }
                {
                    std::array<int, 3> node = {m_cellSize[0], node_j, node_k}; // left side
                    auto search = m_activeNodes.find(node);
                    
                    if(search == m_activeNodes.end())
                        throw std::logic_error("expected to find particle at right side");
                    else
                        m_rightHandleIndices.push_back(search->second);
                }
            }
#endif
    }

    void initializeDeformation()
    {
        // No need to apply any deformation; this example is driven by moving handles
    }

    virtual void addExternalForce(std::vector<Vector3>& f) override {
        Vector3 aGravity( 0., -9.81, 0. );

        for (int p = 0; p < m_particleX.size(); p++)
            f[p] += m_particleMass[p] * aGravity;
    }

    void clearConstrainedParticles(std::vector<Vector3>& x) override
    { 
        for(const auto v: m_leftHandleIndices)
            x[v] = Vector3::Zero();
        for(const auto v: m_rightHandleIndices)
            x[v] = Vector3::Zero();       
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
    simulationMesh.m_cellSize = { 10, 10, 10 };
    simulationMesh.m_gridDX = 0.1;
    simulationMesh.m_nFrames = 10;
    simulationMesh.m_subSteps = 1;
    simulationMesh.m_frameDt = 0.04;

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

