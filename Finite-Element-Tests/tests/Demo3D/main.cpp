#include <iostream>

#include "AnimatedTetrahedronMesh.h"

#include <Eigen/Dense>

#include <map>

template<class T>
struct LatticeMesh : public AnimatedTetrahedonMesh<T>
{
    using Base = AnimatedTetrahedonMesh<T>;

    // from AnimatedTetrahedonMesh
    using Base::m_meshElements;
    using Base::m_particleX;
    using Base::initializeUSD;
    using Base::initializeTopology;
    using Base::initializeParticles;
    using Vector3 = typename Base::Vector3;

    std::array<int, 3> m_cellSize; // dimensions in grid cells
    int m_radius; // radius of sphere in grid cells
    T m_gridDX;

    std::vector<std::array<int, 3>> m_activeCells; // Marks the "active" cells in the lattice
    std::map<std::array<int, 3>, int> m_activeNodes; // Maps the "active" nodes to their particle index

    void initialize()
    {
        initializeUSD("Demo3D.usda");

        // Activate cells within a sphere of radius m_radius (in cells)

        for(int cell_i = 0; cell_i < m_cellSize[0]; cell_i++)
        for(int cell_j = 0; cell_j < m_cellSize[1]; cell_j++)
        for(int cell_k = 0; cell_k < m_cellSize[1]; cell_k++){

            int r = (cell_i - m_cellSize[0]/2) * (cell_i - m_cellSize[0]/2) +
                    (cell_j - m_cellSize[1]/2) * (cell_j - m_cellSize[1]/2) +
                    (cell_k - m_cellSize[2]/2) * (cell_k - m_cellSize[2]/2);

            if(r <= m_radius * m_radius)
                m_activeCells.push_back(std::array<int, 3>{cell_i, cell_j, cell_k});

        }

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

        // Check particle indexing in mesh

        for(const auto& element: m_meshElements)
            for(const auto vertex: element)
                if(vertex < 0 || vertex >= m_particleX.size())
                    throw std::logic_error("mismatch between mesh vertex and particle array");
    }
};

int main(int argc, char *argv[])
{
    LatticeMesh<float> simulationMesh;
    simulationMesh.m_cellSize = { 20, 20, 20 };
    simulationMesh.m_radius = 8;
    simulationMesh.m_gridDX = 0.1;

    // Initialize the simulation example
    simulationMesh.initialize();

    // Output the initial shape of the mesh
    simulationMesh.writeFrame(0);

    // Write the entire timeline to USD
    simulationMesh.writeUSD();

    return 0;
}

