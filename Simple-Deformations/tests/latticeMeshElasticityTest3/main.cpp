
#include "pxr/pxr.h"


#include "pxr/usd/sdf/layer.h"
#include "pxr/usd/sdf/path.h"
#include "pxr/usd/usd/stage.h"
#include "pxr/usd/usdGeom/mesh.h"
#include "pxr/base/vt/array.h"

#include "pxr/base/gf/range3f.h"

#include <iostream>
#include <cmath>

#include <Eigen/Dense>

PXR_NAMESPACE_USING_DIRECTIVE

template<class T, int d> // d is the dimension of the mesh elements, e.g. 3 for triangles, 4 for quads
struct AnimatedMesh
{
    using Vector3 = Eigen::Matrix< T , 3 , 1>;
    
    SdfLayerRefPtr m_layer;
    UsdStageRefPtr m_stage;
    UsdGeomMesh m_mesh;
    UsdAttribute m_pointsAttribute;

    GfRange3f m_extent;
    int m_lastFrame;
    
    std::vector<std::array<int, d>> m_meshElements;
    std::vector<Vector3> m_particleX;

    AnimatedMesh()
        :m_lastFrame(-1)
    {}
    
    void initializeUSD(const std::string filename)
    {
        // Create the layer to populate.
        m_layer = SdfLayer::CreateNew(filename);

        // Create a UsdStage with that root layer.
        m_stage = UsdStage::Open(m_layer);
    }

    void initializeTopology()
    {
        // Create a mesh for this surface
        m_mesh = UsdGeomMesh::Define(m_stage, SdfPath("/MeshSurface"));        

        // Create appropriate buffers for vertex counts and indices, and populate them
        VtIntArray faceVertexCounts, faceVertexIndices;
        for (const auto& element : m_meshElements) {
            faceVertexCounts.push_back(element.size());
            for (const auto& vertex : element)
                faceVertexIndices.push_back(vertex);
        }
        
        // Now set the attributes
        m_mesh.GetFaceVertexCountsAttr().Set(faceVertexCounts);
        m_mesh.GetFaceVertexIndicesAttr().Set(faceVertexIndices);
    }

    void initializeParticles()
    {
        // Grab the points (Positions) attribute, and indicate it is time-varying
        m_pointsAttribute = m_mesh.GetPointsAttr();
        m_pointsAttribute.SetVariability(SdfVariabilityVarying);
    }

    void writeFrame(const int frame)
    {
        std::cout << "Writing frame " << frame << " ..." << std::endl;
        
        // Check that there are any particles to write at all
        if (m_particleX.empty())
            throw std::logic_error("Empty array of input vertices");

        // Check that frames have been written in sequence
        if(frame != m_lastFrame+1)
            throw std::logic_error("Non-consequtive frame sequence requested in writeFrame()");
        m_lastFrame = frame;

        // Update extent
        for (const auto& pt : m_particleX)
            m_extent.UnionWith(GfVec3f(pt[0],pt[1],pt[2]));

        // Copy particleX into VtVec3fArray for Usd
        VtVec3fArray usdPoints(m_particleX.size());
        for(int p = 0; p < m_particleX.size(); p++)
            usdPoints[p] = GfVec3f(m_particleX[p][0],m_particleX[p][1],m_particleX[p][2]);
        
        // Write the points attribute for the given frame
        m_pointsAttribute.Set(usdPoints, (double) frame);
    }
    
    void writeUSD()
    {
        // Set up the timecode
        m_stage->SetStartTimeCode(0.);
        m_stage->SetEndTimeCode((double) m_lastFrame);

        // Set the effective extent
        VtVec3fArray extentArray(2);
        extentArray[0] = m_extent.GetMin();
        extentArray[1] = m_extent.GetMax();
        m_mesh.GetExtentAttr().Set(extentArray);

        // Save USD file
        m_stage->GetRootLayer()->Save();
        std::cout << "USD file saved!" << std::endl;
    }
};

template<class T>
struct LatticeMesh : public AnimatedMesh<T, 4>
{
    using Base = AnimatedMesh<T, 4>;
    using Vector3 = typename Base::Vector3;
    using Base::m_meshElements;
    using Base::m_particleX;
    using Base::initializeUSD;
    using Base::initializeTopology;
    using Base::initializeParticles;

    std::array<int, 2> m_cellSize; // dimensions in grid cells
    T m_gridDX;
    int m_nFrames;
    int m_subSteps;
    T m_frameDt;

    const int m_pinchRadius;
    const T m_particleMass;
    const T m_stiffnessCoeff;
    const T m_dampingCoeff;

    std::vector<Vector3> m_particleV;

    LatticeMesh()
        :m_pinchRadius(1), m_particleMass(1.0), m_stiffnessCoeff(1.0), m_dampingCoeff(0.2)
    {}

    void initialize()
    {
        initializeUSD("latticeMeshElasticityTest3.usda");

        // Create a Cartesian lattice topology
        for(int cell_i = 0; cell_i < m_cellSize[0]; cell_i++)
        for(int cell_j = 0; cell_j < m_cellSize[1]; cell_j++)
            m_meshElements.emplace_back(
                std::array<int, 4>{
                    gridToParticleID(cell_i  , cell_j  ), 
                    gridToParticleID(cell_i+1, cell_j  ),
                    gridToParticleID(cell_i+1, cell_j+1),
                    gridToParticleID(cell_i  , cell_j+1)
                }
            );
        initializeTopology();

        // Also initialize the associated particles

        for(int node_i = 0; node_i <= m_cellSize[0]; node_i++)
        for(int node_j = 0; node_j <= m_cellSize[1]; node_j++)
            if( std::abs(node_i - m_cellSize[0]/2) <= m_pinchRadius &&
                std::abs(node_j - m_cellSize[1]/2) <= m_pinchRadius )
                m_particleX.emplace_back(m_gridDX * (T)node_i + m_gridDX * 0.1 * (T) (m_cellSize[0]+m_cellSize[1]), m_gridDX * (T)node_j + m_gridDX * 0.1 * (T) (m_cellSize[0]+m_cellSize[1]), T());
            else
                m_particleX.emplace_back(m_gridDX * (T)node_i, m_gridDX * (T)node_j, T());
        initializeParticles();

        // Also resize the velocities to match
        m_particleV.resize(m_particleX.size(), Vector3::Zero());
    }

    void prepareFrame(const int frame)
    {
        // Nothing to do; boundary stays static
    }

    void relaxFreeNodes()
    {
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

    void addForce(std::vector<Vector3>& f)
    {
        for(int node_i = 1; node_i < m_cellSize[0]; node_i++)
        for(int node_j = 1; node_j < m_cellSize[1]; node_j++){

            int pCenter = gridToParticleID(node_i  ,node_j  );
            int pPlusX  = gridToParticleID(node_i+1,node_j  );
            int pMinusX = gridToParticleID(node_i-1,node_j  );
            int pPlusY  = gridToParticleID(node_i  ,node_j+1);
            int pMinusY = gridToParticleID(node_i  ,node_j-1);    

            f[pCenter] -= m_stiffnessCoeff * (m_particleX[pCenter] - m_particleX[pPlusX]);
            f[pCenter] -= m_stiffnessCoeff * (m_particleX[pCenter] - m_particleX[pMinusX]);
            f[pCenter] -= m_stiffnessCoeff * (m_particleX[pCenter] - m_particleX[pPlusY]);
            f[pCenter] -= m_stiffnessCoeff * (m_particleX[pCenter] - m_particleX[pMinusY]);

            f[pCenter] -= m_dampingCoeff * (m_particleV[pCenter] - m_particleV[pPlusX]);
            f[pCenter] -= m_dampingCoeff * (m_particleV[pCenter] - m_particleV[pMinusX]);
            f[pCenter] -= m_dampingCoeff * (m_particleV[pCenter] - m_particleV[pPlusY]);
            f[pCenter] -= m_dampingCoeff * (m_particleV[pCenter] - m_particleV[pMinusY]);
        }        
    }

    void simulateSubstep(const T dt)
    {
        const int nParticles = m_particleX.size();
        std::vector<Vector3> force(nParticles, Vector3::Zero());

        addForce(force);
        for(int p = 0; p < nParticles; p++)
            m_particleX[p] += dt * m_particleV[p];
        for(int p = 0; p < nParticles; p++)
            m_particleV[p] += (dt / m_particleMass) * force[p];
    }

    void simulateFrame(const int frame)
    {
        T stepDt = m_frameDt / (T) m_subSteps;

        for(int step = 1; step <= m_subSteps; step++)
            simulateSubstep(stepDt);
    }

private:
    inline int gridToParticleID(const int i, const int j) { return i * (m_cellSize[1]+1) + j; }
};

int main(int argc, char *argv[])
{
    LatticeMesh<float> simulationMesh;
    simulationMesh.m_cellSize = { 40, 40 };
    simulationMesh.m_gridDX = 0.025;
    simulationMesh.m_nFrames = 400;
    simulationMesh.m_subSteps = 10;
    simulationMesh.m_frameDt = 0.1;

    // Initialize the simulation example
    simulationMesh.initialize();
    simulationMesh.relaxFreeNodes();
    
    // Output the initial shape of the mesh
    simulationMesh.writeFrame(0);

    // Perform the animation, output results at each frame
    for(int frame = 1; frame <= simulationMesh.m_nFrames; frame++){
        simulationMesh.prepareFrame(frame);
        simulationMesh.simulateFrame(frame);
        simulationMesh.writeFrame(frame);
    }

    // Write the entire timeline to USD
    simulationMesh.writeUSD();

    return 0;
}

