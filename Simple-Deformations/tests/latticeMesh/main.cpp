
#include "pxr/pxr.h"


#include "pxr/usd/sdf/layer.h"
#include "pxr/usd/sdf/path.h"
#include "pxr/usd/usd/stage.h"
#include "pxr/usd/usdGeom/mesh.h"
#include "pxr/base/vt/array.h"

#include "pxr/base/gf/range3f.h"

#include <iostream>

PXR_NAMESPACE_USING_DIRECTIVE

template<class T, int d> // d is the dimension of the mesh elements, e.g. 3 for triangles, 4 for quads
struct AnimatedMesh
{
    SdfLayerRefPtr m_layer;
    UsdStageRefPtr m_stage;
    UsdGeomMesh m_mesh;
    UsdAttribute m_pointsAttribute;

    GfRange3f m_extent;
    int m_lastFrame;
    
    std::vector<std::array<int, d>> m_meshElements;
    std::vector<GfVec3f> m_particleX;

public:

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
            m_extent.UnionWith(pt);

        // Copy particleX into VtVec3fArray for Usd
        VtVec3fArray usdPoints;
        usdPoints.assign(m_particleX.begin(), m_particleX.end());
        
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

int main(int argc, char *argv[])
{
    AnimatedMesh<float, 4> simulationMesh;

    // Initialize USD data structures
    simulationMesh.initializeUSD("SimulationOutput.usda");

    // Initialize the topology
    simulationMesh.m_meshElements.emplace_back(std::array<int, 4>{0, 1, 3, 2});
    simulationMesh.m_meshElements.emplace_back(std::array<int, 4>{4, 6, 7, 5});
    simulationMesh.m_meshElements.emplace_back(std::array<int, 4>{0, 4, 5, 1});
    simulationMesh.m_meshElements.emplace_back(std::array<int, 4>{2, 3, 7, 6});
    simulationMesh.m_meshElements.emplace_back(std::array<int, 4>{1, 5, 7, 3});
    simulationMesh.m_meshElements.emplace_back(std::array<int, 4>{0, 2, 6, 4});   
    simulationMesh.initializeTopology();

    // Then, initialize the particles
    for(int i = 0; i <= 1; i++)
    for(int j = 0; j <= 1; j++)
    for(int k = 0; k <= 1; k++)
        simulationMesh.m_particleX.emplace_back((float)i,(float)j,(float)k);
    simulationMesh.initializeParticles();

    // Output the initial shape of the mesh
    simulationMesh.writeFrame(0);

    // Perform the animation, output results at each frame
    for(int frame = 1; frame <= 20; frame++){
        for(auto& vertex: simulationMesh.m_particleX)
            vertex += GfVec3f(0.03125,.0625,.125);
        simulationMesh.writeFrame(frame);
    }

    // Write the entire timeline to USD
    simulationMesh.writeUSD();

    return 0;
}

