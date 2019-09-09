
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
struct LatticeMesh
{
    SdfLayerRefPtr m_layer;
    UsdStageRefPtr m_stage;
    UsdGeomMesh m_mesh;
    GfRange3f m_extent;
    int m_nFrames;
    
    std::vector<std::array<int, d>> m_meshElements;

public:

    LatticeMesh()
        :m_nFrames(-1)
    {}
    
    void initializeUSD()
    {
        // Create the layer to populate.
        m_layer = SdfLayer::CreateNew("latticeMesh.usda");

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

    void writeUSD()
    {
        // Set up the timecode
        m_stage->SetStartTimeCode(0.);
        m_stage->SetEndTimeCode((double) m_nFrames);

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
    LatticeMesh<float, 4> latticeMesh;

    // Initialize USD data structures
    latticeMesh.initializeUSD();

    // Initialize the topology
    latticeMesh.m_meshElements.emplace_back(std::array<int, 4>{0, 1, 3, 2});
    latticeMesh.m_meshElements.emplace_back(std::array<int, 4>{4, 6, 7, 5});
    latticeMesh.m_meshElements.emplace_back(std::array<int, 4>{0, 4, 5, 1});
    latticeMesh.m_meshElements.emplace_back(std::array<int, 4>{2, 3, 7, 6});
    latticeMesh.m_meshElements.emplace_back(std::array<int, 4>{1, 5, 7, 3});
    latticeMesh.m_meshElements.emplace_back(std::array<int, 4>{0, 2, 6, 4});   
    latticeMesh.initializeTopology();

    // Now we'll populate the stage with content from the objStream.
    std::vector<GfVec3f> particleX;
    for(int i = 0; i <= 1; i++)
    for(int j = 0; j <= 1; j++)
    for(int k = 0; k <= 1; k++)
        particleX.emplace_back((float)i,(float)j,(float)k);

    std::cout << "particleX.size() = " << particleX.size() << std::endl;

    if (particleX.empty())
        throw std::logic_error("Empty array of input vertices");

    // Copy particleX into VtVec3fArray for Usd.
    VtVec3fArray usdPoints;
    usdPoints.assign(particleX.begin(), particleX.end());

    // Usd currently requires an extent, somewhat unfortunately.
    latticeMesh.m_nFrames = 20;
    for (const auto& pt : usdPoints) {
        latticeMesh.m_extent.UnionWith(pt);
    }
    // VtVec3fArray extentArray(2);
    // extentArray[0] = latticeMesh.m_extent.GetMin();
    // extentArray[1] = latticeMesh.m_extent.GetMax() + GfVec3f(0.03125,.0625,.125) * (float) latticeMesh.m_nFrames;

    // Populate the mesh vertex data
    UsdAttribute pointsAttribute = latticeMesh.m_mesh.GetPointsAttr();
    pointsAttribute.SetVariability(SdfVariabilityVarying);

    std::vector<double> timeSamples;
    timeSamples.push_back(0.);
    for(int frame = 1; frame <= latticeMesh.m_nFrames; frame++)
        timeSamples.push_back((double)frame);
    for(const auto sample: timeSamples){
        pointsAttribute.Set(usdPoints, sample);
        for(auto& vertex: particleX){
            vertex += GfVec3f(0.03125,.0625,.125);
            latticeMesh.m_extent.UnionWith(vertex);
        }
        usdPoints.assign(particleX.begin(), particleX.end());
    }

    // Set extent.
    // latticeMesh.m_mesh.GetExtentAttr().Set(extentArray);

    latticeMesh.writeUSD();

    return 0;
}

