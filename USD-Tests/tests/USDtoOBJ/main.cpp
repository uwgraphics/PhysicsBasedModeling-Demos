#include "pxr/usd/usdGeom/mesh.h"

#include <iostream>
#include <iomanip>
#include <fstream>

PXR_NAMESPACE_USING_DIRECTIVE

int main(int argc, char *argv[])
{
    UsdStageRefPtr stage = UsdStage::Open("SimulationOutput.usda");
    UsdGeomMesh mesh = UsdGeomMesh::Get(stage, SdfPath("/TriangulatedSurface0"));
    VtIntArray faceVertexIndices;
    mesh.GetFaceVertexIndicesAttr().Get(&faceVertexIndices);
    if(faceVertexIndices.size() % 3) throw std::runtime_error("Vertex indices not a multiple of 3");
    int nTriangles = faceVertexIndices.size() / 3;
    std::cout << "Found " << nTriangles << " triangles" << std::endl;
    UsdAttribute pointsAttribute = mesh.GetPointsAttr();
    int nFrames = pointsAttribute.GetNumTimeSamples(); 
    std::cout << "Found " << nFrames << " frames" << std::endl;

    for(int frame = 0; frame < nFrames; frame++){
        std::stringstream sstr;
        sstr << "frame_" << std::setfill('0') << std::setw(5) << frame << ".obj";
        std::string filename = sstr.rdbuf()->str();
        std::cout << "Writing output file " << filename << " ..." << std::endl;
        std::ofstream ofile(filename);

        VtVec3fArray usdPoints;
        pointsAttribute.Get(&usdPoints, frame);

        for(int v = 0; v < usdPoints.size(); v++)
            ofile << "v " << usdPoints[v][0] << " " << usdPoints[v][1] << " " << usdPoints[v][2] << std::endl;
        for(int t = 0; t < nTriangles; t++)
            ofile << "f " << faceVertexIndices[3*t]+1 << " " << faceVertexIndices[3*t+1]+1 << " " << faceVertexIndices[3*t+2]+1 << std::endl;

        ofile.close();
    }
    
    return 0;
}

