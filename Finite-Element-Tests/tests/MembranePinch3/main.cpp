
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
    using Vector2 = Eigen::Matrix< T , 2 , 1>;
    
    SdfLayerRefPtr m_layer;
    UsdStageRefPtr m_stage;
    UsdGeomMesh m_mesh;
    UsdAttribute m_pointsAttribute;

    GfRange3f m_extent;
    int m_lastFrame;
    
    std::vector<std::array<int, d>> m_meshElements;
    std::vector<Vector2> m_particleX;

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
            m_extent.UnionWith(GfVec3f(pt[0],pt[1],T()));

        // Copy particleX into VtVec3fArray for Usd
        VtVec3fArray usdPoints(m_particleX.size());
        for(int p = 0; p < m_particleX.size(); p++)
            usdPoints[p] = GfVec3f(m_particleX[p][0],m_particleX[p][1],T());
        
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
struct LatticeMesh : public AnimatedMesh<T, 3>
{
    using Base = AnimatedMesh<T, 3>;
    using Base::m_meshElements;
    using Base::m_particleX;
    using Base::initializeUSD;
    using Base::initializeTopology;
    using Base::initializeParticles;
    using Vector2 = typename Base::Vector2;
    using Matrix22 = Eigen::Matrix< T , 2 , 2>;

    std::array<int, 2> m_cellSize; // dimensions in grid cells
    T m_gridDX;
    int m_nFrames;
    int m_subSteps;
    T m_frameDt;

    const int m_pinchRadius;
    const T m_particleDensity;
    const T m_mu;
    const T m_lambda;
    const T m_rayleighCoefficient;
    
    std::vector<T> m_particleM;
    std::vector<Vector2> m_particleV;
    std::vector<Matrix22> m_DmInverse;
    std::vector<T> m_restVolume;
    
    LatticeMesh()
        :m_pinchRadius(1), m_particleDensity(1.e3), m_mu(1.0), m_lambda(1.0), m_rayleighCoefficient(0.1)
    {}

    void initialize()
    {
        initializeUSD("MembranePinch3.usda");

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

        // Initialize rest shape and particle mass (based on constant density)
        m_particleM.resize(m_particleX.size(), T()); // Initialize all particle masses to zero
        for(const auto& element: m_meshElements)
        {
            Matrix22 Dm;
            for(int j = 0; j < 2; j++)
                Dm.col(j) = m_particleX[element[j+1]]-m_particleX[element[0]];
            T restVolume = .5 * Dm.determinant();
            if(restVolume < 0)
                throw std::logic_error("Inverted element");
            m_DmInverse.emplace_back(Dm.inverse());
            m_restVolume.push_back(restVolume);
            T elementMass = m_particleDensity * restVolume;
            for(const int v: element)
                m_particleM[v] += (1./3.) * elementMass;
        }

    }

    void prepareFrame(const int frame)
    {
        // Nothing to do; boundary stays static
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
                
                if(pCenter < 0 || pCenter >= m_particleX.size()) throw std::logic_error("out of bounds");
                if(pPlusX < 0 || pPlusX >= m_particleX.size()) throw std::logic_error("out of bounds");
                if(pMinusX < 0 || pMinusX >= m_particleX.size()) throw std::logic_error("out of bounds");
                if(pPlusY < 0 || pPlusY >= m_particleX.size()) throw std::logic_error("out of bounds");
                if(pMinusY < 0 || pMinusY >= m_particleX.size()) throw std::logic_error("out of bounds");

                m_particleX[pCenter] = .25 * ( m_particleX[pPlusX] + m_particleX[pMinusX] + m_particleX[pPlusY] + m_particleX[pMinusY]);
            }
    }

    void addForce(std::vector<Vector2>& f) const
    {
        for(int e = 0; e < m_meshElements.size(); e++)
        {
            const auto& element = m_meshElements[e];

            // Linear Elasticity
            Matrix22 Ds;
            for(int j = 0; j < 2; j++)
                Ds.col(j) = m_particleX[element[j+1]]-m_particleX[element[0]];
            Matrix22 F = Ds * m_DmInverse[e];

            Matrix22 strain = .5 * (F + F.transpose()) - Matrix22::Identity();
            Matrix22 P = 2. * m_mu * strain + m_lambda * strain.trace() * Matrix22::Identity();

            Matrix22 H = -m_restVolume[e] * P * m_DmInverse[e].transpose();
            
            for(int j = 0; j < 2; j++){
                f[element[j+1]] += H.col(j);
                f[element[0]] -= H.col(j);
            }

            // Linear Damping
            Matrix22 Ds_dot;
            for(int j = 0; j < 2; j++)
                Ds_dot.col(j) = m_particleV[element[j+1]]-m_particleV[element[0]];
            Matrix22 F_dot = Ds_dot * m_DmInverse[e];

            Matrix22 strain_rate = .5 * (F_dot + F_dot.transpose());
            Matrix22 P_damping = m_rayleighCoefficient * (2. * m_mu * strain_rate + m_lambda * strain_rate.trace() * Matrix22::Identity());

            Matrix22 H_damping = -m_restVolume[e] * P_damping * m_DmInverse[e].transpose();
            
            for(int j = 0; j < 2; j++){
                f[element[j+1]] += H_damping.col(j);
                f[element[0]] -= H_damping.col(j);
            }
        }
    }

    void simulateSubstep(const T dt)
    {
        const int nParticles = m_particleX.size();
        std::vector<Vector2> force(nParticles, Vector2::Zero());

        addForce(force);
        for(int p = 0; p < nParticles; p++)
            m_particleX[p] += dt * m_particleV[p];
        for(int p = 0; p < nParticles; p++)
            m_particleV[p] += (dt / m_particleM[p]) * force[p];
    }

    void simulateFrame(const int frame)
    {
        T stepDt = m_frameDt / (T) m_subSteps;

        for(int step = 1; step <= m_subSteps; step++)
            simulateSubstep(stepDt);
    }

private:
    inline int gridToParticleID(const int i, const int j) const { return i * (m_cellSize[1]+1) + j; }
};

int main(int argc, char *argv[])
{
    LatticeMesh<float> simulationMesh;
    simulationMesh.m_cellSize = { 40, 40 };
    simulationMesh.m_gridDX = 0.025;
    simulationMesh.m_nFrames = 200;
    simulationMesh.m_subSteps = 10;
    simulationMesh.m_frameDt = 0.1;

    // Initialize the simulation example
    simulationMesh.initialize();
    simulationMesh.initializeDeformation();

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

