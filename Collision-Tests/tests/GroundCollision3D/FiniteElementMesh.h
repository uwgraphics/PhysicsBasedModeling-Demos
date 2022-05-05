#pragma once

#include "AnimatedTetrahedronMesh.h"

#include "CGVectorWrapper.h"
#include "CGSystemWrapper.h"
#include "ConjugateGradient.h"

#include <Eigen/Dense>

#define USE_LINEAR_ELASTICITY
//#define USE_COROTATED_ELASTICITY

template<class T>
struct FiniteElementMesh : public AnimatedTetrahedronMesh<T>
{
    using Base = AnimatedTetrahedronMesh<T>;
    using Base::m_meshElements;
    using Base::m_particleX;
    using Vector3 = typename Base::Vector3;
    using Matrix33 = Eigen::Matrix< T , 3 , 3>;

    int m_nFrames;
    int m_subSteps;
    T m_frameDt;
    T m_stepDt;
    T m_stepEndTime;

    const T m_density;
    const T m_mu;
    const T m_lambda;
    const T m_rayleighCoefficient;
    const T m_singularValueThreshold;
    
    std::vector<T> m_particleMass;
    std::vector<Vector3> m_particleV;
    std::vector<Matrix33> m_DmInverse;
    std::vector<T> m_restVolume;

    mutable std::vector<Matrix33> Ue,Ve;
    mutable std::vector<Vector3> vSigmae;

    std::vector<int> m_surfaceParticles;
    mutable std::vector<bool> m_collisionActive;
    mutable std::vector<Vector3> m_collisionTarget;
    T m_collisionStiffness;

    FiniteElementMesh(const T density, const T mu, const T lambda, const T rayleighCoefficient)
        :m_density(density), m_mu(mu), m_lambda(lambda), m_rayleighCoefficient(rayleighCoefficient), m_singularValueThreshold(-std::numeric_limits<T>::max()),
        m_collisionStiffness(2.e1)
    {}

    void initializeUndeformedConfiguration()
    {
        // Initialize rest shape and particle mass (based on constant density)
        m_particleMass.resize(m_particleX.size(), T()); // Initialize all particle masses to zero
        for(const auto& element: m_meshElements)
        {
            Matrix33 Dm;
            for(int j = 0; j < 3; j++)
                Dm.col(j) = m_particleX[element[j+1]]-m_particleX[element[0]];
            T restVolume = (.1/.6) * Dm.determinant();
            if(restVolume < 0)
                throw std::logic_error("Inverted element");
            m_DmInverse.emplace_back(Dm.inverse());
            m_restVolume.push_back(restVolume);
            T elementMass = m_density * restVolume;
            for(const int v: element)
                m_particleMass[v] += (1./4.) * elementMass;
        }
    }
    
    void addElasticForce(std::vector<Vector3>& f) const
    {
        Ue.resize(m_meshElements.size());
        Ve.resize(m_meshElements.size());
        vSigmae.resize(m_meshElements.size());

#pragma omp parallel for
        for(int e = 0; e < m_meshElements.size(); e++)
        {
            const auto& element = m_meshElements[e];

            // Compute deformation gradient
            Matrix33 Ds;
            for(int j = 0; j < 3; j++)
                Ds.col(j) = m_particleX[element[j+1]]-m_particleX[element[0]];
            Matrix33 F = Ds * m_DmInverse[e];

            // Compute SVD
            Eigen::JacobiSVD<Matrix33> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
            Matrix33& U = Ue[e];
            Matrix33& V = Ve[e];
            Vector3& vSigma = vSigmae[e];
            U = svd.matrixU();
            V = svd.matrixV();
            vSigma = svd.singularValues();
            if ( U.determinant() < 0. ) {
                if ( V.determinant() < 0. ) {
                    // Both determinants negative, just negate last column on both
                    U.col(2) *= -1.f;
                    V.col(2) *= -1.f;
                }
                else {
                    // Only U has negative determinant, negate last column and second singular value
                    U.col(2) *= -1.f;
                    vSigma[2] = -vSigma[2];
                }
            }
            else
                if ( V.determinant() < 0.) {
                    // Only V has negative determinant, negate last column and second singular value
                    V.col(2) *= -1.f;
                    vSigma[2] = -vSigma[2];
                }

            // Apply thresholding of singular values, and re-constitute F
            for (int v = 0; v < 3; v++)
                vSigma[v] = std::max<T>(m_singularValueThreshold, vSigma[v]);
            Matrix33 Sigma = vSigma.asDiagonal();
        }

        for(int e = 0; e < m_meshElements.size(); e++)
        {
            const auto& element = m_meshElements[e];

            // Use precomputed SVD
            const Matrix33& U = Ue[e];
            const Matrix33& V = Ue[e];
            const Vector3& vSigma = vSigmae[e];
            const Matrix33 F = U * vSigma.asDiagonal() * V.transpose();
            
#ifdef USE_LINEAR_ELASTICITY
            Matrix33 strain = .5 * (F + F.transpose()) - Matrix33::Identity();
            Matrix33 P = 2. * m_mu * strain + m_lambda * strain.trace() * Matrix33::Identity();
#endif


#ifdef USE_COROTATED_ELASTICITY
            Vector3 vStrain = vSigma - Vector3::Ones();
            Vector3 vP = 2. * m_mu * vStrain + m_lambda * vStrain.sum() * Vector3::Ones();
            Matrix33 P = U * vP.asDiagonal() * V.transpose();
#endif

            Matrix33 H = -m_restVolume[e] * P * m_DmInverse[e].transpose();
            
            for(int j = 0; j < 3; j++){
                f[element[j+1]] += H.col(j);
                f[element[0]] -= H.col(j);
            }
        }

        // Detect collisions

        const T groundHeight = -.3;
        
        m_collisionActive.resize( m_particleX.size() );
        for(auto b : m_collisionActive) b = false;
        m_collisionTarget.resize( m_particleX.size() );


        for (const int p: m_surfaceParticles) {
            T phi = m_particleX[p][1] - groundHeight;
            if (phi < 0) {
                m_collisionActive[p] = true;
                m_collisionTarget[p] = Vector3( m_particleX[p][0], groundHeight, m_particleX[p][2] );
            }
            else m_collisionActive[p] = false;
        }

        // Apply collision force

        for (int p = 0; p < m_particleX.size(); p++)
            if ( m_collisionActive[p] ) {
                Vector3 collisionDX = m_particleX[p] - m_collisionTarget[p];
                f[p] -= m_collisionStiffness * collisionDX;
            }
    }

#if 0
    struct DiagonalizedStressDerivative
    {
        Matrix22 A;
        Matrix22 B12;
    };
#endif

    void addProductWithStiffnessMatrix(std::vector<Vector3>& dx, std::vector<Vector3>& df, const T scale) const
    {
        for(int e = 0; e < m_meshElements.size(); e++)
        {
            const auto& element = m_meshElements[e];

            // Use precomputed SVD
            const Matrix33& U = Ue[e];
            const Matrix33& V = Ue[e];
            const Vector3& vSigma = vSigmae[e];
            
            // Compute differential(s)
            Matrix33 dDs;
            for(int j = 0; j < 3; j++)
                dDs.col(j) = dx[element[j+1]]-dx[element[0]];
            Matrix33 dF = dDs * m_DmInverse[e];

#ifdef USE_LINEAR_ELASTICITY
            Matrix33 dstrain = .5 * (dF + dF.transpose());
            Matrix33 dP = 2. * m_mu * dstrain + m_lambda * dstrain.trace() * Matrix33::Identity();
#endif

#ifdef USE_COROTATED_ELASTICITY
            // Implementing stiffness via projective dynamics instead!
            Matrix33 dP = 2. * m_mu * dF;
#endif

            Matrix33 dH = (scale * m_restVolume[e]) * dP * m_DmInverse[e].transpose();
            
            for(int j = 0; j < 3; j++){
                df[element[j+1]] += dH.col(j);
                df[element[0]] -= dH.col(j);
            }
        }

        // Apply collision force differential

        for (int p = 0; p < m_particleX.size(); p++)
            if ( m_collisionActive[p] ) {
                Vector3 collisionDX = m_particleX[p] - m_collisionTarget[p];
                df[p] += scale * m_collisionStiffness * dx[p];
            }
    }

    void simulateSubstep()
    {
        using FEMType = FiniteElementMesh<T>;        

        const int nParticles = m_particleX.size();

        // Save last time step velocities

        std::vector<Vector3> lastV = m_particleV;
        
        // Construct initial guess for next-timestep
        //   Velocities -> Same as last timestep
        //   Positions -> Using Forward Euler
        
        for(int p = 0; p < nParticles; p++)
            m_particleX[p] += m_stepDt * m_particleV[p];

        // Overwrite boundary conditions with desired values

        setBoundaryConditions();
        
        // Solve for everything else using Conjugate Gradients

        std::vector<Vector3> dx(nParticles, Vector3::Zero());
        std::vector<Vector3> rhs(nParticles, Vector3::Zero());
        std::vector<Vector3> q(nParticles, Vector3::Zero());
        std::vector<Vector3> s(nParticles, Vector3::Zero());
        std::vector<Vector3> r(nParticles, Vector3::Zero());
        CGVectorWrapper<Vector3> dxWrapper(dx);
        CGVectorWrapper<Vector3> rhsWrapper(rhs);
        CGVectorWrapper<Vector3> qWrapper(q);
        CGVectorWrapper<Vector3> sWrapper(s);
        CGVectorWrapper<Vector3> rWrapper(r);
        CGSystemWrapper<Vector3, FEMType> systemWrapper(*this);
        
        addElasticForce(rhs);
        addExternalForce(rhs);
        for(int p = 0; p < nParticles; p++)
            rhs[p] += (m_particleMass[p] / m_stepDt) * (lastV[p] - m_particleV[p]);
        addProductWithStiffnessMatrix(m_particleV, rhs, -m_rayleighCoefficient);
        clearConstrainedParticles(rhs);

        ConjugateGradient<T>::Solve(systemWrapper,
            dxWrapper, rhsWrapper, qWrapper, sWrapper, rWrapper,
            1e-3, 50);

        // Apply corrections to positions and velocities

        const T oneOverDt = T(1.) / m_stepDt;
        for(int p = 0; p < nParticles; p++){
            m_particleX[p] += dx[p];
            m_particleV[p] += oneOverDt * dx[p];
        }
    }

    void simulateFrame(const int frame)
    {
        m_stepDt = m_frameDt / (T) m_subSteps;

        for(int step = 1; step <= m_subSteps; step++){
            m_stepEndTime = m_frameDt * (T) (frame-1) + m_stepDt * (T) step;
            simulateSubstep();
        }
    }

    virtual void addExternalForce(std::vector<Vector3>& f) {}
    virtual void clearConstrainedParticles(std::vector<Vector3>& x) {}
    virtual void setBoundaryConditions() {}
};

