#pragma once

#include "AnimatedMesh.h"

#include "CGVectorWrapper.h"
#include "CGSystemWrapper.h"
#include "ConjugateGradient.h"

#include <Eigen/Dense>

//#define USE_LINEAR_ELASTICITY
//#define USE_ST_VENANT_KIRCHHOFF
#define USE_COROTATED_ELASTICITY

template<class T>
struct FiniteElementMesh : public AnimatedMesh<T, 3>
{
    using Base = AnimatedMesh<T, 3>;
    using Base::m_meshElements;
    using Base::m_particleX;
    using Vector2 = typename Base::Vector2;
    using Matrix22 = Eigen::Matrix< T , 2 , 2>;

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
    std::vector<Matrix22> m_DmInverse;
    std::vector<T> m_restVolume;
    
    FiniteElementMesh(const T density, const T mu, const T lambda, const T rayleighCoefficient)
        :m_density(density), m_mu(mu), m_lambda(lambda), m_rayleighCoefficient(rayleighCoefficient), m_singularValueThreshold(.2f)
    {}

    void initializeUndeformedConfiguration()
    {
        // Initialize rest shape and particle mass (based on constant density)
        m_particleMass.resize(m_particleX.size(), T()); // Initialize all particle masses to zero
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
            T elementMass = m_density * restVolume;
            for(const int v: element)
                m_particleMass[v] += (1./3.) * elementMass;
        }
    }
    
    void addElasticForce(std::vector<Vector2>& f) const
    {
        for(int e = 0; e < m_meshElements.size(); e++)
        {
            const auto& element = m_meshElements[e];

            // Compute deformation gradient
            Matrix22 Ds;
            for(int j = 0; j < 2; j++)
                Ds.col(j) = m_particleX[element[j+1]]-m_particleX[element[0]];
            Matrix22 F = Ds * m_DmInverse[e];

            // Compute SVD
            Eigen::JacobiSVD<Matrix22> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
            Matrix22 U = svd.matrixU();
            Matrix22 V = svd.matrixV();
            Vector2 vSigma = svd.singularValues();
            if ( U.determinant() < 0. ) {
                if ( V.determinant() < 0. ) {
                    // Both determinants negative, just negate 2nd column on both
                    U.col(1) *= -1.f;
                    V.col(1) *= -1.f;
                }
                else {
                    // Only U has negative determinant, negate 2nd column and second singular value
                    U.col(1) *= -1.f;
                    vSigma[1] = -vSigma[1];
                }
            }
            else
                if ( V.determinant() < 0.) {
                    // Only V has negative determinant, negate 2nd column and second singular value
                    V.col(1) *= -1.f;
                    vSigma[1] = -vSigma[1];
                }
            if ( (F-U*vSigma.asDiagonal()*V.transpose()).norm() > 1e-5 )
                throw std::logic_error("SVD error");

            // Apply thresholding of singular values, and re-constitute F
            for (int v = 0; v < 2; v++)
                vSigma[v] = std::max<T>(m_singularValueThreshold, vSigma[v]);
            Matrix22 Sigma = vSigma.asDiagonal();
            F = U * Sigma * V.transpose();
            
#ifdef USE_LINEAR_ELASTICITY
            Matrix22 strain = .5 * (F + F.transpose()) - Matrix22::Identity();
            Matrix22 P = 2. * m_mu * strain + m_lambda * strain.trace() * Matrix22::Identity();
#endif

#ifdef USE_ST_VENANT_KIRCHHOFF
            Matrix22 E = .5 * ( F.transpose() * F - Matrix22::Identity());
            Matrix22 P = F * (2. * m_mu * E + m_lambda * E.trace() * Matrix22::Identity());
#endif

#ifdef USE_COROTATED_ELASTICITY
            Vector2 vStrain = vSigma - Vector2::Ones();
            Vector2 vP = 2. * m_mu * vStrain + m_lambda * vStrain.sum() * Vector2::Ones();
            Matrix22 P = U * vP.asDiagonal() * V.transpose();
#endif

            Matrix22 H = -m_restVolume[e] * P * m_DmInverse[e].transpose();
            
            for(int j = 0; j < 2; j++){
                f[element[j+1]] += H.col(j);
                f[element[0]] -= H.col(j);
            }
        }
    }

    void addProductWithStiffnessMatrix(std::vector<Vector2>& dx, std::vector<Vector2>& df, const T scale) const
    {
        for(int e = 0; e < m_meshElements.size(); e++)
        {
            const auto& element = m_meshElements[e];

            // Compute deformation gradient
            Matrix22 Ds;
            for(int j = 0; j < 2; j++)
                Ds.col(j) = m_particleX[element[j+1]]-m_particleX[element[0]];
            Matrix22 F = Ds * m_DmInverse[e];

            // Compute SVD
            Eigen::JacobiSVD<Matrix22> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
            Matrix22 U = svd.matrixU();
            Matrix22 V = svd.matrixV();
            Vector2 vSigma = svd.singularValues();
            if ( U.determinant() < 0. ) {
                if ( V.determinant() < 0. ) {
                    // Both determinants negative, just negate 2nd column on both
                    U.col(1) *= -1.f;
                    V.col(1) *= -1.f;
                }
                else {
                    // Only U has negative determinant, negate 2nd column and second singular value
                    U.col(1) *= -1.f;
                    vSigma[1] = -vSigma[1];
                }
            }
            else
                if ( V.determinant() < 0.) {
                    // Only V has negative determinant, negate 2nd column and second singular value
                    V.col(1) *= -1.f;
                    vSigma[1] = -vSigma[1];
                }
            if ( (F-U*vSigma.asDiagonal()*V.transpose()).norm() > 1e-5 )
                throw std::logic_error("SVD error");

            // Apply thresholding of singular values, and re-constitute F
            for (int v = 0; v < 2; v++)
                vSigma[v] = std::max<T>(m_singularValueThreshold, vSigma[v]);
            Matrix22 Sigma = vSigma.asDiagonal();
            F = U * Sigma * V.transpose();
            
            // Compute differential(s)
            Matrix22 dDs;
            for(int j = 0; j < 2; j++)
                dDs.col(j) = dx[element[j+1]]-dx[element[0]];
            Matrix22 dF = dDs * m_DmInverse[e];

#ifdef USE_LINEAR_ELASTICITY
            Matrix22 dstrain = .5 * (dF + dF.transpose());
            Matrix22 dP = scale * (2. * m_mu * dstrain + m_lambda * dstrain.trace() * Matrix22::Identity());
#endif

#ifdef USE_ST_VENANT_KIRCHHOFF
            Matrix22 E = .5 * ( F.transpose() * F - Matrix22::Identity());
            Matrix22 dE = .5 * ( dF.transpose() * F + F.transpose() * dF);
            Matrix22 dP = dF * (2. * m_mu *  E + m_lambda *  E.trace() * Matrix22::Identity()) +
                           F * (2. * m_mu * dE + m_lambda * dE.trace() * Matrix22::Identity());

#endif

#ifdef USE_COROTATED_ELASTICITY
            // Construct diagonalized dP/dF tensor
            Matrix22 A, B12;
            A(0, 0) = A(1, 1) = 2. * m_mu + m_lambda;
            A(0, 1) = A(1, 0) = m_lambda;
            Vector2 vStrain = vSigma - Vector2::Ones();
            T q = std::max<T>( m_lambda * vStrain.sum() - 2. * m_mu, 0. ); // Positive definiteness fix
            B12(0, 0) = B12(1, 1) = m_mu + q;
            B12(0, 1) = B12(1, 0) = m_mu - q;

            // Apply tensor (with required rotations)
            Matrix22 dF_hat = U.transpose() * dF * V;

            Vector2 vdP_A = A * Vector2( dF_hat(0, 0), dF_hat(1, 1) );
            Vector2 vdP_B12 = B12 * Vector2( dF_hat(0, 1), dF_hat(1, 0) );
            Matrix22 dP_hat;
            dP_hat(0, 0) = vdP_A[0];
            dP_hat(1, 1) = vdP_A[1];
            dP_hat(0, 1) = vdP_B12[0];
            dP_hat(1, 0) = vdP_B12[1];

            Matrix22 dP = U * dP_hat * V.transpose();
#endif
            
            Matrix22 dH = m_restVolume[e] * dP * m_DmInverse[e].transpose();
            
            for(int j = 0; j < 2; j++){
                df[element[j+1]] += dH.col(j);
                df[element[0]] -= dH.col(j);
            }
        }
    }

    void simulateSubstep()
    {
        using FEMType = FiniteElementMesh<T>;        

        const int nParticles = m_particleX.size();

        // Construct initial guess for next-timestep
        //   Positions -> Same as last timestep
        
        // Overwrite boundary conditions with desired values

        setBoundaryConditions();
        
        // Solve for everything else using Conjugate Gradients

        std::vector<Vector2> dx(nParticles, Vector2::Zero());
        std::vector<Vector2> rhs(nParticles, Vector2::Zero());
        std::vector<Vector2> q(nParticles, Vector2::Zero());
        std::vector<Vector2> s(nParticles, Vector2::Zero());
        std::vector<Vector2> r(nParticles, Vector2::Zero());
        CGVectorWrapper<Vector2> dxWrapper(dx);
        CGVectorWrapper<Vector2> rhsWrapper(rhs);
        CGVectorWrapper<Vector2> qWrapper(q);
        CGVectorWrapper<Vector2> sWrapper(s);
        CGVectorWrapper<Vector2> rWrapper(r);
        CGSystemWrapper<Vector2, FEMType> systemWrapper(*this);
        
        addElasticForce(rhs);
        clearConstrainedParticles(rhs);

        ConjugateGradient<T>::Solve(systemWrapper,
            dxWrapper, rhsWrapper, qWrapper, sWrapper, rWrapper,
            1e-4, 50);

        // Apply corrections to positions and velocities

        for(int p = 0; p < nParticles; p++)
            m_particleX[p] += dx[p];
    }

    void simulateFrame(const int frame)
    {
        m_stepDt = m_frameDt / (T) m_subSteps;

        for(int step = 1; step <= m_subSteps; step++){
            m_stepEndTime = m_frameDt * (T) (frame-1) + m_stepDt * (T) step;
            simulateSubstep();
        }
    }

    virtual void clearConstrainedParticles(std::vector<Vector2>& x) {}
    virtual void setBoundaryConditions() {}
};

