#pragma once

#include "AnimatedMesh.h"

#include <Eigen/Dense>

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

    const T m_density;
    const T m_mu;
    const T m_lambda;
    const T m_rayleighCoefficient;

    std::vector<T> m_particleMass;
    std::vector<Vector2> m_particleV;
    std::vector<Matrix22> m_DmInverse;
    std::vector<T> m_restVolume;
    
    FiniteElementMesh(const T density, const T mu, const T lambda, const T rayleighCoefficient)
        :m_density(density), m_mu(mu), m_lambda(lambda), m_rayleighCoefficient(rayleighCoefficient)
    {}

    void initializeUndeformedConfiguration()
    {
        // Initialize rest shape and particle mass (based on constant density)
        m_particleMass.resize(m_particleX.size(), T()); // Initialize all particle masses to zero
        for(const auto& element: m_meshElements)
        {
            Matrix22 Dm;
            for(int j = 0; j < 2; j++)
                Dm.col(j) = m_particleV[element[j+1]]-m_particleV[element[0]];
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
            m_particleV[p] += (dt / m_particleMass[p]) * force[p];
    }

    void simulateFrame(const int frame)
    {
        T stepDt = m_frameDt / (T) m_subSteps;

        for(int step = 1; step <= m_subSteps; step++)
            simulateSubstep(stepDt);
    }
};

