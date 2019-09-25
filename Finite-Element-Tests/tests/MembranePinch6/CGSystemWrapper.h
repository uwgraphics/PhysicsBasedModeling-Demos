template<class AttributeType, // AttributeType is the type of each particle's property,
                              // e.g., Eigen::Vector2f for typical positions/velocities
         class FEMType>       // FEMType is the class that defines forces (e.g. FiniteElementMesh)
struct CGSystemWrapper
{
    using T = typename AttributeType::Scalar;                  // Underlying scalar type, e.g. float
    using AggregateType = typename std::vector<AttributeType>; // Aggregate of all particles' properties,
                                                               // this is the "unwrapped" version of this object
    using VectorWrapper = CGVectorWrapper<AttributeType>;

    FEMType& m_finiteElements;

    CGSystemWrapper(FEMType& finiteElements)
        :m_finiteElements(finiteElements)
    {}

    void Multiply(const VectorWrapper& x, VectorWrapper& y) const
    {
        const T oneOverDtSquared = T(1.) / (m_finiteElements.m_stepDt * m_finiteElements.m_stepDt);

        for(int i = 0; i < x.m_data.size(); i++)
            y.m_data[i] = m_finiteElements.m_particleMass[i] * oneOverDtSquared * x.m_data[i];
        m_finiteElements.addProductWithStiffnessMatrix(
            x.m_data,
            y.m_data,
            T(1.) + m_finiteElements.m_rayleighCoefficient / m_finiteElements.m_stepDt
        );
    }

    T Convergence_Norm(const VectorWrapper& x) const
    {
        double result = 0.;
        for(const auto& v: x.m_data)
            result = std::max<double>(result, v.squaredNorm());
        return (T) std::sqrt(result);
    }

    T Inner_Product(const VectorWrapper& x, const VectorWrapper& y) const
    {
        double result = 0.;
        for(int i = 0; i < x.m_data.size(); i++)
            result += (double) x.m_data[i].dot(y.m_data[i]);
        return (T) result;
    }
};
