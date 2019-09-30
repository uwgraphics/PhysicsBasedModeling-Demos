template<class AttributeType> // AttributeType is the type of each particle's property,
                              // e.g., Eigen::Vector2f for typical positions/velocities
struct CGVectorWrapper
{
    using T = typename AttributeType::Scalar;                  // Underlying scalar type, e.g. float
    using AggregateType = typename std::vector<AttributeType>; // Aggregate of all particles' properties,
                                                               // this is the "unwrapped" version of this object
    AggregateType& m_data;

    CGVectorWrapper(AggregateType& data)
        :m_data(data)
    {}

    CGVectorWrapper& operator = (const CGVectorWrapper& v)
    {
        m_data = v.m_data;
        return *this;
    }

    CGVectorWrapper& operator -= (const CGVectorWrapper& v)
    {
        for(int i = 0; i < m_data.size(); i++)
            m_data[i] -= v.m_data[i];
        return *this;
    }

    // replaces current vector with c * x + y
    void Saxpy(const T c, const CGVectorWrapper& x, const CGVectorWrapper& y)
    {
        for(int i = 0; i < m_data.size(); i++)
            m_data[i] = c * x.m_data[i] + y.m_data[i];
    }
};
