template<class T>
struct ConjugateGradient
{
    template<class SystemType, class VectorType>
    static bool Solve(const SystemType& system, VectorType& x, const VectorType& b,
        VectorType& q, VectorType& s, VectorType& r,
        const T tolerance,const int max_iterations)
    {
        // NOTE: you should never try to make copies of VECTOR_T's inside here as they could be indirect.
        static const T small_number=std::numeric_limits<T>::epsilon();
        T rho_old=(T)FLT_MAX;
        T convergence_norm=0;
        int iterations;
        for(iterations=0;;iterations++){
            bool restart=!iterations;
            if(restart){
                std::cout<<"restarting cg"<<std::endl;
                r=b;system.Multiply(x,q);r-=q;}
            // stopping conditions
            convergence_norm=system.Convergence_Norm(r);
            std::cout<<convergence_norm<<std::endl;
            if(convergence_norm<=tolerance) return true;
            if(iterations==max_iterations) break;
            // actual iteration
            T rho=(T)system.Inner_Product(r,r);
            if(restart) s=r;
            else s.Copy(rho/rho_old,s,r);
            system.Multiply(s,q);
            T s_dot_q=(T)system.Inner_Product(s,q);
            if(s_dot_q<=0) std::cout<<"CG: matrix appears indefinite or singular, s_dot_q/s_dot_s="<<s_dot_q/(T)system.Inner_Product(s,s)<<std::endl;
            T alpha=s_dot_q?rho/s_dot_q:(T)FLT_MAX;
            x.Copy(alpha,s,x);
            r.Copy(-alpha,q,r);
            rho_old=rho;
        }

        std::cout<<"cg not converged after "<<max_iterations<<" iterations, error = "<<convergence_norm<<std::endl;
        return false;
    }
};
