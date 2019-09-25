template<class T>
struct ConjugateGradient
{
    template<class SystemType, class VectorType>
    static bool Solve(const SystemType& system, VectorType& x, const VectorType& b,
        VectorType& q, VectorType& p, VectorType& r,
        const T tolerance,const int max_iterations)
    {
        // NOTE: you should never try to make copies of VECTOR_T's inside here as they could be indirect.
        T rhoOld=std::numeric_limits<T>::max();
        T convergence_norm=0;
        r = b;
        system.Multiply(x, q);
        r -= q;
        for(int iterations=0;;iterations++){

            // check for convergence
            convergence_norm = system.Convergence_Norm(r);
            std::cout << convergence_norm << std::endl;
            if(convergence_norm <= tolerance) return true;
            if(iterations == max_iterations) break;

            // actual iteration
            T rho = system.Inner_Product(r,r);
            if(iterations == 0)
                p = r;
            else
                p.Saxpy(rho/rhoOld,p,r);
            system.Multiply(p, q);
            T p_dot_q = system.Inner_Product(p,q);
            if(p_dot_q<=0)
                std::cout << "CG: matrix appears indefinite or singular, p_dot_q/p_dot_p="
                          << p_dot_q / (T)system.Inner_Product(p,p) << std::endl;
            T alpha = rho/p_dot_q;
            x.Saxpy(alpha, p, x);
            r.Saxpy(-alpha, q, r);
            rhoOld = rho;
            
        }

        std::cout<<"cg not converged after "<<max_iterations<<" iterations, error = "<<convergence_norm<<std::endl;
        return false;
    }
};
