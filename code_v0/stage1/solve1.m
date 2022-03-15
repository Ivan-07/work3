function [p,w,z] = solve1(N_T,L,M,rho,u,a,x,y,G,z_pre,h,sigma,lambda1,lambda2,lambda3,gamma_c,g,S_0)


solve_choice = 2;
e = 1e-7;
if solve_choice == 1
    cvx_begin quiet
        variable p nonnegative
        variable w(N_T) complex
        variable z complex
        expression sum1(M,L)
        expression sum2(M,L)
        for m=1:M
            u_m_L = u(:,:,m);
            a_m_L = a(:,:,m);
            for l=1:L
                sum1(m,l) = square_abs(x(m,l)-p*u_m_L(:,l)'*a_m_L(:,l)+rho(1)*lambda1(m,l));
                sum2(m,l) = square_abs(y(m,l)-u_m_L(:,l)'*G*w+rho(2)*lambda2(m,l));
            end
        end
        minimize(L*sum_square_abs(w)+L*p^2+1/(2*rho(1))*sum(sum(sum1))+1/(2*rho(2))*sum(sum(sum2))+...
            1/(2*rho(3))*square_abs(z-h'*w+rho(3)*lambda3))
        subject to
            gamma_c*(p^2*sum_square_abs(g'*S_0)+L*sigma)+L*square_abs(z_pre)-2*L*real(conj(z_pre)*z) <= 0;
            
    cvx_end
    w1 = w; z1 = z; p1 = p;
    
elseif solve_choice == 2
    tmp1 = zeros(N_T,N_T);
    tmp2 = zeros(N_T,1);
    tmp3 = 0;
    tmp4 = 0;
    for m = 1:M
        u_m_L = u(:,:,m);
        a_m_L = a(:,:,m);
        for l=1:L
            tmp = u_m_L(:,l)'*G;
            tmp1 = tmp1+tmp'*tmp;
            tmp2 = tmp2+(y(m,l)+rho(2)*lambda2(m,l))*tmp';
            tmp_ = u_m_L(:,l)'*a_m_L(:,l);
            tmp3 = tmp3+abs(tmp_)^2;
            tmp4 = tmp4+real((x(m,l)+rho(1)*lambda1(m,l))*tmp_');
        end
    end


    w = (1/rho(2)*tmp1+2*L*eye(N_T))^(-1)*(1/rho(2)*tmp2);
    z = h'*w-rho(3)*lambda3;
    p = max((2*L+1/rho(1)*tmp3)^(-1)*(1/rho(1)*tmp4),0);


    if solve1_judge(gamma_c,p,g,S_0,sigma,L,z_pre,z) >= e
        lam_lb = 0;lam_ub = 100;
        while 1
            lam = lam_ub;
            p = max((2*L+1/rho(1)*tmp3+lam*2*gamma_c*norm(g'*S_0)^2)^(-1)*(1/rho(1)*tmp4),0);
            w = (1/rho(2)*tmp1+2*L*eye(N_T))^(-1)*(1/rho(2)*tmp2+2*lam*L*z_pre*h);
            z = h'*w-rho(3)*lambda3+2*rho(3)*lam*L*z_pre;
            if solve1_judge(gamma_c,p,g,S_0,sigma,L,z_pre,z) <= -e
                break;
            else
                lam_lb = lam;
                lam_ub = lam*10;
            end
        end
        while 1
            lam = (lam_lb+lam_ub)/2;
            p = max((2*L+1/rho(1)*tmp3+lam*2*gamma_c*norm(g'*S_0)^2)^(-1)*(1/rho(1)*tmp4),0);
            w = (1/rho(2)*tmp1+2*L*eye(N_T))^(-1)*(1/rho(2)*tmp2+2*lam*L*z_pre*h);
            z = h'*w-rho(3)*lambda3+2*rho(3)*lam*L*z_pre;
            judge = solve1_judge(gamma_c,p,g,S_0,sigma,L,z_pre,z);
            if abs(judge) <= e
                break;
            elseif judge >= 0
                lam_lb = lam;
            else
                lam_ub = lam;
            end
        end
    end
    
%     [p p1]
%     [z z1]
%     [w w1]

    
    
end







%end