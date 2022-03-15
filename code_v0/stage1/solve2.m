function [u,x,y] = solve2(M_R,L,M,gamma_r,sigma,G,w,x_pre,a,rho,lambda2,p,lambda1,gamma_c)

solve_choice = 1;
e = 1e-8;
objValue1 = zeros(M,1);
objValue2 = zeros(M,1);
if solve_choice == 1
    u = zeros(M_R,L,M);
    x = zeros(M,L);
    y = zeros(M,L);
    for m=1:M
        %m
        a_m_L = a(:,:,m);
        
        cvx_begin quiet
        variable u_m(M_R,L) complex
        variable x_m(L) complex
        variable y_m(L) complex
        dual variable mu_dual
        expression sum1(L)
        expression sum2(L)
        expression sum3(L)
        expression sum4(L)
        for l=1:L
            sum1(l) = square_abs(x_m(l)-p*u_m(:,l)'*a_m_L(:,l)+rho(1)*lambda1(m,l));
            sum2(l) = square_abs(y_m(l)-u_m(:,l)'*G*w+rho(2)*lambda2(m,l));
            sum3(l) = sum_square_abs(u_m(:,l));
            sum4(l) = 2*real(conj(x_pre(m,l))*x_m(l));
        end
        minimize(sum(sum(sum1))+sum(sum(sum2)))
        subject to
            gamma_r*(sigma*sum(sum3)+sum_square_abs(y_m))+sum_square_abs(x_pre(m,:))-sum(sum4) <= 0;
            for i = 1:L
                u_m(1,l) == 1;
            end
        cvx_end
        
        u(:,:,m) = u_m;
        x(m,:) = x_m.';
        y(m,:) = y_m.';
        objValue1(m) = solve2_obj(L,a_m_L,rho,u_m,y_m,lambda1,x_m,lambda2,G,w,p,m);
        %[m mu_dual]
    end
    %u1 = u; x1 = x; y1 = y;
    
elseif solve_choice == 2
    u = zeros(M_R,L,M);
    x = zeros(M,L);
    y = zeros(M,L);
    for m=1:M
        u_m = zeros(M_R,L);
        x_m = zeros(1,L);
        y_m = zeros(1,L);
        a_m_L = a(:,:,m);
        x_pre_m = x_pre(m,:);
        
        for l=1:L
            u_m(:,l) = (gamma_r*sigma*eye(M_R)+gamma_r*G*w*w'*G')^(-1)*...
                (p*x_pre(m,l)'*a_m_L(:,l)+gamma_r*rho(2)*lambda2(m,l)*G*w);
            x_m(l) = p*u_m(:,l)'*a_m_L(:,l)-rho(1)*lambda1(m,l);
            y_m(l) = u_m(:,l)'*G*w-rho(2)*lambda2(m,l);
        end
        mu = 0;
        if solve2_judge(L,gamma_r,sigma,u_m,y_m,x_pre_m,x_m) >= e
            mu_lb = 0;mu_ub = 100;
            while 1
                mu = mu_ub;
                for l=1:L
                    u_m(:,l) = (mu*gamma_r*sigma*eye(M_R)+(mu*gamma_r/(1+mu*gamma_r))*G*w*w'*G')^(-1)*...
                        ((mu*x_pre_m(l))'*p*a_m_L(:,l)+mu*gamma_r/(1+mu*gamma_r)*rho(2)*lambda2(m,l)*G*w);
                    x_m(l) = p*u_m(:,l)'*a_m_L(:,l)-rho(1)*lambda1(m,l)+mu*x_pre_m(l);
                    y_m(l) = (u_m(:,l)'*G*w-rho(2)*lambda2(m,l))/(1+mu*gamma_r);
                end
                if solve2_judge(L,gamma_r,sigma,u_m,y_m,x_pre_m,x_m) <= -e
                    break;
                else
                    mu_lb = mu;
                    mu_ub = mu*10;
                end
            end
            
            while 1
                mu = (mu_ub+mu_lb)/2;
                for l=1:L
                    u_m(:,l) = (mu*gamma_r*sigma*eye(M_R)+(mu*gamma_r/(1+mu*gamma_r))*G*w*w'*G')^(-1)*...
                        ((mu*x_pre_m(l))'*p*a_m_L(:,l)+mu*gamma_r/(1+mu*gamma_r)*rho(2)*lambda2(m,l)*G*w);
                    x_m(l) = p*u_m(:,l)'*a_m_L(:,l)-rho(1)*lambda1(m,l)+mu*x_pre_m(l);
                    y_m(l) = (u_m(:,l)'*G*w-rho(2)*lambda2(m,l))/(1+mu*gamma_r);
                end
                judge = solve2_judge(L,gamma_r,sigma,u_m,y_m,x_pre_m,x_m);
                if abs(judge) <= e
                    break;
                elseif judge >= 0
                    mu_lb = mu;
                else
                    mu_ub = mu;
                end
            end
        end
        u(:,:,m) = u_m;
        x(m,:) = x_m;
        y(m,:) = y_m;
        objValue2(m) = solve2_obj(L,a_m_L,rho,u_m,y_m,lambda1,x_m,lambda2,G,w,p,m);
%         norm(x(m,:))^2/(norm(y(m,:))^2+sigma*norm(u(:,:,m)^2))/gamma_r
%         [m mu]

    end
%     fprintf("----------")
%     [objValue1,objValue2]
%     [sum(abs(u-u1),'all') / sum(abs(u),'all') sum(abs(u-u1),'all') / sum(abs(u1),'all')]
%     [sum(abs(x-x1),'all') / sum(abs(x),'all') sum(abs(x-x1),'all') / sum(abs(x1),'all')]
%     [sum(abs(y-y1),'all') / sum(abs(y),'all') sum(abs(y-y1),'all') / sum(abs(y1),'all')]
end
%end

