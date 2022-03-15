function [v,x,y] = solve1(L,rho,u,a,G,x_pre,h,sigma,lambda1,lambda2,gamma_c,g,S,w)

a_ = 0; b_ = 0;

for l=1:L
    a_ = a_+abs(u(:,l)'*G*w)^2+sigma*norm(u(:,l))^2;
    b_ = b_+u(:,l)'*a(:,:,l)*S(:,l);
end

v = b_/a_;

solve_choice = 1;
e = 1e-7;

if solve_choice == 1
    cvx_begin quiet
        variable x complex
        variable y(1,L) complex
        minimize(square_abs(x-h'*w+rho(1)*lambda1)+sum_square_abs(y-g'*S+rho(2)*lambda2))
        subject to
            gamma_c*(sum_square_abs(y)+L*sigma)+L*square_abs(x_pre)-2*L*real(conj(x_pre)*x)<=0
    cvx_end


elseif solve_choice == 2
    x = h'*w-rho(1)*lambda1;
    y = g'*S-rho(2)*lambda2;

    if solve1_judge(gamma_c,y,x,sigma,L,x_pre) >= e
        mu_lb = 0; mu_ub = 100;
        while 1
            mu = mu_ub;
            x = h'*w-rho(1)*lambda1+mu*L*x_pre;
            y = 1/(1+mu*gamma_c)*(g'*S-rho(2)*lambda2);
            if solve1_judge(gamma_c,y,x,sigma,L,x_pre) <= -e
                break;
            else
                mu_lb = mu;
                mu_ub = mu*10;
            end
        end
        while 1
            mu = (mu_lb+mu_ub)/2;
            x = h'*w-rho(1)*lambda1+mu*L*x_pre;
            y = 1/(1+mu*gamma_c)*(g'*S-rho(2)*lambda2);
            judge = solve1_judge(gamma_c,y,x,sigma,L,x_pre);
            if abs(judge) <= e
                break;
            elseif judge >= 0
                mu_lb = mu;
            else
                mu_ub = mu;
            end
        end
    end

end