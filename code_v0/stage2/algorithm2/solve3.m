function w = solve3(L,u,G,v,rho,h,x,lambda1,E_c,N_T)

% L=8;
% u=rand(8,L)+1i*rand(8,L);
% G=rand(8,8)+1i*rand(8,8);
% v=rand+rand*1i;
% rho=rand(1,2);
% lambda1=rand;
% E_c=1;
% N_T=8;
% h=rand(8,1);
% x=rand+rand*1i;
solve_choice = 1;

if solve_choice == 1
    cvx_begin quiet
        variable w(N_T) complex
        expression sum1(L)
        for l=1:L
            sum1(l) = square_abs(u(:,l)'*G*w);
        end
        minimize(square_abs(v)*sum(sum1)+0.5/rho(1)*square_abs(x-h'*w+rho(1)*lambda1))
        subject to 
            sum_square_abs(w)-E_c<=0;
    cvx_end

%     w1 = w;
elseif solve_choice == 2
    e = 1e-8;
    tmp = 0;
    for l=1:L
        tmp = tmp+(u(:,l)'*G)'*(u(:,l)'*G);
    end

    w = (abs(v)^2*tmp+0.5/rho(1)*h*h')^(-1)*(0.5/rho(1)*(x+rho(1)*lambda1)*h);

    if norm(w)^2-E_c >= e
        lam_lb = 0;lam_ub = 100;
        while 1
            lam = lam_ub;
            w = (abs(v)^2*tmp+0.5/rho(1)*h*h'+lam*eye(N_T))^(-1)*(0.5/rho(1)*(x+rho(1)*lambda1)*h);
            if norm(w)^2-E_c <= -e
                break;
            else
                lam_lb = lam;
                lam_ub = lam*10;
            end
        end
        while 1
            lam = (lam_lb+lam_ub)/2;
            w = (abs(v)^2*tmp+0.5/rho(1)*h*h'+lam*eye(N_T))^(-1)*(0.5/rho(1)*(x+rho(1)*lambda1)*h);
            judge = norm(w)^2-E_c;
            if abs(judge) <= e
                break;
            elseif judge >= 0
                lam_lb = lam;
            else
                lam_ub = lam;
            end
%             [abs(judge)]
        end
    end

%     [sum(abs(w-w1),'all') / sum(abs(w),'all') sum(abs(w-w1),'all') / sum(abs(w1),'all')]
end


