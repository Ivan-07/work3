function S = solve4(M_T,L,v,u,a,E_r,y,g,rho,lambda2)
% clc;clear all
% 
% M_T = 8;M_R=8;
% L=8;
% v=rand()+1i*rand();
% u=rand(M_R,L);
% a=rand(M_R,M_T,L);
% g=1*(rand(M_T,1)+1i*rand(M_T,1));
% E_r=1;
% y=rand(1,L)+1i*rand(1,L);
% lambda2 = rand(1,L);
% rho=[0.1 0.1];

cvx_begin quiet
    variable S(M_T,L)
    expression sum1(L)
    expression sum2(L)
    expression sum3(L)
    for l=1:L
        sum1(l) = u(:,l)'*a(:,:,l)*S(:,l);
        sum2(l) = square_abs(y(l)-g'*S(:,l)+rho(2)*lambda2(l));
        sum3(l) = sum_square_abs(S(:,l));
    end
    minimize(-2*real(conj(v)*sum(sum1))+0.5/rho(2)*sum(sum2))
    subject to
        sum(sum3)-E_r<=0;
cvx_end

% S1=S;
% 
% e = 1e-8;
% S = zeros(M_T,L);
% for l=1:L
%     S(:,l) = (-0.5/rho(2)*g*g')^(-1)*(v*(u(:,l)'*a(:,:,l))'-0.5/rho(2)*(y(l)+rho(2)*lambda2(l))*g)
% end
% 
% if norm(S)^2-E_r >= e
%     lam_lb = 0;lam_ub = 100;
%     while 1
%         lam = lam_ub;
%         for l=1:L
%             S(:,l) = (-0.5/rho(2)*g*g'+lam*eye(M_T))^(-1)*(v*(u(:,l)'*a(:,:,l))'-0.5/rho(2)*(y(l)+rho(2)*lambda2(l))*g);
%         end
%         if norm(S)^2-E_r <= -e
%             break;
%         else
%             lam_lb = lam;
%             lam_ub = lam*10;
%         end
%     end
%     while 1
%         lam = (lam_lb+lam_ub)/2;
%         for l=1:L
%             S(:,l) = (-0.5/rho(2)*g*g'+lam*eye(M_T))^(-1)*(v*(u(:,l)'*a(:,:,l))'-0.5/rho(2)*(y(l)+rho(2)*lambda2(l))*g);
%         end
%         judge = norm(S)^2-E_r;
%         if abs(judge) <= e
%             break;
%         elseif judge >= 0
%             lam_ub = lam;
%         else
%             lam_lb = lam;
%         end
%         abs(judge)
%     end
% end   
% [sum(abs(S-S1),'all') / sum(abs(S),'all') sum(abs(S-S1),'all') / sum(abs(S1),'all')]


