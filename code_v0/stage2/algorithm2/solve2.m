function u = solve2(M_R,L,v,G,w,sigma,a,S)
% 
% M_R = 8;
% L=8;
% v=rand+rand*1i;
% G=rand(8,8)+1i*rand(8,8);
% w=rand(8,1)+1i*rand(8,1);
% sigma=rand;
% a = rand(8,8,8)+1i*rand(8,8,8);
% S=rand(8,8)+1i*rand(8,8);
solve_choice = 1;

if solve_choice == 1
    cvx_begin quiet
        variable u(M_R,L) complex
        expression sum1(L)
        expression sum2(L)
        expression sum3(L)
        for l=1:L
            sum1(l) = square_abs(u(:,l)'*G*w);
            sum2(l) = u(:,l)'*a(:,:,l)*S(:,l);
            sum3(l) = sum_square_abs(u(:,l));
        end
        minimize(square_abs(v)*(sum(sum1)+sigma*sum(sum3))-2*real(conj(v)*sum(sum2)))
    cvx_end
%     u1=u;
elseif solve_choice == 2
    u = zeros(M_R,L);
    for l=1:L
        u(:,l) = (abs(v)^2*(G*w)*(G*w)'+abs(v)^2*sigma*eye(M_R))^(-1)*conj(v)*a(:,:,l)*S(:,l);
    end

%     [sum(abs(u-u1),'all') / sum(abs(u),'all') sum(abs(u-u1),'all') / sum(abs(u1),'all')]
end


