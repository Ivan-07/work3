function echo_comm = get_commSignal(PRT,G,w)

x = 1/sqrt(2)*(rand(1,PRT)+1i*rand(1,PRT));
echo_comm = G*w*x;

end