function h_indicator = get_indicator(x,y,w,h,L,S,g)

h_indicator = 0;

for l=1:L
    h_indicator = max(h_indicator,abs(y(l)-g'*S(:,l)));
end

h_indicator = max(h_indicator,abs(x-h'*w));

end