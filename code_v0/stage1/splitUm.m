function returnVec = mySplit(vec)

unitSize = length(vec)/L;
returnVec = zeros(unitSize,L);
for l=1:L
    returnVec(:,l) = u_m(unitSize*(l-1)+1:unitSize*l);
end

end