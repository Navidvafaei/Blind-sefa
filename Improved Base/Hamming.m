function [hw]=Hamming(inp)
hw=zeros(size(inp));
for i=1:size(inp)
hw(i,1)=sum(de2bi(inp(i)));
end
