function [mat] = number2matrix(n,nodes)

binarylength = nodes^2;
binrep = dec2bin(n,binarylength);
vec = str2num(binrep(:));
mat = reshape(vec,[],nodes);
