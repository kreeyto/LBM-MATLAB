clc; clearvars; close all

A = [1 4 5 6 8 7];

aGpu = gpuArray(A);

for i = 1:length(aGpu)
    aGpu(i) = aGpu(i) + 1;
    disp(aGpu(i))
end