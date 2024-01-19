function parcompute()

%1: Using parfor: Phân chia dữ liệu ra nhiều nhóm và tính toán
parfor i=1:10
    a(i) = i;
end
% Control the number of workers
parfor (A = 1:10, 20)
AA=1;
end
%2: Using GPU: Chuyển Data từ CPU sang GPU và tính toán
gpuDevice; % Check GPU
g = gpuArray(1:10);
for i=1:10
    a(i) = g(i);
end
end