YM = zeros(14,7,14);
load('data.mat')
[N_l, N_p, N_y] = size(YM);

for i = 1:N_y
    for j = 1:N_p
        YM(1:YN(i,j),j,i) = 1;
    end
end
YM_r = circshift(rot90(YM,2),1,2);

%% 按照旋转反转构造各种情况
Y_all = zeros(N_l, N_p, N_y, 2, N_p, N_l);
% Y_all(:,:,:,1,1,1) = YM;

for i = 1:N_l
    for j = 1:N_p
        Y_all(i:end,:,:,1,j,i) = circshift(YM(1:end-i+1,:,:), j-1, 2);
        
        Y_all(1:i,:,:,2,j,i) = circshift(YM_r(end-i+1:end,:,:), j-1, 2);
    end
end

%% 构造前系数矩阵
A = zeros(N_l*N_p + N_y + N_l, N_y*2*N_p*N_l);
b = ones(N_l*N_p + N_y + N_l, 1);
c = zeros(N_y*2*N_p*N_l,1);

for i = 1:N_l
    for j = 1:N_p
        A((j-1)*N_l + i,:) = reshape(Y_all(i,j,:,:,:,:),[],1);
    end
end

for i = 1:N_y
    A(N_l*N_p + i, i:N_y:end) = 1;
end

for i = 1:N_l
    A(N_l*N_p + N_y + i, (N_y*2*N_p*(i-1)+1):(N_y*2*N_p*i) ) = 1;
end

%% 利用内部整数规划求解器求解
x = intlinprog(c,1:N_y*2*N_p*N_l,[],[],A,b,zeros(N_y*2*N_p*N_l,1),ones(N_y*2*N_p*N_l,1));

%% 处理最终结果
p = reshape(x, N_y,2,N_p,N_l);
p_out = zeros(N_y,3);
for i = 1:N_y
    for j = 1:N_p
        for k = 1:N_l
            if round(p(i,1,j,k))==1
                p_out(i,1) = 1;
                p_out(i,2) = j;
                p_out(i,3) = k;
            end
            if round(p(i,2,j,k))==1
                p_out(i,1) = 2;
                p_out(i,2) = j;
                p_out(i,3) = k;
            end
        end
    end
end

