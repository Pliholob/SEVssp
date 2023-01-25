clc
clear
tic
% 
N   = 256;
M =10e+3; 
ID_frame = 1;
u = randn(N,1);
u = u/sqrt(sum(u.^2));
SSp_DB = table(ID_frame, {u});
SSp_DB.Properties.VariableNames{'Var2'} = 'u';
% plot(SSp_DB.u{1,:});
for ID_frame = 2:M 
    u = randn(N,1);
    u = u - mean(u);
    u = u/sqrt(sum(u.^2));

    ST_u.u = {u};
    ST_u.ID_frame = ID_frame;
    SSp_DB = [SSp_DB;struct2table(ST_u)];

    disp(ID_frame)
end