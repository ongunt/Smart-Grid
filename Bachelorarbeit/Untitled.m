clear
clc
syms s
[A,B,C,D,x0,states] =power_analyze('modelwithoutprim');
solve=det(s*eye(size(A,1))-A)




