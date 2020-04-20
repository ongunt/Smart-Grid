function R=rouths(C)
n=length(C);

%R=zeros(n,ceil(n/2));
R(1,:)=C(1:2:n);
R(2,1:length(2:2:n))=C(2:2:n);
for i=3:n
 for j=1:ceil(n/2)-1
 R(i,j)=(R(i-1,1)*R(i-2,j+1)-R(i-2,1)*R(i-1,j+1))/R(i-1,1);
 end
end
