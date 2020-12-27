function Y=Ybus_new(casedata)
result=loadcase(casedata); 
Bus=result.bus; 
Branch=result.branch; 
%% Ybus calculation
z= length(Bus(:,1));
Y=zeros(z,z);

for i=1:z
    for j=1:z
        if i==j
            a=find(Branch(:,1)==i);
            b=find(Branch(:,2)==i);
            for k=1:length(a)
                Y(i,j) = Y(i,j)+ 1/(Branch(a(k),3)+Branch(a(k),4)*1i)+Branch(a(k),5)*1i/2;
            end
            for kk=1:length(b)
                Y(i,j) = Y(i,j)+ 1/(Branch(b(kk),3)+Branch(b(kk),4)*1i)+Branch(b(kk),5)*1i/2;
            end         
        else  
             c = find(Branch(:,1)==i);
             d = find(Branch(:,2)==j);
            for m=1:length(c)
                for n=1:length(d)
                    if c(m)==d(n)  
                    Y(i,j) = Y(i,j)-1/(Branch(c(m),3)+Branch(c(m),4)*1i);
                    end
                end 
            end
        end 
         Y(j,i)=Y(i,j);
     end
end

%% system data Xd and R




end 