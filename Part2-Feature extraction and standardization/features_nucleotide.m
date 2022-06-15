%计算核苷酸的频率、平均位置和中心二阶矩，共12维特征
function features_nucleotide(fasta_file_path)
fid=fopen('interaction_ID_nucleotide.txt','wb');
[Header, Se]=fastaread(fasta_file_path);
for j=1:length(Se)
    S=Se(j);
    W=S{1};
    N=length(W);
    n=zeros(1,4);
    s=0;
    for i=1:N 
        if W(i)=='A' 
        n(1)=n(1)+1;
        s(1,n(1))=i;
        end
        if W(i)=='C' 
        n(2)=n(2)+1;
        s(2,n(2))=i;
        end
        if W(i)=='G' 
        n(3)=n(3)+1;
        s(3,n(3))=i;
        end
        if W(i)=='T' 
        n(4)=n(4)+1;
        s(4,n(4))=i;
        end
    end
    for k=1:4
        if n(k)>0
        u(k)=sum(s(k,:))/n(k);
        end
       if n(k)==0
       u(k)=0;
       end
       t=0;
       for i=1:n(k)
       t=t+(s(k,i)-u(k))^2;
       end
       if n(k)>0
       D(k,1)=t/(n(k)*N);
       end
       if n(k)==0
       D(k,1)=0;
       end
    end
    for k=1:4
       d(k)=D(k,1);
    end
    ID=Header(j);
    fprintf(fid,'%s\n',ID{1})
    fprintf(fid,'\r\n')
    T=[n,u,d];
    C(j,:)=T;
    j=j+1;
end
csvwrite('nature_vector_nucleotide.csv',C)
fclose(fid)


