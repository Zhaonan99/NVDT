%��������������Ƶ�ʡ�ƽ��λ�úͶ������ľ࣬��48άά����
function features_dinucleotide(fasta_file_path)
fid=fopen('interaction_ID_2bases.txt','wb') %�������ı��ĵ����洢������ID��
[Header, Se]=fastaread(fasta_file_path); %Header�洢���������ƣ�Se�洢ÿ�������ʵľ�������

for j=1:length(Se)
    S=Se(j);%ѡ����j�������ʵľ���������У����ַ�����ʽ
    W=S{1};%ͬS�����ı���ʽ
    N=length(W);%�˵����ʵĳ���
%*****************************************************************����2bases��Ƶ��
    n1=zeros(1,4);%Ŀ�������ڴ��2bases-A*��Ƶ��
    n2=zeros(1,4);%Ŀ�������ڴ��2bases-C*��Ƶ��
    n3=zeros(1,4);%Ŀ�������ڴ��2bases-G*��Ƶ��
    n4=zeros(1,4);%Ŀ�������ڴ��2bases-T*��Ƶ��   
    s1=zeros(1,max(n1));%Ŀ�������ڴ��2bases-A*��λ��
    s2=zeros(1,max(n2));%Ŀ�������ڴ��2bases-C*��λ��
    s3=zeros(1,max(n3));%Ŀ�������ڴ��2bases-G*��λ��
    s4=zeros(1,max(n4));%Ŀ�������ڴ��2bases-T*��λ��
    for i=1:N-1
        if W(i)=='A'
            if W(i+1)=='A'
                n1(1)=n1(1)+1;
                s1(1,n1(1))=i+(i+1);
            end
            if W(i+1)=='C' 
                n1(2)=n1(2)+1;
                s1(2,n1(2))=i+(i+1);
            end
            if W(i+1)=='G'
                n1(3)=n1(3)+1;
                s1(3,n1(3))=i+(i+1);
            end
            if W(i+1)=='T'
                n1(4)=n1(4)+1;
                s1(4,n1(4))=i+(i+1);
            end
        end
        
        
        if W(i)=='C'
            if W(i+1)=='A'
                n2(1)=n2(1)+1;
                s2(1,n2(1))=i+(i+1);
            end
            if W(i+1)=='C' 
                n2(2)=n2(2)+1;
                s2(2,n2(2))=i+(i+1);
            end
            if W(i+1)=='G'
                n2(3)=n2(3)+1;
                s2(3,n2(3))=i+(i+1);
            end
            if W(i+1)=='T'
                n2(4)=n2(4)+1;
                s2(4,n2(4))=i+(i+1);
            end
        end
        if W(i)=='G'
            if W(i+1)=='A'
                n3(1)=n3(1)+1;
                s3(1,n3(1))=i+(i+1);
            end
            if W(i+1)=='C' 
                n3(2)=n3(2)+1;
                s3(2,n3(2))=i+(i+1);
            end
            if W(i+1)=='G'
                n3(3)=n3(3)+1;
                s3(3,n3(3))=i+(i+1);
            end
            if W(i+1)=='T'
                n3(4)=n3(4)+1;
                s3(4,n3(4))=i+(i+1);
            end
        end
        if W(i)=='T'
            if W(i+1)=='A'
                n4(1)=n4(1)+1;
                s4(1,n4(1))=i+(i+1);
            end
            if W(i+1)=='C' 
                n4(2)=n4(2)+1;
                s4(2,n4(2))=i+(i+1);
            end
            if W(i+1)=='G'
                n4(3)=n4(3)+1;
                s4(3,n4(3))=i+(i+1);
            end
            if W(i+1)=='T'
                n4(4)=n4(4)+1;
                s4(4,n4(4))=i+(i+1);
            end
        end      
    end
    
%*******************************����2bases��ƽ��λ��
%n����n1��¼����ĸΪA�ģ�n2��¼����ĸΪC�ģ�n3��¼����ĸΪG�ģ�n4��¼����ĸΪT��
%��������n1=[nAA,nAC,nAG,nAT],...,n4=[nTA,nTC,nTG,nTT]
%s����s1��¼����ĸΪA�ģ�s2��¼����ĸΪC�ģ�s3��¼����ĸΪG�ģ�s4��¼����ĸΪT��
%������ʽ��¼ÿ�����ϴ�����ʽ����������Ե�λ�ã�
%s1=[nAA1,nAA2��nAA3,..nAA80;nAC1,nAC2,nAC3,...nAC90;...;nAT1,nAT2,nAT3,...,nAT78]
%u����u1��¼����ĸΪA�ģ�u2��¼����ĸΪC�ģ�u3��¼����ĸΪG�ģ�u4��¼����ĸΪT��
%��������u1=[uAA,uAC,uAG,uAT],...,u4=[uTA,uTC,uTG,uTT]
    u1=zeros(1,4);%Ŀ�������ڴ��2bases��Ƶ��
    u2=zeros(1,4);
    u3=zeros(1,4);
    u4=zeros(1,4);   

    
    for k=1:4
        %*****************
        t1=0;
        if n1(k)>0 %������������ֵ����͸�������0
            u1(k)=sum(s1(k,:))/n1(k); %��¼ƽ��λ��
        else
            u1(k)=0;
        end
        for i=1:n1(k)
            t1=t1+(s1(k,i)-u1(k))^2;
        end
        if n1(k)>0
            D1(k,1)=t1/(n1(k)*(N-1));%�Ķ�һ
        else
            n1(k)==0;
            D1(k,1)=0;
        end     

        %***************
        t2=0;        
        if n2(k)>0
            u2(k)=sum(s2(k,:))/n2(k);
        else
            u2(k)=0;
        end
        for i=1:n2(k)
            t2=t2+(s2(k,i)-u2(k))^2;
        end
        if n2(k)>0
            D2(k,1)=t2/(n2(k)*(N-1));%�Ķ���
        else
            n2(k)==0
            D2(k,1)=0;
        end        
        %***************
        t3=0;        
        if n3(k)>0
            u3(k)=sum(s3(k,:))/n3(k);
        else
            u3(k)=0;
        end
        for i=1:n3(k)
            t3=t3+(s3(k,i)-u3(k))^2;
        end
        if n3(k)>0
            D3(k,1)=t3/(n3(k)*(N-1));%�Ķ���
        else
            n3(k)==0
            D3(k,1)=0;
        end
        %*****************
        t4=0;        
        if n4(k)>0
            u4(k)=sum(s4(k,:))/n4(k);
        else
            u4(k)=0;
        end
        for i=1:n4(k)
            t4=t4+(s4(k,i)-u4(k))^2;
        end
        if n4(k)>0
            D4(k,1)=t4/(n4(k)*(N-1));%�Ķ���
        else
            n4(k)==0
            D4(k,1)=0;
        end
    end
    

    for k=1:4
       d1(k)=D1(k,1);
    end
    for k=1:4
       d2(k)=D2(k,1);
    end
    for k=1:4
       d3(k)=D3(k,1);
    end
    for k=1:4
       d4(k)=D4(k,1);
    end
   
    
    ID=Header(j); %�����ʵ�����
    fprintf(fid,'%s\n',ID{1});
    fprintf(fid,'\r\n');
    n=[n1,n2,n3,n4];
    u=[u1,u2,u3,u4];
    d=[d1,d2,d3,d4];

    T=[n,u,d];
    C(j,:)=T;
    j=j+1;
end
csvwrite('nature_vector_2bases.csv',C);
fclose(fid)



    
    
