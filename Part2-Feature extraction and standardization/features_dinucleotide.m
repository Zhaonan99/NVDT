%计算二联核苷酸的频率、平均位置和二阶中心距，共48维维特征
function features_dinucleotide(fasta_file_path)
fid=fopen('interaction_ID_2bases.txt','wb') %建立新文本文档来存储蛋白质ID号
[Header, Se]=fastaread(fasta_file_path); %Header存储蛋白质名称，Se存储每个蛋白质的具体序列

for j=1:length(Se)
    S=Se(j);%选定第j个蛋白质的具体基因序列，是字符串格式
    W=S{1};%同S，是文本格式
    N=length(W);%此蛋白质的长度
%*****************************************************************计算2bases的频率
    n1=zeros(1,4);%目的是用于存放2bases-A*的频率
    n2=zeros(1,4);%目的是用于存放2bases-C*的频率
    n3=zeros(1,4);%目的是用于存放2bases-G*的频率
    n4=zeros(1,4);%目的是用于存放2bases-T*的频率   
    s1=zeros(1,max(n1));%目的是用于存放2bases-A*的位置
    s2=zeros(1,max(n2));%目的是用于存放2bases-C*的位置
    s3=zeros(1,max(n3));%目的是用于存放2bases-G*的位置
    s4=zeros(1,max(n4));%目的是用于存放2bases-T*的位置
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
    
%*******************************计算2bases的平均位置
%n——n1记录首字母为A的，n2记录首字母为C的，n3记录首字母为G的，n4记录首字母为T的
%列向量：n1=[nAA,nAC,nAG,nAT],...,n4=[nTA,nTC,nTG,nTT]
%s——s1记录首字母为A的，s2记录首字母为C的，s3记录首字母为G的，s4记录首字母为T的
%矩阵形式记录每个符合此类形式二联核苷酸对的位置：
%s1=[nAA1,nAA2，nAA3,..nAA80;nAC1,nAC2,nAC3,...nAC90;...;nAT1,nAT2,nAT3,...,nAT78]
%u——u1记录首字母为A的，u2记录首字母为C的，u3记录首字母为G的，u4记录首字母为T的
%列向量：u1=[uAA,uAC,uAG,uAT],...,u4=[uTA,uTC,uTG,uTT]
    u1=zeros(1,4);%目的是用于存放2bases的频率
    u2=zeros(1,4);
    u3=zeros(1,4);
    u4=zeros(1,4);   

    
    for k=1:4
        %*****************
        t1=0;
        if n1(k)>0 %二联核苷酸出现的类型个数大于0
            u1(k)=sum(s1(k,:))/n1(k); %记录平均位置
        else
            u1(k)=0;
        end
        for i=1:n1(k)
            t1=t1+(s1(k,i)-u1(k))^2;
        end
        if n1(k)>0
            D1(k,1)=t1/(n1(k)*(N-1));%改动一
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
            D2(k,1)=t2/(n2(k)*(N-1));%改动二
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
            D3(k,1)=t3/(n3(k)*(N-1));%改动三
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
            D4(k,1)=t4/(n4(k)*(N-1));%改动四
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
   
    
    ID=Header(j); %蛋白质的名称
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



    
    
