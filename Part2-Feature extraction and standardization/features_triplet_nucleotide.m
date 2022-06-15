%计算三联核苷酸的频率、平均位置和二阶中心距，共64*3=192维特征
function features_triplet_nucleotide(fasta_file_path)
fid=fopen('interaction_ID_3bases.txt','wb'); %建立新文本文档来存储蛋白质ID号
[Header, Se]=fastaread(fasta_file_path); %Header存储蛋白质名称，Se存储每个蛋白质的具体序列

for j=1:length(Se)
    S=Se(j);%选定第j个蛋白质的具体基因序列，是字符串格式
    W=S{1};%同S，是文本格式
    N=length(W);%此蛋白质的长度
%%*****************************************************************计算3base的频率
    %n——目的是用于存放3bases的频率
    %列向量：n=[nAAA,nAAC,nAAG,nAAT,nACA,nACC,nACG,nACT,...,nTTA,nTTC,nTTG,nTTT]
    %s——目的是存放每种类型3bases每一次出现的位置
    %矩阵形式64*？，每一行是一种类型的密码子
    %u——目的是存放每种类型3bases的平均位置
    %列向量：u=[uAAA,uAAC,uAAG,uAAT,uACA,uACC,uACG,uACT,...,uTTA,uTTC,uTTG,uTTT]
 n=zeros(1,64); 
 s=zeros(64,N-2);
  
    for i=1:1:N-2
        %%密码子为A**
        if W(i)=='A'
            
            if W(i+1)=='A'
                if W(i+2)=='A'
                    n(1)=n(1)+1;
                    s(1,n(1))=i+(i+1)+(i+2);
                end
                if W(i+2)=='C'
                    n(2)=n(2)+1;
                    s(2,n(2))=i+(i+1)+(i+2);
                end
                if W(i+2)=='G'
                    n(3)=n(3)+1;
                    s(3,n(3))=i+(i+1)+(i+2);
                end
                if W(i+2)=='T'
                    n(4)=n(4)+1;
                    s(4,n(4))=i+(i+1)+(i+2);
                end
            end
            
            if W(i+1)=='C'
                if W(i+2)=='A'
                    n(5)=n(5)+1;
                    s(5,n(5))=i+(i+1)+(i+2);
                end
                if W(i+2)=='C'
                    n(6)=n(6)+1;
                    s(6,n(6))=i+(i+1)+(i+2);
                end
                if W(i+2)=='G'
                    n(7)=n(7)+1;
                    s(7,n(7))=i+(i+1)+(i+2);
                end
                if W(i+2)=='T'
                    n(8)=n(8)+1;
                    s(8,n(8))=i+(i+1)+(i+2);
                end
            end
            
            if W(i+1)=='G'
                if W(i+2)=='A'
                    n(9)=n(9)+1;
                    s(9,n(9))=i+(i+1)+(i+2);
                end
                if W(i+2)=='C'
                    n(10)=n(10)+1;
                    s(10,n(10))=i+(i+1)+(i+2);
                end
                if W(i+2)=='G'
                    n(11)=n(11)+1;
                    s(11,n(11))=i+(i+1)+(i+2);
                end
                if W(i+2)=='T'
                    n(12)=n(12)+1;
                    s(12,n(12))=i+(i+1)+(i+2);
                end
            end
            
            if W(i+1)=='T'
                if W(i+2)=='A'
                    n(13)=n(13)+1;
                    s(13,n(13))=i+(i+1)+(i+2);
                end
                if W(i+2)=='C'
                    n(14)=n(14)+1;
                    s(14,n(14))=i+(i+1)+(i+2);
                end
                if W(i+2)=='G'
                    n(15)=n(15)+1;
                    s(15,n(15))=i+(i+1)+(i+2);
                end
                if W(i+2)=='T'
                    n(16)=n(16)+1;
                    s(16,n(16))=i+(i+1)+(i+2);
                end
            end
            
        end
        
        
        %%密码子为C**
        if W(i)=='C'
            
            if W(i+1)=='A'
                if W(i+2)=='A'
                    n(17)=n(17)+1;
                    s(17,n(17))=i+(i+1)+(i+2);
                end
                if W(i+2)=='C'
                    n(18)=n(18)+1;
                    s(18,n(18))=i+(i+1)+(i+2);
                end
                if W(i+2)=='G'
                    n(19)=n(19)+1;
                    s(19,n(19))=i+(i+1)+(i+2);
                end
                if W(i+2)=='T'
                    n(20)=n(20)+1;
                    s(20,n(20))=i+(i+1)+(i+2);
                end
            end
            
            if W(i+1)=='C'
                if W(i+2)=='A'
                    n(21)=n(21)+1;
                    s(21,n(21))=i+(i+1)+(i+2);
                end
                if W(i+2)=='C'
                    n(22)=n(22)+1;
                    s(22,n(22))=i+(i+1)+(i+2);
                end
                if W(i+2)=='G'
                    n(23)=n(23)+1;
                    s(23,n(23))=i+(i+1)+(i+2);
                end
                if W(i+2)=='T'
                    n(24)=n(24)+1;
                    s(24,n(24))=i+(i+1)+(i+2);
                end
            end
            
            if W(i+1)=='G'
                if W(i+2)=='A'
                    n(25)=n(25)+1;
                    s(25,n(25))=i+(i+1)+(i+2);
                end
                if W(i+2)=='C'
                    n(26)=n(26)+1;
                    s(26,n(26))=i+(i+1)+(i+2);
                end
                if W(i+2)=='G'
                    n(27)=n(27)+1;
                    s(27,n(27))=i+(i+1)+(i+2);
                end
                if W(i+2)=='T'
                    n(28)=n(28)+1;
                    s(28,n(28))=i+(i+1)+(i+2);
                end
            end
            
            if W(i+1)=='T'
                if W(i+2)=='A'
                    n(29)=n(29)+1;
                    s(29,n(29))=i+(i+1)+(i+2);
                end
                if W(i+2)=='C'
                    n(30)=n(30)+1;
                    s(30,n(30))=i+(i+1)+(i+2);
                end
                if W(i+2)=='G'
                    n(31)=n(31)+1;
                    s(31,n(31))=i+(i+1)+(i+2);
                end
                if W(i+2)=='T'
                    n(32)=n(32)+1;
                    s(32,n(32))=i+(i+1)+(i+2);
                end
            end            
        end
        
        
        %%密码子为G**
        if W(i)=='G'
            
            if W(i+1)=='A'
                if W(i+2)=='A'
                    n(33)=n(33)+1;
                    s(33,n(33))=i+(i+1)+(i+2);
                end
                if W(i+2)=='C'
                    n(34)=n(34)+1;
                    s(34,n(34))=i+(i+1)+(i+2);
                end
                if W(i+2)=='G'
                    n(35)=n(35)+1;
                    s(35,n(35))=i+(i+1)+(i+2);
                end
                if W(i+2)=='T'
                    n(36)=n(36)+1;
                    s(36,n(36))=i+(i+1)+(i+2);
                end
            end
            
            if W(i+1)=='C'
                if W(i+2)=='A'
                    n(37)=n(37)+1;
                    s(37,n(37))=i+(i+1)+(i+2);
                end
                if W(i+2)=='C'
                    n(38)=n(38)+1;
                    s(38,n(38))=i+(i+1)+(i+2);
                end
                if W(i+2)=='G'
                    n(39)=n(39)+1;
                    s(39,n(39))=i+(i+1)+(i+2);
                end
                if W(i+2)=='T'
                    n(40)=n(40)+1;
                    s(40,n(40))=i+(i+1)+(i+2);
                end
            end
            
            if W(i+1)=='G'
                if W(i+2)=='A'
                    n(41)=n(41)+1;
                    s(41,n(41))=i+(i+1)+(i+2);
                end
                if W(i+2)=='C'
                    n(42)=n(42)+1;
                    s(42,n(42))=i+(i+1)+(i+2);
                end
                if W(i+2)=='G'
                    n(43)=n(43)+1;
                    s(43,n(43))=i+(i+1)+(i+2);
                end
                if W(i+2)=='T'
                    n(44)=n(44)+1;
                    s(44,n(44))=i+(i+1)+(i+2);
                end
            end
            
            if W(i+1)=='T'
                if W(i+2)=='A'
                    n(45)=n(45)+1;
                    s(45,n(45))=i+(i+1)+(i+2);
                end
                if W(i+2)=='C'
                    n(46)=n(46)+1;
                    s(46,n(46))=i+(i+1)+(i+2);
                end
                if W(i+2)=='G'
                    n(47)=n(47)+1;
                    s(47,n(47))=i+(i+1)+(i+2);
                end
                if W(i+2)=='T'
                    n(48)=n(48)+1;
                    s(48,n(48))=i+(i+1)+(i+2);
                end
            end            
        end
        
        %%密码子为T**
        if W(i)=='T'
            
            if W(i+1)=='A'
                if W(i+2)=='A'
                    n(49)=n(49)+1;
                    s(49,n(49))=i+(i+1)+(i+2);
                end
                if W(i+2)=='C'
                    n(50)=n(50)+1;
                    s(50,n(50))=i+(i+1)+(i+2);
                end
                if W(i+2)=='G'
                    n(51)=n(51)+1;
                    s(51,n(51))=i+(i+1)+(i+2);
                end
                if W(i+2)=='T'
                    n(52)=n(52)+1;
                    s(52,n(52))=i+(i+1)+(i+2);
                end
            end
            
            if W(i+1)=='C'
                if W(i+2)=='A'
                    n(53)=n(53)+1;
                    s(53,n(53))=i+(i+1)+(i+2);
                end
                if W(i+2)=='C'
                    n(54)=n(54)+1;
                    s(54,n(54))=i+(i+1)+(i+2);
                end
                if W(i+2)=='G'
                    n(55)=n(55)+1;
                    s(55,n(55))=i+(i+1)+(i+2);
                end
                if W(i+2)=='T'
                    n(56)=n(56)+1;
                    s(56,n(56))=i+(i+1)+(i+2);
                end
            end
            
            if W(i+1)=='G'
                if W(i+2)=='A'
                    n(57)=n(57)+1;
                    s(57,n(57))=i+(i+1)+(i+2);
                end
                if W(i+2)=='C'
                    n(58)=n(58)+1;
                    s(58,n(58))=i+(i+1)+(i+2);
                end
                if W(i+2)=='G'
                    n(59)=n(59)+1;
                    s(59,n(59))=i+(i+1)+(i+2);
                end
                if W(i+2)=='T'
                    n(60)=n(60)+1;
                    s(60,n(60))=i+(i+1)+(i+2);
                end
            end
            
            if W(i+1)=='T'
                if W(i+2)=='A'
                    n(61)=n(61)+1;
                    s(61,n(61))=i+(i+1)+(i+2);
                end
                if W(i+2)=='C'
                    n(62)=n(62)+1;
                    s(62,n(62))=i+(i+1)+(i+2);
                end
                if W(i+2)=='G'
                    n(63)=n(63)+1;
                    s(63,n(63))=i+(i+1)+(i+2);
                end
                if W(i+2)=='T'
                    n(64)=n(64)+1;
                    s(64,n(64))=i+(i+1)+(i+2);
                end
            end            
        end
    end

%%*******************************计算密码子的平均位置和二阶中心距
    u=zeros(1,64);
    D=zeros(1,64);
    for k=1:64
        t=0;
        if n(k)>0 %密码子出现的类型个数大于0
            u(k)=sum(s(k,:))/n(k); %记录平均位置
        else
            u(k)=0;
        end
        for i=1:n(k)
            t=t+(s(k,i)-u(k))^2;
        end
        if n(k)>0
            D(k,1)=t/(n(k)*(N-2));
        else
            n(k)==0
            D(k,1)=0;
        end     
    end
    for k=1:64
       d(k)=D(k,1);
    end
    ID=Header(j); %蛋白质的名称
    fprintf(fid,'%s\n',ID{1});
    fprintf(fid,'\r\n');
    T=[n,u,d];
    C(j,:)=T;
    j=j+1;
end
csvwrite('nature_vector_3bases.csv',C);
fclose(fid)
      
            
            
            




    
    
