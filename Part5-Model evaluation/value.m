function T=value(X,Y)
[number_data,dim]=size(X);
N=number_data;
TP=length(find(X-Y==0&X==1));
TN=length(find(X-Y==0&X==0));
FP=length(find((X-Y)==1));
FN=length(find((X-Y)==-1));
sen=TP/(TP+FN);
spe=TN/(TN+FP); 
pre=TP/(TP+FP);
acc=(TP+TN)/(TP+TN+FP+FN);
mcc=(TP*TN-FP*FN)/((TP+FN)*(TN+FP)*(TP+FP)*(TN+FN))^(1/2);
f=(2*pre*sen)/(pre+sen);
N1=N/2;

T=[TP,TN,FP,FN,acc,pre,sen,spe,mcc,f];






