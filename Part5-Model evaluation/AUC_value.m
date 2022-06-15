clc
clear

dec1 = load('1pre.txt');
predict1 = dec1(:,2);
value = load('true.txt');
truth = value;

%%
x = 1.0;y = 1.0;  
pos_num = sum(truth==1);neg_num = sum(truth==0);  
x_step = 1.0/neg_num;y_step = 1.0/pos_num;  
[predict1,index1] = sort(predict1);  
truth= truth(index1);  
for i=1:length(truth)  
    if truth(i) == 1  
        y = y - y_step;  
    else  
        x = x - x_step;  
    end  
    X(i)=x;  
    Y(i)=y;  
end  
figure;
XX = X;YY = Y;
plot(X,Y,'r-','LineWidth',0.6,'MarkerSize',1);  
xlim([0 1]);ylim([0 1]);
xlabel('1-Specificity');ylabel('Sensitivity'); 
auc1 = -trapz(X,Y)

legend('AUC = 0.9991')

hold off

