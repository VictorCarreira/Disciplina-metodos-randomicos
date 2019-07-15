clear
clc
d1=dlmread('fit.txt');

tam=size (d1)
np=tam(1,1)

d2(1,1)=0
d2(1,2)=0

d2(2,1)=100
d2(2,2)=0


%figure(1)


plot(d1(:,1),d1(:,2),'-k')
hold on
plot(d1(:,1),d1(:,3),'*b')
hold on
plot(d2(:,1),d2(:,2),'-r')

axis([0 100 -250 250])
