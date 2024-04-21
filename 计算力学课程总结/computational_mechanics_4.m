syms x a b c d e f; 
%%% 基函数
Phi1(x)=(x+1-exp(x));
Phi2(x)=(x+1-exp(x))^2;
Phi3(x)=(x+1-exp(x))^3;
Phi4(x)=(x+1-exp(x))^4;
Phi5(x)=(x+1-exp(x))^5;
Phi6(x)=(x+1-exp(x))^6;
%%% 近似函数
omega2(x)=a*Phi1(x)+b*Phi2(x);
omega3(x)=a*Phi1(x)+b*Phi2(x)+c*Phi3(x);
omega4(x)=a*Phi1(x)+b*Phi2(x)+c*Phi3(x)+d*Phi4(x);
omega5(x)=a*Phi1(x)+b*Phi2(x)+c*Phi3(x)+d*Phi4(x)+e*Phi5(x);
omega6(x)=a*Phi1(x)+b*Phi2(x)+c*Phi3(x)+d*Phi4(x)+e*Phi5(x)+f*Phi6(x);

ddomega2(x)=diff(omega2(x),2);
ddomega3(x)=diff(omega3(x),2);
ddomega4(x)=diff(omega4(x),2);
ddomega5(x)=diff(omega5(x),2);
ddomega6(x)=diff(omega6(x),2);

ddphi1(x)=diff(Phi1(x),2);
ddphi2(x)=diff(Phi2(x),2);
ddphi3(x)=diff(Phi3(x),2);
ddphi4(x)=diff(Phi4(x),2);
ddphi5(x)=diff(Phi5(x),2);
ddphi6(x)=diff(Phi6(x),2);

%%%
re1(x)=ddomega2(x)*ddphi1(x)-Phi1(x);
re2(x)=ddomega2(x)*ddphi2(x)-Phi2(x);

result1=int(re1,0,1);
result2=int(re2,0,1);

eqns2=[result1==0,result2==0];
vars2=[a,b];
[sola2,solb2]=solve(eqns2,vars2);

%%%
re1(x)=ddomega3(x)*ddphi1(x)-Phi1(x);
re2(x)=ddomega3(x)*ddphi2(x)-Phi2(x);
re3(x)=ddomega3(x)*ddphi3(x)-Phi3(x);

result1=int(re1,0,1);
result2=int(re2,0,1);
result3=int(re3,0,1);

eqns3=[result1==0,result2==0,result3==0];
vars3=[a,b,c];
[sola3,solb3,solc3]=solve(eqns3,vars3);

%%%
re1(x)=ddomega4(x)*ddphi1(x)-Phi1(x);
re2(x)=ddomega4(x)*ddphi2(x)-Phi2(x);
re3(x)=ddomega4(x)*ddphi3(x)-Phi3(x);
re4(x)=ddomega4(x)*ddphi4(x)-Phi4(x);

result1=int(re1,0,1);
result2=int(re2,0,1);
result3=int(re3,0,1);
result4=int(re4,0,1);

eqns4=[result1==0,result2==0,result3==0,result4==0];
vars4=[a,b,c,d];
[sola4,solb4,solc4,sold4]=solve(eqns4,vars4);

%%%
re1(x)=ddomega5(x)*ddphi1(x)-Phi1(x);
re2(x)=ddomega5(x)*ddphi2(x)-Phi2(x);
re3(x)=ddomega5(x)*ddphi3(x)-Phi3(x);
re4(x)=ddomega5(x)*ddphi4(x)-Phi4(x);
re5(x)=ddomega5(x)*ddphi5(x)-Phi5(x);

result1=int(re1,0,1);
result2=int(re2,0,1);
result3=int(re3,0,1);
result4=int(re4,0,1);
result5=int(re5,0,1);

eqns5=[result1==0,result2==0,result3==0,result4==0,result5==0];
vars5=[a,b,c,d,e];
[sola5,solb5,solc5,sold5,sole5]=solve(eqns5,vars5);

%%%
re1(x)=ddomega6(x)*ddphi1(x)-Phi1(x);
re2(x)=ddomega6(x)*ddphi2(x)-Phi2(x);
re3(x)=ddomega6(x)*ddphi3(x)-Phi3(x);
re4(x)=ddomega6(x)*ddphi4(x)-Phi4(x);
re5(x)=ddomega6(x)*ddphi5(x)-Phi5(x);
re6(x)=ddomega6(x)*ddphi6(x)-Phi6(x);

result1=int(re1,0,1);
result2=int(re2,0,1);
result3=int(re3,0,1);
result4=int(re4,0,1);
result5=int(re5,0,1);
result6=int(re6,0,1);

eqns6=[result1==0,result2==0,result3==0,result4==0,result5==0,result6==0];
vars6=[a,b,c,d,e,f];
[sola6,solb6,solc6,sold6,sole6,solf6]=solve(eqns6,vars6);

%%%
resomega2(x)=sola2*Phi1(x)+solb2*Phi2(x);
resomega3(x)=sola3*Phi1(x)+solb3*Phi2(x)+solc3*Phi3(x);
resomega4(x)=sola4*Phi1(x)+solb4*Phi2(x)+solc4*Phi3(x)+sold4*Phi4(x);
resomega5(x)=sola5*Phi1(x)+solb5*Phi2(x)+solc5*Phi3(x)+sold5*Phi4(x)+sole5*Phi5(x);
resomega6(x)=sola6*Phi1(x)+solb6*Phi2(x)+solc6*Phi3(x)+sold6*Phi4(x)+sole6*Phi5(x)+solf6*Phi6(x);

fin2=double(sola2*(2-exp(1))+solb2*(2-exp(1))^2);
fin3=double(sola3*(2-exp(1))+solb3*(2-exp(1))^2+solc3*(2-exp(1))^3);
fin4=double(sola4*(2-exp(1))+solb4*(2-exp(1))^2+solc4*(2-exp(1))^3+sold4*(2-exp(1))^4);
fin5=double(sola5*(2-exp(1))+solb5*(2-exp(1))^2+solc5*(2-exp(1))^3+sold5*(2-exp(1))^4+sole5*(2-exp(1))^5);
fin6=double(sola6*(2-exp(1))+solb6*(2-exp(1))^2+solc6*(2-exp(1))^3+sold6*(2-exp(1))^4+sole6*(2-exp(1))^5+solf6*(2-exp(1))^6);

x = 0:1/100:1;
y = x.^2.*(x.^2-4.*x+6)/24;
plot(x,y,x,resomega2(x),x,resomega3(x),x,resomega4(x),x,resomega5(x),x,resomega6(x));

legend('精确解','2','3','4','5','6');
format short
fin2
fin3
fin4
fin5
fin6
