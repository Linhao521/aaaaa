clear all
clc
x=8,y=6;
sx=8,sy=6;
t=0;
h=0.0001;
p=0.93;
tao=0.5;
A1=[-2,1;0,-3];
A2=[-2,1;0,-3];
B1=[-1,0.25;0,-1];
B2=[-1,0.25;0,-1];
K1=[0.06,0;0,0.06];K2=K1;
L1=[0.07,0;0,0.07];L2=L1;
n=1,P=2,m=3,l=4;
ja1=0.14,ja2=4.3;
jb1=0.2,jb2=0.29;
jM=0.7;
jC=0.4;
sja1=0.14,sja2=4.3;
sjb1=0.2,sjb2=0.29;
sjM=0.7;
sjC=0.4;
umax=[5;10];umin=[-10;-5];pl=2;
for n=1:49990
    t(n+1)=t(n)+h;
    M1(n)=0.5*(1+cos(x(n)));
    M2(n)=0.5*(1-cos(x(n)));
    sM1(n)=0.5*(1+cos(sx(n)));
    sM2(n)=0.5*(1-cos(sx(n)));
    d=[1.5*sin(t(n));1.5*cos(t(n))];
    if t(n)<=tao
        v=[8;6],sv=[8;6];
    else
    v=[x(fix((t(n)-tao)/h+1));y(fix((t(n)-tao)/h+1))],sv=[sx(fix((t(n)-tao)/h+1));sy(fix((t(n)-tao)/h+1))];
    end
    w=[x(n);y(n)];
    sw=[sx(n);sy(n)];
    f1(n)=n*M1(n)*norm(w,2)*norm(w,2);
    f2(n)=n*M2(n)*norm(w,2)*norm(w,2);
    g1(n)=P*M1(n)*norm(w,1)*norm(v,1);
    g2(n)=P*M2(n)*norm(w,1)*norm(v,1);
    hhh(n)=m*norm(w,1);
    rrr(n)=l*norm(w,1);
    sf1(n)=n*sM1(n)*norm(sw,2)*norm(sw,2);
    sf2(n)=n*sM2(n)*norm(sw,2)*norm(sw,2);
    sg1(n)=P*sM1(n)*norm(sw,1)*norm(sv,1);
    sg2(n)=P*sM2(n)*norm(sw,1)*norm(sv,1);
    shhh(n)=m*norm(sw,1);
    srrr(n)=l*norm(sw,1);
    G=(t(n+1)-t(1:n)).^(p-1);
    ff1=G.*f1/gamma(p);
    ff2=G.*f2/gamma(p);
    gg1=G.*g1/gamma(p);
    gg2=G.*g2/gamma(p);
    hh1=G.*hhh/gamma(p);
    rr1=G.*rrr/gamma(p);
    ja1(n+1)=ja1(1)+h*sum(ff1);
    ja2(n+1)=ja2(1)+h*sum(ff2);
    jb1(n+1)=jb1(1)+h*sum(gg1);
    jb2(n+1)=jb2(1)+h*sum(gg2);
    jM(n+1)=jM(1)+h*sum(hh1);
    jC(n+1)=jC(1)+h*sum(rr1);
    sff1=G.*sf1/gamma(p);
    sff2=G.*sf2/gamma(p);
    sgg1=G.*sg1/gamma(p);
    sgg2=G.*sg2/gamma(p);
    shh1=G.*shhh/gamma(p);
    srr1=G.*srrr/gamma(p);
    sja1(n+1)=sja1(1)+h*sum(sff1);
    sja2(n+1)=sja2(1)+h*sum(sff2);
    sjb1(n+1)=sjb1(1)+h*sum(sgg1);
    sjb2(n+1)=sjb2(1)+h*sum(sgg2);
    sjM(n+1)=sjM(1)+h*sum(shh1);
    sjC(n+1)=sjC(1)+h*sum(srr1);
    AAA11=-M1(n)*(K1*[x(n);y(n)]+ja1(n)*[x(n);y(n)]+L1*v+jb1(n)*norm(v,1)*atan(10*w)+jM(n)*atan(10*w)+jC(n)*atan(10*w));
    AAA12=-M2(n)*(K2*[x(n);y(n)]+ja2(n)*[x(n);y(n)]+L2*v+jb2(n)*norm(v,1)*atan(10*w)+jM(n)*atan(10*w)+jC(n)*atan(10*w));
    AAA=AAA11+AAA12;
    u1(n)=AAA(1);
    u2(n)=AAA(2);
    sAAA11=-sM1(n)*(K1*[sx(n);sy(n)]+sja1(n)*[sx(n);sy(n)]+L1*sv+sjb1(n)*norm(sv,1)*atan(10*sw)+sjM(n)*atan(10*sw)+sjC(n)*atan(10*sw));
    sAAA12=-sM2(n)*(K2*[sx(n);sy(n)]+sja2(n)*[sx(n);sy(n)]+L2*sv+sjb2(n)*norm(sv,1)*atan(10*sw)+sjM(n)*atan(10*sw)+sjC(n)*atan(10*sw));
    sAAA=sAAA11+sAAA12;
    su1(n)=sAAA(1);
    su2(n)=sAAA(2);
       if u1(n)>=umax(1)
satu1(n)=umax(1);
du1(n)=satu1(n)-u1(n);
elseif u1(n)<=umin(1)
satu1(n)=umin(1);
du1(n)=satu1(n)-u1(n);
else
satu1(n)=pl*u1(n);
du1(n)=satu1(n)-u1(n);
    end
  if u2(n)>=umax(2)
satu2(n)=umax(2);
du2(n)=satu2(n)-u2(n);
elseif u2(n)<=umin(2)
satu2(n)=umin(2);
du2(n)=satu2(n)-u2(n);
else
satu2(n)=pl*u2(n);
du2(n)=satu2(n)-u2(n);
  end
  sat=[satu1(n);satu2(n)];
  ssat=[tanh(su1(n));tanh(su2(n))];
    BBB=M1(n)*(A1*[x(n);y(n)]+B1*v+sat+d)+M2(n)*(A2*[x(n);y(n)]+B2*v+sat+d);
    abcabc1(n)=BBB(1);
    abcabc2(n)=BBB(2);
    fff1=G.*abcabc1/gamma(p);
    fff2=G.*abcabc2/gamma(p);
    x(n+1)=x(1)+h*sum(fff1);
    y(n+1)=y(1)+h*sum(fff2);
    sBBB=sM1(n)*(A1*[sx(n);sy(n)]+B1*sv+ssat+d)+sM2(n)*(A2*[sx(n);sy(n)]+B2*sv+ssat+d);
    sabcabc1(n)=sBBB(1);
    sabcabc2(n)=sBBB(2);
    sfff1=G.*sabcabc1/gamma(p);
    sfff2=G.*sabcabc2/gamma(p);
    sx(n+1)=sx(1)+h*sum(sfff1);
    sy(n+1)=sy(1)+h*sum(sfff2);
    %du1(n)=sat(u1(n))-u1(n);
    %du2(n)=sat(u2(n))-u2(n);
    
    
end
%plot(t,x,'-',t,y,'--','LineWidth',1.8)
%xlabel('Time(second)')
%ylabel('State variables ')
%legend('x_1(t)','x_2(t)')
%legend('boxoff')

%plot(t(1:n),u1,'-',t(1:n),u2,'-.','LineWidth',1.8)
%xlabel('Time(second)')
%ylabel('Control inputs')
%legend('u_1(t)','u_2(t)')
%legend('boxoff')

%plot(t(1:n),satu1,'-',t(1:n),tanh(su1),'-.','LineWidth',1.8)
%xlabel('Time(second)')
%ylabel('sat(u_1)')
%legend('proposed saturation method','exiting saturation method in Ref. [53]')
%legend('boxoff')

%plot(t(1:n),satu2,'-',t(1:n),tanh(su2),'-.','LineWidth',1.8)
%xlabel('Time(second)')
%ylabel('sat(u_2)')
%legend('proposed saturation method','exiting saturation method in Ref. [53]')
%legend('boxoff')

%subplot(2,1,1)
%plot(t,x,'-',t,sx,'--','LineWidth',1.8)
%title('(a)')
%xlabel('Time(second)')
%ylabel('x_1(t) ')
%legend('proposed saturation method','exiting saturation method in Ref. [53]')
%legend('boxoff')
%subplot(2,1,2)
%plot(t,y,'-',t,sy,'--','LineWidth',1.8)
%title('(b)')
%xlabel('Time(second)')
%ylabel('x_2(t) ')
%legend('proposed saturation method','exiting saturation method in Ref. [53]')
%legend('boxoff')


%subplot(2,1,1)
%plot(t,ja1,'-',t,ja2,'--','LineWidth',1.8)
%xlabel('Time(second)')
%ylabel('Design parameters ')
%legend({'$\hat\kappa_1$','$\hat\kappa_2$'},'interpreter','latex')
%legend('boxoff')
%subplot(2,1,2)
%plot(t,jb1,'-',t,jb2,'--','LineWidth',1.8)
%xlabel('Time(second)')
%ylabel('Design parameters ')
%legend({'$\hat\lambda_1$','$\hat\lambda_2$'},'interpreter','latex')
%legend('boxoff')

%plot(t,jM,'-',t,jC,'--','LineWidth',1.8)
%xlabel('Time(second)')
%ylabel('Estimations of a and b ')
%legend({'$\hat a$','$\hat b$'},'interpreter','latex')
%legend('boxoff')

plot(t(1:n),du1,'-',t(1:n),du2,'-.','LineWidth',1.8)
xlabel('Time(second)')
ylabel('The difference \Deltau')
legend('\Deltau_1(t)','\Deltau_2(t)')
legend('boxoff')



 subplot(2,1,1)
plot(t,x,'-',t,sx,'--','LineWidth',1.8)
ylim([-5 10])
title('(a)')
xlabel('Time(second)')
ylabel('x_1(t) ')
legend('proposed saturation method','exiting saturation method in Ref. [53]')
legend('boxoff')
subplot(2,1,2)
plot(t,y,'-',t,sy,'--','LineWidth',1.8)
ylim([-2 6])
title('(b)')
xlabel('Time(second)')
ylabel('x_2(t) ')
legend('proposed saturation method','exiting saturation method in Ref. [53]')
legend('boxoff')


