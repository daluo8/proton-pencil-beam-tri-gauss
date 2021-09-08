clc
clear
close all
%% initialize
SSD=1000;%mm
BeamAngularSpread=1;%mrad
EnergySpread=0.0667;%percent
inisigma=3;%mm
addpath('/media/daluo8/Elements SE/proton data/topasfile/dosedata3D')
[centrax,centray]=deal(80.5,80.5);
[dimx,dimy,dimz]=deal(160,160,350);
energy_spread=2e-3/3;
%% read MC data
range=20;%cm

filename=strcat(num2str(range,'%.1f'),'cm.csv');
m=csvread(filename,8,0);
x=m(:,1);
y=m(:,2);
z=m(:,3);
doseMC=m(:,4);
doseMC=flipud(doseMC)./max(doseMC);
doseMC3D=zeros(dimx,dimy,dimz);
for p=1:length(x)
    ix=x(p)+1;
    iy=y(p)+1;
    iz=z(p)+1;
    doseMC3D(ix,iy,iz)=doseMC(p);
end
%% cal range
alpha=2.2e-3;
p=1.77;
IniEnergy=(range/alpha)^(1/p);
doseMC2D=squeeze(sum(doseMC3D,2));
dose1D=sum(doseMC2D,1);
dose1D_=dose1D/max(dose1D);
for i=0:dimz-1
    if dose1D_(dimz-i)>0.8
        idx=dimz-i;
        break
    end
end
range_1=idx+(1-0.8)/(1-dose1D_(idx+1));
RangeStruggle=sqrt((EnergySpread*IniEnergy/100)^2+(0.012*(range_1/10)^0.935)^2);
%% plot profile
depth=[20,120,150,152];
figure
for i=1:4
    doseslice=doseMC3D(:,:,depth(i));
    doseslice1D=(doseslice(80,:)+doseslice(81,:)+doseslice(:,80)'+doseslice(:,81)')/4;
    doseslice1D=doseslice1D/max(doseslice1D);
    subplot(2,2,i)
    
    semilogy(1:160,doseslice1D,'-')
    title(strcat('depth=',num2str(depth(i)),'mm'))
end
%% get fwhm and sigma
sigma=zeros(1,dimz);
for i=1:dimz
    doseslice=doseMC3D(:,:,i);
    doseslice1D=(doseslice(80,:)+doseslice(81,:)+doseslice(:,80)'+doseslice(:,81)')/4;
    doseslice1D=doseslice1D/max(doseslice1D);
    idx=find(doseslice1D>0.5);
    if idx(1)==1
        idxleft=1;
    else
        idxleft=idx(1)-1+(0.5-doseslice1D(idx(1)-1))/(doseslice1D(idx(1))-doseslice1D(idx(1)-1));
    end
    if idx(end)==160
        idxright=160;
    else
        idxright=idx(end)+(doseslice1D(idx(end))-0.5)/(doseslice1D(idx(end))-doseslice1D(idx(end)+1));
    end
    FWHM=idxright-idxleft;
    sigma(i)=FWHM/2.355;
end
figure
plot(1:dimz,sigma,'-')
xlabel('depth[mm]')
ylabel('sigma[mm]')
lat=1:160;
lat=lat-80;
figure
for i=1:4
    doseslice=doseMC3D(:,:,depth(i));
    doseslice1D=(doseslice(80,:)+doseslice(81,:)+doseslice(:,80)'+doseslice(:,81)')/4;
    doseslice1D=doseslice1D/max(doseslice1D);
    subplot(2,2,i)
    singlegauss=normpdf(1:160,80.5,sigma(depth(i)));
    singlegauss=singlegauss/max(singlegauss);
    semilogy(lat,singlegauss,'-')
    hold on
    scatter(lat,doseslice1D,'filled','^')
    title(strcat('depth=',num2str(depth(i)),'mm'))
%     legend('single Guassian','simulation')
    axis([lat(1) lat(end) 1e-6 1]);
    abs(max(doseslice1D-singlegauss))
end
legend('single Guassian','simulation','Orientation','horizontal',...
    'location',[0.13,0.05,0.74,0.05]);
%% get MCS sigma
angularsigma=zeros(1,dimz);
for t=1:dimz
    angularsigma(t)=(SSD+t)*BeamAngularSpread/1e3;
end
MCSsigma=sqrt(sigma.^2-inisigma^2-angularsigma.^2);
figure
plot(1:dimz,MCSsigma)
figure
plot(1:dimz,angularsigma)
%% get gamma index
[deltaDis,deltaDos]=deal(2,0.02);
GammaRange=range_1+5*RangeStruggle*10;
rmax=3.5*sigma;
doseS=zeros(dimy,dimz);
doseMCslice=squeeze(doseMC3D(80,:,:)+doseMC3D(81,:,:))+squeeze(doseMC3D(:,80,:)+doseMC3D(:,81,:));

doseMCslice=doseMCslice/max(max(doseMCslice));
dosecentra1D=sum(doseMCslice,1);
for i=1:dimz
    dose_lat=normpdf(1:160,80.5,sigma(i));
    doseS(:,i)=dose_lat/sum(dose_lat).*dosecentra1D(i);
end
doseS=doseS/max(max(doseS));
figure
contourf(doseS)
figure
contourf(doseMCslice)
figure
contourf(doseS-doseMCslice)
colorbar
gamma_rateS=getGammaindex(doseS,doseMCslice,deltaDis,deltaDos,GammaRange,rmax)
%% double gauss fit
Dsigma1=sigma;
Dsigma2=3*sigma;
NUCweight=zeros(1,dimz);
MCSsf=zeros(1,dimz);
NUCsf=zeros(1,dimz);
tic
for i=1:dimz
   doseslice=doseMC3D(:,:,i);
   doseslice1D=(doseslice(80,:)+doseslice(81,:)+doseslice(:,80)'+doseslice(:,81)')/4;
   doseslice1D=doseslice1D/max(doseslice1D); 
   fun=@(x)sum((((1-x(1))*normpdf(1:160,80.5,x(2)*Dsigma1(i))+...
       x(1)*normpdf(1:160,80.5,x(3)*Dsigma2(i)))/...
       max((1-x(1))*normpdf(1:160,80.5,x(2)*Dsigma1(i))+...
       x(1)*normpdf(1:160,80.5,x(3)*Dsigma2(i)))-doseslice1D).^2);
   x0=[0.1,1,1];
   vlb=[0.0,0.1,0.5];
   vub=[0.3,5,5];
   A=[];b=[];
   Aeq=[];beq=[];
   options=optimoptions(@fmincon,'MaxFunEvals',10000,'MaxIterations',10000);
   [x,fval,exitflag]=fmincon(fun,x0,A,b,Aeq,beq,vlb,vub,'',options);
   NUCweight(i)=x(1);
   MCSsf(i)=x(2);
   NUCsf(i)=x(3);
end
toc
%% plot fit result
figure
plot(1:dimz,NUCweight,'-')
xlabel('depth[mm]')
ylabel('NUCweight')
figure
plot(1:dimz,Dsigma1.*MCSsf,'-')
hold on
xlabel('depth[mm]')
ylabel('sigma[mm]')
plot(1:dimz,Dsigma2.*NUCsf,'-')
legend('MCSsigma','NUCsigma')
figure
plot(1:dimz,MCSsf)
hold on
plot(1:dimz,NUCsf,'--')
%% plot double gaussian profile
figure
for i=1:length(depth)
    doseslice=doseMC3D(:,:,depth(i));
    doseslice1D=(doseslice(80,:)+doseslice(81,:)+doseslice(:,80)'+doseslice(:,81)')/4;
    doseslice1D=doseslice1D/max(doseslice1D);
    subplot(2,2,i)
    doublegauss=(1-NUCweight(depth(i)))*normpdf(1:160,80.5,MCSsf(depth(i))*Dsigma1(depth(i)))+...
        NUCweight(depth(i))*normpdf(1:160,80.5,NUCsf(depth(i))*Dsigma2(depth(i)));
    normterm=max(doublegauss);
    doublegauss=doublegauss/normterm;
    semilogy(lat,doublegauss,'-')
    hold on
    scatter(lat,doseslice1D,'filled','^')
    
    gauss1st=(1-NUCweight(depth(i)))*normpdf(1:160,80.5,MCSsf(depth(i))*Dsigma1(depth(i)))/normterm;
    gauss2nd=NUCweight(depth(i))*normpdf(1:160,80.5,NUCsf(depth(i))*Dsigma2(depth(i)))/normterm;
    plot(lat,gauss1st,'--')
    plot(lat,gauss2nd,'--')
%     legend('double Guassian','simulation','1stgauss','2ndgauss')
    title(strcat('depth=',num2str(depth(i)),'mm'))
    axis([lat(1) lat(end) 1e-6 1]);
    abs(max(doseslice1D-doublegauss))
end
legend('double Guassian','simulation','1stgauss','2ndgauss','Orientation','horizontal',...
    'location',[0.13,0.05,0.74,0.05]);
%% get gamma index
[deltaDis,deltaDos]=deal(2,0.02);
doseD=zeros(dimy,dimz);
for i=1:dimz
    gauss1st=(1-NUCweight(i))*normpdf(1:160,80.5,MCSsf(i)*Dsigma1(i));
    gauss2nd=NUCweight(i)*normpdf(1:160,80.5,NUCsf(i)*Dsigma2(i));
    dose_lat=gauss1st+gauss2nd;
%     dose_lat=dose_lat/max(dose_lat);
    doseD(:,i)=dose_lat/sum(dose_lat).*dosecentra1D(i);
end
doseD=doseD/max(max(doseD));
figure
contourf(doseD)
figure
contourf(doseMCslice)
figure
contourf(doseD-doseMCslice)
colorbar
gamma_rateD=getGammaindex(doseD,doseMCslice,deltaDis,deltaDos,GammaRange,rmax)
%% tri gaussian fit
Tsigma1=sigma;
Tsigma2=sigma*3;
Tsigma3=sigma*5;
weight2=zeros(1,dimz);
weight3=zeros(1,dimz);
sf1=zeros(1,dimz);
sf2=zeros(1,dimz);
sf3=zeros(1,dimz);
% sigma fit
tic
for i=1:dimz
   doseslice=doseMC3D(:,:,i);
   doseslice1D=(doseslice(80,:)+doseslice(81,:)+doseslice(:,80)'+doseslice(:,81)')/4;
   doseslice1D=doseslice1D/max(doseslice1D); 
%    fun=@(x)sum(abs(log10((normpdf(1:160,80.5,x(3)*Tsigma1(i))+...
%        x(1)*0.1*normpdf(1:160,80.5,x(4)*Tsigma2(i))+...
%        x(2)*0.1*normpdf(1:160,80.5,x(5)*Tsigma3(i)))/max...
%        (normpdf(1:160,80.5,x(3)*Tsigma1(i))+...
%        x(1)*0.1*normpdf(1:160,80.5,x(4)*Tsigma2(i))+...
%        x(2)*0.1*normpdf(1:160,80.5,x(5)*Tsigma3(i))))-log10(doseslice1D)));
   fun=@(x)sum(abs(log10(((1-exp(x(1))-exp(x(2)))*normpdf(1:160,80.5,x(3)*Tsigma1(i))+...
       exp(x(1))*normpdf(1:160,80.5,x(4)*Tsigma2(i))+...
       exp(x(2))*normpdf(1:160,80.5,x(5)*Tsigma3(i)))/max...
       ((1-exp(x(1))-exp(x(2)))*normpdf(1:160,80.5,x(3)*Tsigma1(i))+...
       exp(x(1))*normpdf(1:160,80.5,x(4)*Tsigma2(i))+...
       exp(x(2))*normpdf(1:160,80.5,x(5)*Tsigma3(i))))-log10(doseslice1D)));
   x0=[-2,-2,1,1,1];
   vlb=[-20,-20,0.5,0.5,0.1];
   vub=[-1,-1,2,10,10];
   A=[];b=[];
   Aeq=[];beq=[];
   options=optimoptions(@fmincon,'MaxFunEvals',10000,'MaxIterations',10000);
   [x,fval,exitflag]=fmincon(fun,x0,A,b,Aeq,beq,vlb,vub,'',options);
   weight2(i)=x(1);
   weight3(i)=x(2);
   sf1(i)=x(3);
   sf2(i)=x(4);
   sf3(i)=x(5);
end
toc
%% plot fit result
figure
semilogy(1:dimz,exp(weight2),'-')
hold on
plot(1:dimz,exp(weight3),'-')
xlabel('depth[mm]')
ylabel('weight')
legend('weight 2nd','weight 3rd')
figure
plot(1:dimz,Tsigma1.*sf1,'-')
hold on
xlabel('depth[mm]')
ylabel('sigma[mm]')
plot(1:dimz,Tsigma2.*sf2,'--')
plot(1:dimz,Tsigma3.*sf3,'-')
legend('sigma 1st','sigma 2nd','sigma 3rd')
%% plot triple gaussian profile
figure
for i=1:length(depth)
    doseslice=doseMC3D(:,:,depth(i));
    doseslice1D=(doseslice(80,:)+doseslice(81,:)+doseslice(:,80)'+doseslice(:,81)')/4;
    doseslice1D=doseslice1D/max(doseslice1D);
    subplot(2,2,i)
    triplegauss=(1-exp(weight2(depth(i)))-exp(weight3(depth(i))))*normpdf(1:160,80.5,sf1(depth(i))*Tsigma1(depth(i)))+...
        exp(weight2(depth(i)))*normpdf(1:160,80.5,sf2(depth(i))*Tsigma2(depth(i)))+...
        exp(weight3(depth(i)))*normpdf(1:160,80.5,sf3(depth(i))*Tsigma3(depth(i)));
    normterm=max(triplegauss);
    triplegauss=triplegauss/normterm;
    semilogy(1:160,triplegauss,'-')
    hold on
    scatter(1:160,doseslice1D,'filled','^')
    
    gauss1st=(1-exp(weight2(depth(i)))-exp(weight3(depth(i))))*normpdf(1:160,80.5,sf1(depth(i))*Tsigma1(depth(i)))/normterm;
    gauss2nd=exp(weight2(depth(i)))*normpdf(1:160,80.5,sf2(depth(i))*Tsigma2(depth(i)))/normterm;
    gauss3rd=exp(weight3(depth(i)))*normpdf(1:160,80.5,sf3(depth(i))*Tsigma3(depth(i)))/normterm;
    plot(1:160,gauss1st,'--')
    plot(1:160,gauss2nd,'--')
    plot(1:160,gauss3rd,'--')
    title(strcat('depth=',num2str(depth(i)),'mm'))
%     legend('triple Guassian','simulation','1stgauss','2ndgauss','3rdgauss')
    axis([1 160 1e-6 1]);
    abs(max(doseslice1D-triplegauss))
end
legend('triple Guassian','simulation','1stgauss','2ndgauss','3rdgauss','Orientation','horizontal',...
    'location',[0.13,0.05,0.74,0.05]);
%% get gamma index
[deltaDis,deltaDos]=deal(1,0.01);
doseT=zeros(160,dimz);
for i=1:dimz
    gauss1st=(1-exp(weight2(i))-exp(weight3(i)))*normpdf(1:160,80.5,sf1(i)*Tsigma1(i));
    gauss2nd=exp(weight2(i))*normpdf(1:160,80.5,sf2(i)*Tsigma2(i));
    gauss3rd=exp(weight3(i))*normpdf(1:160,80.5,sf3(i)*Tsigma3(i));
    dose_lat=gauss1st+gauss2nd+gauss3rd;
%     dose_lat=dose_lat/sum(dose_lat);
    doseT(:,i)=dose_lat/sum(dose_lat).*dosecentra1D(i);
end
doseT=doseT/max(max(doseT));
figure
contourf(doseT-doseMCslice)
colorbar
gamma_rateT=getGammaindex(doseT,doseMCslice,deltaDis,deltaDos,GammaRange,rmax)
%%
slicez=2;
doseMC1D=doseMCslice(:,slicez);
doseT1D=doseT(:,slicez);
figure
semilogy(1:160,doseMC1D,'-')
hold on
semilogy(1:160,doseT1D,'--')
legend('MC sim','triGauss')
max(abs(doseT1D-doseMC1D))
%% define the phantom
BoneZTrans=100;%mm
phantom=ones(dimy,dimz);
phantom(:,BoneZTrans-15+1:BoneZTrans+15)=phantom(:,BoneZTrans-15+1:BoneZTrans+15)*1.95;
figure
contourf(phantom)
colorbar
%% cal dose to phantom
doseTEQ=zeros(size(phantom));
depEQ=0;
% depdose=zeros(1,dimz);
for i=1:dimz
    Rho=phantom(80,i);
    depEQ=depEQ+Rho;
    if depEQ>dimz
        break
    end
    depdose=interp1(1:dimz,dosecentra1D,depEQ);
    weight2nd=interp1(1:dimz,weight2,depEQ);
    weight3rd=interp1(1:dimz,weight3,depEQ);
    sf1EQ=interp1(1:dimz,sf1,depEQ);
    sf2EQ=interp1(1:dimz,sf2,depEQ);
    sf3EQ=interp1(1:dimz,sf3,depEQ);
    Tsigma1EQ=interp1(1:dimz,Tsigma1,depEQ);
    Tsigma2EQ=interp1(1:dimz,Tsigma2,depEQ);
    Tsigma3EQ=interp1(1:dimz,Tsigma3,depEQ);
    lat_arr=((1:dimy)-80.5)+80.5;
    gauss1st=(1-exp(weight2nd)-exp(weight3rd))*normpdf(lat_arr,80.5,sf1EQ*Tsigma1EQ);
    gauss2nd=exp(weight2nd)*normpdf(lat_arr,80.5,sf2EQ*Tsigma2EQ);
    gauss3rd=exp(weight3rd)*normpdf(lat_arr,80.5,sf3EQ*Tsigma3EQ);
    dose_lat=gauss1st+gauss2nd+gauss3rd;
%     dose_lat=dose_lat/sum(dose_lat);
    doseTEQ(:,i)=dose_lat/sum(dose_lat).*depdose;
end
figure
doseTEQ=doseTEQ/max(max(doseTEQ));
contourf(doseTEQ)
colorbar
dose1DEQ=sum(doseTEQ,1);
dose1DEQ=dose1DEQ/max(dose1DEQ);
figure
plot(1:dimz,dose1DEQ)
%%
