%% clear and close Command Window : 

clear;
close all;
clc;

%% Input : 

n0i = input('Insert n0i Values:  \n');
MuPLi = input('Insert MuPLi Values : \n');
MuFi = input('Insert MuFi Values : \n');
Ispi = input('Insert Ispi Values : \n');
incl = input('Insert inpli (i) Value : \n');
Ha = input('Insert HA Value : \n');
Hp = input('Insert Hp Value : \n');
Phi = input('Insert Phi Value: \n');
Mpl = input('Insert Mpl Value : \n');

%% Algorithm : 


Grid = input('Insert Grid of MuPLi : \n');
Gridd = input('Insert Grid Of Ispi : \n');
GridMuf = input('Insert Grid Of MuFi : \n');

delup = input('Insert Delta UP Value : \n');
deldown = input('Insert Delta Down Value : \n');
delmufup = input('Insert Delta MuFi UP : \n');
delmufdown = input('Insert Delta MuFi Down : \n');
delIspup = input('Insert Delta Ispi Up : \n');
delIspdown = input('Insert Delta Ispi Down : \n');


L = length(MuPLi);
Mu_pl = prod(MuPLi);
M0 = [];
M0 = [M0 Mpl/Mu_pl];

for i=1:L
    M0(i+1) = MuPLi(i).*M0(i);
end


Epoch = 0;
Itaration = 0;
while true 
Itaration = Itaration + 1;
V_av=0;
g0=9.81;
for i=1:L
    V_av=V_av-(Ispi(i)*log(MuFi(i))*g0*0.75)*0.001;
end

Mu=398620;
RE=6378;
WE=7.2685e-05;

ra=Ha+RE;
rp=Hp+RE;
a=(ra+rp)/2;
     Vc=sqrt((2*Mu/rp)-(Mu/a));
     Phi=deg2rad(Phi);
     incl=deg2rad(incl);
     V_Earth=RE*WE*cos(Phi);
     Vf=sqrt((abs(Vc))^2+(abs(V_Earth))^2-2*(abs(Vc))*(abs(V_Earth))*cos(incl));
     


 if (V_av >= Vf) 
     break;
 else
     Logic=true;
     Epoch=Epoch+1;
      for i=1:L
         MuFi(i)=MuFi(i)-(delmufup(i)-delmufdown(i))/GridMuf;
         if(MuFi(i)<delmufdown(i))
             Logic=false;
         end
      end
  if   (Logic)
      continue;
  else

         for i=1:L
         MuPLi(i)=MuPLi(i)+(delup(i)-deldown(i))/Grid;
         if (MuPLi(i)>delup(i))
           disp('PLEASE Change ISP Range : \n');
              for j=1:L
                 Ispi(j)=Ispi(j)+(delIspup(j)-delIspdown(j))/Gridd;
                 if (Ispi(j)>delIspup(j))
                 disp('Please increace The Number Of Stage');
                 end
              end 
         
         end
         
         end 
         
         
   mupltot=1;
   for i=1:length(MuPLi)
       mupltot=mupltot*MuPLi(i);
   end

   M0(1)=Mpl/mupltot;

   for i=1:L
   M0(i+1)=MuPLi(i).*M0(i);
   end
   
   
  end
 
     
 end
 
 
      
end
%% Results : 
disp(V_av);
[m,n] = size(MuFi);
thi =zeros(m,n);
msti = zeros(m,n);
mp = zeros(m,n);
tbi = zeros(m,n);
mdoti = zeros(m,n);
     for i=1:length(Ispi)
          mp(i)=M0(i)*(1-MuFi(i));
          msti(i)=MuFi(i)*M0(i)-M0(i+1);
          thi(i)=n0i(i)*M0(i);
          tbi(i)=Ispi(i)*(1-MuFi(i))/n0i(i);
          mdoti(i)=thi(i)/Ispi(i);
     end

fprintf('MP =%f \n',mp);
fprintf('msti =%f \n',msti);
fprintf('thi =%f \n',thi);
fprintf('tbi =%f \n',tbi);
fprintf('mdoti =%f \n',mdoti);




