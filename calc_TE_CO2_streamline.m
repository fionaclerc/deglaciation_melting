function [CS,CL,fn,P,l_int,dFdP_int,dFdl_int] = calc_TE_CO2_streamline(DC,dFm,Fcrit0,Dmat,XX,YY,t_ind,XY,XCO2,DF)
% input:
% DC (table of partition coefficients; rows are trace elements; columns describe partitioning with mineral phases:
%      [ cpx, spinel, opx, olivine, garnet] are columns [8 9 10 11 14], respectively,
%      column 13 is initial concentration in mantle)
% dFm (step size, of melt fraction)
% Fcrit0 (retained melt fraction, 0 is fractional limit, 1 is batch limit)
% Dmat (structure containing numerical model output, with fields "P" (pressure)
        %        and "F" (melt fraction) as 2-D grids; can be indexed according to timestep; see t_ind)
% XX (grid of horizontal coordinate, generated using meshgrid, same size as pressure, melt fraction fields.)
% YY (corresponding grid of vertical coordinate)
% t_ind (index of desired timestep for Dmat structure)
% XY (first column is x, second column is y, for location of 1-D streamline along which calculation performed)
% XCO2 (source concentration of CO2)
% DF (same as Dmat only contains dFdP field; derivative od melt fraction w.r.t. pressure)

% output: 1-D columns/profiles for each element.
% CS: concentration in solid mantle/residue
% CL: concentration in melt
% fn: melt fraction (degree of melting)
% P: pressure
% l_int: distance along profile
% dFdP_int: dFdP along profile
% dFdl_int: dFdl along profile



dP=.1;
V = [];

MON91 = 0;
%comp = 0; % set to 1 to compare to shaw figures
inpt = 1;
Pvar = 1;
DLa = [.0536 0 .0004 .0001];
ll = sqrt((XY(:,1) - XY(1,1)).^2 + (XY(:,2) - XY(1,2)).^2);

dFdP_int = zeros([1 1e5]);
dFdl_int = zeros([1 1e5]);

nstep = 7;


dFdP=[];
if inpt==1

    Pmat = interp2(XX,YY,Dmat(t_ind).P,XY(:,1),XY(:,2))*1e-9;
    P = Pmat;
    l_int = zeros([1 1e5]);

    dFdP = -interp2(XX,YY,DF.dFdP,XY(:,1),XY(:,2))*1e9;%1/GPa
    
    if dFdP(1)~=0
       fprintf('\n');
       warning('dFdP at base is not equal to zero.');
       fprintf('\n');
       dFdP(1)=0;
    end
    
    dFdP0=dFdP;
  
    dFdl= streamline_gradient(XX,YY,Dmat(t_ind).F,XY);
end

if Pvar ==1
    P = P(1);
end

gt_cutoff = 3;
DC.data(29,1:4) = DC.data(3,1:4);% barium
DC.data(29,14) = DC.data(29,9);
DC.data(29,13) = XCO2;

Fn = zeros([1 1e5]);
Fn(1) = dFdP(1);


W0=1;

dFdP_sil = .05;dFdP_car = .005;


dF_sil = dFdP_sil.*dP;
dF_car = dFdP_car.*dP;


if inpt==0
    Finit = dF_car;Fn(1) = Finit;
else
    Finit = Fn(1);
end



init = 100;


if MON91
    X0 = [.076,.115, .211,.598];
else
    X0 = [.18,.16,.05,.61];
end


c2s = 1;

D=[];X=[];Cl=[];Cl_avg = [];Cs = [];
n=1;


YMS(1,n) = 0.0016*(P(n)*10)^2 + 0.0006*(P(n)*10) + 0.5; %cpx
YMS(3,n) = -0.0022*(P(n)*10)^2 - 0.0223*(P(n)*10) + 0.9365;%opx
YMS(4,n) = 0.0015*(P(n)*10)^2 - 0.008*(P(n)*10) - 0.3414;%ol
YMS(2,n) = 0.0024*(P(n)*10) + 0.0738;%sp

YMS(:,n) = YMS(:,n)./sum(YMS(:,n)); 
X(1,:) = X0;

D(1,:) = sum(DC.data(:,8:11)'.*meshgrid(X(1,:),[1:length(DC.data)])')';


%%
WS(n) = (1-Fn(n)).*W0;
Wl(n) = Fn(n).*W0;
WL(n) = Wl(n);

Fcrit(n) = Fcrit0;

L(n) = Fn(1).*W0; 

if Fn(n)>0
    Q(n) = min((Fcrit(n)./(1-Fcrit(n))).*WS(n)./(Wl(n)),1);

    WLE(n) = (1-Q(n)).*WL(n);
    WLR(n) = Q(n).*WL(n);

    cln(n,:) = DC.data(:,13)'./(D(n,:).*WS(n)+Fn(1).*W0);
    WR(n) = (1-Fn(n)+Q(n).*Fn(n)).*W0;
    Lr(n) = Q(n).*Fn(n).*W0;
    crn(n,:) = cln(n,:).*(D(n,:).*WS(n)+Lr(n))./WR(n);

    CL(n,:) = DC.data(:,13)'./(D(n,:).*WS(n)+WL(n));    
    CR(n,:) = CL(n,:).*(D(n,:).*WS(n)+WLR(n))./WR(n);
    CLA(n,:) = CL(1:n,:).*meshgrid(WLE(1:n),[1:29])'./(WLE(1:n));
else
    Q(n) = min((Fcrit(n)./(1-Fcrit(n))).*WS(n)./(Wl(n)),1);

    WLE(n) = 0;
    WLR(n) = 0;
    cln(n,:) = DC.data(:,13)'.*0;
    WR(n) = W0;
    Lr(n) = 0;
    crn(n,:) = cln(n,:);

    CL(n,:) = DC.data(:,13)'./(D(n,:).*WS(n)+WL(n));
    CR(n,:) = DC.data(:,13);
    CLA(n,:) = CL(1:n,:).*meshgrid(WLE(1:n),[1:29])'./(WLE(1:n));
end



fn(n) = Fn(1);
use_dF_incr = 1;

while l_int(n)<=ll(end)

    n = n+1;
        dFdP_int(n) = interp1(Pmat,dFdP,P(n-1));
        dFdl_int(n) = interp1(ll,dFdl,l_int(n-1));
        if dFdl_int(n)>0 & dFdP_int(n)>0
            if use_dF_incr
                Fn(n) = dFm;
                dl = dFm./dFdl_int(n);
                if dl > 1000
                    dl=1000;
                    Fn(n) = dFdl_int(n).*dl;
                    
                end
                dP = Fn(n)./dFdP_int(n);
                P(n) = P(n-1)-dP;
                
                l_int(n) = l_int(n-1) + dl;
            else
                dP = dFm;% know when you input this that it will be P in GPa
                Fn(n) = dP.*dFdP_int(n);
                P(n) = P(n-1)-dP;
            end
        else
            Fn(n) = 0;
            dl = 1000;
            l_int(n) = l_int(n-1) + dl;
            P(n) = interp1(ll,Pmat,l_int(n));

        end
        


    fn(n) = nansum(Fn(1:n));
    
    if comp==0
    YMS(1,n) = 0.0016*(P(n)*10)^2 + 0.0006*(P(n)*10) + 0.5; %cpx
    YMS(3,n) = -0.0022*(P(n)*10)^2 - 0.0223*(P(n)*10) + 0.9365;%opx
    YMS(4,n) = 0.0015*(P(n)*10)^2 - 0.008*(P(n)*10) - 0.3414;%ol
    
    YMS(2,n) = 0.0024*(P(n)*10) + 0.0738;%sp
    
    
    else
        YMS(1,n) = .5;
        YMS(2,n)  = 0;
        YMS(3,n) = .3;
        YMS(4,n) = .2;
    end
    YMS(:,n) = YMS(:,n)./sum(YMS(:,n));  
    X(n,:) = max((X(n-1,:)-(Fn(n)'.*YMS(:,n)'))/(1-Fn(n)),0);
   
            %silicate melting
    if P(n) > gt_cutoff
        % garnet field
       D(n,:) = sum(DC.data(:,[8 14 10 11])'.*meshgrid(X(n,:),[1:length(DC.data)])')'; 

    else
        % spinel field 
        D(n,:) = sum(DC.data(:,[8:11])'.*meshgrid(X(n,:),[1:length(DC.data)])')';
    end
 
    
    if comp==1
         D(n,:) = sum(meshgrid(DLa,[1:29])'.*meshgrid(X(n,:),[1:length(DC.data)])')';
         
         X(n,:) = (X0 - (n.*YMS(:,n)'.*Fn(1)))./(1-(n.*Fn(1)));
         
         D00(n+1) = nansum(X(n,:).*DLa);
         D(n,:) = (D00(n) - (n.*.0268.*Fn(1)))./(1-(n.*Fn(1)));
    end

    
    
    Wl(n) = Fn(n).*WS(n-1);
    WS(n) = WS(n-1)-Wl(n);
    
    
    WL(n) = nansum([Wl(n) WLR(n-1)]);
    
     Fcrit(n) = Fcrit0;
    
    if Fn(n)>0
     
    Q(n) = min((Fcrit(n-1)./(1-Fcrit(n-1))).*WS(n)./(WL(n)),1);
    if comp==1
    Q(n) = 0.2;
    end
    WLR(n) = Q(n).*WL(n);
   
    WLE(n) = (1-Q(n)).*WL(n);
    
    
    L(n) = (Fn(1).*W0)+Lr(n-1);
    Lr(n) = Q(n).*L(n);
    cln(n,:) = crn(n-1,:).*WR(n-1)./(D(n,:).*WS(n)+L(n));
    CL(n,:) = CR(n-1,:).*WR(n-1)./(D(n,:).*WS(n)+WL(n));
        
    WR(n) = WLR(n)+WS(n);
    crn(n,:) = cln(n-1,:).*(WR(n-1)./WR(n)).*(D(n,:).*WS(n)+Lr(n))./(D(n,:).*WS(n)+L(n));
    CR(n,:) = CL(n,:).*(D(n,:).*WS(n)+WLR(n))./WR(n);
    CLA(n,:) = nansum(CL(1:n,:).*meshgrid(WLE(1:n),[1:29])')./nansum(meshgrid(WLE(1:n),[1:29])');
    else
        
    Q(n) = 0;WLR(n) = Q(n).*WL(n);
        WLE(n) = (1-Q(n)).*WL(n);

     L(n) = (Fn(1).*W0)+Lr(n-1);
    Lr(n) = Q(n).*L(n);    cln(n,:) = crn(n-1,:).*WR(n-1)./(D(n,:).*WS(n)+L(n));
    CL(n,:) = CR(n-1,:).*WR(n-1)./(D(n,:).*WS(n)+WL(n));
    
    WR(n) = WLR(n)+WS(n);
    crn(n,:) = cln(n-1,:).*(WR(n-1)./WR(n)).*(D(n,:).*WS(n)+Lr(n))./(D(n,:).*WS(n)+L(n));
    CR(n,:) = CL(n,:).*(D(n,:).*WS(n)+WLR(n))./WR(n);
    CLA(n,:) = sum(CL(1:n,:).*meshgrid(nanmax(WLE(1:n),0),[1:29])')./nansum(meshgrid(WLE(1:n),[1:29])');
    end
end
dFdP_int= dFdP_int(1:n); 
dFdl_int= dFdl_int(1:n); 
l_int= l_int(1:n); 

CS = D.*CL;
