clear;
clc;
N = 10;
R = 100;


betaT = 0.25+(0.75-0.25)*rand(N,1);  %user preference for T
betaE = 0.25+(0.75-0.25)*rand(N,1);	 %user preference for E
Fl = (0.5+(1.5-0.5)*rand(N,1))*1e9;
%L0 = R*rand(N,1);          %distance to base station
L0 = 0.8*R*ones(N,1);
v = 20*rand(N,1);
%v = 20*ones(N,1);
%theta = 180*rand(N,1);
theta = 180*ones(N,1);

tpro = 1;
C = 1000e6;
f0 = 20e9;
D = 420*1024*8; %bit
B = 20e6;
W = 1e6;
rou = 1;        %user payment
p = 10e-3*ones(N,1);
Tl = C./Fl;
El = 1e-20*Fl*C;
tao = rou*betaT;


m=5;    %% ant number
Alpha=1;  %% Alpha importance of pheromones
Beta=5;  %% Beta importance of Heuristic factor
Rho=0.1; %% Rho decrease factor of pheromones
NC_max=200; %% max iteration times
Q=1;         %%increasement of the pheromones?


%%step 1 initialization?
Tau=ones(N,2);     %Tau pheromones matrix 
NC=1;               %iteration counter
R_best=zeros(NC_max,N);       %best route at each iteration
L_best=inf.*ones(NC_max,1);   %best fitness at each iteration
L_ave=zeros(NC_max,1);        %average fitness at each iteration?


while NC<=NC_max        %stop when reach the iteration times­¢
	%%step 2 put m ants on different varaiables,each variable is a combine
	%%the cities
	Choice = ones(m,N);
	%%step 3 mth ant choose next city by possibility
    for j=1:N     
        for i=1:m
            %calculate the posibility of the unvisited city?
            for k=1:2       %1 represent the local ,2 represent the cloud
                P(k)=(Tau(j,k))^Alpha; 
            end
            P=P/(sum(P));
            %choose the next city based on possibility
            Pcum=cumsum(P);     
            Select=find(Pcum>=rand); 
			if Select(1) == 1
				Choice(i,j) = 0;
			else
				Choice(i,j) = 1;
			end
        end
    end
    %%step 4 record the best?
    L=zeros(m,1);     %init the fitness at this iteration
    for i=1:m
        L(i) = calcufit(Choice(i,:),Fl,betaT,betaE,L0,v,theta,N,tpro,C,f0,D,B,W,R,rou,p);      
    end
    L_best(NC)=max(L);           %choose the best fitness
    pos=find(L==L_best(NC));
    R_best(NC,:)=Choice(pos(1),:); %the best schedule
    L_ave(NC)=mean(L);           %average fitness?
    NC=NC+1;                      %continue the iteration
 
    
    %%step 5 update the pheromones
    Delta_Tau=zeros(N,2);        %init the delta pheromones matrix
    for i=1:m
        for j=1:N
			if Choice(i,j) == 0
				Delta_Tau(j,1)=Delta_Tau(j,1)+Q*L(i)/N;
			else
				Delta_Tau(j,2)=Delta_Tau(j,2)+Q*L(i)/N;
			end
        end
    end
    Tau=(1-Rho).*Tau+Delta_Tau; %consider the decrease of the pheromones?         
end

%%step 7 output the result
Pos=find(L_best==max(L_best)); %best route?
Best_Schedule=R_best(Pos(1),:); %best schedule?
Biggest_Length=L_best(Pos(1)); %best fitness

% enumerate
for i = 1:2^N
    bi = dec2bin(i-1);
    st = num2str(bi);
    for j = 1:length(st)
        meis(i,j) = (st(length(st)-j+1)=='1');
    end
end
% save the optional solution
enuFinalLocation = zeros(N,1); 
enuFinalUserUti = zeros(N,1);
enuFinalSysUti = 0;
[finalE,finalT] = deal(zeros(N,1));
fiSleave = [];


for i = 1:2^N
    location = meis(i,:);
    Sleave = [];
    if sum(location==1)>B/W
        continue;
    end
    %translate location to schedule
    meiS = find(location==1);
    % calculate fi according to schedule
    su = 0; 
    for j=1:length(meiS)
        su = su + sqrt(tao(meiS(j))*Fl(meiS(j)));
    end
    for j=1:length(meiS)
        k = meiS(j);
        fi(k) = sqrt(tao(k)*Fl(k))*f0/su;
    end
    %caculate the transtime 
    for i=1:N
        Ttrans(i) = transtime(L0(i),v(i),theta(i),tpro,D,W,p(i));
    end
    %caculate the user who finally leaves the cover
    for j=1:length(meiS)
        k = meiS(j);
        Ttol(k) = tpro + Ttrans(k) +  C/fi(k);
        Lf(k) = sqrt(L0(k)^2+(v(k)*Ttol(k))^2-2*L0(k)*v(k)*Ttol(k)*cosd(theta(k)));
        if Lf(k)>R
            Sleave(length(Sleave)+1)=k;
        end
    end 
    %calculate the result
    meiUserUti = zeros(N,1);
    meiSysUti = 0;
    for j=1:N
        if(location(j)==0)
            meiUserUti(j) = 0;
        else
            meiUserUti(j) =  betaT(j)*(Tl(j)-Ttrans(j)-C/fi(j))/Tl(j)...
                    + betaE(j)*(El(j)-p(j)*Ttrans(j))/El(j);
            %punish those who leave 
            Ttol(j) = tpro + Ttrans(j) +  C/fi(j);
            Lf(j) = sqrt(L0(j)^2+(v(j)*Ttol(j))^2-2*L0(j)*v(j)*Ttol(j)*cosd(theta(j)));
            if Lf(j)>R
                meiUserUti(j) = meiUserUti(j) - 0.00002*Lf(j)^2;
            end
			if(meiUserUti(j)<0)
				location(j)=0;
				meiUserUti(j) = 0;
			end
        end
        meiSysUti = meiSysUti + meiUserUti(j);
    end
    
    %compare with the temp optimal
    if meiSysUti > enuFinalSysUti
        enuFinalLocation = location; 
        enuFinalUserUti = meiUserUti;
        enuFinalSysUti = meiSysUti;     
        fiSleave =Sleave;
    end
end

function res = calcufit(location,Fl,betaT,betaE,L0,v,theta,N,tpro,C,f0,D,B,W,R,rou,p)
    Tl = C./Fl;
    El = 1e-20*Fl*C;
    tao = rou*betaT;

    %caculate the transtime 
    for i=1:N
        Ttrans(i) = transtime(L0(i),v(i),theta(i),tpro,D,W,p(i));
    end

	if sum(location==1)>B/W
        res = 0;
		return;
	end
	%translate location to schedule
	S = find(location==1);
	% calculate fi according to schedule
	su = 0; 
	if ~isempty(S)
		for j=1:length(S)
			su = su + sqrt(tao(S(j))*Fl(S(j)));
		end
		for j=1:length(S)
			k = S(j);
			fi(k) = sqrt(tao(k)*Fl(k))*f0/su;
		end
	end
	
	%calculate the result
	UserUti = zeros(N,1);
	SysUti = 0;

	for j=1:N
		if(location(j)==0)
			UserUti(j) = 0;
		else
			UserUti(j) =  betaT(j)*(Tl(j)-Ttrans(j)-C/fi(j))/Tl(j)...
				+ betaE(j)*(El(j)-p(j)*Ttrans(j))/El(j);
			%punish those who leave 
			Ttol(j) = tpro + Ttrans(j) +  C/fi(j);
			Lf(j) = sqrt(L0(j)^2+(v(j)*Ttol(j))^2-2*L0(j)*v(j)*Ttol(j)*cosd(theta(j)));
			if Lf(j)>R
				UserUti(j) = UserUti(j) - 0.00002*Lf(j)^2;
			end
			if(UserUti(j)<0)
				location(j)=0;
				UserUti(j) = 0;
			end
		end
		SysUti = SysUti + UserUti(j);
    end
    res = SysUti;
	return;
end
