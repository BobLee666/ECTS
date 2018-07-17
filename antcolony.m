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



%caculate the transtime 
for i=1:N
	Ttrans(i) = transtime(L0(i),v(i),theta(i),tpro,D,W,p(i));
end
% enumerate
i = 2^N;
location = zeros(N,1);
bi = dec2bin(i-1);
st = num2str(bi);
for j = 1:length(st)
    location(j) = (st(length(st)-j+1)=='1');
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
% calculate myf according to schedule
An = zeros(N,1);
Bn = zeros(N,1);

for j=1:length(meiS)
    t = meiS(j);
    An(t) = (rou*betaT(t)*C/Tl(t)) + 2*(tpro+Ttrans(t))*C;
    Bn(t) = C^2;
end
lam = callambda(meiS,An,Bn,f0);
for j=1:length(meiS)
    t = meiS(j);
    strA = num2str(An(t));
    strB = num2str(Bn(t));
    stra = num2str(-lam(t));
    equal = strcat([stra,'+',strA,'/(x^2)+2*',strB,'/(x^3)=0']);
    x = solve(equal,'x');

    y = x(3);
    myf(t) = double(y);
end


m=5;    %% ant number
Alpha=1;  %% Alpha importance of pheromones
Beta=5;  %% Beta importance of Heuristic factor
Rho=0.1; %% Rho decrease factor of pheromones
NC_max=200; %% max iteration times
Q=1;         %%increasement of the pheromones?


%%step 1 initialization?
Tau=ones(N);     %Tau pheromones matrix 
Tabu=zeros(m,N);   %store the route
NC=1;               %iteration counter
R_best=zeros(NC_max,N);       %best route at each iteration
L_best=inf.*ones(NC_max,1);   %best fitness at each iteration
L_ave=zeros(NC_max,1);        %average fitness at each iteration?


while NC<=NC_max        %stop when reach the iteration times­¢
	%%step 2 put m ants on different varaiables,each variable is a combine
	%%the cities
	Randpos = round(rand(m,N));
	Tabu(:,1) = Randpos(:,1);
	Choice = Randpos;
	%%step 3 mth ant choose next city by possibility
    for j=1:n     
        for i=1:m
            visited=Tabu(i,1:(j-1)); %record the visited cities ?
            J=zeros(1,(n-j+1));       %the unvisited city 
            P=J;                      %init the posibility of the unvisited city
            Jc=1;
            for k=1:n
                if isempty(find(visited==k,1))   %add the unvisited cities
                    J(Jc)=k;
                    Jc=Jc+1;                         
                end
            end
            %calculate the posibility of the unvisited city?
            for k=1:length(J)
				pop(1,:) = Choice(i,:);
				pop(2,:) = Choice(i,:);
				pop(2,J(k)) = ~Choice(i,J(k));
				Eta = nij(pop(2,:),Fl,betaT,betaE,L0,v,theta,N,tpro,C,f0,D,B,W,R,rou,p) ...
				-nij(pop(1,:),Fl,betaT,betaE,L0,v,theta,N,tpro,C,f0,D,B,W,R,rou,p);
				if Eta < 0
					Eta = 0;
				end
                P(k)=(Tau(J(k)))^Alpha*(Eta^Beta); %Eta is the diff between the two schedule?
            end
            P=P/(sum(P));
            %choose the next city based on possibility
            Pcum=cumsum(P);     
            Select=find(Pcum>=rand); 
            to_visit=J(Select(1));
            Tabu(i,j)=to_visit;
			Choice(i,J(Select(1))) = ~Choice(i,J(Select(1)));
        end
    end
    %%step 4 record the best?
    L=zeros(m,1);     %init the fitness at this iteration
    for i=1:m
        L(i) = nij(Choice(i,:),Fl,betaT,betaE,L0,v,theta,N,tpro,C,f0,D,B,W,R,rou,p);      ¦»
    end
    L_best(NC)=max(L);           %choose the best fitness
    pos=find(L==L_best(NC));
    R_best(NC,:)=Choice(pos(1),:); %the best schedule
    L_ave(NC)=mean(L);           %average fitness?
    NC=NC+1;                      %continue the iteration
 
    
    %%step 5 update the pheromones
    Delta_Tau=zeros(N);        %init the delta pheromones matrix
    for i=1:m
        for j=1:N
            Delta_Tau(j)=Delta_Tau(j)+Q*Choice(i,j);
        end
    end
    Tau=(1-Rho).*Tau+Delta_Tau; %consider the decrease of the pheromones?
    %%step 5 clear the history route?
    Tabu=zeros(m,n);            
end

%%step 7 output the result
Pos=find(L_best==max(L_best)); %best route?
Best_Schedule=R_best(Pos(1),:); %best schedule?
Biggest_Length=L_best(Pos(1)); %best fitness

function res = nij(location,Fl,betaT,betaE,L0,v,theta,N,tpro,C,f0,D,B,W,R,rou,p)
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
