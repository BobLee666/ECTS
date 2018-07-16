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


m=5;    %% m 蚂蚁个数
Alpha=1;  %% Alpha 表征信息素重要程度的参数
Beta=5;  %% Beta 表征启发式因子重要程度的参数
Rho=0.1; %% Rho 信息素蒸发系数
NC_max=200; %%最大迭代次数
Q=1;         %%信息素增加强度系数


%%第一步：变量初始化
Tau=ones(N);     %Tau为信息素矩阵
Tabu=zeros(m,N);   %存储并记录路径的生成
NC=1;               %迭代计数器，记录迭代次数
R_best=zeros(NC_max,N);       %各代最佳路线
L_best=inf.*ones(NC_max,1);   %各代最佳路线的长度
L_ave=zeros(NC_max,1);        %各代路线的平均长度(长度=适应度)


while NC<=NC_max        %停止条件之一：达到最大迭代次数，停止
	%%第二步：将m只蚂蚁放到不同变量上，每个变量是城市的组合，随机0-1变量矩阵
	Randpos = round(rand(m,N));
	Tabu(:,1) = Randpos(:,1);
	Choice = Randpos;
	%%第三步：m只蚂蚁按概率函数选择下一座城市，完成各自的周游

    for j=1:n     
        for i=1:m
            visited=Tabu(i,1:(j-1)); %记录已访问的城市，避免重复访问
            J=zeros(1,(n-j+1));       %待访问的城市
            P=J;                      %待访问城市的选择概率分布
            Jc=1;
            for k=1:n
                if length(find(visited==k))==0   %开始时置0
                    J(Jc)=k;
                    Jc=Jc+1;                         %访问的城市个数自加1
                end
            end
            %下面计算待选城市的概率分布
            for k=1:length(J)
				pop(1,:) = Choice(i,:);
				pop(2,:) = Choice(i,:);
				pop(2,J(k)) = ~Choice(i,J(k));
				Eta = nij(pop(2,:),Fl,betaT,betaE,L0,v,theta,N,tpro,C,f0,D,B,W,R,rou,p) ...
				-nij(pop(1,:),Fl,betaT,betaE,L0,v,theta,N,tpro,C,f0,D,B,W,R,rou,p);
				if Eta < 0
					Eta = 0;
				end
                P(k)=(Tau(J(k)))^Alpha)*(Eta^Beta); %Eta就是将J(k)位置的变量取反，计算适应度之差
            end
            P=P/(sum(P));
            %按概率原则选取下一个城市
            Pcum=cumsum(P);     %cumsum，元素累加即求和
            Select=find(Pcum>=rand); %若计算的概率大于原来的就选择这条路线
            to_visit=J(Select(1));
            Tabu(i,j)=to_visit;
			Choice(i,J(Select(1))) = ~Choice(i,J(Select(1)));
        end
    end
    %%第四步：记录本次迭代最佳适应度
    L=zeros(m,1);     %开始距离为0，m*1的列向量
    for i=1:m
        L(i) = nij(Choice(i,:),Fl,betaT,betaE,L0,v,theta,N,tpro,C,f0,D,B,W,R,rou,p);      %一轮下来后走过的距离
    end
    L_best(NC)=max(L);           %最佳适应度取最大
    pos=find(L==L_best(NC));
    R_best(NC,:)=Choice(pos(1),:); %此轮迭代后的最佳决策
    L_ave(NC)=mean(L);           %此轮迭代后的平均适应度
    NC=NC+1                      %迭代继续
 
    
    %%第五步：更新信息素
    Delta_Tau=zeros(N);        %开始时信息素为n*1的0矩阵
    for i=1:m
        for j=1:N
            Delta_Tau(j)=Delta_Tau(j)+Q*Choice(i,j);
            %此次循环在设备j上的信息素增量
        end
    end
    Tau=(1-Rho).*Tau+Delta_Tau; %考虑信息素挥发，更新后的信息素
    %%第六步：禁忌表清零
    Tabu=zeros(m,n);             %%直到最大迭代次数
end

%%第七步：输出结果
Pos=find(L_best==max(L_best)); %找到最佳路径（非0为真）
Best_Schedule=R_best(Pos(1),:) %最大迭代次数后最佳决策
Biggest_Length=L_best(Pos(1)) %最大迭代次数后最大适应度

function res = nij(location,Fl,betaT,betaE,L0,v,theta,N,tpro,C,f0,D,B,W,R,rou,p)
    Tl = C./Fl;
    El = 1e-20*Fl*C;
    tao = rou*betaT;

    %caculate the transtime 
    for i=1:N
        Ttrans(i) = transtime(L0(i),v(i),theta(i),tpro,D,W,p(i));
    end

	if sum(location==1)>B/W
		return 0;
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
	objvalue=SysUti;
	return objvalue;
end
