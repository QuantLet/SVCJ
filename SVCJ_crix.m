
clear all
close all
clc
%data = simSVCJdata(4000);
%Y = data(2:end,1);

%%%load data

load CrixData.mat
dateS=CrixData(:,1);
dateD=year(dateS)*10000+month(dateS)*100+day(dateS);

y=price2ret(CrixData(:,2))*sqrt(250);
Y=y;

%Y = y; %WTI
%Y=yB; %Brent

%load wti1.mat
% Y = 100 * ( log(WTI(2:end,2))-log(WTI(1:end-1,2)) );
%Y = 100 * ( log(US_index(2:end))-log(US_index(1:end-1)) );
T = length(Y);
one = ones(T,1);

% Prior distribution hyperparameter values
a=0;
%a = [0 0]';
%A = 25*eye(2);
A=25;%prior for mu
b = [0 0]'; %prior for alpha and beta, 
B = 1*eye(length(b)); %prior for alpha and beta
c = 2.5;   %prior for sigma2_v
C = 0.1; %prior for signma2_v
d = 10;  %prior for mu_v
D = 20; % prior for mu_v
e = 0;
E = 100; %prior for mu_y
f = 5;  %prior for sigma2_y
F = 40;  %prior for sigma2_y
g = 0;
G = 4;
k = 2;  %prior for lambda
K = 30; %prior for lambda
tmp=[];
% Starting values
m = a; msum = 0; m2sum = 0;                 %% m=mu
kappa = 0; kappasum = 0; kappa2sum = 0;     %% kappa=-alpha/beta, the theta in equation 1
alpha = b(1); alphasum = 0; alpha2sum = 0;   %%alpha in equation 3
beta = b(2); betasum = 0, beta2sum = 0;     %% beta in eq.3
s2V = C/(c - 2); s2Vsum = 0; s2V2sum = 0;  %% sigma_v in eq.3
rho = 0; rhosum = 0; rho2sum = 0;           %% the relation between w1 ad w2
mV = D/(d-2); mVsum = 0; mV2sum = 0;        %% mu_v, the param in expoential distr. of Z_v 
                                            %%(jump size in variance
mJ = e; mJsum = 0; mJ2sum = 0;              %%mu_y, the mean of jump size in price Z_y
s2J = F/(f - 2); s2Jsum = 0; s2J2sum = 0;   %% sigma_Y, the variance of jump size in price Z_y 
rhoJ = g; rhoJsum = 0; rhoJ2sum = 0;        %% rho para in the jump size of price
lambda = 0; lambdasum = 0; lambda2sum = 0;  %% jump intensity

V = 0.1*(Y - mean(Y)).^2 + 0.9*var(Y); %initial values for variance_t;
Vsum = 0;Vsum2 = 0;
J = abs(Y) - mean(Y) > 2 * std(Y); %J = data(2:end,3);
Jsum = 0;
XV = exprnd(mV,T,1); % the jump size in price, Z_t^y
XVsum = 0;
X = mvnrnd((mJ + XV*rhoJ), s2J); % the jump process in variance, Z_t^y*dN_t
Xsum = 0;
stdevrho = 0.01;
dfrho = 6.5; 
stdevV = 0.9;
dfV = 4.5;
acceptsumV = zeros(size(V));
acceptsumrho = 0;
acceptsums2V = 0;
Z = ones(T,1);
N =5.e3;   % Specifies the total number of draws
n =1.e3; % Specifies the burn-in period
test = zeros(T,10);%%matrix for params
%rho = -0.5; alpha = 0.015; beta = -0.03; s2V = 0.01; rhoJ = -1; mJ = -1; s2J = 4; mV = 1;
for i = 1:N

    i
    Rho = 1 / (1 - rho^2);   
    V0 = V(1);  
        
    % Draw m(i+1)
    Q = ( Y - X.*J - rho/s2V^0.5.*...
        ( V - [V0;V(1:end-1)]*(1+beta)-alpha - J.*XV ) )./ ([V0;V(1:end-1)] ).^0.5; 
    %W = [1./([V0;V(1:end-1)]).^0.5, [0;Y(1:end-1)]./([V0;V(1:end-1)]).^0.5];
    W = [1./([V0;V(1:end-1)]).^0.5];
    
    As = inv( inv(A) + 1 / (1-rho^2) * W'*W );
    as = As*( inv(A) * a + 1 / ( 1 - rho^2 ) * W'*Q );
    
    m = mvnrnd(as',As)';
    
    if i > n
        msum = msum + m;
        m2sum = m2sum + m.^2;
    end
    
    % (alpha, beta)
    eY = Y - Z*m - X.*J; %expected return
    eV = V - [V0;V(1:end-1)] - XV.*J; %expected variance
    Q = (eV - rho*sqrt(s2V)*eY)./[V0;V(1:end-1)].^0.5;
    W = [1./[V0;V(1:end-1)].^0.5 [V0;V(1:end-1)].^0.5];
    Bs = inv(inv(B) + (Rho/s2V)*W'*W);
    bs = Bs*(inv(B)*b + (Rho/s2V)*(W'*Q));
    temp = mvnrnd(bs,Bs);
    alpha = temp(1); 
    beta = temp(2); 
    kappa = - alpha/beta;
    if i > n
        alphasum = alphasum + alpha;alpha2sum = alpha2sum + alpha^2;
        betasum = betasum + beta; beta2sum = beta2sum + beta^2;
        kappasum = kappasum + kappa; kappa2sum = kappa2sum + kappa^2;
    end
    % s2V 
    cs = c + T;
    Cs = C + sum( ((V - [V0;V(1:end-1)] - alpha - beta*[V0;V(1:end-1)] - XV.*J).^2)./[V0;V(1:end-1)] );
    s2Vprop = iwishrnd(Cs,cs); 
    q = exp(-0.5*sum((V - [V0;V(1:end-1)]*(1+beta) - alpha - J.*XV).^2./(s2Vprop*[V0;V(1:end-1)])-...
            (V - [V0;V(1:end-1)]*(1+beta) - alpha - J.*XV).^2./(s2V*[V0;V(1:end-1)])));  
    p = exp(-0.5*sum((V - [V0;V(1:end-1)]*(1+beta) - alpha - J.*XV - rho*s2Vprop^0.5*(Y-Z*m-J.*X)).^2./...
           ((1-rho^2)*s2Vprop*[V0;V(1:end-1)])-...
           (V - [V0;V(1:end-1)]*(1+beta) - alpha - J.*XV -rho*s2V^0.5*(Y-Z*m-J.*X)).^2./...
           ((1-rho^2)*s2V*[V0;V(1:end-1)])));  
    x = min(p/q,1);
    u = rand(1);
    if x > u
        s2V = s2Vprop; 
        if i > n; acceptsums2V = acceptsums2V + 1;end;
    end
  
    if i > n
        s2Vsum = s2Vsum + s2V;
        s2V2sum = s2V2sum + s2V^2;
    end
    
    % rho
    %Draw a candidate for rho(i+1)
    rhoprop = rho + stdevrho * trnd(dfrho); %draw rhoc from a t distribution with 8 df and std  of 0.2666
    if abs(rhoprop) < 1
        p = (sqrt( 1 - rho^2 )/ sqrt( 1 - rhoprop^2 )).^T * exp( sum( - 1 / ( 2 * ( 1 - rhoprop^2 ) ) *...
           ( Y - Z*m - J.*X - rhoprop / s2V^0.5 * ( V - alpha - [V0;V(1:end-1)] * ( 1 + beta ) - J.*XV ) ).^2./...
             [V0;V(1:end-1)] + 1 / ( 2 * ( 1 - rho^2 ) ) *...
           ( Y - Z*m - J.*X - rho / s2V^0.5 * ( V(1:end) - alpha - [V0;V(1:end-1)] * ( 1 + beta ) - J.*XV ) ).^2./...
             [V0;V(1:end-1)] ) );
        u = rand(1);
        x = min(p,1);
        if x > u
            rho = rhoprop;
            if i > n; acceptsumrho = acceptsumrho +1;end;
        end
    end
    if i > n
        rhosum = rhosum + rho;
        rho2sum = rho2sum + rho^2;    
    end
    
    % mV
    ds = d + 2*T;
    Ds = D + 2*sum(XV);
    mV = iwishrnd(Ds,ds);
    if i > n
        mVsum = mVsum + mV;
        mV2sum = mV2sum + mV^2;
    end
    % mJ
    Es = 1/(T/s2J + 1/E);
    es = Es * (sum( (X - XV*rhoJ)/s2J ) + e/E);
    mJ = normrnd(es,Es^0.5);
    if i > n
        mJsum = mJsum + mJ;
        mJ2sum = mJ2sum + mJ^2;
    end
    
    % s2Y
    fs = f + T;
    Fs = F + sum((X - mJ - rhoJ*XV).^2);
    s2J = iwishrnd(Fs,fs);
    if i > n
        s2Jsum = s2Jsum + s2J;
        s2J2sum = s2J2sum + s2J^2;
    end
    
    % rhoJ
    Gs = inv(sum(XV.^2)/s2J + 1/G);
    gs = Gs * (sum((X - mJ).*XV)/s2J + g/G);
    rhoJ = normrnd(gs,Gs^0.5); 
    if i > n
        rhoJsum = rhoJsum + rhoJ; rhoJ2sum = rhoJ2sum + rhoJ^2;
    end
    
    % lambda
    ks = k + sum(J);
    Ks = K + T - sum(J);
    lambda = betarnd(ks,Ks); 
    if i > n
        lambdasum = lambdasum + lambda;
        lambda2sum = lambda2sum + lambda^2;
    end
    
    % J
    eY1 = Y - Z*m - X;
    eY2 = Y - Z*m;
    eV1 = V - [V0;V(1:end-1)] - alpha - beta*[V0;V(1:end-1)] - XV;
    eV2 = V - [V0;V(1:end-1)] - alpha - beta*[V0;V(1:end-1)];        
    p1 = lambda*exp( -0.5 * ( ((eY1 - (rho/sqrt(s2V))*eV1).^2)./((1-rho^2)*[V0;V(1:end-1)]) + (eV1.^2)./(s2V*[V0;V(1:end-1)]) ) );
    p2 = (1 - lambda) * exp( -0.5 * ( ((eY2 - (rho/sqrt(s2V))*eV2).^2)./((1-rho^2)*[V0;V(1:end-1)]) + (eV2.^2)./(s2V*[V0;V(1:end-1)]) ) );
    p = p1./(p1 + p2);
     tmp=[tmp;[p1 p2 p]];
   
    u = rand(T,1);
    J = double(u < p);
    if i > n; Jsum = Jsum + J; end; 
    
    Jindex = find(J == 1);
    
    % XV
    XV(logical(~J)) = exprnd(mV,T-sum(J),1); 
    if ~isempty(Jindex)
        if Jindex(1) == 1
            t = 1;
            eV = V(1) - V0 - alpha - beta*V0;
            eY = Y(1) - Z(1,:)*m - X(1);
            H = inv( 1 /((1 - rho^2)*s2V*V0) + rhoJ^2/s2J );
            h = H * ((eV-rho*sqrt(s2V)*eY)/((1 - rho^2)*s2V*V0) + rhoJ*(X(1) - mJ)/s2J - 1/mV);
            if h+5*sqrt(H) > 0; XV(1) = normt_rnd(h,H,0,h+5*sqrt(H)); else XV(1) = 0; end;
            if XV(1) == Inf | XV(1) == NaN; XV(1) = 0; end; 
        else
            t = Jindex(1);
            eV = V(t) - V(t-1) - alpha - beta*V(t-1);
            eY = Y(t) - Z(t,:)*m - X(t);
            H = inv( 1 /((1 - rho^2)*s2V*V(t-1)) + rhoJ^2/s2J );
            h = H * ((eV-rho*sqrt(s2V)*eY)/((1 - rho^2)*s2V*V(t-1)) + rhoJ*(X(t) - mJ)/s2J - 1/mV);
            if h+5*sqrt(H) > 0; XV(t) = normt_rnd(h,H,0,h+5*sqrt(H)); else XV(t) = 0; end;
            if XV(t) == Inf | XV(t) == NaN; XV(t) = 0; end; 
        end
        if length(Jindex) > 1
            for t = Jindex(2:end)'
                eV = V(t) - V(t-1) - alpha - beta*V(t-1);
                eY = Y(t) - Z(t,:)*m - X(t);
                H = inv( 1 /((1 - rho^2)*s2V*V(t-1)) + rhoJ^2/s2J );
                h = H * ((eV-rho*sqrt(s2V)*eY)/((1 - rho^2)*s2V*V(t-1)) + rhoJ*(X(t) - mJ)/s2J - 1/mV);
                if h+5*sqrt(H) > 0; XV(t) = normt_rnd(h,H,0,h+5*sqrt(H)); else XV(t) = 0; end;
                if XV(t) == Inf | XV(t) == NaN; XV(t) = 0; end;
            end
        end
    end
    if i > n; XVsum = XVsum + XV; end;
    
    % X
    X(logical(~J)) = normrnd(mJ + rhoJ*XV(logical(~J)),s2J^0.5)';   
    if ~isempty(Jindex)
        if Jindex(1) == 1
            t = 1;
            eV = V(1) - V0 - alpha - beta*V0 - XV(1);
            eY = Y(1) - Z(1,:)*m;
            L = inv(1/((1 - rho^2)*V0) + 1/s2J);
            l = L * ( (eY - (rho/sqrt(s2V))*eV)/((1 - rho^2)*V0) + (mJ + rhoJ*XV(1))/s2J );
            X(1) = normrnd(l,sqrt(L));
        else
            t = Jindex(1);
            eV = V(t) - V(t-1) - alpha - beta*V(t-1) - XV(t);
            eY = Y(t) - Z(t,:)*m;
            L = inv(1/((1 - rho^2)*V(t-1)) + 1/s2J);
            l = L * ( (eY - (rho/sqrt(s2V))*eV)/((1 - rho^2)*V(t-1)) + (mJ + rhoJ*XV(t))/s2J );
            X(t) = normrnd(l,sqrt(L));
        end
        if length(Jindex) > 1
            for t = Jindex(2:end)'
                eV = V(t) - V(t-1) - alpha - beta*V(t-1) - XV(t);
                eY = Y(t) - Z(t,:)*m;
                L = inv(1/((1 - rho^2)*V(t-1)) + 1/s2J);
                l = L * ( (eY - (rho/sqrt(s2V))*eV)/((1 - rho^2)*V(t-1)) + (mJ + rhoJ*XV(t))/s2J );
                X(t) = normrnd(l,sqrt(L));
            end
        end
    end
    if i > n; Xsum = Xsum + X;end
    
    % Draw V
    epsilon = trnd(dfV,T,1); 
    [mv, v] = tstat(dfV);
    epsilon = (stdevV/sqrt(v)) * epsilon;
    if i == floor(n / 2);
        Vindex1 = find(Vsum2 > prctile(Vsum2, 92.5));
        Vindex2 = find(Vsum2 > prctile(Vsum2, 75) & ...
                       Vsum2 < prctile(Vsum2, 92.5));
        Vindex3 = find(Vsum2 < prctile(Vsum2, 25) & ...
                       Vsum2 > prctile(Vsum2, 2.5));
        Vindex4 = find(Vsum2 < prctile(Vsum2, 2.5));
    end
    if i > floor(n / 2) - 1
        epsilon(Vindex1) = 1.35 * epsilon(Vindex1);
        epsilon(Vindex2) = 1.25 * epsilon(Vindex2);
        epsilon(Vindex3) = 0.75 * epsilon(Vindex3);
        epsilon(Vindex4) = 0.65 * epsilon(Vindex4);
    end
    j = 1;
    Vprop = V + epsilon;
    p1 = max(0,exp( -0.5 * ( ( Y(j+1) - Z(j+1,1)*m(1)  - J(j+1)*X(j+1) - rho / s2V^0.5 *(V(j+1) - Vprop(j) - alpha - Vprop(j) * beta - J(j+1)*XV(j+1) ) )^2/( (1 - rho^2) * Vprop(j) ) +...
            ( Y(j) - Z(j,1)*m(1)  - J(j)*X(j) - rho / s2V^0.5 *(Vprop(j) - V0 - alpha - V0 * beta - J(j)*XV(j)))^2/( (1 - rho^2) * V0 ) +...
            ( V(j+1) - Vprop(j) - alpha - Vprop(j) * beta - J(j+1)*XV(j+1) )^2/( s2V * Vprop(j) ) +...
            ( Vprop(j) - V0 - alpha - V0 * beta - J(j)*XV(j))^2/( s2V * V0 ) ) ) / Vprop(j));
    p2 = max(0,exp( -0.5 * ( ( Y(j+1) - Z(j+1,1)*m(1)  - J(j+1)*X(j+1) - rho / s2V^0.5 *(V(j+1) - V(j) - alpha - V(j) * beta - J(j+1)*XV(j+1) ) )^2/( (1 - rho^2) * V(j) ) +...
            ( Y(j) - Z(j,1)*m(1)  - J(j)*X(j) -rho / s2V^0.5 *(V(j) - V0 - alpha - V0 * beta - J(j)*XV(j)) )^2/( (1 - rho^2) * V0 ) +...
            ( V(j+1) - V(j) - alpha - V(j) * beta - J(j+1)*XV(j+1))^2/( s2V * V(j) ) +...
            ( V(j) - V0 - alpha - V0 * beta - J(j)*XV(j))^2/( s2V * V0 ) ) ) / V(j));
    if p2 ~= 0; acceptV = min(p1/p2, 1); elseif p1 > 0; acceptV = 1; else; acceptV = 0; end;
    u = rand(T,1);
    if u(j) < acceptV 
        V(j) = Vprop(j);
        if i > n; acceptsumV(j) = acceptsumV(j) + 1; end;
    end
    
    for j = 2:T-1
        p1 = max(0,exp( -0.5 * ( ( Y(j+1) - Z(j+1,1)*m(1)  - J(j+1)*X(j+1) - rho / s2V^0.5 *(V(j+1) - Vprop(j) - alpha - Vprop(j) * beta - J(j+1)*XV(j+1)) )^2/( (1 - rho^2) * Vprop(j) ) +...
            ( Y(j) - Z(j,1)*m(1)  - J(j)*X(j) - rho / s2V^0.5 *(Vprop(j) - V(j-1) - alpha - V(j-1) * beta - J(j)*XV(j) ) )^2/( (1 - rho^2) * V(j-1) ) +...
            ( V(j+1) - Vprop(j) - alpha - Vprop(j) * beta - J(j+1)*XV(j+1))^2/( s2V * Vprop(j) ) +...
            ( Vprop(j) - V(j-1) - alpha - V(j-1) * beta - J(j)*XV(j))^2/( s2V * V(j-1) ) ) ) / Vprop(j));
        p2 = max(0,exp( -0.5 * ( ( Y(j+1) - Z(j+1,1)*m(1)  - J(j+1)*X(j+1) - rho / s2V^0.5 *(V(j+1) - V(j) - alpha - V(j) * beta - J(j+1)*XV(j+1)) )^2/( (1 - rho^2) * V(j) ) +...
            ( Y(j) - Z(j,1)*m(1) - J(j)*X(j) - rho / s2V^0.5 *(V(j) - V(j-1) - alpha - J(j)*XV(j) - V(j-1) * beta ) )^2/( (1 - rho^2) * V(j-1) ) +...
            ( V(j+1) - V(j) - alpha - V(j) * beta - J(j+1)*XV(j+1))^2/( s2V * V(j) ) +...
            ( V(j) - V(j-1) - alpha - V(j-1) * beta - J(j)*XV(j))^2/( s2V * V(j-1) ) ) ) / V(j));
        if p2 ~= 0; acceptV = min(p1/p2, 1); elseif p1 > 0; acceptV = 1; else; acceptV = 0; end;
       
        if u(j) < acceptV 
            V(j) = Vprop(j);
            if i > n; acceptsumV(j) = acceptsumV(j) + 1; end;
        end
    end
    j = T;
    p1 = max(0,exp( -0.5 * ( ( Y(j) - Z(j,1)*m(1)  - J(j)*X(j) - rho / s2V^0.5 *(Vprop(j) - V(j-1) - alpha - V(j-1) * beta - J(j)*XV(j)) )^2/( (1 - rho^2) * V(j-1) ) +...
            ( Vprop(j) - V(j-1) - alpha - V(j-1) * beta - J(j)*XV(j))^2/( s2V * V(j-1) ) ) ) / Vprop(j)^0.5);
    p2 = max(0,exp( -0.5 * ( ( Y(j) - Z(j,1)*m(1)  - J(j)*X(j) - rho / s2V^0.5 *(V(j) - V(j-1) - alpha - V(j-1) * beta - J(j)*XV(j)) )^2/( (1 - rho^2) * V(j-1) ) +...
            ( V(j) - V(j-1) - alpha - V(j-1) * beta  - J(j)*XV(j))^2/( s2V * V(j-1) ) ) ) / V(j)^0.5);
    if p2 ~= 0; acceptV = min(p1/p2, 1); elseif p1 > 0; acceptV = 1; else; acceptV = 0; end;
       
    if u(j) < acceptV 
            V(j) = Vprop(j);
            if i > n; acceptsumV(j) = acceptsumV(j) + 1; end;
    end
   
    if i > n; Vsum = Vsum + V;end;
    if i > floor(n / 2) - 100 | i < floor(n / 2); Vsum2 = Vsum2 + V; end;
    test(i,:) = [m mJ s2J lambda alpha beta  rho s2V rhoJ mV];

end

disp(['m'])
disp([msum/(N-n) (m2sum/(N-n)-(msum/(N-n)).^2).^0.5])
disp(['mJ'])
disp([mJsum/(N-n) (mJ2sum/(N-n)-(mJsum/(N-n))^2)^0.5])
disp(['VJ'])
disp([s2Jsum/(N-n) (s2J2sum/(N-n)-(s2Jsum/(N-n))^2)^0.5])
disp(['lambda'])
disp([lambdasum/(N-n) (lambda2sum/(N-n)-(lambdasum/(N-n))^2)^0.5])
disp(['alpha'])
disp([alphasum/(N-n) (alpha2sum/(N-n)-(alphasum/(N-n))^2)^0.5])
disp(['beta'])
disp([betasum/(N-n) (beta2sum/(N-n)-(betasum/(N-n))^2)^0.5])
disp(['kappa'])
disp([kappasum/(N-n) (kappa2sum/(N-n)-(kappasum/(N-n))^2)^0.5])
disp(['rho'])
disp([rhosum/(N-n) (rho2sum/(N-n)-(rhosum/(N-n))^2)^0.5])
disp(['s2V'])
disp([s2Vsum/(N-n) (s2V2sum/(N-n)-(s2Vsum/(N-n))^2)^0.5])
disp(['rhoJ'])
disp([rhoJsum/(N-n) (rhoJ2sum/(N-n)-(rhoJsum/(N-n))^2)^0.5])
disp(['mV'])
disp([mVsum/(N-n) (mV2sum/(N-n)-(mVsum/(N-n))^2)^0.5])


Result=[msum/(N-n) (m2sum/(N-n)-(msum/(N-n)).^2).^0.5;
        mJsum/(N-n) (mJ2sum/(N-n)-(mJsum/(N-n))^2)^0.5;
        s2Jsum/(N-n)  (s2J2sum/(N-n)-(s2Jsum/(N-n))^2)^0.5;
        lambdasum/(N-n) (lambda2sum/(N-n)-(lambdasum/(N-n))^2)^0.5;
        alphasum/(N-n) (alpha2sum/(N-n)-(alphasum/(N-n))^2)^0.5;
        betasum/(N-n) (beta2sum/(N-n)-(betasum/(N-n))^2)^0.5;
        rhosum/(N-n) (rho2sum/(N-n)-(rhosum/(N-n))^2)^0.5;
        s2Vsum/(N-n) (s2V2sum/(N-n)-(s2Vsum/(N-n))^2)^0.5;
        rhoJsum/(N-n) (rhoJ2sum/(N-n)-(rhoJsum/(N-n))^2)^0.5;
        mVsum/(N-n) (mV2sum/(N-n)-(mVsum/(N-n))^2)^0.5];

TS=CrixData(2:end,1);

% figure   
% subplot(4,1,1)
% plot(TS,Y)
% datetick('x', 'mmmyy', 'keeplimits', 'keepticks')
% 
% subplot(4,1,2)
% plot(TS,Jsum/(N-n))
% datetick('x', 'mmmyy', 'keeplimits', 'keepticks')
% 
% subplot(4,1,3)
% plot(TS,Xsum/(N-n))
% datetick('x', 'mmmyy', 'keeplimits', 'keepticks')
% 
% subplot(4,1,4)
% plot(TS,Vsum/(N-n))
% datetick('x', 'mmmyy', 'keeplimits', 'keepticks')

eX=Xsum/(N-n);
eXV=XVsum/(N-n);
eJ=round(Jsum/(N-n));
eV=Vsum/(N-n);



% figure
% subplot(2,1,1)
% plot(TS,eX.*eJ)
% datetick('x', 'mmmyy', 'keeplimits', 'keepticks')
% 
% subplot(2,1,2)
% plot(TS,eXV.*eJ)
% datetick('x', 'mmmyy', 'keeplimits', 'keepticks')
% 
 Jump_V=eXV.*eJ;
 Jump_P=eX.*eJ;

%Diff_V=Jump_V-Xsum/(N-n);
%Diff_p=Jump_V-X/(N-n);

Sig=eV.^.5;
ResY=(Y(2:end)-m-Jump_P(2:end))./Sig(1:end-1);




MSE_SVCJ=mse(eV,Y.^2);


Mdl = garch(1,1);
EstMdl = estimate(Mdl,Y);
 
vG = infer(EstMdl,Y);
MSE_G=mse(vG,Y.^2);
Res_G=Y./sqrt(vG);



% figure
% normplot(Y)
% set(gca, 'xlim', [-4 4])
% title('Log Return');
% 
% figure
% normplot(Res_G);
% set(gca, 'xlim', [-4 4])
% title('GARCH');
% 
% 
% figure
% %normplot(ResY)
% qqplot(ResY)
% set(gca, 'xlim', [-4 4])
% title('SVCJ','FontSize',18);

%%%save('testRaw_svcj','test')
%%%save ('allVars_svcj')
