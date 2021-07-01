clear;
close all;
clc;


% Badania współczynnika korelacji Pearsona
% X = rand(1,10000) - 0.5;
% Y = X.^2;
% corr = rn(X,Y);

iter = 10000;
wek = [];

for i=1:1:iter
    X = rand(1,i) - 0.5;
    Y = X.^2;
    corr = rn(X, Y);
    wek(end+1) = corr;
end

plot(wek);
title('Zaleznosc wspolczynnika korelacji od ilosci iteracji $Y_N = X^2_N$','interpreter','latex');
xlabel('N','interpreter','latex');
ylabel('Corr','interpreter','latex');

% Badania generacji z rozkładu Cauchy'ego
% cvar = Cauchy(10000);
% miN = EX(cvar);
% var = sN(cvar, miN);

% iter = 1000;
% wek1 = [];
% wek2 = [];
% for i=1:1:iter
%     cvar = Cauchy(100000);
%     miN = EX(cvar);
%     var = sN(cvar, miN);
%     wek1(end+1) = miN;
%     wek2(end+1) = var;
% end
% 
% subplot(2,1,1);
% plot(wek1,'.');
% xlabel('N','interpreter','latex');
% ylabel('$\widehat{\mu}_N$','interpreter','latex');
% title('Estymacja wartosci oczekiwanej rozkladu Cauchyego','interpreter','latex');
% subplot(2,1,2);
% plot(wek2,'.');
% xlabel('N','interpreter','latex');
% ylabel('$\widehat{s}^2_N$','interpreter','latex');
% title('Estymacja wariancji rozkladu Cauchyego','interpreter','latex');

% Badania błędów empirycznych estymatorów
% N = 5;
% L = 20;
% iter = 10000;
% 
% Err0 = 0;
% Err1 = 0;
% Err2 = 0;
% for i=1:1:iter
%    Err0 = Err0 + EX_Err(L, N, 0, 1);
%    Err1 = Err1 + Var_Err_Ob(L, N, 0, 1);
%    Err2 = Err2 + Var_Err_Nie(L, N, 0, 1);
% end
% Sredni_Blad_Wart = Err0/iter;
% Sredni_Blad_Obciazonego = Err1/iter;
% Sredni_Blad_Nieobciazonego = Err2/iter;

% Wykresy błędu empirycznego
% L1 = 10;
% L2 = 50;
% L3 = 100;
% N = 1:1:1000;
% 
% Err1 = [];
% Err2 = [];
% Err3 = [];
% 
% for i=1:1:length(N)
%    Err1(end+1) = EX_Err(L1,i,0,1);
%    Err2(end+1) = EX_Err(L1,i,0,1);
%    Err3(end+1) = EX_Err(L1,i,0,1);
% end
% 
% figure(1);
% plot(N,Err1);
% hold on;
% plot(N,Err2);
% hold on;
% plot(N,Err3);
% xlabel('N','interpreter','latex');
% ylabel('Err','interpreter','latex');
% title('Blad empiryczny estymatora wartosci oczekiwanej','interpreter','latex');
% legend('L=10','L=50','L=100','interpreter','latex');


% Err1 = [];
% Err2 = [];
% Err3 = [];
% Err4 = [];
% Err5 = [];
% Err6 = [];
% for i=1:1:1000
%     Err1(end+1) = Var_Err_Ob(L1,i,0,1);
%     Err2(end+1) = Var_Err_Nie(L1,i,0,1);
%     Err3(end+1) = Var_Err_Ob(L2,i,0,1);
%     Err4(end+1) = Var_Err_Nie(L2,i,0,1);
%     Err5(end+1) = Var_Err_Ob(L3,i,0,1);
%     Err6(end+1) = Var_Err_Nie(L3,i,0,1);
% end
% 
% subplot(3,1,1);
% plot(N,Err1);
% hold on;
% plot(N,Err2);
% xlabel('N','interpreter','latex');
% ylabel('Err','interpreter','latex');
% title('Blad empiryczny estymatorow wariancji - $L = 10$','interpreter','latex');
% legend('Ob','Nieob','interpreter','latex');
% subplot(3,1,2);
% plot(N,Err3);
% hold on;
% plot(N,Err4);
% xlabel('N','interpreter','latex');
% ylabel('Err','interpreter','latex');
% title('Blad empiryczny estymatorow wariancji - $L = 50$','interpreter','latex');
% legend('Ob','Nieob','interpreter','latex');
% subplot(3,1,3);
% plot(N,Err5);
% hold on;
% plot(N,Err6);
% xlabel('N','interpreter','latex');
% ylabel('Err','interpreter','latex');
% title('Blad empiryczny estymatorow wariancji - $L = 100$','interpreter','latex');
% legend('Ob','Nieob','interpreter','latex');

% iter = 10000;
% wek = [];
% for i=1:1:iter
%     X = normrnd(0,1,[1,i]);
%     miN = EX(X);
%     var = sN(X, miN);
%     wek(end+1) = var;
% end
% 
% wek2 = [];
% for i=1:1:iter
%     X = normrnd(0,2,[1,i]);
%     miN = EX(X);
%     var = sN(X, miN);
%     wek2(end+1) = var;
% end
% 
% subplot(2,1,1);
% plot(wek);
% xlabel('N','interpreter','latex');
% ylabel('$\widehat{S}^2_N$','interpreter','latex')
% title('Nieobciazony estymator wariancji - $\mathcal{N}(0,1)$','interpreter','latex');
% subplot(2,1,2);
% plot(wek2);
% title('Nieobciazony estymator wariancji - $\mathcal{N}(0,2)$','interpreter','latex');
% xlabel('N','interpreter','latex');
% ylabel('$\widehat{S}^2_N$','interpreter','latex')

% iter = 1000;
% N = 100000;
% mu = 0;
% sigma = 2;
% sum = 0;
% for i=1:1:iter
%     Err = Var_Err_Nie(10,N,mu,sigma);
%     sum = sum + Err;
% end
% sum = sum/iter;

% Estymator wartości oczekiwanej
function miN = EX(X)
    miN = 1/length(X) * sum(X);
end

%Estymator wariancji obciążony
function var = sN(X, miN)
    suma = 0;
    for i=1:1:length(X)
       suma = suma + (X(i) - miN)^2; 
    end
    var = 1/length(X) * suma;
end

%Estymator wariancji nieobciążony
function var = SN(X, miN)
    suma = 0;
    for i=1:1:length(X)
       suma = suma + (X(i) - miN)^2; 
    end
    var = 1/(length(X)-1) * suma;
end

%Estymator korelacji i kowariancji
function corr = rn(X1, X2)
    if(length(X1) == length(X2))
        N = length(X1);
        X1_sr = mean(X1);
        X2_sr = mean(X2);
        Licznik = 0;
        for i=1:1:N
           Licznik = Licznik + ((X1(i) - X1_sr)*(X2(i) - X2_sr)); 
        end
        Licznik = 1/N * Licznik;
        Mianownik = 0;
        suma1 = 0;
        suma2 = 0;
        for i=1:1:N
           suma1 = suma1 + (X1(i) - X1_sr)^2; 
           suma2 = suma2 + (X2(i) - X2_sr)^2;
        end
        suma1 = 1/N * suma1;
        suma2 = 1/N * suma2;
        Mianownik = sqrt(suma1 * suma2);
        corr = Licznik/Mianownik;
    else
       return; 
    end
end

%Wyznaczanie błędu empirycznego wartości oczekiwanej
function Err = EX_Err(L, N, mu, sigma)
    suma = 0;
    for i=1:1:L
       T_n = normrnd(mu, sigma, [1,N]);
       suma = suma + (EX(T_n) - mu)^2;
    end
    Err = 1/L * suma;
end

%Wyznaczenie błędu empirycznego wariancji est. obciążonego wariancji
function Err = Var_Err_Ob(L, N, mu, sigma)
    suma = 0;
    for i=1:1:L
       T_n = normrnd(mu, sigma, [1,N]);
       miN = EX(T_n);
       suma = suma + (sN(T_n, miN) - sigma^2)^2;
    end
    Err = 1/L * suma;
end

%Wyznaczenie błędu empirycznego wariancji est. nieobciążonego wariancji
function Err = Var_Err_Nie(L, N, mu, sigma)
    suma = 0;
    for i=1:1:L
       T_n = normrnd(mu, sigma, [1,N]);
       miN = EX(T_n);
       suma = suma + (SN(T_n, miN) - sigma^2)^2;
    end
    Err = 1/L * suma;
end

%Generowanie liczb losowych z rozkładu Cauchyiego
function Cvar = Cauchy(N)
    %cvar = @(N) tan((rand(1,N) - 0.5)*pi);
    cvar = @(N) randn(1,N)./randn(1,N);
    Cvar = cvar(N);
end

%Generator liczb z rozkładu normalnego metodą Box-Muller
function Rozklad = Box_Muller(iter, mu, sigma)
    
    Rozklad = [];
    for i=1:1:iter/2
        u1 = rand(1,1);
        u2 = rand(1,1);
        Pierw = sigma * sqrt(-2*log(u1));
        Rozklad(end+1) = Pierw * cos(2*pi*u2) + mu;
        Rozklad(end+1) = Pierw * sin(2*pi*u2) + mu;
    end
end

