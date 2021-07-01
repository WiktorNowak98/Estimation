clear;
close all;
clc;

Zad_4();
% x = -1:0.01:2;
% [wek, Fn_emp] = Zad_5(x);
% 
% figure(1);
% plot(x,wek);
% hold on;
% plot(x, (Fn_emp.*(1-Fn_emp))/1000);

% X = 0:0.01:.99;
% N = 1000;
% 
% u = [];
% 
% for i=1:1:length(X)
%    [wek, Fn] = Ojoj(N,X,i);
%    u = [u;wek];
% end
% 
% hm = Ojej(u);
% surf(hm);
% title('Prezentacja zbieznosci wedlug rozkladu dystrybuanty empirycznej do rozkladu normalnego - caly przedzial','interpreter','latex');

% jadro_epa = @(x) ((x >= -1) && (x <= 1)) .* (3/4 *(1 - x^2)) + ((x < -1) && (x > 1)) .* 0;
% h = 0.5;
% x = -5:0.01:5;
% 
% 
% subplot(5,1,1);
% histogram(u(5,:),'Normalization','pdf');
% hold on;
% plot(x,hm(5,:));
% title('Prezentacja zbieznosci wedlug rozkladu dystrybuanty empirycznej do rozkladu normalnego - x = 0.05','interpreter','latex');
% 
% subplot(5,1,2);
% histogram(u(20,:),'Normalization','pdf');
% hold on;
% plot(x,hm(20,:));
% title('Prezentacja zbieznosci wedlug rozkladu dystrybuanty empirycznej do rozkladu normalnego - x = 0.20','interpreter','latex');
% 
% subplot(5,1,3);
% histogram(u(50,:),'Normalization','pdf');
% hold on;
% plot(x,hm(50,:));
% title('Prezentacja zbieznosci wedlug rozkladu dystrybuanty empirycznej do rozkladu normalnego - x = 0.50','interpreter','latex');
% 
% subplot(5,1,4);
% histogram(u(75,:),'Normalization','pdf');
% hold on;
% plot(x,hm(75,:));
% title('Prezentacja zbieznosci wedlug rozkladu dystrybuanty empirycznej do rozkladu normalnego - x = 0.75','interpreter','latex');
% 
% subplot(5,1,5);
% histogram(u(95,:),'Normalization','pdf');
% hold on;
% plot(x,hm(95,:));
% title('Prezentacja zbieznosci wedlug rozkladu dystrybuanty empirycznej do rozkladu normalnego - x = 0.95','interpreter','latex');

% ESTYMATOR JĄDROWY GĘSTOŚCI
function Estymacja = Estymator_jadrowy(jadro, h, Xn, x)
    Estymacja = [];
    for i=1:1:length(x)
        suma = 0;
        for j=1:1:length(Xn)
            K = jadro((Xn(j) - x(i))/h);
            suma = suma + K;
        end
        Estymacja(end+1) = 1/(length(Xn)*h) * suma;
    end
end

% Sprawdzenie centralnego twierdzenia granicznego i zbieżności
function [wek,Fn_emp] = Ojoj(N, x, o)
    Fn = Zrob_dystr(x);
    wek = [];
    for i=1:1:N
       X = Zad_1(i);
       Fn_emp = Dyst_Emp(x,X);
       wek(end+1) = Fn_emp(o) - Fn(o);
    end
    wek = sqrt(N) * wek;
end

function hm = Ojej(u)
    jadro_epa = @(x) ((x >= -1) && (x <= 1)) .* (3/4 *(1 - x^2)) + ((x < -1) && (x > 1)) .* 0;
    x = -5:0.01:5;
    hm = [];
    for i=1:1:100
        a = Estymator_jadrowy(jadro_epa,0.5,u(i,:),x);
        hm = [hm;a];
    end
end

% Błąd empiryczny
function Err = Emp_err(L, x, N)
    Fn = Zrob_dystr(x);
    M = length(x);
    sumaL = 0;
    sumaM = 0;
    for i=1:1:L
       X = Zad_1(N);
       Fn_emp = Dyst_Emp(x, X);
       for j=1:1:M
          sumaM = sumaM + (Fn_emp(j) - Fn(j))^2; 
       end
       sumaL = sumaL + sumaM;
       sumaM = 0;
    end
    Err = 1/(L*M) * sumaL;
end

% Estymacja wariancji dystrybuanty
function [V, Fn_emp] = Zad_5(x)
    U = [];
    for i=1:1:1000
       X = Zad_1(1000);
       Fn_emp = Dyst_Emp(x, X);
       U = [U;Fn_emp];
    end
    V = [];
    for i=1:1:length(Fn_emp)
        col = U(:,i)';
        miN = EX(col);
        V(end+1) = SN(U(:,i)',miN);
    end
end

% Generator liczb metodą odwracania dystrybuanty z rozkładu trójkątnego
function Rozklad = Zad_1(N)
    X = rand(1,N);
    Rozklad = [];
    for i=1:1:length(X)
        Rozklad(end+1) = sqrt(X(i));
    end 
end

% Tworzy dystrybuantę rozkładu trójkątnego na zadanym przedziale
function wek = Zrob_dystr(x)
    f = @(x) (x<0) .* (0) + ((x>=0) && (x<=1)) .* (x^2) + (x>1) .* (1);
    wek = [];
    for i=1:1:length(x)
        wek(end+1) = f(x(i));
    end
end

% Dystrybuanta empiryczna, zakres oraz zmienne losowe z rozkł.
function Fn = Dyst_Emp(x, X)
    Fn = zeros(1,length(x));
    for i=1:1:length(x)
        suma = 0;
        for j=1:1:length(X)
            if(X(j) <= x(i))
                suma = suma + 1;
            else
               suma = suma + 0; 
            end
        end
        Fn(i) = 1/length(X) * suma;
    end
end

% Test Kołmogorowa-Smirnowa
function Dn = Kol_smirn(Fn, Fn_emp)
    Dn = zeros(1,length(Fn));
    for i=1:1:length(Fn)
        Dn(i) = abs(Fn_emp(i) - Fn(i));
    end
    Dn = max(Dn);
end

% Sprawdzenie twierdzenia Gliwienko-Cantelli
function Zad_3()
    x = -1:0.01:2;
    Fn = Zrob_dystr(x);
    Dn_wek = [];
    for i=1:1:5000
        X = Zad_1(i);
        Fn_emp = Dyst_Emp(x, X);
        Dn_wek(end+1) = Kol_smirn(Fn, Fn_emp);
    end
    N = 1:1:5000;
    figure(1);
    plot(N, Dn_wek);
    title('Odleglosc Kolmogorowa-Smirnowa w zaleznosci od ilosci probek','interpreter','latex');
    xlabel('N','interpreter','latex');
    ylabel('$D_N$','interpreter','latex');
    %legend('$F(x)$','$\widehat{F}_N(x)$','interpreter','latex');
end

% Odtworzenie i dopasowanie dystrybuanty nieznanego rozkładu
function Dn = Zad_4()
    load("asd.txt");
    X = asd';
    x = -10:0.001:10;
    N = length(X);
    Fn = zeros(1,length(x));
    
    for i=1:1:length(x)
        suma = 0;
        for j=1:1:N
            if(X(j) <= x(i))
                suma = suma + 1;
            else
               suma = suma + 0; 
            end
        end
        Fn(i) = 1/length(X) * suma;
    end
    % N(1,1)
    y1 = (1/2)*(1+erf((x-1)/sqrt(2)));
    % N(0,5)
    y2 = (1/2)*(1+erf(x/(5*sqrt(2))));
    % Cauchy'ego x0=0 i gamma=1
    y3 = 1/pi * atan(x)+1/2;
    
    Dn = zeros(1,3);
    Dn(1) = Kol_smirn(y1, Fn);
    Dn(2) = Kol_smirn(y2, Fn);
    Dn(3) = Kol_smirn(y3, Fn);
    
    figure(1);
    stairs(x,Fn);
    hold on;
    plot(x,y1);
    hold on;
    plot(x,y2);
    hold on;
    plot(x,y3);
    title('Porownanie dystrybuant rozkladow $\mathcal{N}(1,1)$, $\mathcal{N}(0,5)$ i  Cauchyego z wynikiem estymacji','interpreter','latex');
    xlabel('x','interpreter','latex');
    ylabel('F(x)','interpreter','latex');
    legend('$\widehat{F}_N(x)$','$\mathcal{N}(1,1)$','$\mathcal{N}(0,5)$','Cauchyego','interpreter','latex');
end

%Estymator wartości oczekiwanej
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

%%%%%%%%%%%%%%%%%% ZAD 1 %%%%%%%%%%%%%%%%%%%%%%
% x = -0.5:0.01:1.5;
% X = Zad_1(10);
% Fn = Zrob_dystr(x);
% Fn_est = Dyst_Emp(x,X);
% 
% X1 = Zad_1(100);
% Fn_est1 = Dyst_Emp(x,X1);
% 
% X2 = Zad_1(1000);
% Fn_est2 = Dyst_Emp(x,X2);
% 
% subplot(3,1,1);
% plot(x,Fn);
% hold on;
% stairs(x,Fn_est);
% title('Oryginalna dystrybuanta porownana z estymowana dystrybuanta empiryczna N=10','interpreter','latex');
% xlabel('x','interpreter','latex');
% ylabel('$F(x)$','interpreter','latex');
% legend('$F(x)$','$\widehat{F}_N(x)$','interpreter','latex');
% 
% subplot(3,1,2);
% plot(x,Fn);
% hold on;
% stairs(x,Fn_est1);
% title('Oryginalna dystrybuanta porownana z estymowana dystrybuanta empiryczna N=100','interpreter','latex');
% xlabel('x','interpreter','latex');
% ylabel('$F(x)$','interpreter','latex');
% legend('$F(x)$','$\widehat{F}_N(x)$','interpreter','latex');
% 
% subplot(3,1,3);
% plot(x,Fn);
% hold on;
% stairs(x,Fn_est2);
% title('Oryginalna dystrybuanta porownana z estymowana dystrybuanta empiryczna N=1000','interpreter','latex');
% xlabel('x','interpreter','latex');
% ylabel('$F(x)$','interpreter','latex');
% legend('$F(x)$','$\widehat{F}_N(x)$','interpreter','latex');

% N = 1000;
% wek = [];
% wek1 = [];
% wek2 = [];
% x1 = -0.5:0.02:0.5;
% x2 = -1:0.02:0.99;
% x3 = -2:0.02:0.99;
% 
% for i=1:10:N
%     wek(end+1) = Emp_err(10,x1,i);
%     wek1(end+1) = Emp_err(10,x2,i);
%     wek2(end+1) = Emp_err(10,x3,i);
% end
% N = 1:10:1000;
% 
% figure(1);
% plot(N,wek);
% hold on;
% plot(N,wek1);
% hold on;
% plot(N,wek2);
% title('Blad empiryczny estymatora dystrybuanty empirycznej','interpreter','latex');
% xlabel('N','interpreter','latex');
% ylabel('Err','interpreter','latex');
% legend('M=50','M=100','M=150','interpreter','latex');

% x = -2:0.01:3;
% [V,F] = Zad_5(x);
% F_id = Zrob_dystr(x);
% subplot(2,1,1);
% plot(x,V);
% hold on;
% plot(x,(F.*(1-F))/1000);
% title('Estymacja wariancji dystrybuanty empirycznej w funkcji x','interpreter','latex');
% xlabel('x','interpreter','latex');
% ylabel('Var(Fn(x))','interpreter','latex');
% legend('Var($\widehat{F}_N$)','Var = F(x)[1-F(x)]','interpreter','latex');
% subplot(2,1,2);
% plot(F,V);
% hold on;
% plot(F,(F.*(1-F))/1000)
% title('Estymacja wariancji dystrybuanty empirycznej w funkcji F(x)','interpreter','latex');
% xlabel('F(x)','interpreter','latex');
% ylabel('Var(Fn(x))','interpreter','latex');
% legend('Var($\widehat{F}_N$)','Var = F(x)[1-F(x)]','interpreter','latex');
