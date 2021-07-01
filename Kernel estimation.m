clear;
close all;
clc;

%%%%%%%%%%%%%%%%% ZAD 1 i ZAD 2
% a = -10;
% b = 10;
% x = a:0.1:b;
% mu = 1;
% sigma = 1;
% jadro_nor = @(x) 1/(sigma*sqrt(2*pi))*exp((-(x - mu).^2)/(2*sigma^2));
% jadro_box = @(x) (abs(x) <= 0.5) .* (1) + (abs(x) > 0.5) .* (0);
% jadro_epa = @(x) ((x >= -1) && (x <= 1)) .* (3/4 *(1 - x^2)) + ((x < -1) && (x > 1)) .* 0;
% jadro_tri = @(x) ((x >= -1) && (x <= 1)) .* (70/81 * ((1 - abs(x)^3)^3))  + ((x < -1) && (x > 1)) .* 0;
% % h = 0.01;
% %Xn = Rozklad_Normalny(a,b,sigma,mu, 1000);
% % Xn = Box_Muller(10000,1,1);
% % Est = Estymator_jadrowy(jadro_box, h, Xn, x);
% 
% x = -1:0.001:1;
% %Xn = Box_Muller(10000,1,1);
% %Xn = Wzium();
% Xn = Troj();
% h1 = 0.05;
% h2 = 0.3;
% h3 = 0.5;
% h4 = 1;
% e1 = Estymator_jadrowy(jadro_box, h1, Xn, x);
% e2 = Estymator_jadrowy(jadro_box, h2, Xn, x);
% e3 = Estymator_jadrowy(jadro_box, h3, Xn, x);
% e4 = Estymator_jadrowy(jadro_box, h4, Xn, x);
% 
% figure(1);
% plot(x,e1);
% hold on;
% plot(x,e2);
% hold on;
% plot(x,e3);
% hold on;
% plot(x,e4);
% title('Jadrowy estymator gestosci - gestosc trojkatna','interpreter','latex');
% xlabel('x','interpreter','latex');
% ylabel('$\widehat{f}(x)$','interpreter','latex');
% legend('h=0.05','h=0.3','h=0.5','h=1','interpreter','latex');

%%%%%%%%%%%%%%%% ZAD 3
% a = -10;
% b = 10;
% x = a:0.1:b;
% %x = -2:0.01:2;
% jadro_nor = @(x) 1/(sqrt(2*pi))*exp((-(x).^2)/(2));
% jadro_box = @(x) (abs(x) <= 0.5) .* (1) + (abs(x) > 0.5) .* (0);
% jadro_epa = @(x) ((x >= -1) && (x <= 1)) .* (3/4 *(1 - x^2)) + ((x < -1) && (x > 1)) .* 0;
% jadro_tri = @(x) ((x >= -1) && (x <= 1)) .* (70/81 * ((1 - abs(x)^3)^3))  + ((x < -1) && (x > 1)) .* 0;
% h = 0.1;
% Xn = Box_Muller(10000,1,1);
% %Xn = Wzium();
% %Xn = Troj();
% E1 = Estymator_jadrowy(jadro_nor, h, Xn, x);
% E2 = Estymator_jadrowy(jadro_box, h, Xn, x);
% E3 = Estymator_jadrowy(jadro_epa, h, Xn, x);
% E4 = Estymator_jadrowy(jadro_tri, h, Xn, x);
% 
% figure(1);
% plot(x,E1);
% hold on;
% plot(x,E2);
% hold on;
% plot(x,E3);
% hold on;
% plot(x,E4);
% title('Jadrowy estymator gestosci','interpreter','latex');
% xlabel('x','interpreter','latex');
% ylabel('$f(x)$','interpreter','latex');
% legend('Gauss','Prost','Epa','Tri','interpreter','latex');

%%%%%%%%%%%%%%%% ZAD 4
% jadro_epa = @(x) ((x >= -1) && (x <= 1)) .* (3/4 *(1 - x^2)) + ((x < -1) && (x > 1)) .* 0;
% jadro_box = @(x) (abs(x) <= 0.5) .* (1) + (abs(x) > 0.5) .* (0);
% tab = [];
% x = -2:0.04:1.99;
% sigma = 1;
% mu = 1;
% f = @(x) 1/(sigma*sqrt(2*pi))*exp((-(x - mu).^2)/(2*sigma^2));
% y = [];
% for i=1:1:100
%     y(end+1) = f(x(i));
% end
% 
% for i=0.01:0.04:4
%     tab(end+1) = Emp_err(10,jadro_box,i,y,x);
% end
% 
% h_zak = 0.01:0.04:4;
% plot(h_zak,tab);
% title('Wplyw parametru wygladzania na jakosc estymacji - blad empiryczny','interpreter','latex');
% xlabel('h','interpreter','latex');
% ylabel('$Err(h)$','interpreter','latex');

%%%%%%%%%%%%%%%% ZAD DODATKOWE
jadro_epa = @(x) ((x >= -1) && (x <= 1)) .* (3/4 *(1 - x^2)) + ((x < -1) && (x > 1)) .* 0;
jadro_box = @(x) (abs(x) <= 0.5) .* (1) + (abs(x) > 0.5) .* (0);
%x = Box_Muller(1500,0,1);
g = @(x) ((1<=x)&&(x<=2)) .* (1);
x = Generuj_liczby(g,5000,0,3);
setGlobal(x);

options = optimset('Display','iter');
[h_opt, fval] = fminbnd(@Cross_validation,0,3,options);

zak = 0:0.01:3;
fn = Estymator_jadrowy(jadro_box,h_opt,x,zak);

y = [];
for i=1:1:length(zak)
   y(end+1) = g(zak(i)); 
end

plot(zak,y);
histogram(x,'Normalization','pdf');
hold on;
plot(zak,fn);
hold on;
%x = Wzium();
% setGlobal(x);
% J_n = [];
% for h=0.001:0.001:1
%      J_n(end+1) = Cross_validation(h);
% end
% 
% h0 = 0.0001;
% lb = -Inf; 
% ub = Inf;
% options = optimoptions('fmincon', 'MaxFunctionEvaluations', 10E4);
% [h, fval] = fmincon(@Cross_validation, h0,[],[],[],[],lb,ub,[],options);


% figure(1);
% % subplot(2,1,1);
% plot(h_n, J_n);
% hold on;
% plot(h, fval, 'r*');
% title('Przebieg wielkosci $\hat{J}_{h_n}$ w zaleznosci od h, z zaznaczonym optimum','interpreter','latex');
% xlabel('h','interpreter','latex');
% ylabel('$\hat{J}_{h_n}$','interpreter','latex');
% 
% zak = 0:0.001:1;
% fn = Estymator_jadrowy(jadro_epa,h,x,zak);
% 
% figure(1);
% histogram(x,'Normalization','probability');
% hold on;
% plot(zak, fn);
% title('Wyestymowana gestosc prawdopodobienstwa przy optymalnym wspolczynniku wygladzania','interpreter','latex');
% xlabel('x','interpreter','latex');
% ylabel('f(x)','interpreter','latex');
% 
% sigma = 1;
% mu = 0;
% f = @(x) 1/(sigma*sqrt(2*pi))*exp((-(x - mu).^2)/(2*sigma^2));
% f = @(x) ((x>0) && (x<=1/100)) .* (50) + ((x>1/100) && (x<=1)) .* (100/198);
% y = [];
% for i=1:1:length(zak)
%     y(end+1) = f(zak(i));
% end
% figure(2);
% plot(zak,y);
% hold on;
% plot(zak,fn);
% title('Porownanie wyestymowanej gestosci przy optymalnym wspolczynniku wygladzania z oryginalna gestoscia','interpreter','latex');
% xlabel('x','interpreter','latex');
% ylabel('f(x)','interpreter','latex');

function Rozklad = Troj()
    
    fun = @(x)((-1<x) && (x<=0)) .*(x+1)+((0<x) && (x<1)).*(-x+1);
    U1 = 2*(rand(1,10000))-1;
    U2 = 1*rand(1,10000);
    
    Rozklad = [];
    for i=1:1:length(U1)
        if(U2(i) <= fun(U1(i)))
            Rozklad(end+1) = U1(i);
        end
    end
end

function Rozklad = Wzium()
    fun = @(x) ((x>0) && (x<=1/100)) .* (50) + ((x>1/100) && (x<=1)) .* (100/198);
    U1 = rand(1,80000);
    U2 = 50*rand(1,80000);
    
    Rozklad = [];
    for i=1:1:length(U1)
       if(U2(i) <= fun(U1(i)))
          Rozklad(end+1) = U1(i); 
       end
    end
end

% METODA CROSS VALIDATION
function J_hN = Cross_validation(h)
    global Xn;
    x = 0:0.01:3;
    jadro = @(x) ((x >= -1) && (x <= 1)) .* (3/4 *(1 - x^2)) + ((x < -1) && (x > 1)) .* 0;
    
    fn = Estymator_jadrowy(jadro, h, Xn, x);
    Calka_kwadratu_estymatora = trapz(x, fn.^2);
    fn_k = [];
    
    for k=1:1:length(Xn) % ESTYMATOR LEAVE ONE OUT
       suma_wewnetrzna = 0;
       for i=1:1:length(Xn)
          if(i ~= k)
             K = jadro((Xn(i) - Xn(k))/h);
             suma_wewnetrzna = suma_wewnetrzna + K;
          end
       end
       fn_k(end+1) = 1/((length(Xn)-1) * h) * suma_wewnetrzna;
    end
    suma_zewnetrzna = sum(fn_k);
    Wartosc_oczekiwana = 1/(length(Xn)) * suma_zewnetrzna;
    J_hN = Calka_kwadratu_estymatora - 2 * Wartosc_oczekiwana;
end

function setGlobal(x)
    global Xn;
    Xn = x;
end

% BŁĄD EMPIRYCZNY
function Err = Emp_err(L, jadro, h, y,x)
    sumaL = 0;
    sumaM = 0;
    M = length(x); 
    for i=1:1:L
       Xn = Box_Muller(M,1,1);
       Est = Estymator_jadrowy(jadro,h,Xn,x);
       for j=1:1:M
           sumaM = sumaM + (Est(j) - y(j))^2; 
       end
       sumaL = sumaL + sumaM;
       sumaM = 0;
    end
    Err = 1/(L*M) * sumaL;
end

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

% GENERATOR LICZB ROZKŁAD NORMALNY 2
function Rozklad = Rozklad_Normalny(a,b,sigma,mu, iter)
    f = @(x) 1/(sigma*sqrt(2*pi))*exp((-(x - mu).^2)/(2*sigma^2));
    Rozklad = Generuj_liczby(f,iter,a,b);
end

% GENERATOR LICZB ROZKŁAD NORMALNY
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

% GENERATOR Z KAŻDEGO ROZKŁADU
function Rozklad = Generuj_liczby(fun, iter, a, b)
    x = a:0.01:b;
    y = zeros(1,length(x));
    for i=1:1:length(x)
        y(i) = fun(x(i));
    end
    d = ceil(max(y));
    Rozklad = [];
    going = true;
    while(going)
        U1 = a + (b - a) * rand(1,1);
        U2 = d * rand(1,1);
        if(U2 <= fun(U1))
            Rozklad(end+1) = U1;
        end
        if(length(Rozklad) == iter)
            going = false;
        end
    end
end