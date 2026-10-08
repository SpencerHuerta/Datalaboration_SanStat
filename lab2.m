close all;
clc;
clear;
%% Problem 1: Simulering av konfidensintervall
% Parametrar:
n = 25; %Antal matningar
mu = 2; %Vantevardet
sigma = 1; %Standardavvikelsen
alpha = 0.05;
%Simulerar n observationer for varje intervall
x = normrnd(mu, sigma,n,100); %n x 100 matris med varden
%Skattar mu med medelvardet
xbar = mean(x); %vektor med 100 medelvarden.x
%Beraknar de undre och ovre granserna
undre = xbar - norminv(1-alpha/2)*sigma/sqrt(n);
ovre = xbar + norminv(1-alpha/2)*sigma/sqrt(n);
% Problem 1: Simulering av konfidensintervall (forts.)
%Ritar upp alla intervall
figure(1);
hold on
for k=1:100
    if ovre(k) < mu % Rodmarkerar intervall som missar mu
        plot([undre(k) ovre(k)],[k k],'r')
        elseif undre(k) > mu
            plot([undre(k) ovre(k)],[k k],'r')
        else
            plot([undre(k) ovre(k)],[k k],'b')
        end
    end
%b1 och b2 ar bara till for att figuren ska se snygg ut.
b1 = min(xbar - norminv(1 - alpha/2)*sigma/sqrt(n));
b2 = max(xbar + norminv(1 - alpha/2)*sigma/sqrt(n));
axis([b1 b2 0 101]) %Tar bort outnyttjat utrymme i figuren
%Ritar ut det sanna vardet
plot([mu mu],[0 101],'g')
hold off
%% Problem 2: Maximum likelihood/Minsta kvadrat
figure(2);
M = 1e4;
b = 4;
x = raylrnd(b, M, 1);
hist_density(x, 40)
hold on
my_est_ml = sqrt(1 / (2 * M) * sum(x.^2)); % Skriv in din ML-skattning har
my_est_mk = sqrt(2 / pi) * mean(x);        % Skriv in din MK-skattning har
plot(my_est_ml, 0, 'r*','MarkerSize',14)
plot(my_est_mk, 0, 'g*','MarkerSize',14)
plot(b, 0, 'ro')
plot(0:0.1:6, raylpdf(0:0.1:6, my_est_mk), 'r')
hold off

figure(3);
plot(0:0.1:6, raylpdf(0:0.1:6, my_est_ml), 'r')
hold off

%% Problem 3: Konfidensintervall for Rayleighfordelning
figure(4);
load wave_data.mat
subplot(3,1,1), plot(y(1:100))
subplot(3,1,2), plot(y(1:end))
subplot(3,1,3), hist_density(y)

my_est = sqrt(2 / pi) * mean(y); % skattningen av b, mean(y) är skattningen av väntevärdet
n = length(y); % hur många mätvärden det finns
s = std(y); % matlabs inbyggda standardavvikelseberäknare
alpha3 = 0.05;

upper_bound = sqrt(2/pi) * (mean(y) + norminv(1-alpha3/2)*s/sqrt(n));
lower_bound = sqrt(2/pi) * (mean(y) - norminv(1-alpha3/2)*s/sqrt(n));


hold on % Gor sa att ploten halls kvar
plot(lower_bound, 0, 'b*')
plot(upper_bound, 0, 'g*')

plot(0:0.1:6, raylpdf(0:0.1:6, my_est), 'r')
hold off

%% Problem 4: Fordelningar av givna data
figure(5);
load birth.dat
% kategori 20 är rökning, 26 är alkoholvanor
x = birth(birth(:, 20) < 3, 3); %visar vilka som inte röker
y = birth(birth(:, 20) == 3, 3);% visar vilka som röker

subplot(2,2,1), boxplot(x), % plotten visar kvartiler. det inom lådan är 50% av värdena. sen visar den outliers.
axis([0 2 500 5000])
subplot(2,2,2), boxplot(y),
axis([0 2 500 5000])

subplot(2,2,3:4), ksdensity(x), % ksdensity visar typ ett utjämnat histogram av hur vanligt det är
hold on
[fy, ty] = ksdensity(y);
plot(ty, fy, 'r')
hold off
%test med annan kategori

figure(6); % alkoholvanor
x2 = birth(birth(:, 26) < 2, 3); %visar vilka som inte röker
y2 = birth(birth(:, 26) == 2, 3);% visar vilka som röker

subplot(2,2,1), boxplot(x2), % plotten visar kvartiler. det inom lådan är 50% av värdena. sen visar den outliers.
axis([0 2 500 5000])
subplot(2,2,2), boxplot(y2),
axis([0 2 500 5000])

subplot(2,2,3:4), ksdensity(x2), % ksdensity visar typ ett utjämnat histogram av hur vanligt det är
hold on
[fy2, ty2] = ksdensity(y2);
plot(ty2, fy2, 'r')
hold off

figure(7); % lätt moder
x3 = birth(birth(:, 23) < 1, 3); %visar vilka som inte röker
y3 = birth(birth(:, 23) == 1, 3);% visar vilka som röker

subplot(2,2,1), boxplot(x3), % plotten visar kvartiler. det inom lådan är 50% av värdena. sen visar den outliers.
axis([0 2 500 5000])
subplot(2,2,2), boxplot(y3),
axis([0 2 500 5000])

subplot(2,2,3:4), ksdensity(x3), % ksdensity visar typ ett utjämnat histogram av hur vanligt det är
hold on
[fy3, ty3] = ksdensity(y3);
plot(ty3, fy3, 'r')
hold off

%% Problem 5: Test av normalitet
figure(8);
normplot(birth(:,3));

figure(9);
qqplot(birth(:,3));
mu_births= mean(birth(:,3));
s_births = std(birth(:,3),1);
skyvheten_gamma = mean(((birth(:,3) - mu_births)/s_births).^3);
kurtosisen_kappa = mean(((birth(:,3) - mu_births)/s_births).^4);

n = numel(birth(:,3));
Jacque_beras = (n/ 6) * (skyvheten_gamma ^ 2 + (1/4) * (kurtosisen_kappa-3)^2)

nx = numel(x3);
ny = numel(y3);
[hx, px] = jbtest(x3, 0.05);
[hy, py] = jbtest(y3, 0.05);


fprintf('Jarque–Bera: h = %d, p = %.4g\n', hx, px);
fprintf('Jarque–Bera: h = %d, p = %.4g\n', hy, py);
% h = 1 betyder att nollhypotesen förkastas på 5%-nivån och det inte är
% normalfördelat.


%% Problem 6: Konfidensintervall för skillnad mellan väntevärden för födelsevikter

x6 = birth(birth(:, 20) < 3, 3);

y6 = birth(birth(:, 20) == 3, 3);

skillnad = mean(x6) - mean(y6) % alltså mu_x - mu_y

nx = length(x6);
ny = length(y6);

sx = std(x6);
sy = std(y6);

standardfel = sqrt(sx^2/nx + sy^2/ny);

alpha = 0.05;
z = norminv(1 - alpha/2);

undre6 = skillnad - z*standardfel;
ovre6  = skillnad + z*standardfel;




fprintf('Antal: icke-rökande = %d, rökande = %d\n', nx, ny);
fprintf('Skattad skillnad: %.1f gram\n', skillnad);
fprintf('95%% konfidensintervall: [%.1f, %.1f] gram\n', ...
    undre6, ovre6);

if undre6 > 0
    fprintf('Signifikant skillnad: högre medelvikt i gruppen icke-rökande.\n');
elseif ovre6 < 0
    fprintf('Signifikant skillnad: lägre medelvikt i gruppen icke-rökande.\n');
else
    fprintf('Ingen signifikant skillnad på 5%%-nivån.\n');
end

mu_xy = mean(x6) - mean(y6) % alltså mu_x - mu_y


%% Problem 7: enkel linjär regression
clf;
clc;
clear;
load moore.dat

X = ones(length(moore),2);
X(:,2) = moore(:,1);
y = log(moore(:,2));
[beta_hat,n1,n2,n3,stats] = regress(y,X);
fprintf('R2 = %.5f\n',stats(1))
fprintf('Antalet uppskattaddae transistorer år 2025: %.2f\n',exp([1,2025]*beta_hat))
figure(1);
hold on
plot((X*beta_hat),X(:,2),'r')
plot(y,X(:,2),'g*')
hold off
figure(2);
w=y;
res = w-X*beta_hat;
subplot(2,1,1), normplot(res)
subplot(2,1,2), hist(res)
