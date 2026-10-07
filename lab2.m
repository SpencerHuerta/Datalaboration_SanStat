clf;
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
%% Problem 1: Simulering av konfidensintervall (forts.)
%Ritar upp alla intervall
figure(1)
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
clf;
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
