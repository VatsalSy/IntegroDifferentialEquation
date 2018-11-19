%% Testing the IntegroDifferential Equation solver
clc 
clear
close all
%% Case 2: a''[t] == a[t] - \int_0^t a[x]dx With a[0] = 0
tmax = 2*pi;
t = linspace(0,tmax,500);
%% Theoritical Solution (Using MATHEMATICA)
% ythGuess = t;
yth = (1/3) * exp(-t/2).*(exp(3*t/2) - cos(sqrt(3)*t/2) + sqrt(3)*sin(sqrt(3)*t/2));

figure1 = figure('visible','on','WindowState','fullscreen','Color',[1 1 1]);
axes1 = axes('Parent',figure1);
hold(axes1,'on');
% plot(t, ythGuess,'r-','MarkerSize',20,'LineWidth',2,'DisplayName','Theoritical');
plot(t, yth,'r-','MarkerSize',20,'LineWidth',2,'DisplayName','Theoritical');

%% MATLAB Solution
odeOptions = odeset('AbsTol',1e-8,'RelTol',1e-8);
IC = [0; 1];  % IC = [a(0); a'(0)]
[t,guess] = ode45(@model0,t,IC,odeOptions);
guess = [t guess]; 
counter = 1; err = 1e4; Tol = 1e-8;
% [t,a] = ode113(@model, t, IC, odeOptions, guess);
while err > Tol
    [t,a] = ode113(@model, t, IC, odeOptions, guess);
    err = sum((a(:,1)-guess(:,2)).^2);
    fprintf('  %4i: %8.2e\n',[counter err]);
    wt = 0.5;
    guess = [t (1-wt)*guess(:,2) + wt*a(:,1) guess(:,3)];
    counter = counter+1;
end

nskip = 20;
% plot(t(1:nskip:end), guess(1:nskip:end,1),'r*','MarkerSize',20,'LineWidth',2,'DisplayName','Matlab');
plot(t(1:nskip:end), a(1:nskip:end,1),'r*','MarkerSize',20,'LineWidth',2,'DisplayName','Matlab');

legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.45 0.75 0.263939012799944 0.0847176079734219],...
    'Interpreter','latex',...
    'FontSize',30,...
    'EdgeColor',[1 1 1]);
box(axes1,'on');
set(axes1,'FontName','times new roman','FontSize',30,'FontWeight','bold',...
    'LineWidth',3);
axis square
xlim([0.0 tmax])
% ylim([0.0 1.0])
xlabel('\boldmath{$t$}','LineWidth',2,'FontWeight','bold','FontSize',50,...
            'FontName','times new roman',...
            'Interpreter','latex');
ylabel('\boldmath{$a$}','LineWidth',2,'FontWeight','bold','FontSize',50,...
    'FontName','times new roman',...
    'Interpreter','latex'); 

annotation(figure1,'textbox',...
    [0.45 0.53 0.25 0.10],...
    'String','$\ddot{a} = \int_0^ta(x)dx$',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'HorizontalAlignment','center',...
    'FontSize',50,...
    'FitBoxToText','off');

function dadt = model0(~,a)
    % Generating Initial Guess
    dadt = [a(2); 0];
end

function dadt = model(t,a,guess)
% Inside the while loop, getting to the actual solution
% Integro-differential equation
ys = @(s) interp1(guess(:,1), guess(:,2), s);
dadt = [a(2); integral(@(s) ys(s), 0, t,'RelTol',1e-8,'AbsTol',1e-13)];
end
