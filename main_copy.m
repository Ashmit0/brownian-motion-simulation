set(0,'defaulttextInterpreter','latex')

%%
dt = .01 ; % time step
time = (1:100000) * dt   ;  % time array 
N = length( time ) ; % size of time data array\
N2 = N -1 ; % would be frequently used later on 
tau = dt * ( 1 : N2 ) ;  % All the possible lag time values 

% Defining Constants 
n = 1e-3 ; a = 1e-6 ; T = 300 ; kb = 1.38e-23 ; 
coef = 6 * pi * n * a ; D = kb * T / coef ; 
%%

% X direction data generation 

dx = sqrt( 2 * D * dt ) * randn( N ,1 ) ; 
x = cumsum( dx ) ; 

plot( time , dx ) 
xlabel('Time in sec')
ylabel('dx in m')
title('particle disaplcement from $t$ to $t + \mathrm{d}t$ ')
clf ; 
plot( time , x ) ; 
xlabel('Time (sec)', Interpreter='latex')
ylabel('X (m)' , Interpreter= 'latex')
title('X Position vs Time ' , Interpreter= 'latex')
%%
% Y direction data generation 
% We assume X and Y data to be independent ! 

dy = sqrt( 2 * D * dt ) * randn( N ,1 ) ; 
y = cumsum( dy ) ; 

plot( time , dy ) 
xlabel('Time in sec')
ylabel('dy in m')
title('particle disaplcement from $t$ to $t + \mathrm{d}t$ ')
clf ; 
plot( time , y ) ; 
xlabel('Time (sec)', Interpreter='latex')
ylabel('Y (m)' , Interpreter= 'latex')
title('Y Position vs Time ' , Interpreter= 'latex')
%% 
rsq = x.^2 + y.^2 ; % net displacement from the origin
plot( time , rsq )
xlabel("Time (sec)")
ylabel(" $|r^2|$ ($m^2$)")
title("Net Displacement squared from the origin vs time ")
%%
pcount = 20 ; % we repeat the above for 20 particles  
particle = {} ; % empty list for storing the data 
for i = 1:pcount 
    particle{i} = struct(); 
    particle{i}.dx = sqrt( 2 * D * dt ) .*randn( N , 1 ); 
    particle{i}.dy = sqrt( 2 * D * dt ) .*randn( N , 1 ); 
    particle{i}.x = cumsum( particle{i}.dx); 
    particle{i}.y = cumsum( particle{i}.dy);
    particle{i}.rsq = particle{i}.x.^2 + particle{i}.y.^2 ; 
end 
%%
% X direction Displacement for 20 Particles 
clf ; hold on ; 
for i = 1 : pcount 
    plot( time , particle{i}.x )
end 
xlabel('Time (sec)')
ylabel('X displacement ($m$)')
title("Displacement in the X-Direction vs Time for 20 Particles")
hold off ; 
%%
% Autocorrelation of the x data set for particle I 
auto = zeros( N , 1 ) ; 
for i = 1:N  
    auto( i ) = mean( particle{1}.dx( i : end ).* particle{1}.dx( 1 : end - i + 1 ) )  ; 
end 
%%
plot( dt * (0:N-1) , auto / max( auto ) ) 
xlabel('Time Lag ($\tau$)')
ylabel('$ \langle \Delta x( t + \tau) \Delta x(t)  \rangle $ ')
title('Autocorrelation of Particle 1')
grid on 
%%
cross = zeros( N , 1 ) ; 
for i = 1:N  
    cross( i ) = mean( particle{1}.dx( i : end ).* particle{1}.dy( 1 : end - i + 1 ) )  ; 
end 
%%
plot( dt * (0:N-1) , cross / max( cross ) ) 
xlabel('Time Lag ($\tau$)')
ylabel('$ \langle \Delta x( t + \tau) \Delta y(t)  \rangle $ , ($m^2$) ')
title('Crosscorrelation of Particle 1')
grid on 
%%
% Mean Squared Displacement in the X direction 
msdx = zeros( N -1 ,1 ) ;
for i = 1:N-1
    temp =  particle{1}.x( i+1 : end ) - particle{1}.x(1 : end - i  );
    msdx(i) = mean( temp .* temp ) ; 
end
%%
clf ;
loglog( dt * (1: (N-1) ) , msdx )
title("Mean square displacement in the x direction for particle 1 in LogLog scale")
ylabel(' $\langle (x(t+\tau) - x(t))^2 \rangle$ ($m^2$) ' ) 
xlabel('Time Lag, $\tau$ (sec)')
%%
%Ensemble Mean Squared Displacement in the X direction

for i = 1:pcount
    particle{i}.msdx = zeros( N - 1,1 ); 
    particle{i}.msdy = zeros( N - 1,1 );
    particle{i}.msdr = zeros( N - 1,1 );
end
msdx = zeros( N -1 ,1 ); 
msdy = zeros( N -1 ,1 );
msdr = zeros( N -1 ,1 );

for j = 1 : pcount
    for i = 1:N-1
        temp = particle{j}.x( i + 1: end ) - particle{j}.x( 1 : end - i ) ;
        particle{j}.msdx(i) = mean( temp .* temp ); 
        msdx(i) = msdx(i) +  particle{j}.msdx(i);
        temp2 = particle{j}.y( i + 1: end ) - particle{j}.y( 1 : end - i ) ;
        particle{j}.msdy(i) = mean( temp2 .* temp2 );
        msdy(i) = msdy(i) +  particle{j}.msdy(i);
        particle{j}.msdr(i) = particle{j}.msdx(i) + particle{j}.msdy(i);
        msdr(i) = msdr(i) + particle{j}.msdr(i) ; 
    end
end
%%
clf ; 
loglog( tau , msdx / pcount  , LineWidth= 3 , Color= 'k')
legend( "Mean"  ,'AutoUpdate','off')
hold on ; 
for i = 1:pcount 
    plot( tau , particle{i}.msdx ) 
end 
title("Ensemble Averaged Mean square displacement in the x direction for 20 particle in LogLog scale")
ylabel(' $\langle (x(t+\tau) - x(t))^2 \rangle$ ($m^2$) ')
xlabel('Time Lag, $\tau$ (sec) ')
hold off ;
%%
clf ; 
loglog( tau , msdx / pcount  , LineWidth= 3 , Color= 'k')
legend( "Mean"  ,'AutoUpdate','off')
hold on ; 
for i = 1:pcount 
    plot( tau , particle{i}.msdx ) 
end 
title("Ensemble Averaged Mean square displacement in the x direction for 20 particle in LogLog scale")
ylabel(' $\langle (x(t+\tau) - x(t))^2 \rangle$ ($m^2$) ')
xlabel('Time Lag, $\tau$ (sec) ')
grid on ; 
hold off ;

%%

clf ; 
loglog( tau , msdy / pcount  , LineWidth= 3 , Color= 'k')
legend( "Mean"  ,'AutoUpdate','off')
hold on ; 
for i = 1:pcount 
    plot( tau , particle{i}.msdy ) 
end 
title("Ensemble Averaged Mean square displacement in the y direction for 20 particle in LogLog scale")
ylabel(' $\langle (y(t+\tau) - y(t))^2 \rangle$ ($m^2$) ')
xlabel('Time Lag, $\tau$ (sec) ')
grid on ; 
hold off ;

%% 



clf ; 
loglog( tau , msdy / pcount  , LineWidth= 3 , Color= 'k')
legend( "Mean"  ,'AutoUpdate','off')
hold on ; 
for i = 1:pcount 
    plot( tau , particle{i}.msdy ) 
end 
title("Ensemble Averaged Mean square displacement in the y direction for 20 particle in LogLog scale")
ylabel(' $\langle (y(t+\tau) - y(t))^2 \rangle$ ($m^2$) ')
xlabel('Time Lag, $\tau$ (sec) ')
grid on ; 
hold off ;
%%

clf ; 
loglog( tau , msdr / pcount  , LineWidth= 3 , Color= 'k')
legend( "Mean"  ,'AutoUpdate','off')
hold on ; 
for i = 1:pcount 
    plot( tau , particle{i}.msdr ) 
end 
title("Ensemble Averaged Mean square displacement for 20 particle in LogLog scale")
ylabel(' $\langle \Delta r^2  \rangle$ ($m^2$) ')
xlabel('Time Lag, $\tau$ (sec) ')
grid on ; 
hold off ;

%% 
msdx = msdx / pcount ; 
msdy = msdy / pcount ; 
msdr = msdr / pcount ; 

%%
figure1 = figure;
index = find( tau <= 100 , 1 , 'last') ; 
xfit = polyfit( tau( 1:index) , msdx( 1:index) , 1 ) ; 
plot( tau(1:index) , msdx( 1 : index) ) 
hold on ; 
plot( tau( 1 : index ) , xfit(1)*tau( 1:index) + xfit(2))
legend( "msdx" , "Fit")
grid on ; 
title(" Diffusion constant from the MSDX plot")
ylabel("$\langle ( x(t + \tau ) - x(t))^2 \rangle$ in $m^2$")
xlable(" Time Lag $\tau$ in sec")
mD = xfit(1)/2 ; 
er = (mD - D)/D ; 
annotation(figure1,'textbox',...
    [0.177785714285714 0.735714285714286 0.363285714285714 0.126190476190479],...
    'String',{ append('acctual $D$ =' , string(D)), ...
    append( 'measured $D$($\frac{\text{slope}}{2}$) =' , string(mD)),...
    append( 'error =' , string(er) , '$\%$') },...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

%%

mD = xfit(1)/2 ; 
er = (mD - D)/D ; 
%CREATEFIGURE(X1, YMatrix1)
%  X1:  vector of plot x data
%  YMATRIX1:  matrix of plot y data

%  Auto-generated by MATLAB on 03-Apr-2023 19:55:53

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple line objects using matrix input to plot
% plot1 = plot(X1,YMatrix1);
% set(plot1(1),'DisplayName','msdx');
% set(plot1(2),'DisplayName','Fit');

% Create ylabel
ylabel('$\langle ( x(t + \tau ) - x(t))^2 \rangle$ in $m^2$');

% Create title
title(' Diffusion constant from the MSDX plot (rejecting the last decade)');

plot( tau(1:index) , msdx( 1 : index) ) 
hold on ; 
plot( tau( 1 : index ) , xfit(1)*tau( 1:index) + xfit(2))

box(axes1,'on');
grid(axes1,'on');
hold(axes1,'off');
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.735714285714286 0.200595238095239 0.132142857142857 0.0803571428571422],...
    'Interpreter','latex');

% Create textbox
annotation(figure1,'textbox',...
    [0.177785714285714 0.735714285714286 0.363285714285714 0.126190476190479],...
    'String',{append('acctual D =' , string(D)), ...
    append( 'measured D (slope/2) =' , string(mD)),...
    append( 'error =' , string(er) ,'%')},...
    'Interpreter','latex',...
    'FontSize',10,...
    'FitBoxToText','off',...
    'EdgeColor','none');

%% 

xfit = polyfit( tau , msdx , 1 ) ; 
mD = xfit(1)/2 ; 
er =  ( mD - D )*100/D ; 

figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple line objects using matrix input to plot
% plot1 = plot(X1,YMatrix1);
% set(plot1(1),'DisplayName','msdx');
% set(plot1(2),'DisplayName','Fit');

% Create ylabel
ylabel('$\langle ( x(t + \tau ) - x(t))^2 \rangle$ in $m^2$');

% Create title
title(' Diffusion constant from the MSDX plot');

plot( tau , msdx ) 
hold on ; 
plot( tau , xfit(1)*tau + xfit(2))

box(axes1,'on');
grid(axes1,'on');
hold(axes1,'off');
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.735714285714286 0.200595238095239 0.132142857142857 0.0803571428571422],...
    'Interpreter','latex');

% Create textbox
annotation(figure1,'textbox',...
    [0.177785714285714 0.735714285714286 0.363285714285714 0.126190476190479],...
    'String',{append('acctual D =' , string(D)), ...
    append( 'measured D (slope/2) =' , string(mD)),...
    append( 'error =' , string(er) ,'%')},...
    'Interpreter','latex',...
    'FontSize',10,...
    'FitBoxToText','off',...
    'EdgeColor','none');


%% 


index = find( tau >= 100 , 1 , 'last') ; 
ntau = tau( 1 : index ) ; 
nmsdr = msdr( 1 : index ) ; 

xfit = polyfit( ntau , nmsdr, 1 ) ; 
mD = xfit(1)/4 ; 
er =  ( mD - D )*100/D ; 

figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple line objects using matrix input to plot
% plot1 = plot(X1,YMatrix1);
% set(plot1(1),'DisplayName','msdx');
% set(plot1(2),'DisplayName','Fit');

% Create ylabel
ylabel('$\langle ( x(t + \tau ) - x(t))^2 \rangle$ in $m^2$');

% Create title
title(' Diffusion constant from the MSDX plot (rejecting the last decade)');

plot( ntau , nmsdr ) 
hold on ; 
plot( ntau , xfit(1)*ntau + xfit(2))

box(axes1,'on');
grid(axes1,'on');
hold(axes1,'off');
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.735714285714286 0.200595238095239 0.132142857142857 0.0803571428571422],...
    'Interpreter','latex');

% Create textbox
annotation(figure1,'textbox',...
    [0.177785714285714 0.735714285714286 0.363285714285714 0.126190476190479],...
    'String',{append('acctual D =' , string(D)), ...
    append( 'measured D (slope/4) =' , string(mD)),...
    append( 'error =' , string(er) ,'%')},...
    'Interpreter','latex',...
    'FontSize',10,...
    'FitBoxToText','off',...
    'EdgeColor','none');


