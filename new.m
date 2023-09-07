% % % % Random Time Step Simulation!%%%%%%
tmin = 10^( -3 ) ; tmax = 10^( 5 );
dt = -3 + 8*rand( 1 , 10000 ) ; 
%  dt = half( dt ) ; 
dt = 10.^( dt ) ;
time = cumsum( dt ) ; 

N = length( time ) ; 

n = 1e-3 ; a = 1e-6 ; T = 300 ; kb = 1.38e-23 ; 
coef = 6 * pi * n * a ; D = kb * T / coef ; 
set(0,'defaulttextInterpreter','latex')
%%
dx = sqrt( 2* D * dt ) .* randn( 1 , N  ) ; 
dy = sqrt( 2* D * dt ) .* randn( 1 , N  ) ; 
%%
plot( time , dx )
xlabel('Time in sec')
ylabel('dx in m')
title('particle disaplcement from $t$ to $t + \mathrm{d}t$ ')
%%
x = cumsum( dx ) ;  y = cumsum( dy ) ; 
plot( x , y ) 

%%
tp = (1/5) * (-15:25) ; 
count = zeros( length( tp ) , 1 ) ; 
msdx = zeros( length( tp ) - 1 , 1) ; 
msdy = zeros( length( tp) - 1  ,1 ) ;   
%%
% alpha = [] ;
for i = 1:N-1 
    for j = i+1:N
        tau = time( j ) - time( i ) ; 
        index =  floor(5*log10( tau )) + 16 ; 
%         index = round( index ) ;
%         alpha.append = index  ; 
        if index <= length( tp ) -1  
            count( index ) = count( index ) + 1 ; 
            msdx( index ) = msdx( index ) + ( x( j ) - x( i) )^2 ; 
            msdy( index ) = msdy( index ) + ( y( j ) - y( i ))^2 ; 
        end
    end
end
%%
for i = 1:length( tp ) -1 
    msdx( i ) = msdx( i ) / count( i) ; 
    msdy( i ) = msdy( i ) / count( i ) ; 
end
tau = zeros( length( tp ) -1 , 1 ) ; 
for i = 1:length( tp ) -1 
    tau( i ) = ( tp( i ) + tp( i + 1))/2 ; 
end 
tau = 10.^tau  ; 
%%
loglog( tau , msdx )
title('Mean squared displacement in the x direction for a single particle on LogLog scale')
ylabel('$\langle ( x( t + \tau ) - x(t) )^2 \rangle $ in $m^2$')
xlabel('Time Step $\tau$ in sec')
grid on ;
%%
loglog( tau , msdy )
title('Mean squared displacement in the y direction for a single particle on LogLog scale')
ylabel('$\langle ( y( t + \tau ) - y(t) )^2 \rangle $ in $m^2$')
xlabel('Time Step $\tau$ in sec')
grid on ; 
%%
pcount = 20  ; % we repeat the above for 20 particles  
particle = {} ; % empty list for storing the data 
for i = 1:pcount 
    particle{i} = struct(); 
    particle{i}.dx = sqrt( 2 * D * dt ) .* randn( 1 , N  ); 
    particle{i}.dy = sqrt( 2 * D * dt ) .*randn( 1 , N  ); 
    particle{i}.x = cumsum( particle{i}.dx); 
    particle{i}.y = cumsum( particle{i}.dy);
    particle{i}.rsq = particle{i}.x.^2 + particle{i}.y.^2 ; 
end 
%%
clf ; hold on ; 
for i = 1 : pcount 
    plot( time , particle{i}.x )
end 
xlabel('Time (sec)')
ylabel('X displacement ($m$)')
title("Displacement in the X-Direction vs Time for 20 Particles")
hold off ; 
%%


clf ; hold on ; 
for i = 1:pcount
    plot( particle{i}.x , particle{i}.y )
end 
xlabel('X position ($m$)')
ylabel('Y position ($m$)')
title('Trajectory for 20 Brownian Particles')

%%

count = zeros( length( tp ) , 1 ) ; 
msdx = zeros( length( tp ) - 1 , 1) ; 
msdy = zeros( length( tp) - 1  ,1 ) ;   
msdr = zeros( length( tp) - 1  ,1 ) ; 
for i = 1 : pcount 
    particle{i}.msdx = zeros( length( tp ) - 1 , 1) ; 
    particle{i}.msdy = zeros( length( tp ) - 1 , 1) ;
    particle{i}.msdr = zeros( length( tp ) - 1 , 1) ;
end 

for i = 1:N-1 
    for j = i+1:N
        tau = time( j ) - time( i ) ; 
        index =  floor(5*log10( tau )) + 16 ; 
        if index <= length( tp ) -1  
             count( index ) = count( index ) + 1 ; 
            for k = 1:pcount 
                particle{k}.msdx(index) = particle{k}.msdx(index) ...
                    + ( particle{k}.x( j ) - particle{k}.x( i) )^2 ; 
                particle{k}.msdy(index) = particle{k}.msdy(index) + ...
                    ( particle{k}.y( j ) - particle{k}.y( i ))^2;
                particle{k}.msdr( index ) = particle{k}.msdr( index )  ...
                    + ( particle{k}.x( j ) - particle{k}.x( i) )^2 + ( particle{k}.y( j ) - particle{k}.y( i) )^2 ; 
                msdx( index ) =  msdx( index ) + ( particle{k}.x( j ) - particle{k}.x( i) )^2 ; 
                msdy( index ) = msdy( index ) + ( particle{k}.y( j ) - particle{k}.y( i) )^2; 
                msdr( index ) = msdr( index ) +  ...
                    ( particle{k}.x( j ) - particle{k}.x( i) )^2 + ( particle{k}.y( j ) - particle{k}.y( i) )^2;
            end
        end
    end
end

%%
for i = 1:length( tp) -1 
    msdx(i) = msdx(i) / ( count(i) * pcount ) ; 
    msdy(i) = msdy(i) / ( count(i) * pcount ) ; 
    msdr(i) = msdr(i) / ( count(i) * pcount ) ; 
    for j = 1:pcount 
        particle{j}.msdx(i) = particle{j}.msdx(i)/count(i) ;
        particle{j}.msdy(i) = particle{j}.msdy(i)/count(i) ;
        particle{j}.msdr(i) = particle{j}.msdr(i)/count(i) ;
    end
end 
tau = zeros( length( tp ) -1 , 1 ) ; 
for i = 1:length( tp ) -1 
    tau( i ) = ( tp( i ) + tp( i + 1))/2 ; 
end 
tau = 10.^tau  ;
%%
for i = 1:length( tp) -1 
    for j = 1:pcount 
        particle{j}.msdy(i) = particle{j}.msdy(i)/count(i) ;
    end
end 
%%

poly = polyfit( tau , msdx , 1 ) ; 

mD = poly(1)/2 ; 
er = ( mD - D )*100/D ; 
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
title(' Diffusion constant from the linear MSDX plot (rejecting the last decade)');


plot( tau , msdx , LineWidth= 1.4 , Color= 'k'  )
hold on ; 
plot( tau , poly(1)*tau + poly(2) , LineWidth= 1.4 , Color= 'r' )
legend( "Mean" , 'Fit' ,'AutoUpdate','off')
for i = 1:pcount 
    plot( tau , particle{i}.msdx ) 
end 


box(axes1,'on');
grid(axes1,'on');
hold(axes1,'off');
% Create legend


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

loglog( tau , msdx , LineWidth= 2  , Color='k' )
legend( 'Mean'  , 'AutoUpdate' , 'off' ); 
hold on ; 
for i = 1:pcount 
    plot( tau , particle{i}.msdx  ) 
end 
title("Ensemble Averaged Mean square displacement in the x direction for 20 particle in LogLog scale")
ylabel(' $\langle (x(t+\tau) - x(t))^2 \rangle$ ($m^2$) ')
xlabel('Time Lag, $\tau$ (sec) ')
grid on ; 
hold off ;
%%
loglog( tau ,  msdy , LineWidth= 2  , Color='k' )
legend( 'Mean'  , 'AutoUpdate' , 'off' ); 
hold on ; 
for i = 1:pcount 
    plot( tau , particle{i}.msdy  ) 
end 
title("Ensemble Averaged Mean square displacement in the y direction for 20 particle in LogLog scale")
ylabel(' $\langle (y(t+\tau) - y(t))^2 \rangle$ ($m^2$) ')
xlabel('Time Lag, $\tau$ (sec) ')
grid on ; 
hold off ;
%%
% msdr = msdr *pcount ;
loglog( tau ,  msdr , LineWidth= 2  , Color='k' )
legend( 'Mean'  , 'AutoUpdate' , 'off' ); 
hold on ; 
for i = 1:pcount 
    plot( tau , particle{i}.msdr  ) 
end 
title("Ensemble Averaged Mean square displacement for 20 particle in LogLog scale")
ylabel(' $\langle \Delta r^2  \rangle$ ($m^2$) ')
xlabel('Time Lag, $\tau$ (sec) ')
grid on ; 

%%


poly = polyfit( tau , msdr , 1 ) ; 

mD = poly(1)/4 ; 
er = ( mD - D )*100/D ; 
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
title("Ensemble Averaged Mean square displacement for 20 particle in Linear scale")
ylabel(' $\langle \Delta r^2  \rangle$ ($m^2$) ')
xlabel('Time Lag, $\tau$ (sec) ')


plot( tau , msdr , LineWidth= 1.4 , Color= 'k'  )
hold on ; 
plot( tau , poly(1)*tau + poly(2) , LineWidth= 1.4 , Color= 'r' )
legend( "Mean" , 'Fit' ,'AutoUpdate','off')
for i = 1:pcount 
    plot( tau , particle{i}.msdr ) 
end 


box(axes1,'on');
grid(axes1,'on');
hold(axes1,'off');
% Create legend


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


%% 


% autocorelation 

auto = zeros( 1 , 1e4)  ; 
cross = zeros( 1 , 1e4 ) ; 
count = zeros( 1 , 1e4 ) ; 
timelag = 100*( 0: 1e4 - 1 ) ; 


for i = 1:N
    for j = i:N
        tau = ( time(j) - time(i) ); 
        index = ceil( tau / 1000  ) + 1 ; 
        if index < 1e4
            count( index ) = count( index ) + 1 ; 
            auto( index ) = auto( index ) + ...
                particle{1}.dx(i) * particle{1}.dx(j) ; 
            cross( index ) = cross( index ) + ... 
                particle{1}.dx( j )*particle{1}.dy( i ) ; 
        end 
    end 
end

t = 0 ; 

for i = 1:1e4
    if count( i ) >= 0 
        auto( i ) = auto( i ) / count( i ) ; 
        cross( i ) = cross( i ) / count( i) ; 
    else 
        t =1 
    end 
end 

%%

clf ; 
plot( timelag , auto  / max( auto )) 
xlabel("Time Lag ($\tau$) (sec)")
ylabel("$ \langle \Delta x (t + \tau ) \Delta x (t) \rangle $ ($m^2$)")
title("Autoorrelation Plot for particle 1")
grid on ; 
hold off ; 
 

%%


clf ; 
plot( timelag , cross / max( cross ) ) 
xlabel('Time Lag ($\tau$)')
ylabel('$ \langle \Delta x( t + \tau) \Delta y(t)  \rangle $ , ($m^2$) ')
title('Crosscorrelation of Particle 1')
grid on ; 
hold off ; 


 