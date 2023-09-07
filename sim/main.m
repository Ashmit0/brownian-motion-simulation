%%
dt = .01 ; % time step
time = (1:100000) * dt   ;  % time array 
N = length( time ) ; % size of time data array 
%%
% defining constnats 
n = 1e-3 ; a = 1e-6 ; T = 300 ; kb = 1.38e-23 ; 
coef = 6 * pi * n * a ; D = kb * T / coef ; 

dx = sqrt( 2 * D * dt ) .* randn( N , 1 ) ; % Normal displacements samples 
%%
plot( time , dx )
xlabel('Time in sec')
ylabel('dx in m')
title('particle disaplcement from $t$ to $t + \mathrm{d}t$ ')
%%
x = cumsum( dx ); % Net displacement in the x axis
%%
plot( time , x )
xlabel('Time (sec)', Interpreter='latex')
ylabel('X (m)' , Interpreter= 'latex')
title('X Position vs Time ' , Interpreter= 'latex')
%%

set(0,'defaulttextInterpreter','latex')

dy = sqrt( 2 * D * dt ) .* randn( N , 1 ) ; 
y = cumsum( dy ); 


%%
clf ; 
plot( x , y ); 
xlabel("$x$ position (m)")
ylabel("$y$ position (m)")
title("Brownian Particel Trajectory")
%%

rsq = x.^2 + y.^2 ; % net displacement from the origin

%%
plot( time , rsq )
xlabel("Time (sec)")
ylabel(" $|r^2|$ ($m^2$)")
title("Net Displacement squared from the origin vs time ")
%%

pcount = 5 ; % we repeat the above for 20 particles  
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
clf ; hold on ; 
rsq_mean = zeros( N , 1 ); 
for i = 1:pcount 
    plot( time , particle{i}.rsq )
    rsq_mean = rsq_mean + particle{i}.rsq ; 
end
plot( time , rsq_mean / pcount  , 'color' , 'k' , 'LineWidth', 2 )
xlabel('Time (sec)')
ylabel("$|r^2|$")
title("Net displacement squared for 20 particels")
hold off ; 

%%
beta = xcorr( particle{1}.dx , 'coeff'); 
beta = beta( ( length(beta) + 1)/2 : end );
%%
clf ;  
plot( time , beta)
xlabel('Time Lag ($\tau$)')
ylabel('$ \langle \Delta x( t + \tau) \Delta x(t)  \rangle $ ')
title('Autocorrelation of Particle 1')
%%

beta = xcorr( particle{1}.dx , particle{1}.dy , 'coeff'); 
beta = beta( ( length(beta) + 1)/2 : end ); 
%%
clf ; 
plot( time , beta)
xlabel('Time Lag ($\tau$)')
ylabel('$ \langle \Delta x( t + \tau) \Delta y(t)  \rangle $ , ($m^2$) ')
title('Crosscorrelation of Particle 1')


rsq_mean = zeros( N , 1 ); 
for i = 1:pcount 
    rsq_mean = rsq_mean + particle{i}.rsq ; 
end
%%
clf ; hold on ; 

plot( time , rsq_mean / pcount  , 'color' , 'k' , 'LineWidth', 3 )
legend( "Mean" ,'AutoUpdate','off')
for i = 1:pcount 
    plot( time , particle{i}.rsq )
end
xlabel('Time (sec)')
ylabel("$|r^2|$ , ($m^2$)")
title("Net displacement squared for 20 particels")
hold off ; 
%%
msdx = zeros( N -1 ,1 ) ;
for i = 1:N-1
    temp =  particle{1}.x( i+1 : end ) - particle{1}.x(1 : end - i  );
    msdx(i) = mean( temp .* temp ) ; 
end
%%
clf ;
loglog( dt * (1: (N-1) ) , msdx )
title("Mean square displacement in the x direction for particle 1 in LogLog scale")
ylabel(' $\langle (x(t+\tau) - x(t))^2 \rangle$ ($m^2$) ')
xlabel('Time Lag, $\tau$ (sec)')

%%


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
        particle{j}.msdr(i) = mean( temp .* temp + temp2 .* temp2 );
        msdr(i) = msdr(i) + particle{j}.msdr(i) ; 
    end
end
msdx = msdx / pcount ; 

%%
pfit = polyfit( dt *(1:N-1)  , msdx  , 1 ) ; 
array = dt * (1 : N -1 ) ; 
%%
clf;
plot( dt * ( 1:N-1) , msdx  , LineWidth= 3 , Color= 'k' )
hold on 
plot( array  ,  pfit( 1 ) * array  + pfit( 2 ) )
legend( "Mean"  , "Fit" ,'AutoUpdate','off')
for i = 1:pcount
    plot( dt * ( 1:N-1) , particle{i}.msdx )
end

title("Ensemble Averaged Mean square displacement in the x direction for 20 particle in LogLog scale")
ylabel(' $\langle (x(t+\tau) - x(t))^2 \rangle$ ($m^2$) ')
xlabel('Time Lag, $\tau$ (sec) ')
hold off ;

%%
msdy = msdy / pcount ;
%%
pfit = polyfit(dt *(1:N-1),  msdy, 1 ) ; 
array = dt * (1 : N -1) ; 
%%
clf ; 
plot( dt * ( 1:N-1) , msdy  , LineWidth= 3 , Color= 'k' )
hold on 
plot( array , pfit( 1 )* array + pfit(2) , color = 'r' , LineWidth= 3 )
legend( "Mean" , "Linear Fit" ,'AutoUpdate','off')
for i = 1:pcount
    plot( dt * ( 1:N-1 ) , particle{i}.msdy )
end
title("Ensemble Averaged Mean square displacement in the y direction for 20 particle in LogLog scale")
ylabel(' $\langle (y(t+\tau) - y(t))^2 \rangle$ ($m^2$) ')
xlabel('Time Lag, $\tau$ (sec) ')
annotation('textbox',...
    [0.149842271293375 0.660044150110375 0.736908517350157 0.242825607064018],...
    'String',{ append( 'Acctual D = ' , string(D)) , ...
    append('Measured D (Slope)/2 = ' , string( pfit(1)/2)) ... 
    ,append( 'Error = ' , string((pfit(1)/2 - D )*100/D ) , '%' ) },...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);
hold off ;

%%
msdr = msdr / pcount ; 
%%
pfit = polyfit(dt *(1:N-1),  msdr , 1 ) ; 
array = dt * (1 : N -1) ; 
plot( dt*( 1:N-1) , msdr  , LineWidth= 3 , Color= 'k' )
hold on 
plot( array , pfit( 1 )* array + pfit(2) , color = 'r' , LineWidth= 3 )
legend( "Mean" , "Fit" ,'AutoUpdate','off')
for i = 1:pcount
    plot( dt * ( 1:N-1) , particle{i}.msdr )
end
title("Ensemble Averaged Mean square displacement for 20 particle in LogLog scale")
ylabel(' $\langle \Delta r^2  \rangle$ ($m^2$) ')
xlabel('Time Lag, $\tau$ (sec) ')
annotation('textbox',...
    [0.149842271293375 0.660044150110375 0.736908517350157 0.242825607064018],...
    'String',{ append( 'Acctual D = ' , string(D)) , ...
    append('Measured D (Slope)/4 = ' , string( pfit(1)/4)) ... 
    ,append( 'Error = ' , string((pfit(1)/4 - D )*100/D ) , '%' ) },...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);
hold off ; 
%%

loglog( time , rsq_mean /pcount , Color='k' , LineWidth=3)
legend( "Mean" ,'AutoUpdate','off')
hold on 
for i = 1:pcount 
    plot( time , particle{i}.rsq )
end
xlabel('Time (sec)')
ylabel("$|r^2|$ , ($m^2$)")
title("Net displacement squared for 20 particels on LogLog scale ")
%%
clf ; hold on ; 
xsq_mean = zeros( N , 1 ); 
for i = 1:pcount 
    xsq_mean = xsq_mean + particle{i}.x.^2 ; 
end
plot( time , xsq_mean / pcount  , 'color' , 'k' , 'LineWidth', 2 )
legend( "Mean" ,'AutoUpdate','off')
for i = 1:pcount 
    plot( time , particle{i}.x.^2  )
end
xlabel('Time (sec)')
ylabel("$\langle \Delta x^2 \rangle $ ($m^2$)")
title("Net x direction displacement squared for 20 particles")
hold off ; 


clf ;
loglog( time , xsq_mean / pcount  , 'color' , 'k' , 'LineWidth', 2 )
legend( "Mean" ,'AutoUpdate','off')
hold on ; 
for i = 1:pcount 
    plot( time , particle{i}.x.^2  )
end
xlabel('Time (sec)')
ylabel("$\langle \Delta x^2 \rangle $ ($m^2$)")
title("Net x direction displacement squared for 20 particles on LogLog scale")
hold off ; 



clf ; hold on ; 
ysq_mean = zeros( N , 1 ); 
for i = 1:pcount 
    ysq_mean = ysq_mean + particle{i}.y.^2 ; 
end
plot( time , ysq_mean / pcount  , 'color' , 'k' , 'LineWidth', 2 )
legend( "Mean" ,'AutoUpdate','off')
for i = 1:pcount 
    plot( time , particle{i}.y.^2  )
end
xlabel('Time (sec)')
ylabel("$\langle \Delta y^2 \rangle $  ($m^2$)")
title("Net y direction displacement squared for 20 particels")
hold off ;



clf ; 
loglog( time , ysq_mean / pcount  , 'color' , 'k' , 'LineWidth', 2 )
legend( "Mean" ,'AutoUpdate','off')
hold on ; 
for i = 1:pcount 
    plot( time , particle{i}.y.^2  )
end
xlabel('Time (sec)')
ylabel("$\langle \Delta y^2 \rangle $  ($m^2$)")
title("Net y direction displacement squared for 20 particles on LogLog Scale")
hold off ;


%%
beta = xcorr( particle{1}.dx , 'coeff'); 
beta = beta( (length(beta) +1 )/2 : end );
%%

alpha = zeros( N , 1 ) ; 
for i = 1:N 
alpha(i) = mean( particle{1}.dx( 1 : end - i + 1 ) .* particle{1}.dx( i : end )) ;
end
%%

%%
clf; hold on ; 
plot( dt * ( 1: N ) , alpha / max( alpha )  )
plot( dt * ( 1: N ) , beta  )
legend( "My function" , "Xcorr function")
xlabel("Time Lag ($\tau$) (sec)")
ylabel("$ \langle \Delta x (t + \tau ) \Delta x (t) \rangle $ ($m^2$)")
title("Autoorrelation Plot Differance")
hold off ; 
%%

% 
% clf ; hold on ; 
% cor = xcorr( particle{1}.dx , 'coeff'); 
% cor = cor( ( length(cor) + 1)/2 : end ); 
% lis = zeros( N , 1 );
% for i = 1:N 
%     lis(i) = mean( particle{1}.dx( 1 : end - i + 1 ) .* particle{1}.dx( i : end )); 
% end 
% plot( time , lis )
% plot( time , cor)
% legend( 'My Function' , 'xcorr Matlab Function')
% hold off ;



% movie = zeros( N , 1 ); 

% figh = figure('position',[100 100 850 600]); 
% for i = 1:100
%     clf ; hold on ; 
%     for j = 1 : pcount 
%         plot( particle{j}.x( 1 : i ) , particle{j}.y( 1 : i ) )
%     end
%     xlabel('X position ($m$)')
%     ylabel('Y position ($m$)')
%     title('Trajectory for 20 Brownian Particles')
% 
% %     drawnow
% 
%     movie(i) = getframe( figh ); 
% end

% 
% 
% 
% 
% wit = VideoWriter('traj' , 'MPEG-4'); 
% wit.FrameRate = 30 ; 
% open(wit);
% writeVideo( wit , movie ) 
% close(wit)