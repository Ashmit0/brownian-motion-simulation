N = 1000 ; 
disp = randn( 1 , N ); 
plot( disp ) % normal sample with sd = 1 

hist( disp , 25 )

x = cumsum( disp ); 
plot( x )
ylabel("pos")
xlabel("time step")
title("Position in 1D vs time")


% for 2D motion we have 
x = cumsum( randn( N , 1 ) );
y = cumsum( randn( N , 1 ) );
plot( x , y )
xlabel("x position")
ylabel("y position")
title("Position vs time in 2D")

% square displacement of the particle 
dispsq = x.^2 + y.^2 ; 
plot( dispsq )

d = 1e-6 ; 
eta = 1e-3 ; 
kb = 1.38e-23 ; 
t = 293 ; 
D = kb*t / ( 3 * pi * eta * d ); 



tau = .1 ;  % time step 
time = tau * ( 1 : N ) ; % 100 sec time interval with step = .1 
k = sqrt( 2*D*tau ); 
dx = k*randn( N , 1 ); % dispacement vec 
dy = k*randn( N , 1 );
x = cumsum( dx );  % position vec 
y = cumsum( dy );
ddispsq = dx.^2 + dy.^2 ; % dquared change in pos 
dispsq = x.^2 + y.^2 ; % net sq displacement 
plot( x ,y )
title("2D brownian motion")


clf; 
hold on ; 
plot( time , (0 : N -1 )*2*k^2 , 'k'); 
plot( time , dispsq ); 
hold off; 
xlabel("time")
ylabel("sq_disp")
title("disp sq vs time")



simD  = mean( dispsq )/ ( 2 *4* tau ); 

simD/D


% displacement in bulk flow ! 
% 
% clf ; hold on ; 
% plot( time , ( 0 : N -1 )* 2 * k^2 , 'k') ;
% plot( time , dispsq );
% hold off ; 
% xlabel("Time")
% ylabel("disp sq ")
% title("disp sq vs time in bulk flow")

pcount = 10 ; 
N = 50 ; 
tau = .1 ; % time step 
time = ( 0: N -1 ) * tau ; % time array 
particle = {} ; % empty array ! to store the results ! 
for i = 1 : pcount 
    particle{i} = struct(); 
    particle{i}.dx = k * randn( N ,1 ); 
    particle{i}.dy = k * randn( N ,1 ); 
    particle{i}.x = cumsum( particle{i}.dx );
    particle{i}.y = cumsum( particle{i}.dy );
    particle{i}.drsq = particle{i}.dx.^2 + particle{i}.dy.^2;
    particle{i}.rsq = particle{i}.x.^2 + particle{i}.y.^2;
end 


clf ; 
hold on ; 
for i = 1:pcount
    plot( particle{i}.x , particle{i}.y );
end 
xlabel("X pos (m)")
ylabel("Y pos (m)")
title("combined particel track")
hold off ; 

% we now compute the ensamble average of the d sq values over time 

rsq_sum = zeros( N , 1 ); 
for i  = 1:pcount 
    rsq_sum = rsq_sum + particle{i}.rsq; 
end 
ensamble_avg = rsq_sum / pcount ; 

clf ; 
hold on ; 
plot( time , (0:N -1 )* 2 * k^2 , 'b' , 'LineWidth',3);  % theorectical line
plot( time , ensamble_avg , 'k' , 'LineWidth', 3 );  % averaged sq disp for all the particels 
legend( 'theoretical' , 'averaged'); 
for i = 1:pcount 
    plot( time , particle{i}.rsq ); 
end 
xlabel('Time');
ylabel("Disaplacement Squared (m^2)");
title("Disp Sq vs Time ");
hold off ; 

hist( randn( 1e6 , 1  ) ,100)
xlabel("value")
ylabel("freq")
title("histogram of 1e6 normal(0 , 1) samples")


% smapling uncertinty ! 

clf ; hold on ; 
for i = 1: 100 
    for j  = 1: 50 
        y = randn( 1 , i ); 
        m = mean( y ) ;
        plot( i , m , 'o' , 'color' , 'b');
    end 
end


% deviation decreases as inverse squared n ! 


plot( 1:100 , 1./sqrt( 1:100) , 'k' , 'LineWidth', 2 );
plot( 1:100 , -1./sqrt( 1:100) , 'k' , 'LineWidth', 2 );
hold off ; 

xlabel("Pop Size ");
ylabel("Pop Mean "); 
title("Pop mean of normal random distributuion")

dx = randn( 1e6  , 1 ); 
dy = randn( 1e6  , 1 );
drsq = dx.^2 + dy.^2 ; 
mean( drsq )
var( dx ) + var( dy )
clf ; 
hist( drsq  , 100 )
title("chi sq distribution with degree of freadom = 2 , 1e6 sample size ")


clf ; hold on ; 
for i = 1: 100 
    for  j = 1:50 
        dx = randn( 1, i);
        dy = randn( 1 ,i);
        m = mean( dx.^2 + dy.^2 );
        plot( i , m , 'o' , 'color' , 'b');
    end
end

plot( 1:100 , 2 + 2./sqrt( 1:100 ) , 'k' , 'LineWidth', 2 );
plot( 1:100 , 2 - 2./sqrt( 1:100 ) , 'k' , 'LineWidth', 2 );

hold off ; 

xlabel("Population Size (N)")
ylabel("Population Mean ")
title("Population Mean of Chi Squared Distribution with dof = 2 vs Population size")



% Autocorrelation !!! xcorr 

clf ; 
c = xcorr( particle{1}.dx , 'coeff');
alpha = c(( length(c))/2   : length(c)); 

plot( time , alpha )
xlabel( "Lag time ")
ylabel("Autocorrelation")



% Nice !!! 


% crosscorrelation !! 

clf ;
c = xcorr( particle{1}.dx , particle{2}.dx  , 'coeff');
c = c( ( length(c) + 1)/2 : length(c)); 
hold on; 
plot( time , c ); 
m = mean( c ); 
plot( time , ones( N , 1 )*m )
hold off ; 



 c_my = zeros( N ,1 );
 for i = 1:N
     c_my(i) = mean( particle{1}.dx( i : end ) .* particle{1}.dx( 1 : end - i + 1 ) );
 end 

hold on 
plot( time , c_my / max( c_my ))

c = xcorr( particle{1}.dx  , 'coeff');
c = c( ( length(c) + 1 )/2 : end );

plot( time , c )

legend( "my func" , "xcorr func ")

hold off 



