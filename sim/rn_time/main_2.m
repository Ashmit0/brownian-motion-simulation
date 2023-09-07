set(0,'defaulttextInterpreter','latex')

maxtime = 1000 ; 
dt = 3 * randn( 500 , 1 );
dt = abs( dt );
time = cumsum( dt );
index = find( time > 1000  , 1 ); 
time = time( 1 : index -1); 
dt = dt( 1 : index -1 ); 


% defining constnats 
n = 1e-3 ; a = 1e-6 ; T = 300 ; kb = 1.38e-23 ; 
coef = 6 * pi * n * a ; D = kb * T / coef ; 
sigma = sqrt( 2 * D * dt ) ; 
N = length( dt );
dx = zeros( N ,1); 
dy = zeros( N ,1); 
dr = zeros( N ,1);
for i = 1:N 
    dx(i) = sigma(i) * randn ; 
    dy(i) = sigma(i) * randn ; 
end 

x = cumsum( dx ) ; 
y = cumsum( dy ) ; 


clf ; 
plot( time , dx )
ylabel('$ \mathrm{d}x$ ($m$)')
xlabel('Time (sec)')
title('Particle displacement in x-direction at random time steps')


plot( x , y ) 
xlabel("$x$ position (m)")
ylabel("$y$ position (m)")
title("Brownian Particel Trajectory")



rsq = x.*x + y.*y; 
plot( time , rsq )



timelag = 5 * (1:200) ; 
count = zeros( 200 , 1 ); 
msdx = zeros( 200 , 1 ); 
msdy = zeros( 200 , 1 ); 


for i = 1:N-1
    for j = i+1:N 
        k = ceil( (time(j) - time(i))/5 ) ; 
        count( k ) = count( k ) + 1 ; 
        msdx( k ) = msdx(k) +  (x(i) - x(j) )^2 ; 
        msdy( k ) = msdy(k) + ( y(i) - y(j) )^2 ; 
    end
end 

for i = 1:200 
    msdx(i) = msdx(i)/count(i) ; 
    msdy(i) = msdy(i)/count(i) ; 
end

loglog( timelag , msdx )
hold on ; 
plot( timelag , msdy ) 
xlabel('Time Lag ($\tau$) in sec')
ylabel('$ \langle \Delta \alpha ^2 \rangle $')
legend('$\alpha = x $' , '$\alpha = y $' )
title('Msd in the x and y direction')
hold off ; 


particle = {} ; % empty list 
pcount = 20 ; 
for i = 1:20 
    particle{i} = struct() ; 
    particle{i}.dx = zeros( N ,1); 
    particle{i}.dy = zeros( N ,1); 
    for j = 1:N 
        particle{i}.dx(j) = sigma(j) * randn ; 
        particle{i}.dy(j) = sigma(j) * randn ;
    end
    particle{i}.x = cumsum( particle{i}.dx );
    particle{i}.y = cumsum( particle{i}.dy );
    particle{i}.msdx = zeros( 200 , 1 ) ; 
    particle{i}.msdy = zeros( 200 , 1 ) ; 
end 
 


enmsdx = zeros( 200 , 1 ) ; 
enmsdy = zeros( 200 , 1 ) ; 
enmsdr = zeros( 200 , 1 ) ; 
for k = 1: pcount 
    count = zeros(200,1) ;
    for i = 1:N-1 
        for j = i+1:N 
            index =  ceil( ( time(j) - time(i))/5 ) ; 
            particle{k}.msdx(index) = particle{k}.msdx(index) + ...
                ( particle{k}.x(j) - particle{k}.x(i))^2 ;
            count(index) = count( index ) + 1 ; 
            particle{k}.msdy(index) = particle{k}.msdy(index) + ...
                ( particle{k}.y(j) - particle{k}.y(i))^2 ; 
        end
    end
    for i =  1:200 
        if count(i) > 0 
        particle{k}.msdx(i) = particle{k}.msdx(i) / count(i) ; 
        particle{k}.msdy(i) = particle{k}.msdy(i) / count(i) ; 
        end 
    end
    particle{k}.msdr = particle{k}.msdx + particle{k}.msdy ;
    enmsdy = enmsdy + particle{k}.msdy ; 
    enmsdx = enmsdx + particle{k}.msdx ; 
    enmsdr = enmsdr + particle{k}.msdr ; 
end

enmsdx = enmsdx / pcount ; 
enmsdy = enmsdy / pcount ; 
enmsdr = enmsdr / pcount ; 


%%
clf ; 
loglog( timelag , enmsdx , LineWidth= 2 , Color = 'k') 
legend( "Mean" ,'AutoUpdate','off')
hold on ; 
for i = 1:pcount 
    plot( timelag , particle{i}.msdx )
end 
title('Ensembeld MSD in the x direction for 20 particles')
xlabel('Time Lag ($\tau$) in sec')
ylabel('$\langle (x(t+\tau) - x(t))^2 \rangle $ in $m^2$')
%%
clf ; 
loglog( timelag , enmsdy , LineWidth= 2 , Color = 'k') 
legend( "Mean" ,'AutoUpdate','off')
hold on ; 
for i = 1:pcount 
    plot( timelag , particle{i}.msdy )
end 
title('Ensembeld MSD in the y direction for 20 particles')
xlabel('Time Lag ($\tau$) in sec')
ylabel('$\langle (y(t+\tau) - y(t))^2 \rangle $ in $m^2$')
%%
clf ; 
loglog( timelag , enmsdr , LineWidth= 2 , Color = 'k') 
legend( "Mean" ,'AutoUpdate','off')
hold on ; 
for i = 1:pcount 
    plot( timelag , particle{i}.msdr )
end 
title('Ensembeld MSD for 20 particles')
xlabel('Time Lag ($\tau$) in sec')
ylabel('$\langle |\vec{r}(t+\tau) - \vec{r}(t)|^2 \rangle $ in $m^2$')
%%
clf ;
cor = zeros( N ,1 ) ; 
for i = 1:N 
    temp = particle{1}.dx(i : end ) .* particle{1}.dx( 1 : end - i + 1 ) ; 
    cor(i) = mean( temp .* temp ) ; 
end
plot( time , cor /max( cor ) ) ; 
title('Autocorrelation for the x direction')
xlabel('Time Lag ($ \tau $) in sec')
ylabel('$\langle \Delta x(t + \tau ) \Delta x(t)\rangle $ in $m^2$')
%%
clf ;
cor = zeros( N ,1 ) ; 
for i = 1:N 
    temp = particle{1}.dx(i : end ) .* particle{1}.dy( 1 : end - i + 1 ) ; 
    cor(i) = mean( temp .* temp ) ; 
end
plot( time , cor /max( cor ) ) ; 
title('Crosscorrelation between the x and y direction')
xlabel('Time Lag ($ \tau $) in sec')
ylabel('$\langle \Delta x(t + \tau ) \Delta y(t)\rangle $ in $m^2$')
%%
