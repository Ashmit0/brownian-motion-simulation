function[ index ] = func( time , i , j  )
    if time( ceil( ( i + j )/2 )  ) > 1000 
        if time( floor( ( i + j )/2 )  ) <= 1000
            index = floor( ( i + j )/2 ) ;
        else
            func( time , i , floor( ( i + j )/2 )  );
        end 
    else 
        func( time , ceil( ( i + j )/2 )  , j );
    end 
end