try
    disp('NBP of water is:')
    CoolProp.PropsSI('T','P',101325,'Q',0,'Water')
    exit
catch err
	error('err')
end