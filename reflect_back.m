%% Mirror Dual Field back 

function dGdV_xr = reflect_back(dGdV)

    two_r = size(dGdV,2);
   % dGdV_xr = zeros(size(dGdV,1),round(two_r/2),size(dGdV,3));
    
    dGdV_xr = dGdV(:,(round(two_r/2):end));
    dGdV_xr(:,2:end) = dGdV_xr(:,2:end) + dGdV(:,(floor(two_r/2):-1:1));

end 