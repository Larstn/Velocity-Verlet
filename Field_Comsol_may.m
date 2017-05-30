Nx = 50;
Ny = 50;
Nz = 50;

E_x = zeros(Nx, Ny, Nz);
E_y = E_x; 
E_z = E_x; 

cnt = 0;
for ii = 1:50
    for jj = 1:50
        for kk = 1:50
            cnt = cnt + 1; 
            E_x(kk,jj,ii) = GeoMayEall(cnt,4); 
            E_y(kk,jj,ii) = GeoMayEall(cnt,5); 
            E_z(kk,jj,ii) = GeoMayEall(cnt,6);


        end 
    end 
end 
%%
figure(1001)
imagesc(E_x(:,:,26)')

figure(1002)
imagesc(E_y(:,:,26)')

figure(1003)
imagesc(E_z(:,:,26)')

%%