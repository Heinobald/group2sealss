density0 = sw_dens0(SAL, TEMP);
density = sw_dens(SAL, TEMP, Z);

o_pres = ncread('ct34-2447-08_prof.nc', 'PRES');
o_psal = ncread('ct34-2447-08_prof.nc', 'PSAL');

figure(1);
subplot(2,2,1);
histogram(density);
subplot(2,2,2);
histogram(density0);
subplot(2,2,3);
plot(o_psal(:,466), -o_pres(:,466));
subplot(2,2,4);
scatter(reshape(SAL.',1,[]), reshape(TEMP.',1,[]))
xlim([33.5,35])

figure(2);
subplot(1,2,1);
newlat = reshape(repmat(LAT,57,1),1,[])
scatter(SAL(3,:), TEMP(3,:), 2, LAT);
title('colors encode latitude')
xlim([33.5,35])
subplot(1,2,2);
newtime = reshape(repmat(DATE,57,1),1,[])
scatter(reshape(SAL.',1,[]), reshape(TEMP.',1,[]), 2, newtime);
xlim([33.5,35])
title('colors encode time')
    

%scatter(TEMP, SAL);
%contourf(density0, 20);
%caxis([1027,1030]) 
%imshow(density0);