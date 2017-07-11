function DrawContoursAndMassCenters(atestImage, aPsi, a_sz_i)

   imagesc(uint8(atestImage)); axis('off'); %axis('image');
   hold on;
   %contour(aPsi,[0 0],'r');
   contour(aPsi,'LineWidth',3,'LineColor',[1 0 0],'LevelList',0);
   hold off;
   drawnow;
end