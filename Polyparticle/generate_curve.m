%Generates a moving curve sequence and saves as a series of tiff files

x = -1000:1000;

j = 100;

figure;

for i = 0 : 0.005 : 0.2
    
    plot(x,x.^(3+i),'k-','linewidth',3)
    ylim([-4*10^9 4*10^9])
    axis off
    
    print(['D:\dkenwright\Tracking\test_images\moving_curve3\' num2str(j) '.tif'], '-dtiffn')
    j=j+1;
end