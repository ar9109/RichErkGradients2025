function print_png(f,outfilename)

set(f, 'PaperUnits','centimeters');
set(f, 'Units','centimeters');
pos=get(f,'Position');
set(f, 'PaperSize', [pos(3) pos(4)]);
set(f, 'PaperPositionMode', 'manual');
set(f, 'PaperPosition',[0 0 pos(3) pos(4)]);
print('-dpng',outfilename);

end

