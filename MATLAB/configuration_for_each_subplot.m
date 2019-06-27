colorbar;
caxis([0, 0.45]);
xticks(1:7);
xticklabels({'1', '2', '3', '4', '5', '6', '7'});
xlabel('Threshold 2');
yticks(2:2:10);
yticklabels({'6', '8',  '10',  '12', '14'});
ylabel('Threshold 1');
set(gca, 'YDir', 'normal');