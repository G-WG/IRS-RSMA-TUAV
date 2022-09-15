function saveFiguresMultipleFormats(namePlot, filePath)

%% saveFiguresMultipleFormats(namePlot, filePath)
% namePlot = 'surfWSR_TUAV-0m_Cityscenery_RSMA_IC_AMBF_MRT';
% filePath = '.\IRS-RSMA-TUAV\plots\WSR_TUAV-ALtitude\v2\';
fileExt = {'.png', '.fig', '.eps'};
for fE = 1:length(fileExt)
    saveas(gcf, [filePath, namePlot, fileExt{fE}])
end

end