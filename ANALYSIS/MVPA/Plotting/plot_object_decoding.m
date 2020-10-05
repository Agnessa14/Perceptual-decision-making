function plot_object_decoding(taskType,subjects)
%PLOT_OBJECT_DECODING
%
%
%
%
%Get the task name
if taskType == 1
    taskName = 'categorization';  
elseif taskType == 2
    taskName = 'fixation';
end

%Loop: create a matrix containing the decoding accuracy matrices of each subject
for subject = subjects
    if subject < 10
        subNum = ['0',num2str(subject)];
    else
        subNum = num2str(subject);
    end
    cd(fullfile('/home/agnek95/SMST/PDM/RESULTS/',subNum)); %change directory to the current subject
    load(['mvnn_svm_decoding_accuracy_pdm_',taskName,'.mat']);
    
    %Random figure number
    figure(abs(round(randn*10)))
    
    %plot
    avg = squeeze(nanmean(nanmean(DA_end,1),2));
    plot(avg,'LineWidth',4);
    set(gcf, 'Position', get(0, 'Screensize'));
    
    %Define the domain
    xlim([0 201]);
    
    %Label of x
    xlabel('Timepoints');
    
    %Label of y
    ylabel('Decoding accuracy');
    
    %Title of plot
    if taskType == 1 
        title(['Classification accuracy of individual scenes per timepoint in a categorization task, subject ',subNum]);
    elseif taskType == 2
        title(['Classification accuracy of individual scenes per timepoint in a distracted task, subject ',subNum]);
    end
    
    % %Legend if needed
    legend('SVM accuracy');
    
    %Grid
    grid on
    
    %Save the plot
    saveas(gcf,['svm_object_pdm_',taskName]); %save as matlab figure
    saveas(gcf,['svm_object_pdm_',taskName,'.svg']); %save as svg
    close(gcf);
    
end

end