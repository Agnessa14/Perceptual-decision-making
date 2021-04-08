function tests = test_helper_functions
addpath(genpath('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/'));
tests = functiontests(localfunctions);
end

function test_get_subject_name(testCase)
actual_solution = get_subject_name(5);
expected_solution = '05';
verifyEqual(testCase,actual_solution,expected_solution);
end

function test_get_task_name(testCase)
actual_solution = get_task_name(1);
expected_solution = 'categorization';
verifyEqual(testCase,actual_solution,expected_solution);
end

% function test_create_data_matrix(testCase)
%initialize variables
% numConditions = 60;
% minNumTrials = 25;
% load('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/OTHER/Testing/triggers.mat'); %randomly generated triggers
% load('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/OTHER/Testing/data_random.mat'); %randomly generated data
% load('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/OTHER/Testing/data_organized.mat'); %data organized by conditions 
% 
% actual_solution = create_data_matrix(numConditions,triggers,minNumTrials,data_random);
% expected_solution = data_organized;
% verifyEqual(testCase,actual_solution,expected_solution);

%%%
% numConditions = 60;
% numElectrodes = 63;
% numTimepoints = 200;
% numAllTrials = 30;
% minNumTrials = 25;
% trialsMat_randomized = NaN(numConditions,numAllTrials);
% data_random = rand(numConditions*minNumTrials,numElectrodes,numTimepoints);
% data_organized = NaN(numConditions,minNumTrials,numElectrodes,numTimepoints);
% triggers = repmat(1:numConditions,1,numAllTrials);
% triggers = triggers(randperm(size(triggers,2)));
% save('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/OTHER/Testing/triggers','triggers');
% for c = 1:numConditions
%     trialsMat = find(triggers==c);
%     trialsMat_randomized(c,:) = trialsMat(randperm(numel(trialsMat))); %randomize
%     save('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/OTHER/Testing/trialsMat_randomized','trialsMat_randomized');
%     if numel(trialsMat_randomized(c,:))<minNumTrials
%         continue;
%     else
%         trialsMat_minimum = trialsMat_randomized(c,1:minNumTrials);
%         data_organized(c,:,:,:) = data_random(trialsMat_minimum,:,:); 
%     end
% end   
% save('/home/agnek95/SMST/PDM_PILOT_2/ANALYSIS_Full_experiment/OTHER/Testing/data_organized','data_organized');
%     
% end
