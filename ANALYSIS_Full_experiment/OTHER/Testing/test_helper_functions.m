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

