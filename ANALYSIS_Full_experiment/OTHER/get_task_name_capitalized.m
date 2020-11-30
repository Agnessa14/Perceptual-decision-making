function task_name = get_task_name_capitalized(task_num)
%TASK_NAME Transform task number into a string task name. 
%
%Input: Task number: 1 for categorization or 2 for fixation
%
%Output: Task name: 'categorization' or 'fixation'
%
if task_num == 1
    task_name = 'Categorization';
elseif task_num == 2
    task_name = 'Fixation';
end

end