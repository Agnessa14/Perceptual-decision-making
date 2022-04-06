function task_name = get_task_name(task_num)
%TASK_NAME Transform task number into a string task name. 
%
%Input: Task number: 1 for categorization or 2 for fixation
%
%Output: Task name: 'categorization' or 'fixation'
%
if task_num == 1
    task_name = 'categorization';
elseif task_num == 2
    task_name = 'fixation';
end

end