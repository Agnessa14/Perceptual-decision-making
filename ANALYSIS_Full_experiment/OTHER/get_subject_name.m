function subname = get_subject_name(subject)
%GET_SUBJECT_NAME Transform the subject ID into a usable format.
%
%Input: subject ID (integer); e.g., 3
%
%Output: subject name (string); e.g, '03'

if subject < 10
    subname = ['0',num2str(subject)];
else
    subname = num2str(subject);
end 

end